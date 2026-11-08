//! The optimizer

use super::*;

pub struct ConformerInput {
    pub model_index: usize,
    /// Report label for the model (Reduce2 prints -1 in compat mode).
    pub report_model: i64,
    /// (alt, atoms of that conformer in GetAtomsForConformer order)
    pub alts: Vec<(String, Vec<u32>)>,
    /// Atom i_seqs in the model (for restoring between alternates).
    pub model_atoms: Vec<u32>,
}

/// Movers and final states of the last run (for validation tools).
#[derive(Default)]
pub struct DebugMovers {
    pub movers: Vec<Mover>,
    pub final_cfg: Vec<(usize, i64)>,
    pub final_score: Vec<f64>,
}

pub fn optimize(w: &mut World, models: &[ConformerInput], rotatable_h: &[u32], p: &OptParams, n_models: usize) -> OptOutput {
    optimize_debug(w, models, rotatable_h, p, n_models, None)
}

pub fn optimize_debug(
    w: &mut World,
    models: &[ConformerInput],
    rotatable_h: &[u32],
    p: &OptParams,
    n_models: usize,
    mut debug: Option<&mut DebugMovers>,
) -> OptOutput {
    let mut out = OptOutput::default();
    let flip_states = parse_flip_states(&p.flip_states, p.compat);
    let mut all_deletes: FxHashSet<u32> = FxHashSet::default();
    for mi in models {
        // only a later alternate starts again from these
        let initial_pos: Vec<(u32, Vec3, AtomInfo)> = if mi.alts.len() > 1 {
            mi.model_atoms.iter().map(|&a| (a, w.pos[a as usize], w.info[a as usize])).collect()
        } else {
            Vec::new()
        };
        let mut first = true;
        let mut model_deletes: FxHashSet<u32> = FxHashSet::default();
        for (alt, conf_atoms) in &mi.alts {
            if !first {
                for &(a, pos, inf) in &initial_pos {
                    let al = w.altloc(a);
                    if al.is_empty() || al == " " || al == alt {
                        w.pos[a as usize] = pos;
                        w.info[a as usize] = inf;
                    }
                }
            }
            first = false;
            let old_bonds = if p.compat { Vec::new() } else { keep_bonds_within(w, conf_atoms) };
            let dels = run_one(w, p, mi.report_model, alt, conf_atoms.clone(), rotatable_h, &flip_states, n_models, &mut out, debug.as_deref_mut());
            if p.compat {
                model_deletes = dels;
            } else {
                // Keep earlier alternates' decisions for their own atoms; shared
                // atoms follow the last alternate processed.
                model_deletes.retain(|&a| {
                    let al = w.altloc(a);
                    !(al.is_empty() || al == " " || al == alt)
                });
                model_deletes.extend(dels);
            }
            w.clear_phantoms();
            for (a, bonds) in old_bonds {
                w.bonded[a as usize] = bonds;
            }
        }
        if p.compat {
            all_deletes = model_deletes;
        } else {
            all_deletes.extend(model_deletes);
        }
    }
    if out.num_calculated > 0 {
        out.info += &vcheck(
            p.verbosity,
            1,
            &format!(
                "Calculated : cached atom scores: {} : {}; fraction calculated {:.2}\n",
                out.num_calculated,
                out.num_cached,
                out.num_calculated as f64 / (out.num_calculated + out.num_cached) as f64
            ),
        );
    }
    let mut d: Vec<u32> = all_deletes.into_iter().collect();
    d.sort_unstable();
    out.hydrogens_to_delete = d;
    out
}

/// Limit the bonded-neighbor lists to one conformer and return the full
/// lists: a neighbor outside the conformer is dropped when the conformer has
/// its own alternate of that atom. Reduce2 keeps every alternate's neighbors,
/// so a CB shared by two OG alternates seems to have six bonds and its OH gets
/// no rotator (nor do the methyls of an alternate threonine or valine, or the
/// ring of a histidine split at CG). A neighbor with no alternate in the
/// conformer stays (an HG kept in the main conformation of a serine whose only
/// OG is labeled A). Returns the lists it changed, with their old contents.
fn keep_bonds_within(w: &mut World, conf_atoms: &[u32]) -> Vec<(u32, Vec<u32>)> {
    let mut inside = vec![false; w.len()];
    for &a in conf_atoms {
        inside[a as usize] = true;
    }
    // a neighbor can be replaced only in a residue group with atoms outside
    // the conformer
    let split: FxHashSet<u32> = w.labels.iter().zip(&inside).filter(|x| !*x.1).map(|x| x.0.rg).collect();
    let mut present: FxHashSet<(u32, &str)> = FxHashSet::default();
    for &a in conf_atoms {
        let l = &w.labels[a as usize];
        if split.contains(&l.rg) {
            present.insert((l.rg, l.name.as_str()));
        }
    }
    let replaced = |m: u32| !inside[m as usize] && present.contains(&(w.labels[m as usize].rg, w.labels[m as usize].name.as_str()));
    let within: Vec<(u32, Vec<u32>)> = w
        .bonded
        .iter()
        .enumerate()
        .filter(|(_, l)| l.iter().any(|&m| replaced(m)))
        .map(|(a, l)| (a as u32, l.iter().copied().filter(|&m| !replaced(m)).collect()))
        .collect();
    within.into_iter().map(|(a, l)| (a, std::mem::replace(&mut w.bonded[a as usize], l))).collect()
}

/// Phantom hydrogens for a water oxygen (`getPhantomHydrogensFor`).
fn phantom_hydrogens_for(
    w: &World,
    atom: u32,
    q: &NeighborIndex,
    removed: &FxHashSet<u32>,
    min_occ: f64,
    radius: f64,
    dist: f64,
) -> Vec<Vec3> {
    let mut nearby = q.neighbors(w, w.pos[atom as usize], 0.001, 4.0, removed);
    nearby.sort_unstable();
    let mut cands: Vec<(u32, f64)> = Vec::new();
    let model = w.labels[atom as usize].model_index;
    for a in nearby {
        if !w.compatible(atom, a) {
            continue;
        }
        if w.labels[a as usize].model_index != model {
            continue;
        }
        let overlap = w.pos[atom as usize].dist(w.pos[a as usize]) - (radius + w.info[a as usize].vdw_radius + dist);
        if overlap <= -0.1 && w.occ[a as usize] >= min_occ && w.elem(a) != "H" {
            let mut skip = false;
            if is_aromatic_acceptor(w.resname(a), w.name(a)) {
                for c in cands.iter_mut() {
                    if is_aromatic_acceptor(w.resname(c.0), w.name(c.0)) && w.labels[a as usize].ag == w.labels[c.0 as usize].ag {
                        if overlap < c.1 {
                            *c = (a, overlap);
                        }
                        skip = true;
                        break;
                    }
                }
            }
            if !skip {
                cands.push((a, overlap));
            }
        }
    }
    let mut ret = Vec::new();
    const BEST_HBOND_OVERLAP: f64 = 0.6;
    for (c, ov) in cands {
        let d = dist + (-dist).max((0.0f64).min(ov + BEST_HBOND_OVERLAP));
        let off = w.pos[c as usize] - w.pos[atom as usize];
        let l = off.length();
        if l == 0.0 {
            continue;
        }
        ret.push(w.pos[atom as usize] + off / l * d);
    }
    ret
}

#[allow(clippy::too_many_arguments)]
fn run_one(
    w: &mut World,
    p: &OptParams,
    report_model: i64,
    alt: &str,
    mut atoms: Vec<u32>,
    rotatable_h: &[u32],
    flip_states: &[FlipMoverState],
    n_models: usize,
    out: &mut OptOutput,
    debug: Option<&mut DebugMovers>,
) -> FxHashSet<u32> {
    let v = p.verbosity;
    let mut info = String::new();
    info += &vcheck(v, 1, &format!("Running Reduce optimization on model index {}, alternate '{}'\n", report_model, alt));
    info += &vcheck(v, 1, &format!("  bondedNeighborDepth = {}\n", p.bonded_neighbor_depth));
    info += &vcheck(v, 1, &format!("  probeRadius = {}\n", py_float(p.probe.probe_radius)));
    info += &vcheck(v, 1, &format!("  useNeutronDistances = {}\n", if p.use_neutron_distances && !p.compat { "True" } else { "False" }));
    info += &vcheck(v, 1, &format!("  probeDensity = {}\n", py_float(p.probe.density)));
    info += &vcheck(v, 1, &format!("  minOccupancy = {}\n", py_float(p.min_occupancy)));
    info += &vcheck(v, 1, &format!("  preferenceMagnitude = {}\n", py_float(p.preference_magnitude)));
    let mut tm = Timer::new(v);
    info += &tm.report("construct spatial query");

    let mut max_vdw: f64 = 1.0;
    for &a in &atoms {
        max_vdw = max_vdw.max(w.info[a as usize].vdw_radius);
    }

    // ---------------- water phantom hydrogens
    let neutron = p.use_neutron_distances && !p.compat;
    let (ph_rad, ph_dist) = if neutron { (1.0, 0.98) } else { (1.05, 0.84) };
    let mut removed: FxHashSet<u32> = FxHashSet::default();
    let mut nidx = NeighborIndex::new(w, &atoms, p.compat);
    // before the low-quality waters are dropped: they keep a hydrogen out too
    let clearance = if p.compat { None } else { movers::Clearance::new(w, &atoms, p.rotor_clearance) };
    let is_bad_water = |w: &World, a: u32| -> Option<bool> {
        if w.elem(a) == "O" && w.labels[a as usize].is_water {
            Some(!(w.occ[a as usize] >= 0.66 && w.b[a as usize] < 40.0))
        } else {
            None
        }
    };
    let mut waters_to_delete = Vec::new();
    if !p.compat {
        // Fixed: low-quality waters are dropped before phantoms are aimed.
        for &a in &atoms {
            if is_bad_water(w, a) == Some(true) {
                removed.insert(a);
            }
        }
    }
    let mut phantoms: Vec<(Vec3, u32)> = Vec::new();
    for &a in &atoms {
        match is_bad_water(w, a) {
            Some(false) => {
                let bonded = !w.bonded[a as usize].is_empty();
                if bonded {
                    info += &vcheck(v, 3, &format!("Not adding phantom Hydrogens on {} due to bonded neighbors\n", w.res_name_and_id(a)));
                }
                w.info[a as usize].is_donor = false;
                w.info[a as usize].is_acceptor = true;
                if bonded && !p.compat {
                    continue;
                }
                let newp = phantom_hydrogens_for(w, a, &nidx, &removed, p.min_occupancy, ph_rad, ph_dist);
                if !newp.is_empty() {
                    info += &vcheck(v, 3, &format!("Added {} phantom Hydrogens on {}\n", newp.len(), w.res_name_and_id(a)));
                    for ph in newp {
                        phantoms.push((ph, a));
                    }
                }
            }
            Some(true) => {
                let l = &w.labels[a as usize];
                info += &vcheck(
                    v,
                    3,
                    &format!(
                        "Ignoring {} {} {} {} with occupancy {} and B factor {}\n",
                        l.name, l.resname, l.resseq, l.chain, py_float(w.occ[a as usize]), py_float(w.b[a as usize])
                    ),
                );
                waters_to_delete.push(a);
            }
            None => {}
        }
    }
    if !waters_to_delete.is_empty() {
        info += &vcheck(v, 1, &format!("Ignored {} waters due to occupancy or B factor\n", waters_to_delete.len()));
        let set: FxHashSet<u32> = waters_to_delete.iter().copied().collect();
        atoms.retain(|a| !set.contains(a));
        for &a in &waters_to_delete {
            nidx.remove(w, a);
        }
        removed.extend(set);
    }
    let first_phantom = atoms.len();
    if !phantoms.is_empty() {
        let orig = atoms.len();
        for (pp, parent) in &phantoms {
            let alt_c = AtomInfo { vdw_radius: ph_rad, is_acceptor: false, is_donor: true, is_dummy_hydrogen: true, is_ion: false, charge: 0, alt: b' ' };
            let id = w.add_phantom(*parent, *pp, alt_c);
            atoms.push(id);
            nidx.add(w, id);
        }
        info += &vcheck(v, 1, &format!("Added {} phantom Hydrogens on waters", phantoms.len()));
        info += &vcheck(v, 1, &format!(" (Old total {}, new total {})\n", orig, atoms.len()));
    }
    info += &tm.report("place water phantom Hydrogens");
    fixup_explicit_donors(w, &atoms);
    info += &tm.report("fixup explicit doners");

    // ---------------- Mover placement
    if !p.compat {
        // include the phantoms (they count as potential touches); no atom has
        // moved, and the removed waters are in `removed`
        nidx.add_late(w, &atoms[first_phantom..]);
    }
    let mut placement = Placement {
        movers: Vec::new(),
        mover_info: Vec::new(),
        delete_atoms: Vec::new(),
        info: String::new(),
        amide_flips: Vec::new(),
        his_flips: Vec::new(),
        docked_alpha_carbons: Vec::new(),
    };
    place_movers(w, p, &atoms, rotatable_h, flip_states, alt, max_vdw, &mut nidx, &removed, n_models, clearance.as_ref(), &mut placement);
    info += &placement.info;
    let movers = placement.movers;
    let mut mover_info = placement.mover_info;
    out.amide_flips.extend(placement.amide_flips);
    out.his_flips.extend(placement.his_flips);
    out.docked_alpha_carbons.extend(placement.docked_alpha_carbons);
    info += &vcheck(v, 1, &format!("Inserted {} Movers\n", movers.len()));
    info += &vcheck(v, 1, &format!("Marked {} atoms for deletion\n", placement.delete_atoms.len()));
    info += &tm.report("place movers");
    let placement_deletes: FxHashSet<u32> = placement.delete_atoms.iter().copied().collect();

    // Initialize Movers to coarse state 0 (positions, info, deletions)
    let mut state_deleted: FxHashMap<u32, bool> = FxHashMap::default();
    for mv in &movers {
        apply_coarse(w, mv, 0, &mut state_deleted);
    }
    info += &tm.report("initialize Movers");

    // ---------------- interaction graph
    let nm = movers.len();
    let pr = p.probe.probe_radius;
    let mut atom_movers: FxHashMap<u32, SmallVec<[u32; 4]>> = FxHashMap::default();
    for (mi, mv) in movers.iter().enumerate() {
        for &a in &mv.atoms {
            let e = atom_movers.entry(a).or_default();
            if !e.contains(&(mi as u32)) {
                e.push(mi as u32);
            }
        }
    }
    let mut bbox: Vec<(Vec3, Vec3)> = Vec::with_capacity(nm);
    for mv in &movers {
        let mut lo = v3(1e10, 1e10, 1e10);
        let mut hi = v3(-1e10, -1e10, -1e10);
        for ps in mv.all_poses() {
            for (k, pt) in ps.iter().enumerate() {
                let r = pr + w.info[mv.atoms[k] as usize].vdw_radius;
                lo = v3(lo.x.min(pt.x - r), lo.y.min(pt.y - r), lo.z.min(pt.z - r));
                hi = v3(hi.x.max(pt.x + r), hi.y.max(pt.y + r), hi.z.max(pt.z + r));
            }
        }
        bbox.push((lo, hi));
    }
    // per Mover atom, a sphere holding it in every pose (None when no pose
    // has it); two atoms whose spheres lie farther apart than their reach
    // meet in no pose
    let pose_spheres: Vec<Vec<Option<(Vec3, f64)>>> = movers
        .iter()
        .map(|mv| {
            (0..mv.atoms.len())
                .map(|k| {
                    let mut lo = v3(f64::INFINITY, f64::INFINITY, f64::INFINITY);
                    let mut hi = -lo;
                    for q in mv.all_poses().filter_map(|ps| ps.get(k)) {
                        lo = v3(lo.x.min(q.x), lo.y.min(q.y), lo.z.min(q.z));
                        hi = v3(hi.x.max(q.x), hi.y.max(q.y), hi.z.max(q.z));
                    }
                    if lo.x > hi.x {
                        return None;
                    }
                    let c = (lo + hi) * 0.5;
                    Some((c, mv.all_poses().filter_map(|ps| ps.get(k)).map(|q| q.dist(c)).fold(0.0, f64::max)))
                })
                .collect()
        })
        .collect();
    let mut edges: Vec<(usize, usize)> = Vec::new();
    {
        // candidate pairs by bounding-box overlap using a sweep on x
        let mut order: Vec<usize> = (0..nm).collect();
        order.sort_by(|&a, &b| bbox[a].0.x.partial_cmp(&bbox[b].0.x).unwrap());
        let mut cand: Vec<(usize, usize)> = Vec::new();
        for (ii, &i) in order.iter().enumerate() {
            for &j in &order[ii + 1..] {
                if bbox[j].0.x > bbox[i].1.x {
                    break;
                }
                let (a, b) = (&bbox[i], &bbox[j]);
                if a.0.y <= b.1.y && a.1.y >= b.0.y && a.0.z <= b.1.z && a.1.z >= b.0.z {
                    cand.push((i.min(j), i.max(j)));
                }
            }
        }
        cand.sort_unstable();
        let results: Vec<(usize, usize, Vec<(u32, u32)>)> = crate::par::map_collect(&cand, |&(i, j)| {
                let mut entries: Vec<(u32, u32)> = Vec::new();
                let mi = &movers[i];
                let mj = &movers[j];
                let ri: Vec<f64> = mi.atoms.iter().map(|&a| w.info[a as usize].vdw_radius).collect();
                let rj: Vec<f64> = mj.atoms.iter().map(|&a| w.info[a as usize].vdw_radius).collect();
                let mut hit_i = vec![false; mi.atoms.len()];
                let mut hit_j = vec![false; mj.atoms.len()];
                // whether any poses bring the two atoms within reach; a pair
                // already counted for both atoms, or whose pose spheres keep
                // it out of reach, needs no search
                for a1 in 0..mi.atoms.len() {
                    let lim1 = 2.0 * pr + ri[a1];
                    for a2 in 0..mj.atoms.len() {
                        if hit_i[a1] && hit_j[a2] {
                            continue;
                        }
                        let lim = lim1 + rj[a2];
                        let (Some((c1, r1)), Some((c2, r2))) = (pose_spheres[i][a1], pose_spheres[j][a2]) else { continue };
                        let apart = r1 + r2 + lim + POSE_SPHERE_SLACK;
                        if c1.dist_sq(c2) > apart * apart {
                            continue;
                        }
                        let hit = mi.all_poses().filter_map(|p1| p1.get(a1)).any(|x1| {
                            mj.all_poses().filter_map(|p2| p2.get(a2)).any(|x2| x1.dist_sq(*x2) <= lim * lim)
                        });
                        if hit {
                            hit_i[a1] = true;
                            hit_j[a2] = true;
                        }
                    }
                }
                for (k, &h) in hit_i.iter().enumerate() {
                    if h {
                        entries.push((mi.atoms[k], j as u32));
                    }
                }
                for (k, &h) in hit_j.iter().enumerate() {
                    if h {
                        entries.push((mj.atoms[k], i as u32));
                    }
                }
                (i, j, entries)
            });
        for (i, j, entries) in results {
            if !entries.is_empty() {
                edges.push((i, j));
                for (a, m) in entries {
                    let e = atom_movers.entry(a).or_default();
                    if !e.contains(&m) {
                        e.push(m);
                    }
                }
            }
        }
    }
    // connected components ordered by first Mover index, members in order
    let mut parent: Vec<usize> = (0..nm).collect();
    fn find(p: &mut Vec<usize>, x: usize) -> usize {
        let mut r = x;
        while p[r] != r {
            r = p[r];
        }
        let mut y = x;
        while p[y] != r {
            let n = p[y];
            p[y] = r;
            y = n;
        }
        r
    }
    for &(i, j) in &edges {
        let a = find(&mut parent, i);
        let b = find(&mut parent, j);
        if a != b {
            parent[a.max(b)] = a.min(b);
        }
    }
    let mut comp_index: FxHashMap<usize, usize> = FxHashMap::default();
    let mut components: Vec<Vec<u32>> = Vec::new();
    for i in 0..nm {
        let r = find(&mut parent, i);
        let ci = *comp_index.entry(r).or_insert_with(|| {
            components.push(Vec::new());
            components.len() - 1
        });
        components[ci].push(i as u32);
    }
    let singles: Vec<usize> = (0..components.len()).filter(|&c| components[c].len() == 1).collect();
    let groups: Vec<usize> = (0..components.len()).filter(|&c| components[c].len() > 1).collect();
    let max_len = components.iter().map(|c| c.len()).max().unwrap_or(0);
    info += &vcheck(
        v,
        1,
        &format!(
            "Found {} Cliques ({} are singletons); largest Clique size = {}\n",
            components.len(),
            singles.len(),
            max_len
        ),
    );
    info += &tm.report("compute interaction graph");

    // ---------------- excluded atoms and dots
    let mut mover_atoms: Vec<u32> = Vec::new();
    {
        let mut seen = FxHashSet::default();
        for mv in &movers {
            for &a in &mv.atoms {
                if seen.insert(a) {
                    mover_atoms.push(a);
                }
            }
        }
    }
    let depth = p.bonded_neighbor_depth;
    let exclude: FxHashMap<u32, Vec<u32>> =
        crate::par::map_collect(&mover_atoms, |&a| (a, atoms_within_n_bonds(w, a, pr, depth, 3)));
    info += &tm.report("determine excluded atoms");
    let mut cache = DotSphereCache::new(p.probe.density);
    let mut dots: FxHashMap<u32, Arc<Vec<DotPair>>> = FxHashMap::default();
    let mut dot_reach: FxHashMap<u32, f64> = FxHashMap::default();
    let mut dots_soa: FxHashMap<u32, Arc<DotsSoa>> = FxHashMap::default();
    let mut soa_by_radius: FxHashMap<u64, Arc<DotsSoa>> = FxHashMap::default();
    let mut by_radius: FxHashMap<u64, Arc<Vec<DotPair>>> = FxHashMap::default();
    for &a in &mover_atoms {
        let r = w.info[a as usize].vdw_radius;
        let dp = by_radius
            .entry(r.to_bits())
            .or_insert_with(|| {
                let s = cache.get(r);
                Arc::new(
                    s.iter()
                        .map(|&d| {
                            let l = d.length();
                            DotPair { d, probe: if l > 0.0 { d / l * (l + pr) } else { Vec3::ZERO } }
                        })
                        .collect(),
                )
            })
            .clone();
        dot_reach.insert(a, dp.iter().map(|x| x.d.length()).fold(0.0, f64::max));
        let soa = soa_by_radius.entry(r.to_bits()).or_insert_with(|| Arc::new(DotsSoa::new(&dp))).clone();
        dots_soa.insert(a, soa);
        dots.insert(a, dp);
    }
    info += &tm.report("construct dot scorer");

    // static targets: run atoms (incl. phantoms) that do not move with a Mover
    let mut dyn_of = DynOf(vec![(NONE, 0); w.len()]);
    for (mi, mv) in movers.iter().enumerate() {
        for slot in 0..mv.n_moved {
            dyn_of.insert(mv.atoms[slot], (mi as u32, slot as u16));
        }
    }
    let mut static_pts: Vec<(u32, Vec3)> = Vec::new();
    for &a in &atoms {
        if dyn_of.contains_key(&a) {
            continue;
        }
        if !p.compat && placement_deletes.contains(&a) {
            continue;
        }
        static_pts.push((a, w.pos[a as usize]));
    }
    let grid = if p.compat {
        let present: FxHashSet<u32> = static_pts.iter().map(|x| x.0).collect();
        nidx.static_grid(w, |a| present.contains(&a)).unwrap()
    } else {
        SpatialGrid::new(&static_pts, 2.0)
    };
    let snapshot = if p.compat { Some(w.info.clone()) } else { None };
    let mut local = vec![NONE; nm];
    for comp in &components {
        for (k, &m) in comp.iter().enumerate() {
            local[m as usize] = k as u32;
        }
    }
    info += &tm.report("construct OptimizerC");
    // static grid entries near each Mover atom over all of its positions;
    // fixed mode only: compat's grid files atoms by another position than it
    // measures (Reduce2's query misses some), so a filtered superset differs
    let mut near: FxHashMap<u32, (Vec3, f64, u32, u32)> = FxHashMap::default();
    let mut near_entries: Vec<(Vec3, AtomInfo)> = Vec::new();
    for mv in movers.iter().filter(|_| !p.compat) {
        for slot in 0..mv.n_moved {
            // the atom's positions over all coarse and fine poses
            let poses = || mv.coarse_pos.iter().filter_map(move |c| c.get(slot)).chain(mv.fine.iter().filter_map(move |f| f.get(slot)));
            let (sum, n) = poses().fold((Vec3::ZERO, 0usize), |(acc, n), &q| (acc + q, n + 1));
            if n == 0 {
                continue;
            }
            let a = mv.atoms[slot];
            let c = sum / n as f64;
            // the root of the largest square is the largest distance
            let spread = poses().map(|q| q.dist_sq(c)).fold(0.0, f64::max).sqrt();
            let r = max_vdw + w.info[a as usize].vdw_radius + 2.0 * p.probe.probe_radius + spread + NEAR_MARGIN;
            let excl: &[u32] = exclude.get(&a).map(|v| v.as_slice()).unwrap_or(&[]);
            let lo = near_entries.len() as u32;
            grid.for_each_within(c, 0.0, r, |b, pb, _| {
                if !excl.contains(&b) && w.occ[b as usize].abs() >= p.min_occupancy {
                    near_entries.push((pb, w.info[b as usize]));
                }
            });
            near.insert(a, (c, r, lo, near_entries.len() as u32));
        }
    }

    let ctx = Ctx {
        w: &*w,
        movers: &movers,
        p,
        scorer: DotScorer::new(&p.probe),
        snapshot,
        grid,
        dyn_of,
        atom_movers,
        exclude,
        dots,
        dots_soa,
        dot_reach,
        near,
        near_entries,
        max_vdw,
        local,
    };

    // ---------------- optimize every clique (in parallel)
    let comp_results: Vec<CliqueResult> = if p.compat {
        let mut ident_ctx = ctx;
        ident_ctx.local = (0..nm as u32).collect();
        compat_optimize(&ident_ctx, &components, &edges, v)
    } else {
        crate::par::map_collect(&components, |comp| optimize_clique(&ctx, comp))
    };

    // initial scores into mover info
    for (ci, comp) in components.iter().enumerate().filter(|_| p.report) {
        for (k, &m) in comp.iter().enumerate() {
            mover_info[m as usize] += &format!(" Initial score: {:.2}", comp_results[ci].initial[k]);
        }
    }
    for &ci in singles.iter().filter(|_| p.report) {
        info += &comp_results[ci].coarse_info;
        info += &vcheck(v, 1, &format!("Singleton optimized with score {:.2}\n", comp_results[ci].coarse_best));
    }
    info += &tm.report("optimize singletons (coarse)");
    for &ci in groups.iter().filter(|_| p.report) {
        info += &comp_results[ci].coarse_info;
        info += &vcheck(v, 1, &format!("Clique optimized with score {:.2}\n", comp_results[ci].coarse_best));
    }
    info += &tm.report("optimize cliques (coarse)");
    info += &vcheck(v, 1, "Fine optimization on all Movers\n");
    // fine-phase messages in global Mover order
    let mut fine_msgs: Vec<(u32, String)> = Vec::new();
    for r in comp_results.iter().filter(|_| p.report) {
        fine_msgs.extend(r.fine_info.iter().cloned());
    }
    fine_msgs.sort_by_key(|x| x.0);
    for (_, s) in fine_msgs {
        info += &s;
    }
    info += &tm.report("optimize all Movers (fine)");

    // ---------------- report
    let mut final_cfg: Vec<Cfg> = vec![Cfg::Coarse(0); nm];
    let mut high: Vec<f64> = vec![0.0; nm];
    let mut pose_suffix: Vec<String> = vec![String::new(); nm];
    for (ci, comp) in components.iter().enumerate() {
        for (k, &m) in comp.iter().enumerate() {
            final_cfg[m as usize] = comp_results[ci].cfg[k];
            high[m as usize] = comp_results[ci].high[k];
            if let Some(a) = comp_results[ci].annot.get(k) {
                pose_suffix[m as usize] = a.clone();
            }
        }
    }
    if p.report {
        let describe = |m: usize| -> String {
            let mv = &movers[m];
            let cfg = final_cfg[m];
            let fine = match cfg {
                Cfg::Fine(_, f) => Some(f as usize),
                Cfg::Coarse(_) => None,
            };
            let mut d = mv.pose_description(cfg.coarse(), fine, !p.skip_bond_fixup);
            d += &pose_suffix[m];
            format!("  {} final score: {:.2} pose {}\n", mover_info[m], high[m], d)
        };
        info += &vcheck(v, 1, &format!("BEGIN REPORT: Model {} Alt '{}':\n", report_model, alt));
        let mut sorted_groups = groups.clone();
        sorted_groups.sort_by(|&a, &b| components[b].len().cmp(&components[a].len()));
        for &ci in &sorted_groups {
            let comp = &components[ci];
            info += &vcheck(v, 1, &format!(" Set of {} Movers:", comp.len()));
            let mut initial = 0.0;
            let mut fin = 0.0;
            for (k, &m) in comp.iter().enumerate() {
                if p.compat {
                    // The original re-parses the printed (rounded) initial score.
                    initial += format!("{:.2}", comp_results[ci].initial[k]).parse::<f64>().unwrap();
                } else {
                    initial += comp_results[ci].initial[k];
                }
                fin += high[m as usize];
            }
            info += &vcheck(v, 1, &format!(" Totals: initial score {:.2}, final score {:.2}\n", initial, fin));
            for &m in comp {
                info += &vcheck(v, 1, &describe(m as usize));
            }
        }
        info += &vcheck(v, 1, " Singleton Movers:\n");
        for &ci in &singles {
            let m = components[ci][0] as usize;
            info += &vcheck(v, 1, &describe(m));
        }
        info += &vcheck(v, 1, "END REPORT\n");
    }

    // ---------------- apply final states and fix-ups to the world
    let mut deletes: FxHashSet<u32> = placement_deletes.clone();
    for (m, mv) in movers.iter().enumerate() {
        let cfg = final_cfg[m];
        let c = cfg.coarse();
        for slot in 0..mv.n_moved {
            let a = mv.atoms[slot];
            w.pos[a as usize] = match cfg {
                Cfg::Coarse(c) => mv.coarse_pos[c as usize][slot],
                Cfg::Fine(c, f) => mv.fine.pos(c as usize, f as usize)[slot],
            };
        }
        for (slot, inf) in in_state(&mv.coarse_info, c).iter().enumerate() {
            w.info[mv.atoms[slot] as usize] = *inf;
        }
        for (slot, &d) in in_state(&mv.coarse_del, c).iter().enumerate() {
            if d {
                deletes.insert(mv.atoms[slot]);
            } else {
                deletes.remove(&mv.atoms[slot]);
            }
        }
    }
    let _ = state_deleted;
    if !p.skip_bond_fixup {
        info += &vcheck(v, 1, "FixUp on all Movers\n");
        for (m, mv) in movers.iter().enumerate() {
            let loc = final_cfg[m].coarse();
            if p.report {
                info += &vcheck(v, 3, &format!("FixUp on {} coarse location {}\n", mover_info[m], loc));
            }
            for (i, &pp) in in_state(&mv.fixup_pos, loc).iter().enumerate() {
                w.pos[mv.atoms[i] as usize] = pp;
            }
            if !in_state(&mv.fixup_pos, loc).is_empty() {
                out.docked_alpha_carbons.extend(mv.alpha_carbon);
            }
            for (i, inf) in in_state(&mv.fixup_info, loc).iter().enumerate() {
                w.info[mv.atoms[i] as usize] = *inf;
            }
            for (i, &d) in in_state(&mv.fixup_del, loc).iter().enumerate() {
                if d {
                    deletes.insert(mv.atoms[i]);
                } else {
                    deletes.remove(&mv.atoms[i]);
                }
            }
        }
        info += &tm.report("fix up Movers");
    }
    for (m, mv) in movers.iter().enumerate() {
        let c = final_cfg[m].coarse();
        match mv.kind {
            MoverKind::AmideFlip if c == 1 => out.amide_flips.push(FlippedMoverInfo { alt: alt.into(), base_atom: mv.atoms[0] }),
            MoverKind::HisFlip { .. } if c == 4 => out.his_flips.push(FlippedMoverInfo { alt: alt.into(), base_atom: mv.atoms[0] }),
            _ => {}
        }
    }
    for r in &comp_results {
        out.num_calculated += r.calculated;
        out.num_cached += r.cached;
    }
    if let Some(dbg) = debug {
        dbg.final_cfg = final_cfg
            .iter()
            .map(|c| match c {
                Cfg::Coarse(c) => (*c as usize, -1),
                Cfg::Fine(c, f) => (*c as usize, *f as i64),
            })
            .collect();
        dbg.final_score = high.clone();
        dbg.movers = movers.clone();
    }
    out.info += &info;
    deletes
}

fn apply_coarse(w: &mut World, mv: &Mover, c: usize, deleted: &mut FxHashMap<u32, bool>) {
    for slot in 0..mv.n_moved {
        w.pos[mv.atoms[slot] as usize] = mv.coarse_pos[c][slot];
    }
    for (slot, inf) in in_state(&mv.coarse_info, c).iter().enumerate() {
        w.info[mv.atoms[slot] as usize] = *inf;
    }
    for (slot, &d) in in_state(&mv.coarse_del, c).iter().enumerate() {
        deleted.insert(mv.atoms[slot], d);
    }
}

pub(super) struct CliqueResult {
    pub(super) initial: Vec<f64>,
    pub(super) cfg: Vec<Cfg>,
    pub(super) high: Vec<f64>,
    pub(super) coarse_best: f64,
    pub(super) coarse_info: String,
    pub(super) fine_info: Vec<(u32, String)>,
    pub(super) annot: Vec<String>,
    pub(super) calculated: usize,
    pub(super) cached: usize,
}

/// Score tables of one moving atom, decomposed by dot. A dot's result depends
/// on another Mover only if one of that Mover's atoms, in some state, can trim
/// the dot (an excluded atom covering it) or be a target for it. Dots that no
/// other Mover reaches are scored once per own state; the rest are grouped by
/// the set of Movers that reach them, and each group gets a table over just
/// those Movers. The per-dot evaluation is the same fold over the same targets
/// in the same order as `score_atom`, so the entries are the same sums, added
/// in a different order. Returns the factors and the numbers of dot groups
/// scored and of full atom scores they stand for.
/// Slack on the radius of a Mover atom's precomputed static neighborhood.
pub(super) const NEAR_MARGIN: f64 = 1e-3;

/// Slack on the distance beyond which an excluded atom cannot cover a dot
/// (far above rounding error).
pub(super) const EXCLUSION_MARGIN: f64 = 0.01;

/// Largest number of joint states of the other Movers one group of an atom's
/// dots is scored over, and largest table variable elimination may build. A
/// clique past either is optimized by coordinate ascent instead.
const DOT_GROUP_STATES_LIMIT: usize = 1 << 17;
/// Largest number of dot scorings (dots times joint states) one atom's factors
/// may take; past it the clique is optimized by coordinate ascent.
const ATOM_DOT_WORK_LIMIT: usize = 1 << 23;
const ELIMINATION_TABLE_LIMIT: usize = 1 << 22;
/// Slack on the pose spheres of the interaction graph (far above rounding
/// error).
const POSE_SPHERE_SLACK: f64 = 1e-6;

/// The atom's score factors, or None when a group of its dots depends on more
/// joint states than `DOT_GROUP_STATES_LIMIT`. When the atom keeps its place
/// and info in every own state and its own Mover offers it no targets, its
/// dots are sorted and scored once, and each own state only drops the dots
/// inside the own Mover's excluded atoms in that state.
fn atom_factors_dotwise(
    ctx: &Ctx,
    comp: &[u32],
    doms: &[usize],
    i: usize,
    slot: usize,
    abandoned: &std::sync::atomic::AtomicBool,
) -> Option<(Vec<Factor>, usize, usize)> {
    thread_local! {
        static BUFS: std::cell::RefCell<DotwiseBufs> = std::cell::RefCell::new(DotwiseBufs::default());
    }
    BUFS.with_borrow_mut(|bufs| atom_factors_dotwise_in(ctx, comp, doms, i, slot, abandoned, bufs))
}

/// The table of dependency set `dep` among `tables`, added as `size` zeros
/// when new.
fn table_for<'t>(tables: &'t mut Vec<(Vec<usize>, Vec<f64>)>, dep: &[usize], size: usize) -> &'t mut Vec<f64> {
    let k = match tables.iter().position(|t| t.0 == dep) {
        Some(k) => k,
        None => {
            tables.push((dep.to_vec(), vec![0.0; size]));
            tables.len() - 1
        }
    };
    &mut tables[k].1
}

/// Scratch space of `atom_factors_dotwise`, kept per thread between calls.
/// Each call clears what it uses, so nothing carries over.
#[derive(Default)]
struct DotwiseBufs {
    seq: Vec<(u32, Option<usize>)>,
    others: Vec<usize>,
    live: Vec<(usize, Vec3, AtomInfo)>,
    excl_fixed: Vec<(Vec3, f64)>,
    excl_dyn: Vec<Vec<Vec<(Vec3, f64)>>>,
    static_t: Vec<Target>,
    static_p: Vec<PTarget>,
    own_p: Vec<PTarget>,
    lane_list: Vec<PTarget>,
    seg_targets: Vec<Vec<Vec<Target>>>,
    seg_prepared: Vec<Vec<Vec<PTarget>>>,
    spheres: Vec<DepSphere>,
    cbufs: ClassifyBufs,
    own_t: Vec<Target>,
    trims: Vec<Vec<(Vec3, f64)>>,
    base_adds: Vec<(DotAdd, u64)>,
    adds: Vec<DotAdd>,
}

/// Shape `v` as one empty list per state for each entry of `dims`, keeping
/// the space of lists already there.
fn shape_lists<T>(v: &mut Vec<Vec<Vec<T>>>, dims: impl ExactSizeIterator<Item = usize>) {
    v.resize_with(dims.len(), Vec::new);
    for (per_state, d) in v.iter_mut().zip(dims) {
        per_state.resize_with(d, Vec::new);
        per_state.iter_mut().for_each(Vec::clear);
    }
}

#[allow(clippy::too_many_arguments)]
fn atom_factors_dotwise_in(
    ctx: &Ctx,
    comp: &[u32],
    doms: &[usize],
    i: usize,
    slot: usize,
    abandoned: &std::sync::atomic::AtomicBool,
    bufs: &mut DotwiseBufs,
) -> Option<(Vec<Factor>, usize, usize)> {
    use std::sync::atomic::Ordering::Relaxed;
    let DotwiseBufs {
        seq,
        others,
        live,
        excl_fixed,
        excl_dyn,
        static_t,
        static_p,
        own_p,
        lane_list,
        seg_targets,
        seg_prepared,
        spheres,
        cbufs,
        own_t,
        trims,
        base_adds,
        adds,
    } = bufs;
    if abandoned.load(Relaxed) {
        return None;
    }
    let p = ctx.p;
    let w = ctx.w;
    let m = comp[i];
    let mv = &ctx.movers[m as usize];
    let a = mv.atoms[slot];
    if w.occ[a as usize] < p.min_occupancy {
        return Some((Vec::new(), 0, 0));
    }
    let sc = &ctx.scorer;
    let pr = p.probe.probe_radius;
    let density = p.probe.density;
    let min_occ = p.min_occupancy;
    let r_live = w.info[a as usize].vdw_radius;
    let nearby = ctx.max_vdw + r_live + 2.0 * pr;
    let excl: &[u32] = ctx.exclude.get(&a).map(|v| v.as_slice()).unwrap_or(&[]);
    let soa: &DotsSoa = &ctx.dots_soa[&a];
    let [dx, dy, dz] = &soa.d;
    // Movers of the clique that can touch the atom, in `atom_movers` order;
    // `None` marks the atom's own Mover
    seq.clear();
    if let Some(ms) = ctx.atom_movers.get(&a) {
        for &om in ms {
            let li = ctx.local[om as usize];
            if li == NONE || comp.get(li as usize) != Some(&om) {
                continue;
            }
            seq.push((om, if li as usize == i { None } else { Some(li as usize) }));
        }
    }
    others.clear();
    others.extend(seq.iter().filter_map(|x| x.1));
    others.sort_unstable();
    others.dedup();
    let (seq, others) = (&seq[..], &others[..]);
    if others.len() > 64 {
        // the dot sorting keeps the reaching Movers as bits of a u64; an atom
        // within reach of more Movers settles the clique for coordinate ascent
        abandoned.store(true, Relaxed);
        return None;
    }
    let kidx = |o: usize| others.iter().position(|&x| x == o).unwrap();
    let dep_of = |mask: u64| -> Vec<usize> { others.iter().enumerate().filter(|&(k, _)| mask >> k & 1 != 0).map(|(_, &o)| o).collect() };
    let joint_others = others.iter().fold(1usize, |acc, &o| acc.saturating_mul(doms[o]));
    // tables by dependency set (local indices of other Movers, sorted)
    let mut tables: Vec<(Vec<usize>, Vec<f64>)> = Vec::new();
    let scope_of = |dep: &[usize]| -> Vec<usize> {
        let mut sc: Vec<usize> = dep.to_vec();
        sc.push(i);
        sc.sort_unstable();
        sc
    };
    let mut st = CliqueState { movers: comp, cfg: vec![Cfg::Coarse(0); comp.len()] };
    let (mut evaluated, mut stand_for) = (0usize, 0usize);
    let mut dot_work = 0usize;
    // the own states in which the atom scores, with its position and info
    live.clear();
    for own in 0..doms[i] {
        let own_cfg = Cfg::Coarse(own as u16);
        if ctx.mover_atom_deleted(m, slot, own_cfg) {
            continue;
        }
        st.cfg[i] = own_cfg;
        let (pa, _, ia) = ctx.atom_state(&st, a);
        if !(ia.is_ion && sc.ignore_ions) {
            live.push((own, pa, ia));
        }
    }
    // Mover `om`'s atoms in state `cfg` that are targets of the atom at `pa`
    let mover_targets = |om: u32, cfg: Cfg, pa: Vec3, ia: &AtomInfo, v: &mut Vec<Target>| {
        let omv = &ctx.movers[om as usize];
        v.clear();
        for oslot in 0..omv.n_moved {
            let b = omv.atoms[oslot];
            if b == a || ctx.mover_atom_deleted(om, oslot, cfg) {
                continue;
            }
            let pb = ctx.mover_atom_pos(om, oslot, cfg);
            let d = pb.dist(pa);
            if d < 1e-5 || d > nearby || excl.contains(&b) || w.occ[b as usize].abs() < min_occ {
                continue;
            }
            let ib = ctx.mover_atom_info(om, oslot, cfg, b);
            if d <= ia.vdw_radius + ib.vdw_radius + 2.0 * pr {
                v.push(Target { pos: pb, info: ib });
            }
        }
    };
    // the clique Mover an excluded atom moves with: (local index, Mover, slot)
    let owner_of = |e: u32| {
        ctx.dyn_of.get(&e).and_then(|&(em, eslot)| {
            let li = ctx.local[em as usize];
            if li != NONE && comp.get(li as usize) == Some(&em) { Some((li as usize, em, eslot as usize)) } else { None }
        })
    };
    // one farther than the dots reach plus its radius cannot cover a dot
    let reach = ctx.dot_reach.get(&a).copied().unwrap_or(f64::INFINITY) + EXCLUSION_MARGIN;
    let covers = |pa: Vec3, pe: Vec3, re: f64| pe.dist_sq(pa) < (reach + re) * (reach + re);
    // excluded atoms that can cover a dot at `pa`: fixed ones (static, Movers
    // outside the clique, and the own Mover unless `own_apart`) and, per
    // other Mover and state, the dynamic ones
    let fill_excl = |st: &CliqueState, pa: Vec3, own_apart: bool, fixed: &mut Vec<(Vec3, f64)>, dynm: &mut Vec<Vec<Vec<(Vec3, f64)>>>| {
        fixed.clear();
        dynm.iter_mut().flatten().for_each(Vec::clear);
        for &e in excl {
            match owner_of(e) {
                Some((o, em, eslot)) if o != i => {
                    let k = kidx(o);
                    for sidx in 0..doms[o] {
                        let cfg = Cfg::Coarse(sidx as u16);
                        if ctx.mover_atom_deleted(em, eslot, cfg) {
                            continue;
                        }
                        let pe = ctx.mover_atom_pos(em, eslot, cfg);
                        let ie = ctx.mover_atom_info(em, eslot, cfg, e);
                        if covers(pa, pe, ie.vdw_radius) {
                            dynm[k][sidx].push((pe, ie.vdw_radius));
                        }
                    }
                }
                Some(_) if own_apart => {}
                _ => {
                    let (pe, deleted, ie) = ctx.atom_state(st, e);
                    if !deleted && covers(pa, pe, ie.vdw_radius) {
                        fixed.push((pe, ie.vdw_radius));
                    }
                }
            }
        }
    };
    // the own Mover's excluded atoms that can cover a dot at `pa`
    let own_trims = |st: &CliqueState, pa: Vec3, out: &mut Vec<(Vec3, f64)>| {
        out.clear();
        for &e in excl {
            if matches!(owner_of(e), Some((o, _, _)) if o == i) {
                let (pe, deleted, ie) = ctx.atom_state(st, e);
                if !deleted && covers(pa, pe, ie.vdw_radius) {
                    out.push((pe, ie.vdw_radius));
                }
            }
        }
    };
    // scratch, of the same shape for every state
    let states = |o: Option<usize>| o.map_or(1, |o| doms[o]);
    shape_lists(excl_dyn, others.iter().map(|&o| doms[o]));
    shape_lists(seg_targets, seq.iter().map(|&(_, o)| states(o)));
    shape_lists(seg_prepared, seq.iter().map(|&(_, o)| states(o)));
    // static targets, Mover targets (own Mover in `own_cfg`, others per
    // state), all prepared, and the own Mover's prepared targets in `own_p`
    let fill_targets = |pa: Vec3,
                            ia: &AtomInfo,
                            own_cfg: Cfg,
                            static_t: &mut Vec<Target>,
                            static_p: &mut Vec<PTarget>,
                            own_p: &mut Vec<PTarget>,
                            seg_targets: &mut Vec<Vec<Vec<Target>>>,
                            seg_prepared: &mut Vec<Vec<Vec<PTarget>>>| {
        // static targets (same filters and order as score_atom)
        static_t.clear();
        ctx.static_targets(a, pa, ia, nearby, excl, static_t);
        for (&(om, o), per_state) in seq.iter().zip(seg_targets.iter_mut()) {
            match o {
                None => mover_targets(om, own_cfg, pa, ia, &mut per_state[0]),
                Some(_) => {
                    for (sidx, v) in per_state.iter_mut().enumerate() {
                        mover_targets(om, Cfg::Coarse(sidx as u16), pa, ia, v)
                    }
                }
            }
        }
        // the folds use prepared targets; the reach test keeps the raw lists
        prepare_targets(sc, ia, pa, static_t, pr, static_p);
        for (ts, vs) in seg_targets.iter().flatten().zip(seg_prepared.iter_mut().flatten()) {
            prepare_targets(sc, ia, pa, ts, pr, vs);
        }
        own_p.clear();
        for (q, &(_, so)) in seq.iter().enumerate() {
            if so.is_none() {
                own_p.extend_from_slice(&seg_prepared[q][0]);
            }
        }
    };
    // spheres through which each other Mover reaches a dot
    let fill_spheres = |excl_dyn: &[Vec<Vec<(Vec3, f64)>>], seg_targets: &[Vec<Vec<Target>>], out: &mut Vec<DepSphere>| {
        out.clear();
        for (k, &o) in others.iter().enumerate() {
            for &(c, r) in excl_dyn[k].iter().flatten() {
                out.push(DepSphere { k: k as u32, probe: false, c, r, r2: r * r });
            }
            for (q, &(_, so)) in seq.iter().enumerate() {
                if so == Some(o) {
                    for b in seg_targets[q].iter().flatten() {
                        let r = b.info.vdw_radius + pr;
                        out.push(DepSphere { k: k as u32, probe: true, c: b.pos, r, r2: r * r });
                    }
                }
            }
        }
    };
    // the shape of a group's table: scope, dims, strides, position of the own
    // Mover in the scope and of each dependency
    struct Shape {
        size: usize,
        strides: Vec<usize>,
        own_pos: usize,
        dep_pos: Vec<usize>,
        ddims: Vec<usize>,
        ks: Vec<usize>,
        ncombo: usize,
    }
    let shape_of = |dep: &[usize]| {
        let scope = scope_of(dep);
        let dims: Vec<usize> = scope.iter().map(|&x| doms[x]).collect();
        let mut strides = vec![0usize; scope.len()];
        let mut acc_s = 1usize;
        for q in (0..scope.len()).rev() {
            strides[q] = acc_s;
            acc_s = acc_s.saturating_mul(dims[q]);
        }
        let ddims: Vec<usize> = dep.iter().map(|&o| doms[o]).collect();
        Shape {
            size: dims.iter().product(),
            own_pos: scope.iter().position(|&x| x == i).unwrap(),
            dep_pos: dep.iter().map(|o| scope.iter().position(|x| x == o).unwrap()).collect(),
            ks: dep.iter().map(|&o| kidx(o)).collect(),
            ncombo: ddims.iter().fold(1usize, |a, &b| a.saturating_mul(b)),
            ddims,
            strides,
        }
    };
    // the fold of a grouped dot over the Mover targets with the Movers of
    // `dep` in states `sv`, or None when one of them trims it
    let group_dot = |pa: Vec3, ia: &AtomInfo, k: usize, start: DotFold, dep: &[usize], sh: &Shape, sv: &[usize], excl_dyn: &[Vec<Vec<(Vec3, f64)>>], seg_prepared: &[Vec<Vec<PTarget>>]| {
        let dp = soa.dot(k);
        let dot_abs = pa + dp.d;
        let probe_loc = pa + dp.probe;
        for (q, &kk) in sh.ks.iter().enumerate() {
            for &(pe, re) in &excl_dyn[kk][sv[q]] {
                if dot_abs.dist_sq(pe) < re * re {
                    return None;
                }
            }
        }
        let mut fold = start;
        for (qq, &(_, so)) in seq.iter().enumerate() {
            let ts: &[PTarget] = match so {
                None => &seg_prepared[qq][0],
                Some(o) => match dep.iter().position(|&x| x == o) {
                    Some(q) => &seg_prepared[qq][sv[q]],
                    None => continue,
                },
            };
            for b in ts {
                fold.step(sc, ia, dot_abs, probe_loc, b);
            }
        }
        Some(fold.finish(sc, ia, dot_abs))
    };
    let decode = |mut r: usize, ddims: &[usize], sv: &mut [usize]| {
        for q in (0..ddims.len()).rev() {
            sv[q] = r % ddims[q];
            r /= ddims[q];
        }
    };
    let shared = live.len() > 1 && live.len() <= 64 && {
        let (_, p0, i0) = live[0];
        live.iter().all(|&(_, q, inf)| inf == i0 && q.dist_sq(p0) <= SHARE_TOL * SHARE_TOL)
            && seq.iter().all(|&(om, so)| {
                so.is_some()
                    || live.iter().all(|&(own, pa, ia)| {
                        mover_targets(om, Cfg::Coarse(own as u16), pa, &ia, own_t);
                        own_t.is_empty()
                    })
            })
    };
    if shared {
        // one evaluation at the first live state's place
        let (own0, pa, ia) = live[0];
        st.cfg[i] = Cfg::Coarse(own0 as u16);
        fill_excl(&st, pa, true, excl_fixed, excl_dyn);
        fill_targets(pa, &ia, st.cfg[i], static_t, static_p, own_p, seg_targets, seg_prepared);
        trims.resize_with(live.len(), Vec::new);
        for (&(own, pa_s, _), v) in live.iter().zip(trims.iter_mut()) {
            st.cfg[i] = Cfg::Coarse(own as u16);
            own_trims(&st, pa_s, v);
        }
        let trims = &trims[..live.len()];
        // the live states (bits) in which dot `k` stays out of the own trims
        let mask_of = |k: usize| -> u64 {
            let mut mask = 0u64;
            for (q, (&(_, pa_s, _), tr)) in live.iter().zip(trims).enumerate() {
                let d = v3(pa_s.x + dx[k], pa_s.y + dy[k], pa_s.z + dz[k]);
                if !tr.iter().any(|&(pe, re)| d.dist_sq(pe) < re * re) {
                    mask |= 1 << q;
                }
            }
            mask
        };
        let total = |r: &ScoreDotsResult| r.bump / density + r.hbond / density + r.attract / density;
        base_adds.clear();
        let mut sgroups: Vec<(u64, Vec<(u32, DotFold, u64)>)> = Vec::new();
        if others.is_empty() {
            adds.clear();
            fold_dot_lanes(sc, &ia, pa, soa, &excl_fixed, &static_p, &mut cbufs.lanes, |k, r| {
                if let Some((kind, val)) = dot_addend(sc, r) {
                    adds.push(DotAdd { k: k as u32, kind, val });
                }
            });
            let t = table_for(&mut tables, &[], doms[i]);
            for (&(own, pa_s, _), tr) in live.iter().zip(trims) {
                t[own] += total(&sum_addends(adds, soa, pa_s, tr));
                evaluated += 1;
            }
        } else {
            fill_spheres(&excl_dyn, &seg_targets, spheres);
            lane_list.clear();
            lane_list.extend_from_slice(&static_p);
            lane_list.extend_from_slice(&own_p);
            let mut last = usize::MAX;
            classify_dots(
                sc,
                &ia,
                pa,
                soa,
                &excl_fixed,
                &lane_list,
                static_p.len(),
                &spheres,
                cbufs,
                |k, r| {
                    if let Some((kind, val)) = dot_addend(sc, r) {
                        base_adds.push((DotAdd { k: k as u32, kind, val }, mask_of(k)));
                    }
                },
                |dep, k, fold| {
                    if sgroups.get(last).is_none_or(|g| g.0 != dep) {
                        last = sgroups.iter().position(|g| g.0 == dep).unwrap_or_else(|| {
                            sgroups.push((dep, Vec::new()));
                            sgroups.len() - 1
                        });
                    }
                    sgroups[last].1.push((k as u32, fold, mask_of(k)));
                },
            );
            let t = table_for(&mut tables, &[], doms[i]);
            for (q, &(own, _, _)) in live.iter().enumerate() {
                let mut r = ScoreDotsResult::default();
                for &(ad, mask) in base_adds.iter() {
                    if mask >> q & 1 != 0 {
                        add_to(&mut r, ad.kind, ad.val);
                    }
                }
                t[own] += total(&r);
                evaluated += 1;
            }
            // each group exists in the states that keep one of its dots;
            // count its work as the per-state search would
            let mut work: Vec<(Vec<usize>, Vec<(u32, DotFold, u64)>, Vec<usize>, Shape)> = Vec::new();
            for (mask, members) in sgroups {
                let dep = dep_of(mask);
                let counts: Vec<usize> = (0..live.len()).map(|q| members.iter().filter(|x| x.2 >> q & 1 != 0).count()).collect();
                if counts.iter().all(|&c| c == 0) {
                    continue;
                }
                let sh = shape_of(&dep);
                for &c in counts.iter().filter(|&&c| c > 0) {
                    dot_work = dot_work.saturating_add(sh.ncombo.saturating_mul(c));
                    evaluated += sh.ncombo;
                }
                if sh.ncombo > DOT_GROUP_STATES_LIMIT {
                    dot_work = usize::MAX;
                }
                work.push((dep, members, counts, sh));
            }
            if dot_work > ATOM_DOT_WORK_LIMIT || abandoned.load(Relaxed) {
                // one atom past the budget settles it for the whole clique
                abandoned.store(true, Relaxed);
                return None;
            }
            let mut accs = vec![ScoreDotsResult::default(); live.len()];
            for (dep, members, counts, sh) in work {
                let mut sv = vec![0usize; dep.len()];
                let t = table_for(&mut tables, &dep, sh.size);
                for combo in 0..sh.ncombo {
                    decode(combo, &sh.ddims, &mut sv);
                    accs.iter_mut().for_each(|r| *r = ScoreDotsResult::default());
                    for &(k, start, mask) in &members {
                        if mask == 0 {
                            continue;
                        }
                        let Some(r) = group_dot(pa, &ia, k as usize, start, &dep, &sh, &sv, &excl_dyn, &seg_prepared) else { continue };
                        if let Some((kind, val)) = dot_addend(sc, r) {
                            let mut bits = mask;
                            while bits != 0 {
                                add_to(&mut accs[bits.trailing_zeros() as usize], kind, val);
                                bits &= bits - 1;
                            }
                        }
                    }
                    let idx_dep: usize = sv.iter().zip(&sh.dep_pos).map(|(&s, &pos)| s * sh.strides[pos]).sum();
                    for (q, &(own, _, _)) in live.iter().enumerate() {
                        if counts[q] > 0 {
                            t[own * sh.strides[sh.own_pos] + idx_dep] = total(&accs[q]);
                        }
                    }
                }
            }
        }
        stand_for = stand_for.saturating_add(joint_others.saturating_mul(live.len()));
    } else {
        let mut groups: Vec<(u64, Vec<(u32, DotFold)>)> = Vec::new();
        for &(own, pa, ia) in live.iter() {
            let own_cfg = Cfg::Coarse(own as u16);
            st.cfg[i] = own_cfg;
            fill_excl(&st, pa, false, excl_fixed, excl_dyn);
            fill_targets(pa, &ia, own_cfg, static_t, static_p, own_p, seg_targets, seg_prepared);
            let mut base = ScoreDotsResult::default();
            lane_list.clear();
            lane_list.extend_from_slice(&static_p);
            lane_list.extend_from_slice(&own_p);
            if others.is_empty() {
                // no other Mover reaches a dot: each dot folds the static
                // targets and then its own Mover's
                score_dot_lanes(sc, &ia, pa, soa, &excl_fixed, &lane_list, &mut cbufs.lanes, &mut base);
            } else {
                fill_spheres(&excl_dyn, &seg_targets, spheres);
                let mut last = usize::MAX;
                classify_dots(
                    sc,
                    &ia,
                    pa,
                    soa,
                    &excl_fixed,
                    &lane_list,
                    static_p.len(),
                    &spheres,
                    cbufs,
                    |_, r| accumulate_dot(sc, r, &mut base),
                    |dep, k, fold| {
                        if groups.get(last).is_none_or(|g| g.0 != dep) {
                            last = groups.iter().position(|g| g.0 == dep).unwrap_or_else(|| {
                                groups.push((dep, Vec::new()));
                                groups.len() - 1
                            });
                        }
                        groups[last].1.push((k as u32, fold));
                    },
                );
            }
            // the dots that no other Mover reaches
            {
                let t = table_for(&mut tables, &[], doms[i]);
                t[own] += base.bump / density + base.hbond / density + base.attract / density;
                evaluated += 1;
            }
            // the others, per group, over the states of the Movers that reach them
            for (mask, dots) in groups.drain(..) {
                let dep = dep_of(mask);
                let sh = shape_of(&dep);
                dot_work = dot_work.saturating_add(sh.ncombo.saturating_mul(dots.len()));
                if sh.ncombo > DOT_GROUP_STATES_LIMIT || dot_work > ATOM_DOT_WORK_LIMIT || abandoned.load(Relaxed) {
                    // one atom past the budget settles it for the whole clique
                    abandoned.store(true, Relaxed);
                    return None;
                }
                let mut sv = vec![0usize; dep.len()];
                let t = table_for(&mut tables, &dep, sh.size);
                for combo in 0..sh.ncombo {
                    decode(combo, &sh.ddims, &mut sv);
                    let mut acc = ScoreDotsResult::default();
                    for &(k, start) in &dots {
                        if let Some(r) = group_dot(pa, &ia, k as usize, start, &dep, &sh, &sv, &excl_dyn, &seg_prepared) {
                            accumulate_dot(sc, r, &mut acc);
                        }
                    }
                    let idx_dep: usize = sv.iter().zip(&sh.dep_pos).map(|(&s, &pos)| s * sh.strides[pos]).sum();
                    t[own * sh.strides[sh.own_pos] + idx_dep] = acc.bump / density + acc.hbond / density + acc.attract / density;
                }
                evaluated += sh.ncombo;
            }
            stand_for = stand_for.saturating_add(joint_others);
        }
    }
    let mut out = tables;
    out.sort_by(|x, y| x.0.cmp(&y.0));
    let factors = out
        .into_iter()
        .map(|(dep, table)| {
            let scope = scope_of(&dep);
            let dims = scope.iter().map(|&x| doms[x]).collect();
            Factor { scope, dims, table }
        })
        .collect();
    Some((factors, evaluated, stand_for.saturating_sub(evaluated)))
}

/// The score table of one moving atom over the states of its own Mover and
/// of the other Movers in the clique that can touch it. The other Movers'
/// states are grouped into classes that place the same atoms (with the same
/// flags) within reach of the atom, and each combination of classes is scored
/// once. Returns the factor and the numbers of scored and reused entries.
fn atom_factor(ctx: &Ctx, comp: &[u32], doms: &[usize], i: usize, slot: usize) -> (Factor, usize, usize) {
    let p = ctx.p;
    let m = comp[i];
    let mv = &ctx.movers[m as usize];
    let a = mv.atoms[slot];
    let mut scope: Vec<usize> = vec![i];
    if let Some(ms) = ctx.atom_movers.get(&a) {
        for &om in ms {
            let li = ctx.local[om as usize];
            if li != NONE && comp.get(li as usize) == Some(&om) && !scope.contains(&(li as usize)) {
                scope.push(li as usize);
            }
        }
    }
    scope.sort_unstable();
    let dims: Vec<usize> = scope.iter().map(|&s| doms[s]).collect();
    let size = Factor::size(&dims);
    let mut table = vec![0.0; size];
    let own_pos = scope.iter().position(|&s| s == i).unwrap();
    let others: Vec<usize> = scope.iter().copied().filter(|&s| s != i).collect();
    // strides of each scope position in the table
    let mut strides = vec![0usize; scope.len()];
    {
        let mut acc = 1usize;
        for q in (0..scope.len()).rev() {
            strides[q] = acc;
            acc *= dims[q];
        }
    }
    let other_pos: Vec<usize> = (0..scope.len()).filter(|&q| q != own_pos).collect();
    let r_a = ctx.w.info[a as usize].vdw_radius;
    let reach_extra = r_a + 2.0 * p.probe.probe_radius;
    let mut st = CliqueState { movers: comp, cfg: vec![Cfg::Coarse(0); comp.len()] };
    let mut buf = ScoreBuf::default();
    let (mut calculated, mut cached) = (0usize, 0usize);
    for own in 0..doms[i] {
        if ctx.mover_atom_deleted(m, slot, Cfg::Coarse(own as u16)) {
            continue; // stays 0
        }
        let pa = mv.coarse_pos[own][slot];
        // class id per other Mover state, and one representative state per class
        let mut classes: Vec<Vec<u32>> = Vec::with_capacity(others.len());
        let mut reps: Vec<Vec<usize>> = Vec::with_capacity(others.len());
        for &o in &others {
            let om = comp[o];
            let omv = &ctx.movers[om as usize];
            let mut sigs: Vec<Vec<(u16, u64, u64, u64, bool, u8)>> = Vec::new();
            let mut cls = Vec::with_capacity(doms[o]);
            let mut rep = Vec::new();
            for s in 0..doms[o] {
                let cfg = Cfg::Coarse(s as u16);
                let mut sig = Vec::new();
                for oslot in 0..omv.n_moved {
                    let b = omv.atoms[oslot];
                    let pb = omv.coarse_pos[s][oslot];
                    let rb = ctx.w.info[b as usize].vdw_radius;
                    if pb.dist(pa) <= reach_extra + rb + 1e-9 {
                        let del = ctx.mover_atom_deleted(om, oslot, cfg);
                        let inf = ctx.mover_atom_info(om, oslot, cfg, b);
                        let flags = (inf.is_acceptor as u8) | ((inf.is_donor as u8) << 1);
                        sig.push((oslot as u16, pb.x.to_bits(), pb.y.to_bits(), pb.z.to_bits(), del, flags));
                    }
                }
                let id = match sigs.iter().position(|x| *x == sig) {
                    Some(id) => id,
                    None => {
                        sigs.push(sig);
                        rep.push(s);
                        sigs.len() - 1
                    }
                };
                cls.push(id as u32);
            }
            classes.push(cls);
            reps.push(rep);
        }
        // score every combination of classes once
        let ucount: Vec<usize> = reps.iter().map(|r| r.len()).collect();
        let ncombo: usize = ucount.iter().product();
        let mut vals = vec![0.0f64; ncombo];
        let mut cv = vec![0usize; others.len()];
        for combo in 0..ncombo {
            let mut r = combo;
            for kk in (0..others.len()).rev() {
                cv[kk] = r % ucount[kk];
                r /= ucount[kk];
            }
            st.cfg[i] = Cfg::Coarse(own as u16);
            for (kk, &o) in others.iter().enumerate() {
                st.cfg[o] = Cfg::Coarse(reps[kk][cv[kk]] as u16);
            }
            vals[combo] = ctx.score_atom(&st, a, &mut buf).total();
        }
        calculated += ncombo;
        // fill the table for every assignment of the other Movers
        let ocount: usize = others.iter().map(|&o| doms[o]).product();
        cached += ocount.saturating_sub(ncombo);
        let base = own * strides[own_pos];
        let mut ovals = vec![0usize; others.len()];
        for _ in 0..ocount {
            let mut idx = base;
            let mut combo = 0usize;
            for kk in 0..others.len() {
                idx += ovals[kk] * strides[other_pos[kk]];
                combo = combo * ucount[kk] + classes[kk][ovals[kk]] as usize;
            }
            table[idx] = vals[combo];
            // next assignment (last other Mover varies fastest)
            let mut kk = others.len();
            while kk > 0 {
                kk -= 1;
                ovals[kk] += 1;
                if ovals[kk] < doms[others[kk]] {
                    break;
                }
                ovals[kk] = 0;
            }
        }
    }
    (Factor { scope, dims, table }, calculated, cached)
}

/// Largest table `variable_elimination` builds for these factor scopes (the
/// same greedy order, without the tables).
fn largest_elimination_table(nvars: usize, doms: &[usize], mut scopes: Vec<Vec<usize>>) -> usize {
    let mut eliminated = vec![false; nvars];
    let mut largest = 0usize;
    for _ in 0..nvars {
        let union_with = |v: usize, scopes: &[Vec<usize>]| -> Vec<usize> {
            let mut uni: Vec<usize> = vec![v];
            for sc in scopes.iter().filter(|sc| sc.contains(&v)) {
                for &x in sc {
                    if !uni.contains(&x) {
                        uni.push(x);
                    }
                }
            }
            uni
        };
        let cost = |uni: &[usize]| uni.iter().fold(1usize, |a, &x| a.saturating_mul(doms[x]));
        let mut best: Option<(usize, usize)> = None;
        for v in (0..nvars).filter(|&v| !eliminated[v]) {
            let c = cost(&union_with(v, &scopes));
            if best.map(|b| c < b.1).unwrap_or(true) {
                best = Some((v, c));
            }
        }
        let Some((v, c)) = best else { break };
        largest = largest.max(c);
        eliminated[v] = true;
        let mut rest: Vec<usize> = union_with(v, &scopes).into_iter().filter(|&x| x != v).collect();
        rest.sort_unstable();
        scopes.retain(|sc| !sc.contains(&v));
        if !rest.is_empty() {
            scopes.push(rest);
        }
    }
    largest
}

/// The Movers of the clique that can touch each Mover's atoms (local indices).
fn touching_movers(ctx: &Ctx, comp: &[u32]) -> Vec<Vec<usize>> {
    let k = comp.len();
    let mut touching: Vec<Vec<usize>> = vec![Vec::new(); k];
    for (i, &m) in comp.iter().enumerate() {
        let mv = &ctx.movers[m as usize];
        for &a in &mv.atoms[..mv.n_moved] {
            for &om in ctx.atom_movers.get(&a).map(|v| v.as_slice()).unwrap_or(&[]) {
                let j = ctx.local[om as usize];
                if j == NONE || j as usize == i || comp.get(j as usize) != Some(&om) {
                    continue;
                }
                let j = j as usize;
                if !touching[i].contains(&j) {
                    touching[i].push(j);
                }
                if !touching[j].contains(&i) {
                    touching[j].push(i);
                }
            }
        }
    }
    touching
}

/// Block coordinate ascent over the coarse states of a clique too dense for
/// exact search: each Mover, then each pair of Movers that touch, takes its
/// best states given the others, until nothing changes. A block is scored by
/// its Movers' preferences and atoms and by the atoms of other Movers it can
/// reach. Returns the states and the number of atom scores computed.
fn ascend_clique(ctx: &Ctx, comp: &[u32], doms: &[usize], buf: &mut ScoreBuf) -> (Vec<usize>, usize) {
    let k = comp.len();
    // (other Mover, slot) of the atoms each Mover can touch
    let mut reach: Vec<Vec<(usize, usize)>> = vec![Vec::new(); k];
    for (j, &mj) in comp.iter().enumerate() {
        let mv = &ctx.movers[mj as usize];
        for slot in 0..mv.n_moved {
            for &om in ctx.atom_movers.get(&mv.atoms[slot]).map(|v| v.as_slice()).unwrap_or(&[]) {
                let i = ctx.local[om as usize];
                if i != NONE && i as usize != j && comp.get(i as usize) == Some(&om) && !reach[i as usize].contains(&(j, slot)) {
                    reach[i as usize].push((j, slot));
                }
            }
        }
    }
    let mut pairs: Vec<(usize, usize)> = Vec::new();
    for i in 0..k {
        for &(j, _) in &reach[i] {
            let pair = (i.min(j), i.max(j));
            if !pairs.contains(&pair) {
                pairs.push(pair);
            }
        }
    }
    pairs.sort_unstable();
    let mut st = CliqueState { movers: comp, cfg: vec![Cfg::Coarse(0); k] };
    let mut calculated = 0usize;
    let mut block_score = |st: &CliqueState, block: &[usize]| -> f64 {
        let mut s = 0.0;
        for &i in block {
            s += ctx.pref(comp[i], st.cfg[i]) + ctx.score_mover_atoms(st, comp[i], buf, true).0;
            calculated += ctx.movers[comp[i] as usize].n_moved;
        }
        let mut seen: Vec<(usize, usize)> = Vec::new();
        for &i in block {
            for &(j, slot) in &reach[i] {
                if block.contains(&j) || seen.contains(&(j, slot)) {
                    continue;
                }
                seen.push((j, slot));
                if !ctx.mover_atom_deleted(comp[j], slot, st.cfg[j]) {
                    s += ctx.score_atom(st, ctx.movers[comp[j] as usize].atoms[slot], buf).total();
                    calculated += 1;
                }
            }
        }
        s
    };
    for _ in 0..100 {
        let mut changed = false;
        for i in 0..k {
            let cur = st.cfg[i].coarse();
            let mut best = (block_score(&st, &[i]), cur);
            for c in (0..doms[i]).filter(|&c| c != cur) {
                st.cfg[i] = Cfg::Coarse(c as u16);
                let s = block_score(&st, &[i]);
                if s > best.0 {
                    best = (s, c);
                }
            }
            st.cfg[i] = Cfg::Coarse(best.1 as u16);
            changed |= best.1 != cur;
        }
        if changed {
            continue;
        }
        for &(i, j) in &pairs {
            let cur = (st.cfg[i].coarse(), st.cfg[j].coarse());
            let mut best = (block_score(&st, &[i, j]), cur);
            for ci in 0..doms[i] {
                for cj in 0..doms[j] {
                    if (ci, cj) == cur {
                        continue;
                    }
                    st.cfg[i] = Cfg::Coarse(ci as u16);
                    st.cfg[j] = Cfg::Coarse(cj as u16);
                    let s = block_score(&st, &[i, j]);
                    if s > best.0 + 1e-9 {
                        best = (s, (ci, cj));
                    }
                }
            }
            st.cfg[i] = Cfg::Coarse(best.1 .0 as u16);
            st.cfg[j] = Cfg::Coarse(best.1 .1 as u16);
            changed |= best.1 != cur;
        }
        if !changed {
            break;
        }
    }
    (st.cfg.iter().map(|c| c.coarse()).collect(), calculated)
}

fn optimize_clique(ctx: &Ctx, comp: &[u32]) -> CliqueResult {
    let p = ctx.p;
    let v = p.verbosity;
    let k = comp.len();
    let mut st = CliqueState { movers: comp, cfg: vec![Cfg::Coarse(0); k] };
    // scratch space kept per thread from clique to clique; every use clears it
    thread_local! {
        static CLIQUE_BUF: std::cell::Cell<ScoreBuf> = std::cell::Cell::new(ScoreBuf::default());
    }
    let mut buf = CLIQUE_BUF.take();
    let mut calculated = 0usize;
    let mut cached = 0usize;

    // initial scores (all at coarse 0), for the report; a lone Mover's come
    // from its atom factors below
    let mut initial = if p.report { vec![0.0; k] } else { Vec::new() };
    if k > 1 && p.report {
        for (i, &m) in comp.iter().enumerate() {
            let (s, _) = ctx.score_mover_atoms(&st, m, &mut buf, true);
            calculated += ctx.movers[m as usize].n_moved;
            initial[i] = ctx.pref(m, Cfg::Coarse(0)) + s;
        }
    }

    // ---------- coarse optimization
    let doms: Vec<usize> = comp.iter().map(|&m| ctx.movers[m as usize].num_coarse()).collect();
    // the factors of variable elimination (a lone Mover needs none)
    let mut factors: Vec<Factor> = Vec::new();
    for (i, &m) in comp.iter().enumerate().filter(|_| k > 1) {
        factors.push(Factor {
            scope: vec![i],
            dims: vec![doms[i]],
            table: (0..doms[i]).map(|c| ctx.pref(m, Cfg::Coarse(c as u16))).collect(),
        });
    }
    // a clique whose Mover contact graph alone needs too large a table is not
    // searched exactly (the atom factors only add to it)
    let dense = k > 1 && {
        let touching = touching_movers(ctx, comp);
        let pairs: Vec<Vec<usize>> = (0..k).flat_map(|i| touching[i].iter().filter(move |&&j| j > i).map(move |&j| vec![i, j])).collect();
        largest_elimination_table(k, &doms, pairs) > ELIMINATION_TABLE_LIMIT
    };
    // per-atom factors (kept for final per-Mover scores), built in parallel
    let tasks: Vec<(usize, usize)> = if dense { Vec::new() } else {
        comp.iter().enumerate().flat_map(|(i, &m)| (0..ctx.movers[m as usize].n_moved).map(move |slot| (i, slot))).collect()
    };
    let dotwise = !debug_switches().atomwise;
    let abandoned = std::sync::atomic::AtomicBool::new(false);
    let built: Vec<Option<(usize, usize, Vec<Factor>, usize, usize)>> = crate::par::map_collect(&tasks, |&(i, slot)| {
            if dotwise {
                let (f, c, h) = atom_factors_dotwise(ctx, comp, &doms, i, slot, &abandoned)?;
                Some((i, slot, f, c, h))
            } else {
                let (f, c, h) = atom_factor(ctx, comp, &doms, i, slot);
                Some((i, slot, vec![f], c, h))
            }
        });
    let mut atom_factors: Vec<(usize, Factor)> = Vec::with_capacity(built.len());
    let mut exact = !dense && built.iter().all(|b| b.is_some());
    // a lone Mover's atom scores per coarse state: its atoms' values in slot
    // order, as score_mover_atoms sums them (one factor per atom, no other
    // Movers)
    let mut lone = vec![0.0; if k == 1 { doms[0] } else { 0 }];
    for (i, slot, fs, c, h) in built.into_iter().flatten() {
        // the partial work of an abandoned search depends on thread timing
        if exact {
            calculated += c;
            cached += h;
        }
        for f in fs {
            for (s, l) in lone.iter_mut().enumerate() {
                if !ctx.mover_atom_deleted(comp[0], slot, Cfg::Coarse(s as u16)) {
                    *l += f.table[s];
                }
            }
            if k > 1 {
                factors.push(Factor { scope: f.scope.clone(), dims: f.dims.clone(), table: f.table.clone() });
            }
            atom_factors.push((i, f));
        }
    }
    if k == 1 && p.report {
        initial[0] = if exact {
            ctx.pref(comp[0], Cfg::Coarse(0)) + lone[0]
        } else {
            calculated += ctx.movers[comp[0] as usize].n_moved;
            ctx.pref(comp[0], Cfg::Coarse(0)) + ctx.score_mover_atoms(&st, comp[0], &mut buf, true).0
        };
    }
    if exact && k > 1 {
        let scopes: Vec<Vec<usize>> = factors.iter().map(|f| f.scope.clone()).collect();
        exact = largest_elimination_table(k, &doms, scopes) <= ELIMINATION_TABLE_LIMIT;
    }
    let dbg = debug_switches().debug_atom.as_deref().map(|t| {
        let a0 = ctx.movers[comp[0] as usize].atoms[0] as usize;
        let l = &ctx.w.labels[a0];
        format!("{} {} {} {}", l.chain, l.resname, l.resseq, l.name) == t
    }) == Some(true);
    let best = if !exact {
        let (best, c) = ascend_clique(ctx, comp, &doms, &mut buf);
        calculated += c;
        best
    } else if k == 1 {
        // singleton: pick the first maximum in state order
        let mut bs = 0usize;
        let mut bv = f64::NEG_INFINITY;
        for c in 0..doms[0] {
            let mut s = ctx.pref(comp[0], Cfg::Coarse(c as u16));
            let mut atoms_sum = 0.0;
            for (_, f) in &atom_factors {
                atoms_sum += f.table[c];
            }
            s += atoms_sum;
            if dbg {
                let mv = &ctx.movers[comp[0] as usize];
                eprintln!(
                    "DEBUG state {} score {:.17} atoms {:?} pos {:?}",
                    c,
                    s,
                    atom_factors.iter().map(|(_, f)| f.table[c]).collect::<Vec<_>>(),
                    mv.coarse_pos[c].iter().map(|p| format!("{:.17} {:.17} {:.17}", p.x, p.y, p.z)).collect::<Vec<_>>()
                );
            }
            if c == 0 || s > bv {
                bv = s;
                bs = c;
            }
        }
        vec![bs]
    } else {
        variable_elimination(k, &doms, factors)
    };
    for i in 0..k {
        st.cfg[i] = Cfg::Coarse(best[i] as u16);
    }
    let mut high = vec![0.0; k];
    for i in 0..k {
        high[i] = ctx.pref(comp[i], st.cfg[i]);
    }
    if exact {
        for (i, f) in &atom_factors {
            high[*i] += f.table[table_index(&f.scope, &f.dims, &best)];
        }
    } else if p.report || p.compat {
        // fixed mode reads these scores only in the report
        for i in 0..k {
            high[i] += ctx.score_mover_atoms(&st, comp[i], &mut buf, true).0;
            calculated += ctx.movers[comp[i] as usize].n_moved;
        }
    }
    let coarse_best: f64 = high.iter().sum();
    let mut coarse_info = String::new();
    if !exact && p.report {
        let _ = writeln!(coarse_info, "   Clique of {} Movers is too dense for exact search; optimized by coordinate ascent", k);
    }
    if v >= 3 && p.report {
        for i in 0..k {
            let _ = write!(
                coarse_info,
                "   Setting {}Mover{} to coarse orientation {}, max score = {:.2}\n",
                if k == 1 { "single " } else { "" },
                if k == 1 { "" } else { " in clique" },
                best[i],
                high[i]
            );
        }
    }

    // ---------- fine optimization (global Mover order within the clique)
    let mut fine_info: Vec<(u32, String)> = Vec::new();
    for i in 0..k {
        let m = comp[i];
        let mv = &ctx.movers[m as usize];
        let c = st.cfg[i].coarse();
        let initial_score = if k == 1 && exact { lone[c] } else { ctx.score_mover_atoms(&st, m, &mut buf, true).0 };
        calculated += mv.n_moved;
        let nf = mv.fine.count();
        let mut msg = String::new();
        if nf > 0 {
            let mut best_f = 0usize;
            let mut best_s = f64::NEG_INFINITY;
            // the fine poses lie close together: each atom scores only the
            // dots that can meet a target in one of them
            let shared: Vec<Option<FineShared>> = (0..mv.n_moved).map(|slot| ctx.fine_shared(&st, i, c, slot, &mut buf)).collect();
            for f in 0..nf {
                st.cfg[i] = Cfg::Fine(c as u16, f as u16);
                let (s, _) = ctx.score_mover_atoms_in(&st, m, &mut buf, true, &shared, f);
                calculated += mv.n_moved;
                let s = s + ctx.pref(m, st.cfg[i]);
                if f == 0 || s > best_s {
                    best_s = s;
                    best_f = f;
                }
            }
            let compare = if p.compat { high[i] } else { initial_score + ctx.pref(m, Cfg::Coarse(c as u16)) };
            if best_s > compare {
                if v >= 3 && p.report {
                    msg += &format!(
                        "   Setting Mover to fine orientation {}, max score = {:.2} (coarse score {:.2})\n",
                        best_f, best_s, initial_score
                    );
                }
                st.cfg[i] = Cfg::Fine(c as u16, best_f as u16);
                high[i] = best_s;
            } else {
                st.cfg[i] = Cfg::Coarse(c as u16);
                if v >= 3 && p.report {
                    msg += "   Leaving Mover at coarse orientation\n";
                }
                if !p.compat {
                    high[i] = compare;
                }
            }
        } else if !p.compat {
            high[i] = initial_score + ctx.pref(m, Cfg::Coarse(c as u16));
        }
        if p.report {
            fine_info.push((m, msg));
        }
    }
    if !p.compat && k > 1 && p.report {
        // Report each Mover's score in the final configuration of its clique
        // (a lone Mover's is the score its fine search ended on).
        for i in 0..k {
            let (s, _) = ctx.score_mover_atoms(&st, comp[i], &mut buf, true);
            calculated += ctx.movers[comp[i] as usize].n_moved;
            high[i] = s + ctx.pref(comp[i], st.cfg[i]);
        }
    }

    // ---------- flip annotations (BothClash / Uncertain)
    let mut annot = if p.report { vec![String::new(); k] } else { Vec::new() };
    for i in (0..k).filter(|_| p.report) {
        let m = comp[i];
        let mv = &ctx.movers[m as usize];
        if !mv.kind.is_flip() {
            annot[i] = " .".into();
            continue;
        }
        if let MoverKind::HisFlip { enabled, .. } = mv.kind {
            if enabled != 3 && !p.compat {
                annot[i] = " .".into();
                continue;
            }
        }
        let n = mv.num_coarse();
        let fin = st.cfg[i].coarse();
        let other = (fin + n / 2) % n;
        let saved = st.cfg[i];
        st.cfg[i] = Cfg::Coarse(other as u16);
        let (os, ob) = ctx.score_mover_atoms(&st, m, &mut buf, !p.compat);
        let os = os + ctx.pref(m, st.cfg[i]);
        st.cfg[i] = Cfg::Coarse(fin as u16);
        let (fs, fb) = ctx.score_mover_atoms(&st, m, &mut buf, !p.compat);
        let fs = fs + ctx.pref(m, st.cfg[i]);
        st.cfg[i] = saved;
        let desc = mv.pose_description(fin, None, !p.skip_bond_fixup);
        annot[i] = if ob && fb {
            " BothClash".into()
        } else if desc.contains("Unflipped") && os > fs && os - fs <= p.non_flip_preference {
            " Uncertain".into()
        } else {
            " .".into()
        };
    }
    CLIQUE_BUF.set(buf);

    CliqueResult {
        initial,
        cfg: st.cfg.clone(),
        high,
        coarse_best,
        coarse_info,
        fine_info,
        annot,
        calculated,
        cached,
    }
}
