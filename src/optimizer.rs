//! Placement and optimization of Movers (port of mmtbx/reduce/Optimizers.py,
//! Optimizers.cpp and InteractionGraph.*).
//!
//! The objective is exactly Reduce2's: for every Mover, the preference energy
//! of its state plus the Probe dot score of each of its (non-deleted) atoms
//! against everything else. Cliques of interacting Movers are optimized
//! exactly with variable elimination over cached per-atom score tables
//! instead of Reduce2's exhaustive/vertex-cut search, and independent cliques
//! are processed in parallel.

use crate::geom::*;
use crate::movers::{self, in_state, Mover, MoverKind, SingleHOptions};
use crate::probe::*;
use crate::resclass::ResClass;
use crate::world::World;
use rustc_hash::{FxHashMap, FxHashSet};
use smallvec::SmallVec;
use std::fmt::Write as _;
use std::sync::Arc;
use std::time::Instant;

/// Behavior switches: `compat` reproduces Reduce2 quirks that are bugs.
#[derive(Clone, Debug)]
pub struct OptParams {
    pub probe: ProbeParams,
    pub add_flip_movers: bool,
    pub alt_id: Option<String>,
    pub bonded_neighbor_depth: usize,
    pub use_neutron_distances: bool,
    pub min_occupancy: f64,
    pub preference_magnitude: f64,
    pub non_flip_preference: f64,
    /// Fixed mode: preference for a hydroxyl hydrogen on a planar atom (a
    /// phenol, an enol, a carboxylic acid) to lie in its plane.
    pub planar_hydroxyl_preference: f64,
    /// Fixed mode: further preference for an acid's hydroxyl hydrogen to be
    /// syn to the carbonyl oxygen rather than anti.
    pub acid_syn_preference: f64,
    pub skip_bond_fixup: bool,
    pub flip_states: String,
    pub verbosity: i32,
    pub compat: bool,
    /// Build the report text. Without it the optimizer skips the scores and
    /// strings that only the report shows; the chosen states, deletions and
    /// fix-ups stay the same.
    pub report: bool,
}

impl Default for OptParams {
    fn default() -> Self {
        OptParams {
            probe: ProbeParams::reduce2_defaults(),
            // Reduce2 defaults to False; flips are on unless asked otherwise.
            add_flip_movers: true,
            alt_id: None,
            bonded_neighbor_depth: 4,
            use_neutron_distances: false,
            min_occupancy: 0.02,
            preference_magnitude: 1.0,
            non_flip_preference: 0.5,
            planar_hydroxyl_preference: 1.0,
            acid_syn_preference: 1.0,
            skip_bond_fixup: false,
            flip_states: String::new(),
            verbosity: 2,
            compat: false,
            report: true,
        }
    }
}

#[derive(Clone, Debug)]
pub struct FlipMoverState {
    pub mover_type: String,
    pub model_num: i64,
    pub alt_id: String,
    pub any_alt: bool,
    pub chain: String,
    pub res_name: String,
    pub res_id: i32,
    pub icode: String,
    pub flipped: bool,
    pub fixed_up: bool,
}

pub fn parse_flip_states(s: &str, compat: bool) -> Vec<FlipMoverState> {
    let mut ret = Vec::new();
    for state in s.split(',') {
        let words: Vec<&str> = state.split_whitespace().collect();
        if words.len() != 6 && words.len() != 7 {
            continue;
        }
        let model_num: i64 = words[0].parse().unwrap_or(1);
        let mut alt = words[1].to_ascii_lowercase();
        let any_alt = alt == ".";
        if any_alt {
            alt.clear();
        }
        let chain = if compat { words[2].to_ascii_uppercase() } else { words[2].to_string() };
        let res_name = words[3].to_ascii_uppercase();
        let t = if res_name == "HIS" { "HisFlip" } else { "AmideFlip" };
        let rid = words[4];
        let (res_id, icode) = match rid.parse::<i32>() {
            Ok(v) => (v, String::new()),
            Err(_) => (rid[..rid.len() - 1].parse().unwrap_or(0), rid[rid.len() - 1..].to_string()),
        };
        let flipped = words[5] == "Flipped";
        let fixed_up = words.len() == 7 && words[6] == "AnglesAdjusted";
        ret.push(FlipMoverState {
            mover_type: t.into(),
            // compat: the original adds one here and one more when comparing.
            model_num: if compat { model_num + 1 } else { model_num },
            alt_id: alt,
            any_alt,
            chain,
            res_name,
            res_id,
            icode,
            flipped,
            fixed_up,
        });
    }
    ret
}

fn find_flip_state<'a>(w: &World, a: u32, states: &'a [FlipMoverState], compat: bool, n_models: usize) -> Option<&'a FlipMoverState> {
    let l = &w.labels[a as usize];
    for fs in states {
        let model_ok = if compat {
            // The original compares a string model id against an integer, so
            // only files without MODEL records ever match.
            l.model_id.is_empty()
        } else {
            (l.model_index as i64 + 1 == fs.model_num) || (n_models == 1 && fs.model_num == 1)
        };
        let alt = l.altloc.trim();
        let alt_ok = if compat {
            alt.is_empty() || alt.to_ascii_lowercase() == fs.alt_id.to_ascii_lowercase()
        } else {
            fs.any_alt || alt.is_empty() || alt.to_ascii_lowercase() == fs.alt_id
        };
        if model_ok && *l.chain == *fs.chain && l.resseq == fs.res_id && alt_ok && l.icode.trim() == fs.icode.trim() {
            return Some(fs);
        }
    }
    None
}

/// Flip performed (by Mover or lock-down), for reporting/Flipkins.
#[derive(Clone, Debug)]
pub struct FlippedMoverInfo {
    pub alt: String,
    pub base_atom: u32,
}

#[derive(Default)]
pub struct OptOutput {
    pub info: String,
    pub warnings: String,
    pub hydrogens_to_delete: Vec<u32>,
    pub amide_flips: Vec<FlippedMoverInfo>,
    pub his_flips: Vec<FlippedMoverInfo>,
    pub num_calculated: usize,
    pub num_cached: usize,
}

fn vcheck(verbosity: i32, level: i32, msg: &str) -> String {
    if verbosity >= level { format!("{}{}", " ".repeat(level as usize), msg) } else { String::new() }
}

struct Timer {
    t: Instant,
    verbosity: i32,
}
impl Timer {
    fn new(v: i32) -> Self {
        Timer { t: Instant::now(), verbosity: v }
    }
    fn report(&mut self, msg: &str) -> String {
        let d = self.t.elapsed().as_secs_f64();
        self.t = Instant::now();
        vcheck(self.verbosity, 2, &format!("Time to {}: {:.3}\n", msg, d))
    }
}

/// Python `str(float)`.
pub fn py_float(v: f64) -> String {
    let s = format!("{:?}", v);
    s
}

/// `getAtomsWithinNBonds` from probe Helpers.
pub fn atoms_within_n_bonds(w: &World, atom: u32, probe_rad: f64, n: usize, non_h_n: usize) -> Vec<u32> {
    let a_loc = w.pos[atom as usize];
    let a_rad = w.info[atom as usize].vdw_radius;
    let atom_is_h = w.is_h(atom);
    let n = if w.info[atom as usize].is_dummy_hydrogen { 1 } else { n };
    let mut set: Vec<u32> = vec![atom];
    for i in 0..n {
        // the atoms found before this round
        for k in 0..set.len() {
            let a = set[k];
            for &nb in &w.bonded[a as usize] {
                if i < non_h_n || atom_is_h || w.is_h(nb) {
                    if w.compatible(atom, nb) {
                        let d = w.pos[nb as usize].dist(a_loc);
                        if d <= a_rad + w.info[nb as usize].vdw_radius + 2.0 * probe_rad && !set.contains(&nb) {
                            set.push(nb);
                        }
                    }
                }
            }
        }
    }
    set.retain(|&x| x != atom);
    set
}

/// `fixupExplicitDonors`.
fn fixup_explicit_donors(w: &mut World, atoms: &[u32]) {
    for &a in atoms {
        if w.is_h(a) {
            for k in 0..w.bonded[a as usize].len() {
                let n = w.bonded[a as usize][k];
                let e = w.elem(n);
                if e == "N" || e == "O" || e == "S" {
                    w.info[a as usize].is_donor = true;
                    w.info[n as usize].is_donor = false;
                }
            }
        } else {
            let e = w.elem(a);
            if e == "N" || e == "O" || e == "S" {
                w.info[a as usize].is_donor = false;
            }
        }
    }
}

/// Aromatic-ring atoms treated as acceptors (`AtomTypes.IsAromaticAcceptor`).
pub fn is_aromatic_acceptor(resname: &str, atom_name: &str) -> bool {
    crate::atomtypes::is_aromatic_acceptor(resname, atom_name)
}

/// Neighbor search used while placing Movers.
///
/// Fixed mode finds atoms at their current positions. Compat mode emulates
/// Reduce2's `SpatialQuery`, which files each atom in the grid cell of the
/// position it had when last added: atoms that Movers or methyl staggering
/// move without removing and re-adding them stay filed in their old cell, so
/// queries can miss them.
struct NeighborIndex {
    loose: Option<SpatialGrid>,
    /// Fixed mode: the atoms added after `loose` was built (phantoms).
    late: Option<SpatialGrid>,
    margin: f64,
    reg: Option<RegGrid>,
}

struct RegGrid {
    geom: SpatialGrid,
    cells: Vec<SmallVec<[u32; 8]>>,
}

impl NeighborIndex {
    fn new(w: &World, atoms: &[u32], compat: bool) -> Self {
        if compat {
            let (lo, bin, dims) = reduce2_grid_geometry(atoms.iter().map(|&a| w.pos[a as usize]));
            let geom = SpatialGrid::with_geometry(lo, bin, dims, &[]);
            let mut cells = vec![SmallVec::new(); dims[0] * dims[1] * dims[2]];
            for &a in atoms {
                let c = geom.cell_of(w.pos[a as usize]);
                if !cells[c].contains(&a) {
                    cells[c].push(a);
                }
            }
            NeighborIndex { loose: None, late: None, margin: 0.0, reg: Some(RegGrid { geom, cells }) }
        } else {
            let pts: Vec<(u32, Vec3)> = atoms.iter().map(|&a| (a, w.pos[a as usize])).collect();
            NeighborIndex { loose: Some(SpatialGrid::new(&pts, 3.0)), late: None, margin: 2.5, reg: None }
        }
    }

    /// Fixed mode: add atoms `ids` to the index, as a rebuild over all atoms
    /// would. No atom may have moved since the index was built.
    fn add_late(&mut self, w: &World, ids: &[u32]) {
        if self.loose.is_some() && !ids.is_empty() {
            let pts: Vec<(u32, Vec3)> = ids.iter().map(|&a| (a, w.pos[a as usize])).collect();
            self.late = Some(SpatialGrid::new(&pts, 3.0));
        }
    }

    /// `SpatialQuery::remove(a)` at the atom's current position.
    fn remove(&mut self, w: &World, a: u32) {
        if let Some(r) = &mut self.reg {
            let c = r.geom.cell_of(w.pos[a as usize]);
            r.cells[c].retain(|x| *x != a);
        }
    }
    /// `SpatialQuery::add(a)` at the atom's current position.
    fn add(&mut self, w: &World, a: u32) {
        if let Some(r) = &mut self.reg {
            let c = r.geom.cell_of(w.pos[a as usize]);
            if !r.cells[c].contains(&a) {
                r.cells[c].push(a);
            }
        } else if let Some(g) = &self.loose {
            // atoms added later (phantoms) are not in the loose grid; callers
            // keep them out of placement queries, as the original does not
            // search for them there either except for ion checks.
            let _ = g;
        }
    }

    fn neighbors(&self, w: &World, p: Vec3, min_d: f64, max_d: f64, removed: &FxHashSet<u32>) -> Vec<u32> {
        let mut out = Vec::new();
        let (min2, max2) = (min_d * min_d, max_d * max_d);
        if let Some(r) = &self.reg {
            let [(x0, x1), (y0, y1), (z0, z1)] = r.geom.cell_ranges(p, max_d);
            let dims = r.geom.dims();
            for z in z0..=z1 {
                for y in y0..=y1 {
                    for x in x0..=x1 {
                        for &id in &r.cells[x + dims[0] * (y + dims[1] * z)] {
                            let d2 = w.pos[id as usize].dist_sq(p);
                            if d2 >= min2 && d2 <= max2 {
                                out.push(id);
                            }
                        }
                    }
                }
            }
        } else {
            // the distance test first: it is cheaper than the set lookup and
            // rejects most of the grid's candidates
            for g in self.loose.iter().chain(&self.late) {
                g.for_each_within(p, 0.0, max_d + self.margin, |id, _, _| {
                    let d2 = w.pos[id as usize].dist_sq(p);
                    if d2 >= min2 && d2 <= max2 && !removed.contains(&id) {
                        out.push(id);
                    }
                });
            }
        }
        out.sort_unstable();
        out.dedup();
        out
    }

    /// Grid of the static scoring targets. Compat mode files them in their
    /// registered cells with Reduce2's geometry.
    fn static_grid(&self, w: &World, include: impl Fn(u32) -> bool) -> Option<SpatialGrid> {
        let r = self.reg.as_ref()?;
        let mut pts: Vec<(u32, Vec3, Vec3)> = Vec::new();
        for (ci, cell) in r.cells.iter().enumerate() {
            if cell.is_empty() {
                continue;
            }
            let center = r.geom.cell_center(ci);
            for &a in cell {
                if include(a) {
                    pts.push((a, center, w.pos[a as usize]));
                }
            }
        }
        let g = &r.geom;
        Some(SpatialGrid::with_geometry(g.lower_bound(), g.bins(), g.dims(), &pts))
    }
}

// ----------------------------------------------------------------------------
// Scoring context

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Cfg {
    Coarse(u16),
    Fine(u16, u16),
}
impl Cfg {
    fn coarse(self) -> usize {
        match self {
            Cfg::Coarse(c) | Cfg::Fine(c, _) => c as usize,
        }
    }
}

#[derive(Clone, Copy)]
struct DotPair {
    d: Vec3,
    probe: Vec3,
}

const NONE: u32 = u32::MAX;

struct Ctx<'a> {
    w: &'a World,
    movers: &'a [Mover],
    p: &'a OptParams,
    scorer: DotScorer,
    /// Scorer view of atom info in compat mode (Reduce2's DotScorer keeps a copy).
    snapshot: Option<Vec<AtomInfo>>,
    grid: SpatialGrid,
    /// (mover, slot) for atoms that move with a Mover (slot < n_moved).
    dyn_of: DynOf,
    atom_movers: FxHashMap<u32, SmallVec<[u32; 4]>>,
    exclude: FxHashMap<u32, Vec<u32>>,
    dots: FxHashMap<u32, Arc<Vec<DotPair>>>,
    /// The same dots as coordinate arrays.
    dots_soa: FxHashMap<u32, Arc<DotsSoa>>,
    /// Largest distance of an atom's dots from its center.
    dot_reach: FxHashMap<u32, f64>,
    /// Per Mover atom: the static atoms near any of its positions that can
    /// be its targets (not excluded, occupied), in grid order with their
    /// info (center, radius covered, range of `near_entries`), filtered per
    /// state instead of querying the grid for each.
    near: FxHashMap<u32, (Vec3, f64, u32, u32)>,
    near_entries: Vec<(Vec3, AtomInfo)>,
    max_vdw: f64,
    /// Local index of each mover inside its clique.
    local: Vec<u32>,
}

/// (mover, slot) of each atom that moves with a Mover, by atom index.
struct DynOf(Vec<(u32, u16)>);

impl DynOf {
    #[inline]
    fn get(&self, a: &u32) -> Option<&(u32, u16)> {
        self.0.get(*a as usize).filter(|e| e.0 != NONE)
    }
    #[inline]
    fn contains_key(&self, a: &u32) -> bool {
        self.get(a).is_some()
    }
    fn insert(&mut self, a: u32, v: (u32, u16)) {
        if self.0.len() <= a as usize {
            self.0.resize(a as usize + 1, (NONE, 0));
        }
        self.0[a as usize] = v;
    }
}

/// The configuration of the Movers of one clique.
struct CliqueState<'s> {
    movers: &'s [u32],
    cfg: Vec<Cfg>,
}

impl<'a> Ctx<'a> {
    #[inline]
    fn cfg_of(&self, st: &CliqueState, m: u32) -> Cfg {
        st.cfg[self.local[m as usize] as usize]
    }

    #[inline]
    fn mover_atom_pos(&self, m: u32, slot: usize, cfg: Cfg) -> Vec3 {
        let mv = &self.movers[m as usize];
        match cfg {
            Cfg::Coarse(c) => mv.coarse_pos[c as usize][slot],
            Cfg::Fine(c, f) => mv.fine.pos(c as usize, f as usize)[slot],
        }
    }

    #[inline]
    fn mover_atom_deleted(&self, m: u32, slot: usize, cfg: Cfg) -> bool {
        let d = in_state(&self.movers[m as usize].coarse_del, cfg.coarse());
        slot < d.len() && d[slot]
    }

    #[inline]
    fn mover_atom_info(&self, m: u32, slot: usize, cfg: Cfg, a: u32) -> AtomInfo {
        if let Some(s) = &self.snapshot {
            return s[a as usize];
        }
        let inf = in_state(&self.movers[m as usize].coarse_info, cfg.coarse());
        if slot < inf.len() { inf[slot] } else { self.w.info[a as usize] }
    }

    #[inline]
    fn static_info(&self, a: u32) -> AtomInfo {
        match &self.snapshot {
            Some(s) => s[a as usize],
            None => self.w.info[a as usize],
        }
    }

    /// Position, deletion and info of any atom under the clique state.
    #[inline]
    fn atom_state(&self, st: &CliqueState, a: u32) -> (Vec3, bool, AtomInfo) {
        match self.dyn_of.get(&a) {
            Some(&(m, slot)) if self.local[m as usize] != NONE && st.movers.get(self.local[m as usize] as usize) == Some(&m) => {
                let cfg = self.cfg_of(st, m);
                (
                    self.mover_atom_pos(m, slot as usize, cfg),
                    self.mover_atom_deleted(m, slot as usize, cfg),
                    self.mover_atom_info(m, slot as usize, cfg, a),
                )
            }
            _ => (self.w.pos[a as usize], false, self.static_info(a)),
        }
    }

    /// The static targets of atom `a` at `pa` (`for_each_within` over the grid
    /// with the target filters), from its cached neighborhood when that
    /// covers the position.
    #[inline]
    fn static_targets(&self, a: u32, pa: Vec3, ia: &AtomInfo, nearby: f64, excl: &[u32], out: &mut Vec<Target>) {
        let pr = self.p.probe.probe_radius;
        if let Some(&(c, r, lo, hi)) = self.near.get(&a) {
            if pa.dist(c) + nearby + NEAR_MARGIN <= r {
                let (min2, max2) = (1e-5 * 1e-5, nearby * nearby);
                for &(pb, ib) in &self.near_entries[lo as usize..hi as usize] {
                    let d2 = pb.dist_sq(pa);
                    if d2 >= min2 && d2 <= max2 && d2.sqrt() <= ia.vdw_radius + ib.vdw_radius + 2.0 * pr {
                        out.push(Target { pos: pb, info: ib });
                    }
                }
                return;
            }
        }
        let (w, min_occ) = (self.w, self.p.min_occupancy);
        self.grid.for_each_within(pa, 1e-5, nearby, |b, pb, d2| {
            if excl.contains(&b) || w.occ[b as usize].abs() < min_occ {
                return;
            }
            let ib = self.static_info(b);
            if d2.sqrt() <= ia.vdw_radius + ib.vdw_radius + 2.0 * pr {
                out.push(Target { pos: pb, info: ib });
            }
        });
    }

    /// One evaluation of the dots of Mover atom `slot` of the Mover at local
    /// index `i` for all of its fine poses in coarse state `c`, with the rest
    /// of the clique as in `st`. None unless the atom keeps its place and info
    /// in every pose and its own Mover offers it no targets.
    fn fine_shared(&self, st: &CliqueState, i: usize, c: usize, slot: usize, buf: &mut ScoreBuf) -> Option<FineShared> {
        let (w, p) = (self.w, self.p);
        let m = st.movers[i];
        let mv = &self.movers[m as usize];
        let a = mv.atoms[slot];
        let nf = mv.fine.count();
        let coarse = Cfg::Coarse(c as u16);
        if nf < 2 || p.compat || w.occ[a as usize] < p.min_occupancy || self.mover_atom_deleted(m, slot, coarse) {
            return None;
        }
        let ia = self.mover_atom_info(m, slot, coarse, a);
        if ia.is_ion && self.scorer.ignore_ions {
            return None;
        }
        let pa0 = mv.fine.pos(c, 0)[slot];
        if (0..nf).any(|f| mv.fine.pos(c, f)[slot].dist_sq(pa0) > SHARE_TOL * SHARE_TOL) {
            return None;
        }
        let pa: Vec<Vec3> = (0..nf).map(|f| mv.fine.pos(c, f)[slot]).collect();
        let pr = p.probe.probe_radius;
        let nearby = self.max_vdw + w.info[a as usize].vdw_radius + 2.0 * pr;
        let excl: &[u32] = self.exclude.get(&a).map(|v| v.as_slice()).unwrap_or(&[]);
        let min_occ = p.min_occupancy;
        // the own Mover's other atoms must be no target in any pose
        for oslot in 0..mv.n_moved {
            let b = mv.atoms[oslot];
            if b == a || self.mover_atom_deleted(m, oslot, coarse) || excl.contains(&b) || w.occ[b as usize].abs() < min_occ {
                continue;
            }
            let ib = self.mover_atom_info(m, oslot, coarse, b);
            if (0..nf).any(|f| {
                let d = mv.fine.pos(c, f)[oslot].dist(pa[f]);
                (1e-5..=nearby).contains(&d) && d <= ia.vdw_radius + ib.vdw_radius + 2.0 * pr
            }) {
                return None;
            }
        }
        // excluded atoms: the own Mover's per pose, the others once
        let reach = self.dot_reach.get(&a).copied().unwrap_or(f64::INFINITY) + EXCLUSION_MARGIN;
        let covers = |pa: Vec3, pe: Vec3, re: f64| pe.dist_sq(pa) < (reach + re) * (reach + re);
        buf.excl.clear();
        let mut trims = vec![Vec::new(); nf];
        for &e in excl {
            match self.dyn_of.get(&e) {
                Some(&(em, eslot)) if em == m => {
                    if self.mover_atom_deleted(m, eslot as usize, coarse) {
                        continue;
                    }
                    let re = self.mover_atom_info(m, eslot as usize, coarse, e).vdw_radius;
                    for (f, tr) in trims.iter_mut().enumerate() {
                        let pe = mv.fine.pos(c, f)[eslot as usize];
                        if covers(pa[f], pe, re) {
                            tr.push((pe, re));
                        }
                    }
                }
                _ => {
                    let (pe, deleted, ie) = self.atom_state(st, e);
                    if !deleted && covers(pa[0], pe, ie.vdw_radius) {
                        buf.excl.push((pe, ie.vdw_radius));
                    }
                }
            }
        }
        // targets: static, then the other clique Movers' (as in score_atom_ext)
        buf.targets.clear();
        self.static_targets(a, pa[0], &ia, nearby, excl, &mut buf.targets);
        for &om in self.atom_movers.get(&a).map(|v| v.as_slice()).unwrap_or(&[]) {
            let li = self.local[om as usize];
            if om == m || li == NONE || st.movers.get(li as usize) != Some(&om) {
                continue;
            }
            let cfg = st.cfg[li as usize];
            let omv = &self.movers[om as usize];
            for oslot in 0..omv.n_moved {
                let b = omv.atoms[oslot];
                if b == a || self.mover_atom_deleted(om, oslot, cfg) {
                    continue;
                }
                let pb = self.mover_atom_pos(om, oslot, cfg);
                let d = pb.dist(pa[0]);
                if d < 1e-5 || d > nearby || excl.contains(&b) || w.occ[b as usize].abs() < min_occ {
                    continue;
                }
                let ib = self.mover_atom_info(om, oslot, cfg, b);
                if d <= ia.vdw_radius + ib.vdw_radius + 2.0 * pr {
                    buf.targets.push(Target { pos: pb, info: ib });
                }
            }
        }
        let s = &self.scorer;
        prepare_targets(s, &ia, pa[0], &buf.targets, pr, &mut buf.prepared);
        let mut adds = Vec::new();
        fold_dot_lanes(s, &ia, pa[0], &self.dots_soa[&a], &buf.excl, &buf.prepared, &mut buf.lanes, |k, r| {
            if let Some((kind, val)) = dot_addend(s, r) {
                adds.push(DotAdd { k: k as u32, kind, val });
            }
        });
        Some(FineShared { adds, pa, trims })
    }

    /// Score one Mover atom under the clique state (`OptimizerC::scoreAtom`).
    fn score_atom(&self, st: &CliqueState, a: u32, buf: &mut ScoreBuf) -> ScoreDotsResult {
        self.score_atom_ext(st, a, buf, None)
    }

    /// `score_atom` with Reduce2's trimmed-dot cache keyed by (atom, location
    /// index) when `trim` is given (compat mode).
    fn score_atom_ext(
        &self,
        st: &CliqueState,
        a: u32,
        buf: &mut ScoreBuf,
        trim: Option<(&std::cell::RefCell<FxHashMap<(u32, u32), Vec<DotPair>>>, u32)>,
    ) -> ScoreDotsResult {
        let w = self.w;
        if w.occ[a as usize] < self.p.min_occupancy {
            return ScoreDotsResult::default();
        }
        let (pa, _, ia) = self.atom_state(st, a);
        let pr = self.p.probe.probe_radius;
        let r_live = w.info[a as usize].vdw_radius;
        let nearby = self.max_vdw + r_live + 2.0 * pr;
        let excl: &[u32] = self.exclude.get(&a).map(|v| v.as_slice()).unwrap_or(&[]);
        let min_occ = self.p.min_occupancy;
        buf.targets.clear();
        buf.excl.clear();
        // excluded atoms' current geometry (for trimming); one farther than
        // the dots reach plus its radius cannot cover a dot
        let fixed = !self.p.compat;
        let reach = self.dot_reach.get(&a).copied().unwrap_or(f64::INFINITY) + EXCLUSION_MARGIN;
        for &e in excl {
            let (pe, deleted, ie) = self.atom_state(st, e);
            if fixed && deleted {
                continue;
            }
            if pe.dist_sq(pa) >= (reach + ie.vdw_radius) * (reach + ie.vdw_radius) {
                continue;
            }
            buf.excl.push((pe, ie.vdw_radius));
        }
        // static neighbors
        self.static_targets(a, pa, &ia, nearby, excl, &mut buf.targets);
        // dynamic neighbors from interacting Movers
        if let Some(ms) = self.atom_movers.get(&a) {
            for &m in ms {
                let li = self.local[m as usize];
                if li == NONE || st.movers.get(li as usize) != Some(&m) {
                    continue;
                }
                let cfg = st.cfg[li as usize];
                let mv = &self.movers[m as usize];
                for slot in 0..mv.n_moved {
                    let b = mv.atoms[slot];
                    if b == a || self.mover_atom_deleted(m, slot, cfg) {
                        continue;
                    }
                    let pb = self.mover_atom_pos(m, slot, cfg);
                    let d = pb.dist(pa);
                    if d < 1e-5 || d > nearby {
                        continue;
                    }
                    if excl.contains(&b) || w.occ[b as usize].abs() < min_occ {
                        continue;
                    }
                    let ib = self.mover_atom_info(m, slot, cfg, b);
                    if d <= ia.vdw_radius + ib.vdw_radius + 2.0 * pr {
                        buf.targets.push(Target { pos: pb, info: ib });
                    }
                }
            }
        }
        // trim dots inside excluded atoms
        if let Some((cache, loc)) = trim {
            if let Some(d) = cache.borrow().get(&(a, loc)) {
                return self.scorer_score(pa, &ia, d, &buf.targets, &mut buf.prepared);
            }
        }
        if trim.is_none() {
            // all dots in coordinate arrays: trimming and scoring in one pass
            let s = &self.scorer;
            let mut ret = ScoreDotsResult::default();
            if ia.is_ion && s.ignore_ions {
                return ret;
            }
            prepare_targets(s, &ia, pa, &buf.targets, pr, &mut buf.prepared);
            score_dot_lanes(s, &ia, pa, &self.dots_soa[&a], &buf.excl, &buf.prepared, &mut buf.lanes, &mut ret);
            let density = self.p.probe.density;
            ret.bump /= density;
            ret.hbond /= density;
            ret.attract /= density;
            return ret;
        }
        let dots = &self.dots[&a];
        buf.dots.clear();
        'dot: for dp in dots.iter() {
            let p = pa + dp.d;
            for &(pe, re) in &buf.excl {
                if p.dist_sq(pe) < re * re {
                    continue 'dot;
                }
            }
            buf.dots.push(*dp);
        }
        if let Some((cache, loc)) = trim {
            cache.borrow_mut().insert((a, loc), buf.dots.clone());
        }
        self.scorer_score(pa, &ia, &buf.dots, &buf.targets, &mut buf.prepared)
    }

    fn scorer_score(&self, pa: Vec3, ia: &AtomInfo, dots: &[DotPair], targets: &[Target], prepared: &mut Vec<PTarget>) -> ScoreDotsResult {
        // Same arithmetic as DotScorer::score_dots with precomputed probe offsets.
        let pr = self.p.probe.probe_radius;
        let density = self.p.probe.density;
        let s = &self.scorer;
        let mut ret = ScoreDotsResult::default();
        if ia.is_ion && s.ignore_ions {
            return ret;
        }
        prepare_targets(s, ia, pa, targets, pr, prepared);
        for dp in dots {
            let dot_abs = pa + dp.d;
            let probe_loc = pa + dp.probe;
            accumulate_dot(s, check_dot(s, ia, dot_abs, probe_loc, prepared), &mut ret);
        }
        ret.bump /= density;
        ret.hbond /= density;
        ret.attract /= density;
        ret
    }

    /// Sum of the scores of a Mover's scored atoms under the clique state.
    fn score_mover_atoms(&self, st: &CliqueState, m: u32, buf: &mut ScoreBuf, skip_deleted: bool) -> (f64, bool) {
        self.score_mover_atoms_in(st, m, buf, skip_deleted, &[], 0)
    }

    /// `score_mover_atoms` in fine pose `f`, with each atom scored from its
    /// shared evaluation `shared[slot]` when present.
    fn score_mover_atoms_in(&self, st: &CliqueState, m: u32, buf: &mut ScoreBuf, skip_deleted: bool, shared: &[Option<FineShared>], f: usize) -> (f64, bool) {
        let mv = &self.movers[m as usize];
        let cfg = self.cfg_of(st, m);
        let density = self.p.probe.density;
        let mut tot = 0.0;
        let mut bad = false;
        for slot in 0..mv.n_moved {
            if skip_deleted && self.mover_atom_deleted(m, slot, cfg) {
                continue;
            }
            let a = mv.atoms[slot];
            let r = match shared.get(slot).and_then(Option::as_ref) {
                Some(sh) => {
                    let mut r = sum_addends(&sh.adds, &self.dots_soa[&a], sh.pa[f], &sh.trims[f]);
                    r.bump /= density;
                    r.hbond /= density;
                    r.attract /= density;
                    r
                }
                None => self.score_atom_ext(st, a, buf, None),
            };
            tot += r.total();
            bad |= r.has_bad_bump;
        }
        (tot, bad)
    }

    fn pref(&self, m: u32, cfg: Cfg) -> f64 {
        let mv = &self.movers[m as usize];
        self.p.preference_magnitude
            * match cfg {
                Cfg::Coarse(c) => mv.coarse_pref[c as usize],
                Cfg::Fine(c, f) => mv.fine.pref(c as usize, f as usize),
            }
    }
}

/// A target prepared for one source atom state: the ion and alternate
/// filters applied, and what a dot test needs precomputed (the probe reach and
/// the annular test's surface point and kissing distance), with the same
/// arithmetic as the per-dot test.
#[derive(Clone, Copy)]
struct PTarget {
    pos: Vec3,
    info: AtomInfo,
    prb2: f64,
    surface: Vec3,
    kiss: f64,
}

fn prepare_targets(s: &DotScorer, src: &AtomInfo, src_pos: Vec3, targets: &[Target], probe_radius: f64, out: &mut Vec<PTarget>) {
    out.clear();
    let src_r = src.vdw_radius;
    for b in targets {
        let bi = &b.info;
        if bi.is_ion && s.ignore_ions {
            continue;
        }
        if !compatible_alts(src.alt, bi.alt) {
            continue;
        }
        let prb = bi.vdw_radius + probe_radius;
        let tv = b.pos - src_pos;
        let surface = tv / tv.length() * src_r + src_pos;
        let tr = bi.vdw_radius;
        let kiss = 2.0 * src_r * (tr * probe_radius / ((src_r + tr) * (src_r + probe_radius))).sqrt();
        out.push(PTarget { pos: b.pos, info: *bi, prb2: prb * prb, surface, kiss });
    }
}

/// Running state of `check_dot` over a dot's targets: the target with the
/// smallest gap so far (first one wins ties) and its hydrogen-bond status.
#[derive(Clone, Copy)]
struct DotFold {
    best_gap: f64,
    best: Option<PTarget>,
    is_hbond: bool,
    too_close: bool,
    hb_min: f64,
    cause_dummy: bool,
}

impl DotFold {
    const START: DotFold =
        DotFold { best_gap: 1e100, best: None, is_hbond: false, too_close: false, hb_min: 0.0, cause_dummy: false };

    #[inline]
    fn step(&mut self, s: &DotScorer, src: &AtomInfo, dot_abs: Vec3, probe_loc: Vec3, b: &PTarget) {
        let bi = &b.info;
        if probe_loc.dist_sq(b.pos) > b.prb2 {
            return;
        }
        let gap = dot_abs.dist(b.pos) - bi.vdw_radius;
        if gap < self.best_gap {
            let cs = src.charge as i32;
            let cb = bi.charge as i32;
            let both = cs != 0 && cb != 0;
            let comp = both && cs * cb < 0;
            let could = (src.is_donor && bi.is_acceptor) || (src.is_acceptor && bi.is_donor);
            if could && (!both || comp) {
                self.is_hbond = true;
                self.hb_min = if both { s.max_charged_h_overlap } else { s.max_regular_h_overlap };
                self.too_close = gap < -self.hb_min;
            } else {
                if src.is_dummy_hydrogen || bi.is_dummy_hydrogen {
                    return;
                }
                self.is_hbond = false;
                self.too_close = false;
            }
            self.cause_dummy = bi.is_dummy_hydrogen;
            self.best_gap = gap;
            self.best = Some(*b);
        }
    }

    /// (overlap type, gap, overlap, annular) of the dot, or None without a target.
    #[inline]
    fn finish(&self, s: &DotScorer, src: &AtomInfo, dot_abs: Vec3) -> Option<(OverlapType, f64, f64, bool)> {
        let b = self.best?;
        Some(finish_dot(s, src, self.best_gap, self.is_hbond, self.too_close, self.hb_min, self.cause_dummy, &b, dot_abs))
    }
}

/// The end of `check_dot` for a dot whose nearest target is `b`.
#[allow(clippy::too_many_arguments)]
#[inline(always)]
fn finish_dot(
    s: &DotScorer,
    src: &AtomInfo,
    best_gap: f64,
    is_hbond: bool,
    too_close: bool,
    hb_min: f64,
    cause_dummy: bool,
    b: &PTarget,
    dot_abs: Vec3,
) -> (OverlapType, f64, f64, bool) {
    let mut gap = best_gap;
    let overlap;
    let mut ot;
    if gap > 0.0 {
        overlap = 0.0;
        ot = if s.weak_hbonds && is_hbond { OverlapType::HydrogenBond } else { OverlapType::NoOverlap };
    } else if is_hbond {
        if too_close {
            gap += hb_min;
            overlap = -0.5 * gap;
            ot = OverlapType::Clash;
        } else {
            overlap = -0.5 * gap;
            ot = OverlapType::HydrogenBond;
        }
    } else {
        overlap = -0.5 * gap;
        ot = OverlapType::Clash;
    }
    // annularDots: dot2srcCenter > kissEdge2bullsEye
    let d2s = (b.surface - dot_abs).length();
    let annular = d2s > b.kiss;
    if (src.is_dummy_hydrogen || cause_dummy) && (too_close || ot != OverlapType::HydrogenBond) {
        ot = OverlapType::Ignore;
    }
    (ot, gap, overlap, annular)
}

/// `DotScorer::check_dot` with a precomputed probe location.
#[inline(always)]
fn check_dot(s: &DotScorer, src: &AtomInfo, dot_abs: Vec3, probe_loc: Vec3, targets: &[PTarget]) -> Option<(OverlapType, f64, f64, bool)> {
    let mut f = DotFold::START;
    for b in targets {
        f.step(s, src, dot_abs, probe_loc, b);
    }
    f.finish(s, src, dot_abs)
}

/// Dots per chunk in `score_dot_lanes`.
const LANES: usize = 4;

/// The dot offsets of one dot sphere as coordinate arrays (dot and probe),
/// padded to whole chunks.
struct DotsSoa {
    n: usize,
    d: [Vec<f64>; 3],
    p: [Vec<f64>; 3],
    /// Per chunk: center and radius of a sphere holding its dot offsets, and
    /// of one holding its probe offsets.
    bounds: Vec<(Vec3, f64, Vec3, f64)>,
    /// The probe spheres of `bounds` as coordinate arrays and radii.
    probe_bounds: ChunkSpheres,
    /// The dot spheres of `bounds` as coordinate arrays and radii.
    dot_bounds: ChunkSpheres,
}

/// One bounding sphere per chunk, as coordinate arrays and radii.
struct ChunkSpheres {
    c: [Vec<f64>; 3],
    r: Vec<f64>,
}

impl ChunkSpheres {
    fn new(spheres: impl Iterator<Item = (Vec3, f64)>) -> ChunkSpheres {
        let mut out = ChunkSpheres { c: [Vec::new(), Vec::new(), Vec::new()], r: Vec::new() };
        for (c, r) in spheres {
            out.c[0].push(c.x);
            out.c[1].push(c.y);
            out.c[2].push(c.z);
            out.r.push(r);
        }
        out
    }

    /// Set bit `bit` of word `word` (of `nw` per chunk) in `masks` for each
    /// chunk whose sphere comes within `r` of `q`, a point given relative to
    /// the atom center. The spheres carry `CHUNK_MARGIN`, which keeps this a
    /// superset of the chunks that an exact test of their dots can pass.
    #[inline]
    fn mark(&self, q: Vec3, r: f64, masks: &mut [u64], nw: usize, word: usize, bit: u32) {
        let [cx, cy, cz] = &self.c;
        let n = self.r.len();
        let (cx, cy, cz, cr) = (&cx[..n], &cy[..n], &cz[..n], &self.r[..n]);
        if nw == 1 {
            let masks = &mut masks[..n];
            for c in 0..n {
                let (x, y, z) = (cx[c] - q.x, cy[c] - q.y, cz[c] - q.z);
                let reach = r + cr[c];
                masks[c] |= ((x * x + y * y + z * z <= reach * reach) as u64) << bit;
            }
        } else {
            for c in 0..n {
                let (x, y, z) = (cx[c] - q.x, cy[c] - q.y, cz[c] - q.z);
                let reach = r + cr[c];
                masks[c * nw + word] |= ((x * x + y * y + z * z <= reach * reach) as u64) << bit;
            }
        }
    }
}

/// Bit masks of chunk spheres reached: `nw` words of 64 bits per chunk, bit
/// `i` for item `i`.
#[derive(Default)]
struct ChunkMasks {
    nw: usize,
    bits: Vec<u64>,
}

impl ChunkMasks {
    fn reset(&mut self, nc: usize, items: usize) {
        self.nw = items.div_ceil(64).max(1);
        self.bits.clear();
        self.bits.resize(nc * self.nw, 0);
    }

    #[inline]
    fn mark(&mut self, spheres: &ChunkSpheres, q: Vec3, r: f64, item: usize) {
        spheres.mark(q, r, &mut self.bits, self.nw, item / 64, (item % 64) as u32);
    }

    #[inline(always)]
    fn chunk(&self, c: usize) -> &[u64] {
        &self.bits[c * self.nw..(c + 1) * self.nw]
    }
}

/// Call `f` with each item whose bit is set in `words`, in order.
#[inline(always)]
fn for_bits(words: &[u64], mut f: impl FnMut(usize)) {
    for (wi, &w) in words.iter().enumerate() {
        let mut bits = w;
        while bits != 0 {
            f(wi * 64 + bits.trailing_zeros() as usize);
            bits &= bits - 1;
        }
    }
}

/// The debugging switches read from the environment: REDUCE3_ATOMWISE
/// (score atoms one state at a time) and REDUCE3_DEBUG_ATOM (print the
/// coarse state scores of the singleton Mover of one atom). They are read
/// once.
struct DebugSwitches {
    atomwise: bool,
    debug_atom: Option<String>,
}

fn debug_switches() -> &'static DebugSwitches {
    static SWITCHES: std::sync::OnceLock<DebugSwitches> = std::sync::OnceLock::new();
    SWITCHES.get_or_init(|| DebugSwitches {
        atomwise: std::env::var_os("REDUCE3_ATOMWISE").is_some(),
        debug_atom: std::env::var("REDUCE3_DEBUG_ATOM").ok(),
    })
}

/// Slack on the chunk bounds (far above rounding error).
const CHUNK_MARGIN: f64 = 1e-6;

impl DotsSoa {
    fn new(dots: &[DotPair]) -> DotsSoa {
        let padded = dots.len().div_ceil(LANES) * LANES;
        let col = |f: &dyn Fn(&DotPair) -> f64| {
            let mut v: Vec<f64> = dots.iter().map(f).collect();
            v.resize(padded, 0.0);
            v
        };
        let sphere = |pts: &[Vec3]| {
            let c = pts.iter().fold(Vec3::ZERO, |a, &q| a + q) / pts.len() as f64;
            (c, pts.iter().map(|q| q.dist(c)).fold(0.0, f64::max))
        };
        let bounds: Vec<(Vec3, f64, Vec3, f64)> = dots
            .chunks(LANES)
            .map(|ch| {
                let (cd, rd) = sphere(&ch.iter().map(|x| x.d).collect::<Vec<_>>());
                let (cp, rp) = sphere(&ch.iter().map(|x| x.probe).collect::<Vec<_>>());
                (cd, rd + CHUNK_MARGIN, cp, rp + CHUNK_MARGIN)
            })
            .collect();
        let probe_bounds = ChunkSpheres::new(bounds.iter().map(|b| (b.2, b.3)));
        let dot_bounds = ChunkSpheres::new(bounds.iter().map(|b| (b.0, b.1)));
        DotsSoa {
            n: dots.len(),
            d: [col(&|x| x.d.x), col(&|x| x.d.y), col(&|x| x.d.z)],
            p: [col(&|x| x.probe.x), col(&|x| x.probe.y), col(&|x| x.probe.z)],
            bounds,
            probe_bounds,
            dot_bounds,
        }
    }

    fn dot(&self, k: usize) -> DotPair {
        DotPair { d: v3(self.d[0][k], self.d[1][k], self.d[2][k]), probe: v3(self.p[0][k], self.p[1][k], self.p[2][k]) }
    }
}

/// What `DotFold::step` needs of a target for one source atom.
#[derive(Clone, Copy)]
struct LaneTarget {
    pos: Vec3,
    prb2: f64,
    /// Probe reach (for the chunk cull).
    prb: f64,
    vdw: f64,
    hb: bool,
    hbm: f64,
    /// Fold state bits of the target: hydrogen bond, dummy cause, index.
    base: u64,
}

const LANE_NO_TARGET: u64 = u64::MAX;

/// Scratch space of `score_dot_lanes`, kept between calls.
#[derive(Default)]
struct LaneBufs {
    targets: Vec<LaneTarget>,
    /// Per chunk, the targets that can reach its probes.
    hits: ChunkMasks,
}

/// Fill `hits` with the chunks each of `lt` can reach from the atom at `pa`.
fn target_hits(dots: &DotsSoa, pa: Vec3, lt: &[LaneTarget], hits: &mut ChunkMasks) {
    hits.reset(dots.bounds.len(), lt.len());
    for (ti, t) in lt.iter().enumerate() {
        hits.mark(&dots.probe_bounds, t.pos - pa, t.prb, ti);
    }
}

/// The per-target constants of `DotFold::step` for `src`, in order; a target
/// the step leaves alone for every dot is dropped.
fn lane_targets(s: &DotScorer, src: &AtomInfo, targets: &[PTarget], lt: &mut Vec<LaneTarget>) {
    lt.clear();
    let cs = src.charge as i32;
    for (t, b) in targets.iter().enumerate() {
        let bi = &b.info;
        let cb = bi.charge as i32;
        let both = cs != 0 && cb != 0;
        let comp = both && cs * cb < 0;
        let could = (src.is_donor && bi.is_acceptor) || (src.is_acceptor && bi.is_donor);
        let hb = could && (!both || comp);
        if !hb && (src.is_dummy_hydrogen || bi.is_dummy_hydrogen) {
            continue;
        }
        let hbm = if both { s.max_charged_h_overlap } else { s.max_regular_h_overlap };
        let base = ((t as u64) << 8) | hb as u64 | ((bi.is_dummy_hydrogen as u64) << 2);
        lt.push(LaneTarget { pos: b.pos, prb2: b.prb2, prb: b.prb2.sqrt(), vdw: bi.vdw_radius, hb, hbm, base });
    }
}

/// The dots and probes of one chunk with the atom at `pa`.
struct Chunk {
    a: [[f64; LANES]; 3],
    q: [[f64; LANES]; 3],
}

impl Chunk {
    #[inline(always)]
    fn at(pa: Vec3, dots: &DotsSoa, start: usize) -> Chunk {
        let mut c = Chunk { a: [[0.0; LANES]; 3], q: [[0.0; LANES]; 3] };
        let pc = [pa.x, pa.y, pa.z];
        for ax in 0..3 {
            for l in 0..LANES {
                c.a[ax][l] = pc[ax] + dots.d[ax][start + l];
                c.q[ax][l] = pc[ax] + dots.p[ax][start + l];
            }
        }
        c
    }

    /// Lanes whose dot lies within (strictly) `r2` of `c`.
    #[inline(always)]
    fn dot_within(&self, c: Vec3, r2: f64) -> [bool; LANES] {
        let mut out = [false; LANES];
        for (l, o) in out.iter_mut().enumerate() {
            let (x, y, z) = (self.a[0][l] - c.x, self.a[1][l] - c.y, self.a[2][l] - c.z);
            *o = x * x + y * y + z * z < r2;
        }
        out
    }

    /// Lanes whose probe lies within `r2` of `c`.
    #[inline(always)]
    fn probe_within(&self, c: Vec3, r2: f64) -> [bool; LANES] {
        let mut out = [false; LANES];
        for (l, o) in out.iter_mut().enumerate() {
            let (x, y, z) = (self.q[0][l] - c.x, self.q[1][l] - c.y, self.q[2][l] - c.z);
            *o = x * x + y * y + z * z <= r2;
        }
        out
    }

    /// Live lanes of the first `m`: those outside every atom of `excl`. Only
    /// atoms that reach the chunk's dot sphere (center `dc`, radius `rd`) can
    /// remove one.
    #[inline(always)]
    fn live(&self, m: usize, excl: &[(Vec3, f64)], dc: Vec3, rd: f64) -> [bool; LANES] {
        let mut live = [false; LANES];
        for (l, v) in live.iter_mut().enumerate() {
            *v = l < m;
        }
        for &(pe, re) in excl {
            let reach = re + rd;
            if dc.dist_sq(pe) < reach * reach {
                let inside = self.dot_within(pe, re * re);
                for l in 0..LANES {
                    live[l] &= !inside[l];
                }
            }
        }
        live
    }

    #[inline(always)]
    fn dot(&self, l: usize) -> Vec3 {
        v3(self.a[0][l], self.a[1][l], self.a[2][l])
    }
}

/// `DotFold` of each lane of a chunk, as the target index and state bits of
/// `LaneTarget::base`.
#[derive(Clone, Copy)]
struct LaneFold {
    bg: [f64; LANES],
    st: [u64; LANES],
    hm: [f64; LANES],
}

impl LaneFold {
    const START: LaneFold = LaneFold { bg: [1e100; LANES], st: [LANE_NO_TARGET; LANES], hm: [0.0; LANES] };

    /// `DotFold::step` on the live lanes.
    #[inline(always)]
    fn step(&mut self, c: &Chunk, live: &[bool; LANES], t: &LaneTarget) {
        let mut within = c.probe_within(t.pos, t.prb2);
        for l in 0..LANES {
            within[l] &= live[l];
        }
        if !within.iter().any(|&x| x) {
            return;
        }
        let neg_hbm = -t.hbm;
        for l in 0..LANES {
            let (rx, ry, rz) = (c.a[0][l] - t.pos.x, c.a[1][l] - t.pos.y, c.a[2][l] - t.pos.z);
            let gap = (rx * rx + ry * ry + rz * rz).sqrt() - t.vdw;
            let upd = within[l] & (gap < self.bg[l]);
            let new_state = t.base | (((t.hb & (gap < neg_hbm)) as u64) << 1);
            self.bg[l] = if upd { gap } else { self.bg[l] };
            self.st[l] = if upd { new_state } else { self.st[l] };
            self.hm[l] = if upd & t.hb { t.hbm } else { self.hm[l] };
        }
    }

    /// `DotFold::finish` for lane `l`, which has a target (`targets` as
    /// indexed by `base`).
    #[inline(always)]
    fn finish(&self, l: usize, s: &DotScorer, src: &AtomInfo, targets: &[PTarget], dot_abs: Vec3) -> (OverlapType, f64, f64, bool) {
        let st = self.st[l];
        let b = &targets[(st >> 8) as usize];
        finish_dot(s, src, self.bg[l], st & 1 != 0, st & 2 != 0, self.hm[l], st & 4 != 0, b, dot_abs)
    }

    /// The `DotFold` of lane `l` (`targets` as indexed by `base`).
    #[inline(always)]
    fn fold(&self, l: usize, targets: &[PTarget]) -> DotFold {
        let st = self.st[l];
        if st == LANE_NO_TARGET {
            return DotFold::START;
        }
        DotFold {
            best_gap: self.bg[l],
            best: Some(targets[(st >> 8) as usize]),
            is_hbond: st & 1 != 0,
            too_close: st & 2 != 0,
            hb_min: self.hm[l],
            cause_dummy: st & 4 != 0,
        }
    }
}

/// Score the dots of one atom state against `targets` in order, as
/// `check_dot` and `accumulate_dot` do one dot at a time: dots inside an atom
/// of `excl` are skipped and the rest accumulate into `ret` in dot order.
/// Dots go through in chunks of `LANES`, each test over a chunk at once (the
/// compiler vectorizes the fixed-size loops) with the scalar arithmetic, so
/// the results are the same bits; each dot's fold keeps the target order.
#[allow(clippy::too_many_arguments)]
fn score_dot_lanes(
    s: &DotScorer,
    src: &AtomInfo,
    pa: Vec3,
    dots: &DotsSoa,
    excl: &[(Vec3, f64)],
    targets: &[PTarget],
    lanes: &mut LaneBufs,
    ret: &mut ScoreDotsResult,
) {
    fold_dot_lanes(s, src, pa, dots, excl, targets, lanes, |_, r| accumulate_dot(s, r, ret));
}

/// The chunked fold of `score_dot_lanes`, handing each dot that finds a
/// target to `emit` (dot index and result), in dot order.
#[allow(clippy::too_many_arguments)]
#[inline]
fn fold_dot_lanes(
    s: &DotScorer,
    src: &AtomInfo,
    pa: Vec3,
    dots: &DotsSoa,
    excl: &[(Vec3, f64)],
    targets: &[PTarget],
    lanes: &mut LaneBufs,
    mut emit: impl FnMut(usize, Option<(OverlapType, f64, f64, bool)>),
) {
    let LaneBufs { targets: lt, hits } = lanes;
    lane_targets(s, src, targets, lt);
    if lt.is_empty() {
        return;
    }
    // a chunk with no target in reach of its probes cannot score
    target_hits(dots, pa, lt, hits);
    for (ci, &(cd, rd, _, _)) in dots.bounds.iter().enumerate() {
        let near_t = hits.chunk(ci);
        if near_t.iter().all(|&w| w == 0) {
            continue;
        }
        let start = ci * LANES;
        let m = (dots.n - start).min(LANES);
        let c = Chunk::at(pa, dots, start);
        let live = c.live(m, excl, pa + cd, rd);
        if !live.iter().any(|&x| x) {
            continue;
        }
        let mut f = LaneFold::START;
        for_bits(near_t, |ti| f.step(&c, &live, &lt[ti]));
        // finish and accumulate in dot order
        for l in 0..m {
            if live[l] && f.st[l] != LANE_NO_TARGET {
                emit(start + l, Some(f.finish(l, s, src, targets, c.dot(l))));
            }
        }
    }
}

/// Add one dot's result to the running sums (the body of `score_dots`).
#[inline]
fn accumulate_dot(s: &DotScorer, r: Option<(OverlapType, f64, f64, bool)>, ret: &mut ScoreDotsResult) {
    if let Some((kind, val)) = dot_addend(s, r) {
        add_to(ret, kind, val);
    }
}

/// Which sum of `ScoreDotsResult` a dot adds to.
const ADD_ATTRACT: u8 = 0;
const ADD_BUMP: u8 = 1;
/// A bump that also sets `has_bad_bump`.
const ADD_BAD_BUMP: u8 = 2;
const ADD_HBOND: u8 = 3;

/// The sum one dot's result adds to and the amount, as `accumulate_dot`
/// adds it, or None when it adds nothing.
#[inline]
fn dot_addend(s: &DotScorer, r: Option<(OverlapType, f64, f64, bool)>) -> Option<(u8, f64)> {
    let (ot, gap, overlap, annular) = r?;
    match s.interaction_type(ot, gap, true) {
        InteractionType::Invalid => None,
        InteractionType::WideContact | InteractionType::CloseContact | InteractionType::WeakHydrogenBond => {
            if annular {
                return None;
            }
            let sg = gap / s.gap_scale;
            Some((ADD_ATTRACT, (-sg * sg).exp()))
        }
        InteractionType::SmallOverlap | InteractionType::Bump => Some((ADD_BUMP, -s.bump_weight * overlap)),
        InteractionType::BadBump => Some((ADD_BAD_BUMP, -s.bump_weight * overlap)),
        InteractionType::StandardHydrogenBond => Some((ADD_HBOND, s.hbond_weight * overlap)),
    }
}

#[inline]
fn add_to(ret: &mut ScoreDotsResult, kind: u8, val: f64) {
    match kind {
        ADD_ATTRACT => ret.attract += val,
        ADD_HBOND => ret.hbond += val,
        _ => {
            ret.bump += val;
            ret.has_bad_bump |= kind == ADD_BAD_BUMP;
        }
    }
}

/// A Mover atom's dots evaluated once for all of its fine poses: the addends
/// of the dots that find a target, and per pose the atom's place and the own
/// Mover's excluded atoms that can trim them.
struct FineShared {
    adds: Vec<DotAdd>,
    pa: Vec<Vec3>,
    trims: Vec<Vec<(Vec3, f64)>>,
}

/// One dot's addend to an atom score: dot index, sum and amount.
#[derive(Clone, Copy)]
struct DotAdd {
    k: u32,
    kind: u8,
    val: f64,
}

/// An atom's score at `pa` from its dot addends, leaving out the dots inside
/// an atom of `trim`, added in dot order as `score_dot_lanes` adds them.
fn sum_addends(adds: &[DotAdd], dots: &DotsSoa, pa: Vec3, trim: &[(Vec3, f64)]) -> ScoreDotsResult {
    let mut ret = ScoreDotsResult::default();
    let [dx, dy, dz] = &dots.d;
    for ad in adds {
        let k = ad.k as usize;
        let q = v3(pa.x + dx[k], pa.y + dy[k], pa.z + dz[k]);
        if !trim.iter().any(|&(pe, re)| q.dist_sq(pe) < re * re) {
            add_to(&mut ret, ad.kind, ad.val);
        }
    }
    ret
}

/// Largest distance between an atom's positions in several poses for one
/// evaluation of its dots to serve them all: far above the rounding error of
/// a pose that should leave the atom in place, and far below any distance
/// that moves a score by more than rounding.
const SHARE_TOL: f64 = 1e-9;

/// A sphere through which another Mover of the clique reaches a dot: one of
/// its excluded atoms in some state, tested on the dot (`probe` false, strict
/// test), or one of its targets in some state, tested on the dot's probe.
#[derive(Clone, Copy)]
struct DepSphere {
    /// Bit of the Mover in the atom's list of other Movers.
    k: u32,
    probe: bool,
    c: Vec3,
    r: f64,
    r2: f64,
}

/// Scratch space of `classify_dots`, kept between calls.
#[derive(Default)]
struct ClassifyBufs {
    lanes: LaneBufs,
    /// Per chunk, the spheres that can reach it.
    shits: ChunkMasks,
}

/// Sort the dots of one atom state for `atom_factors_dotwise`, in dot order,
/// with the chunked tests of `fold_dot_lanes`. `targets` holds the static
/// targets (the first `n_static`) and then the own Mover's. A dot inside an
/// atom of `excl` drops out. A dot that no sphere of `spheres` reaches goes
/// to `base` with its result over all of `targets`; any other goes to `group`
/// with the bits (`DepSphere::k`) of the Movers that reach it and its fold
/// over the static targets alone.
#[allow(clippy::too_many_arguments)]
fn classify_dots(
    s: &DotScorer,
    src: &AtomInfo,
    pa: Vec3,
    dots: &DotsSoa,
    excl: &[(Vec3, f64)],
    targets: &[PTarget],
    n_static: usize,
    spheres: &[DepSphere],
    bufs: &mut ClassifyBufs,
    mut base: impl FnMut(usize, Option<(OverlapType, f64, f64, bool)>),
    mut group: impl FnMut(u64, usize, DotFold),
) {
    let ClassifyBufs { lanes: LaneBufs { targets: lt, hits }, shits } = bufs;
    lane_targets(s, src, targets, lt);
    // the chunks each target and each sphere can reach
    target_hits(dots, pa, lt, hits);
    shits.reset(dots.bounds.len(), spheres.len());
    for (si, sp) in spheres.iter().enumerate() {
        let b = if sp.probe { &dots.probe_bounds } else { &dots.dot_bounds };
        shits.mark(b, sp.c - pa, sp.r, si);
    }
    for (ci, &(cd, rd, _, _)) in dots.bounds.iter().enumerate() {
        // a chunk that meets neither a target nor another Mover adds nothing
        let (near_t, near_s) = (hits.chunk(ci), shits.chunk(ci));
        if near_t.iter().chain(near_s).all(|&w| w == 0) {
            continue;
        }
        let start = ci * LANES;
        let m = (dots.n - start).min(LANES);
        let c = Chunk::at(pa, dots, start);
        let live = c.live(m, excl, pa + cd, rd);
        if !live.iter().any(|&x| x) {
            continue;
        }
        let mut dep = [0u64; LANES];
        for_bits(near_s, |si| {
            let sp = &spheres[si];
            let hit = if sp.probe { c.probe_within(sp.c, sp.r2) } else { c.dot_within(sp.c, sp.r2) };
            for l in 0..LANES {
                dep[l] |= ((hit[l] & live[l]) as u64) << sp.k;
            }
        });
        // the static targets, then the own Mover's
        let mut f = LaneFold::START;
        let mut static_fold = f;
        let mut own_started = false;
        for_bits(near_t, |ti| {
            let t = &lt[ti];
            if !own_started && (t.base >> 8) as usize >= n_static {
                static_fold = f;
                own_started = true;
            }
            f.step(&c, &live, t);
        });
        if !own_started {
            static_fold = f;
        }
        for l in 0..m {
            if !live[l] {
                continue;
            }
            if dep[l] == 0 {
                if f.st[l] != LANE_NO_TARGET {
                    base(start + l, Some(f.finish(l, s, src, targets, c.dot(l))));
                }
            } else {
                group(dep[l], start + l, static_fold.fold(l, targets));
            }
        }
    }
}

#[derive(Default)]
struct ScoreBuf {
    targets: Vec<Target>,
    prepared: Vec<PTarget>,
    lanes: LaneBufs,
    excl: Vec<(Vec3, f64)>,
    dots: Vec<DotPair>,
}

// ----------------------------------------------------------------------------
// Exact clique optimization by variable elimination

struct Factor {
    scope: Vec<usize>,
    dims: Vec<usize>,
    table: Vec<f64>,
}

impl Factor {
    fn size(dims: &[usize]) -> usize {
        dims.iter().product()
    }
}

/// Index of an assignment (`vals` indexed by local variable) in a factor table.
#[inline]
fn table_index(scope: &[usize], dims: &[usize], vals: &[usize]) -> usize {
    let mut idx = 0;
    for (k, &v) in scope.iter().enumerate() {
        idx = idx * dims[k] + vals[v];
    }
    idx
}

/// Maximize the sum of factors; returns the best assignment.
fn variable_elimination(nvars: usize, doms: &[usize], mut factors: Vec<Factor>) -> Vec<usize> {
    let mut eliminated = vec![false; nvars];
    // (var, remaining scope, dims, argmax table)
    let mut trace: Vec<(usize, Vec<usize>, Vec<usize>, Vec<u16>)> = Vec::new();
    for _ in 0..nvars {
        // choose variable producing the smallest combined table
        let mut best_v = usize::MAX;
        let mut best_cost = usize::MAX;
        for v in 0..nvars {
            if eliminated[v] {
                continue;
            }
            let mut uni: Vec<usize> = vec![v];
            for f in &factors {
                if f.scope.contains(&v) {
                    for &s in &f.scope {
                        if !uni.contains(&s) {
                            uni.push(s);
                        }
                    }
                }
            }
            let cost: usize = uni.iter().map(|&s| doms[s]).fold(1usize, |a, b| a.saturating_mul(b));
            if cost < best_cost {
                best_cost = cost;
                best_v = v;
            }
        }
        let v = best_v;
        eliminated[v] = true;
        let (with_v, rest): (Vec<Factor>, Vec<Factor>) = factors.into_iter().partition(|f| f.scope.contains(&v));
        factors = rest;
        let mut uni: Vec<usize> = Vec::new();
        for f in &with_v {
            for &s in &f.scope {
                if s != v && !uni.contains(&s) {
                    uni.push(s);
                }
            }
        }
        uni.sort_unstable();
        let udims: Vec<usize> = uni.iter().map(|&s| doms[s]).collect();
        let n = Factor::size(&udims);
        let mut table = vec![f64::NEG_INFINITY; n];
        let mut arg = vec![0u16; n];
        let mut vals = vec![0usize; nvars];
        for idx in 0..n {
            // decode idx into vals for the union scope
            let mut r = idx;
            for k in (0..uni.len()).rev() {
                vals[uni[k]] = r % udims[k];
                r /= udims[k];
            }
            let mut best = f64::NEG_INFINITY;
            let mut bi = 0u16;
            for sv in 0..doms[v] {
                vals[v] = sv;
                let mut s = 0.0;
                for f in &with_v {
                    s += f.table[table_index(&f.scope, &f.dims, &vals)];
                }
                if s > best || sv == 0 {
                    if sv == 0 || s > best {
                        best = s;
                        bi = sv as u16;
                    }
                }
            }
            table[idx] = best;
            arg[idx] = bi;
        }
        trace.push((v, uni.clone(), udims.clone(), arg));
        factors.push(Factor { scope: uni, dims: udims, table });
    }
    let mut assign = vec![0usize; nvars];
    for (v, scope, dims, arg) in trace.into_iter().rev() {
        let idx = table_index(&scope, &dims, &assign);
        assign[v] = arg[idx] as usize;
    }
    assign
}

// ----------------------------------------------------------------------------
// The optimizer

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
    };
    place_movers(w, p, &atoms, rotatable_h, flip_states, alt, max_vdw, &mut nidx, &removed, n_models, &mut placement);
    info += &placement.info;
    let movers = placement.movers;
    let mut mover_info = placement.mover_info;
    out.amide_flips.extend(placement.amide_flips);
    out.his_flips.extend(placement.his_flips);
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

struct CliqueResult {
    initial: Vec<f64>,
    cfg: Vec<Cfg>,
    high: Vec<f64>,
    coarse_best: f64,
    coarse_info: String,
    fine_info: Vec<(u32, String)>,
    annot: Vec<String>,
    calculated: usize,
    cached: usize,
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
const NEAR_MARGIN: f64 = 1e-3;

/// Slack on the distance beyond which an excluded atom cannot cover a dot
/// (far above rounding error).
const EXCLUSION_MARGIN: f64 = 0.01;

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

// ----------------------------------------------------------------------------
// Placement (`Optimizer._PlaceMovers`)

struct Placement {
    movers: Vec<Mover>,
    mover_info: Vec<String>,
    delete_atoms: Vec<u32>,
    info: String,
    amide_flips: Vec<FlippedMoverInfo>,
    his_flips: Vec<FlippedMoverInfo>,
}

#[allow(clippy::too_many_arguments)]
fn place_movers(
    w: &mut World,
    p: &OptParams,
    atoms: &[u32],
    rotatable_h: &[u32],
    flip_states: &[FlipMoverState],
    alt: &str,
    max_vdw: f64,
    q: &mut NeighborIndex,
    removed: &FxHashSet<u32>,
    n_models: usize,
    pl: &mut Placement,
) {
    let v = p.verbosity;
    let clamp = !p.compat;
    let polar_h_radius = 1.05;
    for &a in atoms {
        if a as usize >= w.n_real {
            continue;
        }
        // the tests below on borrowed labels: an atom that passes none adds
        // nothing; the labels stay the same while Movers move atoms, so the
        // messages read them again
        let (is_n, is_c, asx_glx, amide, his) = {
            let (name, res, elem) = (w.name(a), w.resname(a), w.elem(a));
            (
                elem == "N",
                elem == "C",
                (name == "XD2" && res == "ASX") || (name == "XE2" && res == "GLX"),
                (name == "ND2" && res == "ASN") || (name == "NE2" && res == "GLN"),
                name == "NE2" && res == "HIS",
            )
        };
        let four = w.bonded[a as usize].len() == 4;
        if !((four && (is_n || is_c)) || (p.add_flip_movers && (asx_glx || amide)) || his) {
            continue;
        }
        let rid = |w: &World| w.res_name_and_id(a);
        let rep = p.report;

        if is_n && four {
            let num_h = w.bonded[a as usize].iter().filter(|&&n| w.is_h(n)).count();
            if num_h == 3 {
                match movers::nh3_rotator(w, a, !p.compat) {
                    Ok(m) => {
                        pl.movers.push(m);
                        if rep {
                            pl.info += &vcheck(v, 1, &format!("Added MoverNH3Rotator {} to {}\n", pl.movers.len(), rid(w)));
                            pl.mover_info.push(format!("NH3Rotator at {} {}", rid(w), w.name(a)));
                        } else {
                            pl.mover_info.push(String::new());
                        }
                    }
                    Err(e) => pl.info += &vcheck(v, 0, &format!("Could not add MoverNH3Rotator to {}: {}\n", rid(w), e)),
                }
            }
        }
        if is_c && four {
            let mut num_h = 0;
            let mut neighbor = None;
            for &n in &w.bonded[a as usize] {
                if w.is_h(n) {
                    num_h += 1;
                } else {
                    neighbor = Some(n);
                }
            }
            if num_h == 3 {
                let neighbor = neighbor.unwrap();
                if w.bonded[neighbor as usize].len() == 3 {
                    match movers::aromatic_methyl_rotator(w, a, !p.compat) {
                        Ok(m) => {
                            pl.movers.push(m);
                            if rep {
                                pl.info += &vcheck(
                                    v,
                                    1,
                                    &format!("Added MoverAromaticMethylRotator {} to {} {}\n", pl.movers.len(), rid(w), w.name(a)),
                                );
                                pl.mover_info.push(format!("AromaticMethylRotator at {} {}", rid(w), w.name(a)));
                            } else {
                                pl.mover_info.push(String::new());
                            }
                        }
                        Err(e) => {
                            pl.info += &vcheck(v, 0, &format!("Could not add MoverAromaticMethylRotator to {} {}: {}\n", rid(w), w.name(a), e))
                        }
                    }
                } else {
                    match movers::stagger_tetrahedral_methyl(w, a, !p.compat) {
                        Ok(()) => {
                            if rep {
                                pl.info += &vcheck(v, 1, &format!("Used MoverTetrahedralMethylRotator to stagger {} {}\n", rid(w), w.name(a)))
                            }
                        }
                        Err(e) => {
                            pl.info += &vcheck(v, 0, &format!("Could not add MoverTetrahedralMethylRotator to {} {}: {}\n", rid(w), w.name(a), e))
                        }
                    }
                }
            }
        }
        if p.add_flip_movers && asx_glx && rep {
            pl.info += &vcheck(v, 1, &format!("Not attempting to adjust {} {}\n", rid(w), w.name(a)));
        }
        if p.add_flip_movers && amide {
            let mut found_ion = false;
            let mut oxygen = None;
            for &b in &w.bonded[a as usize] {
                if w.elem(b) == "C" {
                    for &b2 in &w.bonded[b as usize] {
                        if w.elem(b2) == "O" {
                            oxygen = Some(b2);
                        }
                    }
                }
            }
            let mut lock_flipped = false;
            if let Some(o) = oxygen {
                if p.compat {
                    // Original: only the unflipped O position, with N's radius.
                    let my_rad = w.info[a as usize].vdw_radius;
                    let max_d = 0.25 + my_rad + max_vdw;
                    for n in q.neighbors(w, w.pos[o as usize], my_rad, max_d, removed) {
                        if w.is_positive_ion(n) {
                            let dist = w.pos[o as usize].dist(w.pos[n as usize]);
                            let expected = my_rad + w.info[n as usize].vdw_radius;
                            if dist >= expected - 0.65 && dist <= expected + 0.25 {
                                found_ion = true;
                            }
                        }
                    }
                } else {
                    // Fixed: test the O in both of its possible positions (its own and
                    // the N's) with the O radius; lock the flip so O faces the ion.
                    let o_rad = w.info[o as usize].vdw_radius;
                    let max_d = 0.25 + o_rad + max_vdw;
                    let ion_near = |pos: Vec3| -> bool {
                        q.neighbors(w, pos, o_rad, max_d, removed).into_iter().any(|n| {
                            w.is_positive_ion(n) && {
                                let dist = pos.dist(w.pos[n as usize]);
                                let expected = o_rad + w.info[n as usize].vdw_radius;
                                dist >= expected - 0.65 && dist <= expected + 0.25
                            }
                        })
                    };
                    if ion_near(w.pos[o as usize]) {
                        found_ion = true;
                    } else if ion_near(w.pos[a as usize]) {
                        found_ion = true;
                        lock_flipped = true;
                    }
                }
            }
            if found_ion && lock_flipped {
                match movers::amide_flip(w, a, "CA", p.non_flip_preference, clamp) {
                    Ok(flip) => {
                        pl.info += &vcheck(v, 1, &format!("Locking MoverAmideFlip on {} flipped (ionic contact)\n", rid(w)));
                        set_state_now(w, &flip, 1, q);
                        apply_fixup(w, &flip, 1, &mut pl.delete_atoms);
                        pl.amide_flips.push(FlippedMoverInfo { alt: alt.into(), base_atom: a });
                    }
                    Err(e) => pl.info += &vcheck(v, 0, &format!("Did not add MoverAmideFlip to {}: {}\n", rid(w), e)),
                }
            } else if !found_ion {
                let fs = find_flip_state(w, a, flip_states, p.compat, n_models);
                match movers::amide_flip(w, a, "CA", p.non_flip_preference, clamp) {
                    Ok(flip) => {
                        if let Some(s) = fs {
                            pl.info += &vcheck(
                                v,
                                1,
                                &format!(
                                    "Setting MoverAmideFlip for {}: flipped = {}, angles adjusted = {}\n",
                                    rid(w),
                                    if s.flipped { "True" } else { "False" },
                                    if s.fixed_up { "True" } else { "False" }
                                ),
                            );
                            let index = if s.flipped { 1 } else { 0 };
                            if s.flipped {
                                pl.amide_flips.push(FlippedMoverInfo { alt: alt.into(), base_atom: a });
                            }
                            set_state_now(w, &flip, index, q);
                            if s.fixed_up {
                                apply_fixup(w, &flip, index, &mut pl.delete_atoms);
                            }
                        } else {
                            pl.movers.push(flip);
                            if rep {
                                pl.info += &vcheck(v, 1, &format!("Added MoverAmideFlip {} to {}\n", pl.movers.len(), rid(w)));
                                pl.mover_info.push(format!("AmideFlip at {} {}", rid(w), w.name(a)));
                            } else {
                                pl.mover_info.push(String::new());
                            }
                        }
                    }
                    Err(e) => pl.info += &vcheck(v, 0, &format!("Did not add MoverAmideFlip to {}: {}\n", rid(w), e)),
                }
            }
        }
        if his {
            let rid = rid(w);
            place_his(w, p, a, &rid, flip_states, alt, max_vdw, q, removed, n_models, pl);
        }
    }

    // single-hydrogen rotators
    let mut in_atoms = vec![false; w.len()];
    for &a in atoms {
        in_atoms[a as usize] = true;
    }
    let opts = SingleHOptions {
        circular_angle_spacing: !p.compat,
        any_partner_valence: !p.compat,
        constrain_planar_hydroxyls: !p.compat,
        planar_preference: if p.compat { 0.0 } else { p.planar_hydroxyl_preference },
        syn_preference: if p.compat { 0.0 } else { p.acid_syn_preference },
    };
    for &h in rotatable_h {
        if !in_atoms.get(h as usize).copied().unwrap_or(false) {
            continue;
        }
        let Some(&neighbor) = w.bonded[h as usize].first() else {
            let (a_name, rid) = (w.name(h), w.res_name_and_id(h));
            pl.info += &vcheck(v, 0, &format!("Could not add MoverSingleHydrogenRotator to {} {}: list index out of range\n", rid, a_name));
            continue;
        };
        if w.bonded[neighbor as usize].len() != 2 {
            continue;
        }
        let pr = p.probe.probe_radius;
        let bonded = atoms_within_n_bonds(w, h, pr, p.bonded_neighbor_depth, 3);
        let nrad = w.info[neighbor as usize].vdw_radius;
        let nearby = q.neighbors(w, w.pos[neighbor as usize], nrad, 4.0, removed);
        let xh = if w.elem(neighbor) == "S" { 1.3 } else { 1.0 };
        let mut touches = Vec::new();
        let mut cands = Vec::new();
        for n in nearby {
            if bonded.contains(&n) {
                continue;
            }
            touches.push(n);
            let d = w.pos[neighbor as usize].dist(w.pos[n as usize]);
            if d <= xh + w.info[n as usize].vdw_radius + polar_h_radius {
                cands.push(n);
            }
        }
        let mut acceptors = Vec::new();
        for c in cands {
            let cn = w.name(c);
            let rn = w.resname(c);
            let flip_partner = (cn == "ND2" && rn == "ASN")
                || (cn == "NE2" && rn == "GLN")
                || (cn == "CE1" && rn == "HIS")
                || (cn == "CD2" && rn == "HIS");
            if w.info[c as usize].is_acceptor || flip_partner {
                acceptors.push(c);
            }
        }
        match movers::single_hydrogen_rotator(w, h, &acceptors, &touches, &opts) {
            Ok(m) => {
                pl.movers.push(m);
                if p.report {
                    let (a_name, rid) = (w.name(h), w.res_name_and_id(h));
                    pl.info += &vcheck(
                        v,
                        1,
                        &format!(
                            "Added MoverSingleHydrogenRotator {} to {} {} with {} potential nearby acceptors\n",
                            pl.movers.len(),
                            rid,
                            a_name,
                            acceptors.len()
                        ),
                    );
                    pl.mover_info.push(format!("SingleHydrogenRotator at {} {}", rid, a_name));
                } else {
                    pl.mover_info.push(String::new());
                }
            }
            Err(e) => {
                let (a_name, rid) = (w.name(h), w.res_name_and_id(h));
                pl.info += &vcheck(v, 0, &format!("Could not add MoverSingleHydrogenRotator to {} {}: {}\n", rid, a_name, e))
            }
        }
    }
}

/// Move atoms into a coarse state right away (placement-time lock-down).
fn set_state_now(w: &mut World, mv: &Mover, c: usize, q: &mut NeighborIndex) {
    // `Optimizer._setMoverState`: remove, move, re-add each moved atom.
    for slot in 0..mv.coarse_pos[c].len() {
        let a = mv.atoms[slot];
        q.remove(w, a);
        w.pos[a as usize] = mv.coarse_pos[c][slot];
        q.add(w, a);
    }
    for (slot, inf) in in_state(&mv.coarse_info, c).iter().enumerate() {
        w.info[mv.atoms[slot] as usize] = *inf;
    }
    for (slot, &d) in in_state(&mv.coarse_del, c).iter().enumerate() {
        if d {
            q.remove(w, mv.atoms[slot]);
        } else {
            q.add(w, mv.atoms[slot]);
        }
    }
}

fn apply_fixup(w: &mut World, mv: &Mover, c: usize, deletes: &mut Vec<u32>) {
    for (i, &pp) in in_state(&mv.fixup_pos, c).iter().enumerate() {
        w.pos[mv.atoms[i] as usize] = pp;
    }
    for (i, inf) in in_state(&mv.fixup_info, c).iter().enumerate() {
        w.info[mv.atoms[i] as usize] = *inf;
    }
    for (i, &d) in in_state(&mv.fixup_del, c).iter().enumerate() {
        if d && !deletes.contains(&mv.atoms[i]) {
            deletes.push(mv.atoms[i]);
        }
    }
}

/// `_HisRingNitrogensWithNonstandardBond`.
fn his_nonstandard(w: &World, ne2: u32) -> Vec<(u32, Option<u32>)> {
    let ag = w.labels[ne2 as usize].ag;
    let mut ret = Vec::new();
    // atoms of the same atom group named ND1/NE2, in hierarchy order
    let lo = (ne2 as usize).saturating_sub(64);
    let hi = (ne2 as usize + 64).min(w.n_real);
    for at in lo..hi {
        let at = at as u32;
        if w.labels[at as usize].ag != ag {
            continue;
        }
        let nm = w.name(at);
        if nm != "ND1" && nm != "NE2" {
            continue;
        }
        let mut hyd = None;
        let mut nonstd = false;
        for &n in &w.bonded[at as usize] {
            if w.is_h(n) {
                hyd = Some(n);
            } else if w.elem(n) != "C" {
                nonstd = true;
            }
        }
        if nonstd {
            ret.push((at, hyd));
        }
    }
    ret
}

#[allow(clippy::too_many_arguments)]
fn place_his(
    w: &mut World,
    p: &OptParams,
    a: u32,
    rid: &str,
    flip_states: &[FlipMoverState],
    alt: &str,
    max_vdw: f64,
    q: &mut NeighborIndex,
    removed: &FxHashSet<u32>,
    n_models: usize,
    pl: &mut Placement,
) {
    let v = p.verbosity;
    let clamp = !p.compat;
    let a_name = w.name(a).to_string();
    let hist = match movers::his_flip(w, a, p.non_flip_preference, 3, true, clamp, p.compat) {
        Ok(h) => h,
        Err(e) => {
            pl.info += &vcheck(v, 0, &format!("Did not add MoverHisFlip to {}: {}\n", rid, e));
            return;
        }
    };
    let cp = &hist.coarse_pos;
    let checks = [cp[0][0], cp[0][4], cp[4][0], cp[4][4]];
    let my_rad = w.info[a as usize].vdw_radius;
    let max_d = 0.25 + my_rad + max_vdw;
    let mut bonded_config: Option<usize> = None;
    for (i, pos) in checks.iter().enumerate() {
        for n in q.neighbors(w, *pos, my_rad, max_d, removed) {
            if w.is_positive_ion(n) {
                let dist = pos.dist(w.pos[n as usize]);
                let expected = my_rad + w.info[n as usize].vdw_radius;
                if dist >= expected - 0.55 && dist <= expected + 0.25 {
                    bonded_config = Some((i / 2) * 4);
                    break;
                }
            }
        }
        if bonded_config.is_some() {
            break;
        }
    }
    if let Some(bc) = bonded_config {
        let coarse_positions = hist.coarse_pos[bc].clone();
        for (i, &at) in hist.atoms.iter().enumerate() {
            if i < hist.fixup_pos[bc].len() {
                w.pos[at as usize] = hist.fixup_pos[bc][i];
            }
            if i < hist.fixup_info[bc].len() {
                w.info[at as usize] = hist.fixup_info[bc][i];
            }
        }
        if bc == 4 {
            // The original records the wrong atom here (its loop variable was
            // reused); record the NE2 like every other flip.
            pl.his_flips.push(FlippedMoverInfo { alt: alt.into(), base_atom: if p.compat { *hist.atoms.last().unwrap() } else { a } });
        }
        let modify = |w: &mut World, nitro: u32, coarse_pos: Vec3, hydro: u32, pl: &mut Placement| {
            let my_rad = w.info[nitro as usize].vdw_radius;
            let max_d = 0.25 + my_rad + max_vdw;
            for n in q.neighbors(w, coarse_pos, my_rad, max_d, removed) {
                if w.is_positive_ion(n) {
                    let dist = coarse_pos.dist(w.pos[n as usize]);
                    let expected = my_rad + w.info[n as usize].vdw_radius;
                    if dist >= expected - 0.55 && dist <= expected + 0.25 {
                        pl.info += &vcheck(
                            v,
                            1,
                            &format!(
                                "Not adding Hydrogen to {}{} and marking as an acceptor (ionic bond to {})\n",
                                rid,
                                w.labels[nitro as usize].raw_name,
                                w.name(n)
                            ),
                        );
                        w.info[nitro as usize].is_acceptor = true;
                        pl.delete_atoms.push(hydro);
                        break;
                    }
                }
            }
        };
        modify(w, hist.atoms[0], coarse_positions[0], hist.atoms[1], pl);
        modify(w, hist.atoms[4], coarse_positions[4], hist.atoms[5], pl);
        pl.info += &vcheck(v, 1, &format!("Set MoverHisFlip on {} to state {}\n", rid, bc));
    } else if p.add_flip_movers && !his_nonstandard(w, a).is_empty() {
        pl.info += &vcheck(v, 1, &format!("Did not add MoverHisFlip to {} (bond outside standard ring connectivity)\n", rid));
        for (nitro, hydro) in his_nonstandard(w, a) {
            w.info[nitro as usize].is_acceptor = true;
            if let Some(h) = hydro {
                pl.delete_atoms.push(h);
                pl.info += &vcheck(
                    v,
                    1,
                    &format!(
                        "Not adding Hydrogen to {}{} and marking as an acceptor (bond outside standard ring connectivity)\n",
                        rid, w.labels[nitro as usize].raw_name
                    ),
                );
            }
        }
    } else if p.add_flip_movers {
        if let Some(s) = find_flip_state(w, a, flip_states, p.compat, n_models) {
            let enabled = if s.flipped { 2 } else { 1 };
            if s.flipped {
                pl.his_flips.push(FlippedMoverInfo { alt: alt.into(), base_atom: a });
            }
            match movers::his_flip(w, a, p.non_flip_preference, enabled, s.fixed_up, clamp, p.compat) {
                Ok(h) => {
                    set_state_now(w, &h, 0, q);
                    pl.movers.push(h);
                    pl.info += &vcheck(v, 1, &format!("Added MoverHisPlace {} to {}\n", pl.movers.len(), rid));
                    pl.mover_info.push(format!("HisPlace at {} {}", rid, a_name));
                }
                Err(e) => pl.info += &vcheck(v, 0, &format!("Did not add MoverHisFlip to {}: {}\n", rid, e)),
            }
        } else {
            pl.movers.push(hist);
            pl.info += &vcheck(v, 1, &format!("Added MoverHisFlip {} to {}\n", pl.movers.len(), rid));
            pl.mover_info.push(format!("HisFlip at {} {}", rid, a_name));
        }
    }
}

#[allow(dead_code)]
fn _unused(_: ResClass) {}

// ----------------------------------------------------------------------------
// Compat mode: a faithful port of Reduce2's OptimizerC (sequential; vertex-cut
// clique search; trimmed-dot and atom-score caches with the same keys and
// evaluation order). Used to verify Reduce3 against Reduce2 exactly.

type TrimCache = std::cell::RefCell<FxHashMap<(u32, u32), Vec<DotPair>>>;

struct CompatOpt<'c, 'a> {
    ctx: &'c Ctx<'a>,
    st: CliqueState<'c>,
    coarse_loc: Vec<usize>,
    high: Vec<f64>,
    trim: TrimCache,
    score_cache: Option<FxHashMap<(u32, SmallVec<[u16; 6]>), f64>>,
    buf: ScoreBuf,
    calculated: usize,
    cached: usize,
    edges: &'c [(usize, usize)],
}

impl<'c, 'a> CompatOpt<'c, 'a> {
    fn mag(&self) -> f64 {
        self.ctx.p.preference_magnitude
    }

    fn score_atom(&mut self, a: u32, loc: u32) -> f64 {
        self.calculated += 1;
        self.ctx.score_atom_ext(&self.st, a, &mut self.buf, Some((&self.trim, loc))).total()
    }

    fn score_atom_cached(&mut self, a: u32, loc: u32) -> f64 {
        let key: SmallVec<[u16; 6]> = match self.ctx.atom_movers.get(&a) {
            Some(ms) => ms.iter().map(|&m| self.coarse_loc[m as usize] as u16).collect(),
            None => SmallVec::new(),
        };
        if let Some(&v) = self.score_cache.as_ref().unwrap().get(&(a, key.clone())) {
            self.cached += 1;
            return v;
        }
        let v = self.score_atom(a, loc);
        self.score_cache.as_mut().unwrap().insert((a, key), v);
        v
    }

    /// `scorePosition(states, index, offset)` for coarse states.
    fn score_position(&mut self, m: usize, index: usize, offset: u32) -> f64 {
        let mv = &self.ctx.movers[m];
        let mut ret = 0.0;
        for slot in 0..mv.n_moved {
            let d = in_state(&mv.coarse_del, index);
            if slot < d.len() && d[slot] {
                continue;
            }
            let a = mv.atoms[slot];
            ret += if self.score_cache.is_some() {
                self.score_atom_cached(a, index as u32 + offset)
            } else {
                self.score_atom(a, index as u32 + offset)
            };
        }
        ret
    }

    fn set(&mut self, m: usize, c: usize) {
        self.st.cfg[m] = Cfg::Coarse(c as u16);
    }

    fn all_states(num: &[usize]) -> Vec<Vec<usize>> {
        // generateAllStates: the first entry varies fastest
        let total: usize = num.iter().product();
        let mut out = Vec::with_capacity(total);
        let mut cur = vec![0usize; num.len()];
        for _ in 0..total {
            out.push(cur.clone());
            for k in 0..num.len() {
                cur[k] += 1;
                if cur[k] < num[k] {
                    break;
                }
                cur[k] = 0;
            }
        }
        out
    }

    fn brute_force(&mut self, movers: &[usize]) -> f64 {
        let num: Vec<usize> = movers.iter().map(|&m| self.ctx.movers[m].num_coarse()).collect();
        let mut best = -1e100;
        let mut best_state: Option<Vec<usize>> = None;
        let mag = self.mag();
        for cur in Self::all_states(&num) {
            for (k, &m) in movers.iter().enumerate() {
                if self.coarse_loc[m] != cur[k] {
                    self.set(m, cur[k]);
                    self.coarse_loc[m] = cur[k];
                }
            }
            let mut score = 0.0;
            for (k, &m) in movers.iter().enumerate() {
                score += mag * self.ctx.movers[m].coarse_pref[cur[k]];
                score += self.score_position(m, cur[k], 0);
            }
            if score > best || best_state.is_none() {
                best = score;
                best_state = Some(cur.clone());
            }
        }
        self.finish(movers, &best_state.unwrap())
    }

    fn finish(&mut self, movers: &[usize], best: &[usize]) -> f64 {
        for (k, &m) in movers.iter().enumerate() {
            self.set(m, best[k]);
            self.coarse_loc[m] = best[k];
        }
        let mag = self.mag();
        let mut ret = 0.0;
        for (k, &m) in movers.iter().enumerate() {
            let mut my = mag * self.ctx.movers[m].coarse_pref[best[k]];
            my += self.score_position(m, best[k], 0);
            self.high[m] = my;
            ret += my;
        }
        ret
    }

    fn components(verts: &[usize], edges: &[(usize, usize)]) -> Vec<Vec<usize>> {
        // connected components numbered by first vertex (boost DFS order)
        let mut comp_of: FxHashMap<usize, usize> = FxHashMap::default();
        let mut comps: Vec<Vec<usize>> = Vec::new();
        for &v in verts {
            if comp_of.contains_key(&v) {
                continue;
            }
            let id = comps.len();
            let mut stack = vec![v];
            comp_of.insert(v, id);
            let mut members = vec![];
            while let Some(x) = stack.pop() {
                members.push(x);
                for &(a, b) in edges {
                    let y = if a == x { b } else if b == x { a } else { continue };
                    if verts.contains(&y) && !comp_of.contains_key(&y) {
                        comp_of.insert(y, id);
                        stack.push(y);
                    }
                }
            }
            comps.push(members);
        }
        // members listed in vertex order
        for c in comps.iter_mut() {
            c.sort_by_key(|x| verts.iter().position(|y| y == x).unwrap());
        }
        comps
    }

    fn n_choose_m(n: usize, m: usize) -> Vec<Vec<usize>> {
        let mut result = Vec::new();
        let mut idx: Vec<usize> = (1..=m).collect();
        loop {
            result.push(idx.iter().map(|i| i - 1).collect());
            let mut i = m as isize - 1;
            while i >= 0 && idx[i as usize] == n - m + i as usize + 1 {
                i -= 1;
            }
            if i < 0 {
                break;
            }
            idx[i as usize] += 1;
            for j in (i as usize + 1)..m {
                idx[j] = idx[j - 1] + 1;
            }
        }
        result
    }

    fn vertex_cut(&mut self, movers: &[usize]) -> f64 {
        if movers.len() <= 2 {
            return self.brute_force(movers);
        }
        let edges: Vec<(usize, usize)> =
            self.edges.iter().copied().filter(|(a, b)| movers.contains(a) && movers.contains(b)).collect();
        // findVertexCut
        let n = movers.len();
        let mut cut: Vec<usize> = Vec::new();
        let mut keep: Vec<usize> = movers.to_vec();
        'search: for size in 0..n {
            for set in Self::n_choose_m(n, size) {
                let k: Vec<usize> = (0..n).filter(|i| !set.contains(i)).map(|i| movers[i]).collect();
                if Self::components(&k, &edges).len() > 1 {
                    cut = set.iter().map(|&i| movers[i]).collect();
                    keep = k;
                    break 'search;
                }
            }
        }
        if cut.is_empty() {
            return self.brute_force(movers);
        }
        let num: Vec<usize> = cut.iter().map(|&m| self.ctx.movers[m].num_coarse()).collect();
        let mut best = -1e100;
        let mut best_state: Option<Vec<usize>> = None;
        let mag = self.mag();
        let comps = Self::components(&keep, &edges);
        for cur in Self::all_states(&num) {
            for (k, &m) in cut.iter().enumerate() {
                if self.coarse_loc[m] != cur[k] {
                    self.set(m, cur[k]);
                    self.coarse_loc[m] = cur[k];
                }
            }
            let mut score = 0.0;
            for comp in &comps {
                score += self.vertex_cut(comp);
            }
            for (k, &m) in cut.iter().enumerate() {
                score += mag * self.ctx.movers[m].coarse_pref[cur[k]];
                score += self.score_position(m, cur[k], 0);
            }
            if score > best || best_state.is_none() {
                best = score;
                best_state = Some(movers.iter().map(|&m| self.coarse_loc[m]).collect());
            }
        }
        self.finish(movers, &best_state.unwrap())
    }
}

/// Run Reduce2's optimizer exactly (compat mode). Returns per-component results.
fn compat_optimize(ctx: &Ctx, components: &[Vec<u32>], edges: &[(usize, usize)], v: i32) -> Vec<CliqueResult> {
    let nm = ctx.movers.len();
    let all: Vec<u32> = (0..nm as u32).collect();
    let mut co = CompatOpt {
        ctx,
        st: CliqueState { movers: &all, cfg: vec![Cfg::Coarse(0); nm] },
        coarse_loc: vec![0; nm],
        high: vec![0.0; nm],
        trim: std::cell::RefCell::new(FxHashMap::default()),
        score_cache: None,
        buf: ScoreBuf::default(),
        calculated: 0,
        cached: 0,
        edges,
    };
    let mag = ctx.p.preference_magnitude;
    // Initialize
    let mut initial = vec![0.0; nm];
    for m in 0..nm {
        let mut s = mag * ctx.movers[m].coarse_pref[0];
        s += co.score_position(m, 0, 0);
        co.high[m] = s;
        initial[m] = s;
    }
    let mut coarse_best = vec![0.0; components.len()];
    // singletons first, then groups (the original's order)
    for (ci, comp) in components.iter().enumerate() {
        if comp.len() != 1 {
            continue;
        }
        let m = comp[0] as usize;
        let n = ctx.movers[m].num_coarse();
        let mut scores: Vec<f64> = ctx.movers[m].coarse_pref.iter().map(|p| p * mag).collect();
        for i in 0..n {
            co.set(m, i);
            scores[i] += co.score_position(m, i, 0);
        }
        let mut max = scores[0];
        let mut mi = 0;
        for i in 1..n {
            if scores[i] > max {
                max = scores[i];
                mi = i;
            }
        }
        co.set(m, mi);
        co.coarse_loc[m] = mi;
        co.high[m] = max;
        coarse_best[ci] = max;
    }
    for (ci, comp) in components.iter().enumerate() {
        if comp.len() == 1 {
            continue;
        }
        co.score_cache = Some(FxHashMap::default());
        let movers: Vec<usize> = comp.iter().map(|&m| m as usize).collect();
        coarse_best[ci] = co.vertex_cut(&movers);
        co.score_cache = None;
    }
    // Fine optimization in global Mover order
    let mut fine_info: Vec<(u32, String)> = Vec::new();
    for m in 0..nm {
        let mv = &ctx.movers[m];
        let c = co.coarse_loc[m];
        let initial_score = co.score_position(m, c, 0);
        let nf = mv.fine.count();
        let mut msg = String::new();
        if nf > 0 {
            let mut scores: Vec<f64> = mv.fine.prefs(c).iter().map(|p| p * mag).collect();
            for i in 0..nf {
                co.st.cfg[m] = Cfg::Fine(c as u16, i as u16);
                let mut s = 0.0;
                for slot in 0..mv.n_moved {
                    s += co.score_atom(mv.atoms[slot], 100000 + i as u32);
                }
                scores[i] += s;
            }
            let mut max = scores[0];
            let mut mi = 0;
            for i in 1..nf {
                if scores[i] > max {
                    max = scores[i];
                    mi = i;
                }
            }
            if max > co.high[m] {
                if v >= 3 {
                    msg += &format!(
                        "   Setting Mover to fine orientation {}, max score = {:.2} (coarse score {:.2})\n",
                        mi, max, initial_score
                    );
                }
                co.st.cfg[m] = Cfg::Fine(c as u16, mi as u16);
                co.high[m] = max;
            } else {
                co.st.cfg[m] = Cfg::Coarse(c as u16);
                if v >= 3 {
                    msg += "   Leaving Mover at coarse orientation\n";
                }
            }
        }
        fine_info.push((m as u32, msg));
    }
    // Flip annotations, in report order does not matter for state (reset after)
    let p = ctx.p;
    let mut annot = vec![String::new(); nm];
    for m in 0..nm {
        let mv = &ctx.movers[m];
        if !mv.kind.is_flip() {
            annot[m] = " .".into();
            continue;
        }
        let n = mv.num_coarse();
        let fin = co.st.cfg[m].coarse();
        let other = (fin + n / 2) % n;
        let saved = co.st.cfg[m];
        co.st.cfg[m] = Cfg::Coarse(other as u16);
        let (os, ob) = ctx.score_mover_atoms(&co.st, m as u32, &mut co.buf, false);
        let os = os + ctx.pref(m as u32, co.st.cfg[m]);
        co.st.cfg[m] = Cfg::Coarse(fin as u16);
        let (fs, fb) = ctx.score_mover_atoms(&co.st, m as u32, &mut co.buf, false);
        let fs = fs + ctx.pref(m as u32, co.st.cfg[m]);
        co.st.cfg[m] = saved;
        let desc = mv.pose_description(fin, None, !p.skip_bond_fixup);
        annot[m] = if ob && fb {
            " BothClash".into()
        } else if desc.contains("Unflipped") && os > fs && os - fs <= p.non_flip_preference {
            " Uncertain".into()
        } else {
            " .".into()
        };
    }
    let mut out = Vec::with_capacity(components.len());
    for (ci, comp) in components.iter().enumerate() {
        out.push(CliqueResult {
            initial: comp.iter().map(|&m| initial[m as usize]).collect(),
            cfg: comp.iter().map(|&m| co.st.cfg[m as usize]).collect(),
            high: comp.iter().map(|&m| co.high[m as usize]).collect(),
            coarse_best: coarse_best[ci],
            coarse_info: String::new(),
            fine_info: comp.iter().map(|&m| fine_info[m as usize].clone()).collect(),
            annot: comp.iter().map(|&m| annot[m as usize].clone()).collect(),
            calculated: if ci == 0 { co.calculated } else { 0 },
            cached: if ci == 0 { co.cached } else { 0 },
        });
    }
    out
}
