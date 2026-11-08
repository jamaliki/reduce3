//! Placement (`Optimizer._PlaceMovers`)

use super::*;

pub(super) struct Placement {
    pub(super) movers: Vec<Mover>,
    pub(super) mover_info: Vec<String>,
    pub(super) delete_atoms: Vec<u32>,
    pub(super) info: String,
    pub(super) amide_flips: Vec<FlippedMoverInfo>,
    pub(super) his_flips: Vec<FlippedMoverInfo>,
    pub(super) docked_alpha_carbons: Vec<u32>,
}

#[allow(clippy::too_many_arguments)]
pub(super) fn place_movers(
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
    clearance: Option<&movers::Clearance>,
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
                match movers::nh3_rotator(w, a, !p.compat, clearance) {
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
                        apply_fixup(w, &flip, 1, &mut pl.delete_atoms, &mut pl.docked_alpha_carbons);
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
                                apply_fixup(w, &flip, index, &mut pl.delete_atoms, &mut pl.docked_alpha_carbons);
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
        clearance,
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

fn apply_fixup(w: &mut World, mv: &Mover, c: usize, deletes: &mut Vec<u32>, docked: &mut Vec<u32>) {
    for (i, &pp) in in_state(&mv.fixup_pos, c).iter().enumerate() {
        w.pos[mv.atoms[i] as usize] = pp;
    }
    if !in_state(&mv.fixup_pos, c).is_empty() {
        docked.extend(mv.alpha_carbon);
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
        if !hist.fixup_pos[bc].is_empty() {
            pl.docked_alpha_carbons.extend(hist.alpha_carbon);
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
