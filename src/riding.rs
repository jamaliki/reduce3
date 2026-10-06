//! Riding hydrogen placement (port of mmtbx.hydrogens connectivity,
//! parameterization and `compute_h_position`).

use crate::geom::*;
use crate::interp::{angle_delta_deg, Interp};
use rustc_hash::{FxHashMap, FxHashSet};
use std::f64::consts::PI;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum HType {
    Flat2,
    TwoNeigbs,
    TwoTetra,
    ThreeNeigbs,
    Alg1a,
    Alg1b,
    Prop,
}

impl HType {
    pub fn name(self) -> &'static str {
        match self {
            HType::Flat2 => "flat_2neigbs",
            HType::TwoNeigbs => "2neigbs",
            HType::TwoTetra => "2tetra",
            HType::ThreeNeigbs => "3neigbs",
            HType::Alg1a => "alg1a",
            HType::Alg1b => "alg1b",
            HType::Prop => "prop",
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct RidingCoef {
    pub htype: HType,
    pub ih: u32,
    pub a0: u32,
    pub a1: u32,
    pub a2: i64,
    pub a3: i64,
    pub a: f64,
    pub b: f64,
    pub h: f64,
    pub n: i32,
    pub disth: f64,
}

#[derive(Clone, Copy, Debug, Default)]
struct NbAngle {
    iseq: u32,
    angle: f64,
}

#[derive(Clone, Debug, Default)]
struct Nb {
    ih: u32,
    a0: Option<(u32, f64)>,
    angle_a1a0a2: Option<f64>,
    angle_a2a0a3: Option<f64>,
    angle_a3a0a1: Option<f64>,
    a: [Option<NbAngle>; 3],
    h: [Option<NbAngle>; 2],
    b1_iseq: Option<u32>,
    b1_ideal: Option<f64>,
    nh: usize,
    nnonh: usize,
    is_in_plane: bool,
    valid: bool,
}

/// Inputs describing the atoms.
pub struct RidingAtoms<'a> {
    pub is_h: &'a [bool],
    pub altloc: &'a [String],
    pub name: &'a [String],
    pub occ: &'a [f64],
    /// Atoms of the same residue group (for alternate-neighbor reduction).
    pub rg_of: &'a [u64],
    /// Expected heavy-neighbor count of an atom from its dictionary (None if unknown).
    pub expected_heavy: &'a dyn Fn(u32) -> Option<usize>,
    /// Take the torsion of an NH2 hydrogen from the dictionary instead of the
    /// periodic image nearest its current position. Reduce2 uses the image,
    /// and the temporary position then decides which H gets which name
    /// (HD21/HD22 come out swapped for a few percent of amides).
    pub dictionary_nh2_torsion: bool,
    /// Orient a one-neighbor group by another atom when its third neighbor lies
    /// on the parent bond axis (a ligand on a crystallographic two-fold, say).
    /// Reduce2 divides by zero there.
    pub reroute_axial_reference: bool,
}

pub struct RidingResult {
    pub coef: Vec<Option<RidingCoef>>,
    /// H that could not be parameterized.
    pub unparameterized: Vec<u32>,
    /// H whose third neighbor lies on the parent bond axis, where Reduce2
    /// divides by zero (left empty when the reference is rerouted).
    pub axial_reference: Vec<u32>,
    pub warnings: Vec<String>,
}

fn reduce_alternates(at: &RidingAtoms, cands: &[u32]) -> Vec<u32> {
    // `process_alternate_neighbors`
    if cands.len() == 1 {
        return vec![cands[0]];
    }
    let mut used: Vec<u32> = Vec::new();
    let mut reduced: Vec<u32> = Vec::new();
    for &i in cands {
        if used.contains(&i) {
            continue;
        }
        let alt_i = &at.altloc[i as usize];
        if let Some(&prev) = reduced.first() {
            if at.altloc[prev as usize] != *alt_i {
                used.push(i);
                continue;
            }
        }
        let mut best: Vec<(u32, f64)> = vec![(i, at.occ[i as usize])];
        for &j in cands {
            if j == i || used.contains(&j) {
                continue;
            }
            if at.name[j as usize] == at.name[i as usize] && at.rg_of[j as usize] == at.rg_of[i as usize] {
                best.push((j, at.occ[j as usize]));
            }
        }
        // max by occupancy; first maximum in insertion order
        let mut bi = best[0];
        for &x in &best[1..] {
            if x.1 > bi.1 {
                bi = x;
            }
        }
        reduced.push(bi.0);
        for (k, _) in best {
            used.push(k);
        }
    }
    reduced
}

/// Build the riding parameterization and new H positions.
pub fn riding(it: &Interp, sites: &mut [Vec3], at: &RidingAtoms, use_ideal_dihedral: bool) -> RidingResult {
    let n = sites.len();
    let mut warnings = Vec::new();
    // bond order: sorted by (i, j)
    let mut bonds: Vec<(u32, u32, f64)> = it.bonds.iter().map(|b| (b.i, b.j, b.ideal)).collect();
    bonds.sort_by(|x, y| (x.0, x.1).cmp(&(y.0, y.1)));
    let mut fsc0: Vec<Vec<u32>> = vec![Vec::new(); n];
    for &(i, j, _) in &bonds {
        fsc0[i as usize].push(j);
        fsc0[j as usize].push(i);
    }
    let mut conn: Vec<Option<Nb>> = vec![None; n];
    let mut double_h: FxHashSet<u32> = FxHashSet::default();
    let mut parents: FxHashSet<u32> = FxHashSet::default();
    // 1. first neighbors
    for &(i, j, ideal) in &bonds {
        let (hi, hj) = (at.is_h[i as usize], at.is_h[j as usize]);
        if !hi && !hj {
            continue;
        }
        if hi && hj {
            continue;
        }
        let (ih, parent) = if hi { (i, j) } else { (j, i) };
        if conn[ih as usize].is_some() {
            continue;
        }
        if fsc0[ih as usize].len() > 1 {
            double_h.insert(ih);
        }
        conn[ih as usize] = Some(Nb { ih, a0: Some((parent, ideal)), valid: true, ..Default::default() });
        parents.insert(parent);
    }
    // 2. slipped H (no bond)
    let slipped: Vec<u32> = (0..n as u32).filter(|&k| at.is_h[k as usize] && conn[k as usize].is_none()).collect();
    // 3. second neighbors (angle proxy order)
    let mut second_raw: Vec<Vec<u32>> = vec![Vec::new(); n];
    let mut angle_dict: FxHashMap<(u32, u32, u32), f64> = FxHashMap::default();
    for ap in &it.angles {
        for &i_test in &ap.i {
            if conn[i_test as usize].is_none() && !at.is_h[i_test as usize] {
                continue;
            }
            let ih = i_test;
            let parent = ap.i[1];
            if !parents.contains(&parent) {
                continue;
            }
            let Some(nb) = conn[ih as usize].as_ref() else { continue };
            let second = match ap.i.iter().find(|&&x| x != ih && x != parent) {
                Some(&s) => s,
                None => continue,
            };
            if nb.a0.unwrap().0 != parent {
                if double_h.contains(&ih) {
                    continue;
                }
                warnings.push(format!("angle and bond restraints conflict for H {}", ih));
                continue;
            }
            second_raw[ih as usize].push(second);
            angle_dict.insert((ih, parent, second), ap.ideal);
        }
    }
    // 4. plane proxies: first H of each plane is "in plane"
    let mut plane_h: FxHashSet<u32> = FxHashSet::default();
    for pl in &it.planes {
        if let Some(&first_h) = pl.iter().find(|&&a| at.is_h[a as usize]) {
            plane_h.insert(first_h);
        }
    }
    // 5. process second neighbors
    let mut a1_atoms: FxHashSet<u32> = FxHashSet::default();
    let mut a0a1: FxHashMap<u32, Vec<u32>> = FxHashMap::default();
    for ih in 0..n {
        let Some(nb) = conn[ih].as_mut() else { continue };
        let parent = nb.a0.unwrap().0;
        let mut reduced: Vec<u32> = Vec::new();
        let mut alts: Vec<u32> = Vec::new();
        for &s in &second_raw[ih] {
            if at.altloc[s as usize].is_empty() {
                reduced.push(s);
            } else {
                alts.push(s);
            }
        }
        if !alts.is_empty() {
            reduced.extend(reduce_alternates(at, &alts));
        }
        let hs: Vec<u32> = reduced.iter().copied().filter(|&x| at.is_h[x as usize]).collect();
        // non-H: Python set order (iteration by hash); positions do not depend on it
        let mut non_h: Vec<u32> = Vec::new();
        for &x in &reduced {
            if !at.is_h[x as usize] && !non_h.contains(&x) {
                non_h.push(x);
            }
        }
        let mut cnt = 0;
        for (k, &x) in non_h.iter().take(3).enumerate() {
            nb.a[k] = Some(NbAngle { iseq: x, angle: angle_dict[&(ih as u32, parent, x)] });
            cnt += 1;
        }
        let mut hc = 0;
        for (k, &x) in hs.iter().take(2).enumerate() {
            nb.h[k] = Some(NbAngle { iseq: x, angle: angle_dict[&(ih as u32, parent, x)] });
            hc += 1;
        }
        nb.nh = hc;
        nb.nnonh = cnt;
        if cnt == 1 {
            let a1 = nb.a[0].unwrap().iseq;
            a1_atoms.insert(a1);
            a0a1.entry(parent).or_default().push(a1);
        }
    }
    // 6. third neighbors from dihedral proxies (const first, then var; last write wins)
    let assign_group = |conn: &mut Vec<Option<Nb>>, ih: u32, third: u32| {
        let (nh, h1, h2) = {
            let nb = conn[ih as usize].as_ref().unwrap();
            (nb.nh, nb.h[0].map(|x| x.iseq), nb.h[1].map(|x| x.iseq))
        };
        if nh > 0 {
            if let Some(h1) = h1 {
                if let Some(o) = conn[h1 as usize].as_mut() {
                    o.b1_iseq = Some(third);
                    o.b1_ideal = None;
                }
            }
            if nh == 2 {
                if let Some(h2) = h2 {
                    if let Some(o) = conn[h2 as usize].as_mut() {
                        o.b1_iseq = Some(third);
                        o.b1_ideal = None;
                    }
                }
            }
        }
    };
    for dp in &it.const_dihedrals {
        for &t in &dp.i {
            if conn[t as usize].is_none() {
                continue;
            }
            let (third, ideal) = if t == dp.i[0] {
                (dp.i[3], dp.ideal)
            } else if t == dp.i[3] {
                (dp.i[0], -dp.ideal)
            } else {
                continue;
            };
            {
                let nb = conn[t as usize].as_mut().unwrap();
                nb.b1_iseq = Some(third);
                nb.b1_ideal = Some(ideal);
            }
            assign_group(&mut conn, t, third);
        }
    }
    'dp: for dp in &it.dihedrals {
        for &t in &dp.i {
            if conn[t as usize].is_none() {
                continue;
            }
            let third = if t == dp.i[0] {
                dp.i[3]
            } else if t == dp.i[3] {
                dp.i[0]
            } else {
                continue;
            };
            let Some(model) = dihedral_rad_cpp(sites[dp.i[0] as usize], sites[dp.i[1] as usize], sites[dp.i[2] as usize], sites[dp.i[3] as usize])
            else {
                continue 'dp;
            };
            let model_deg = model.to_degrees();
            let delta = angle_delta_deg(model_deg, dp.ideal, dp.period);
            {
                let nb = conn[t as usize].as_mut().unwrap();
                nb.b1_iseq = Some(third);
                let nh2 = at.dictionary_nh2_torsion && dp.period == 2 && nb.nh == 1;
                // alg1a places the H at its torsion parameter + 180 degrees
                nb.b1_ideal = Some(if nh2 { dp.ideal - 180.0 } else { model_deg + delta });
            }
            assign_group(&mut conn, t, third);
        }
    }
    // 7. angles at the parent and fallback third neighbors
    let mut parent_angles: FxHashMap<u32, FxHashMap<(u32, u32), f64>> = FxHashMap::default();
    let mut third_raw: FxHashMap<u32, Vec<u32>> = FxHashMap::default();
    for ap in &it.angles {
        let [ix, iy, iz] = ap.i;
        let (hx, hz) = (at.is_h[ix as usize], at.is_h[iz as usize]);
        if parents.contains(&iy) && !hx && !hz {
            parent_angles.entry(iy).or_default().insert((ix, iz), ap.ideal);
        }
        if hx || hz {
            continue;
        }
        if !a1_atoms.contains(&iy) {
            continue;
        }
        let (parent, third) = if parents.contains(&ix) && a0a1.contains_key(&ix) {
            if a0a1[&ix].contains(&iy) {
                (ix, iz)
            } else {
                continue;
            }
        } else if parents.contains(&iz) && a0a1.contains_key(&iz) {
            if a0a1[&iz].contains(&iy) {
                (iz, ix)
            } else {
                continue;
            }
        } else {
            continue;
        };
        let e = third_raw.entry(parent).or_default();
        if !e.contains(&third) {
            e.push(third);
            if parents.contains(&third) {
                third_raw.entry(third).or_default().push(parent);
            }
        }
    }
    let get_pa = |a0: u32, x: u32, z: u32| -> Option<f64> {
        let m = parent_angles.get(&a0)?;
        m.get(&(x, z)).or_else(|| m.get(&(z, x))).copied()
    };
    for ih in 0..n {
        let Some(nb0) = conn[ih].clone() else { continue };
        // assign_a0_angles
        if nb0.nnonh > 1 {
            let a0 = nb0.a0.unwrap().0;
            let a1 = nb0.a[0].unwrap().iseq;
            let a2 = nb0.a[1].unwrap().iseq;
            let ang = get_pa(a0, a1, a2);
            let nb = conn[ih].as_mut().unwrap();
            nb.angle_a1a0a2 = ang;
            if ang.is_none() {
                *nb = Nb { ih: ih as u32, valid: false, ..Default::default() };
                continue;
            }
            if nb0.nnonh == 3 {
                let a3 = nb0.a[2].unwrap().iseq;
                nb.angle_a2a0a3 = get_pa(a0, a2, a3);
                nb.angle_a3a0a1 = get_pa(a0, a3, a1);
                if nb.angle_a2a0a3.is_none() || nb.angle_a3a0a1.is_none() {
                    *nb = Nb { ih: ih as u32, valid: false, ..Default::default() };
                    continue;
                }
            }
        }
        if nb0.nnonh != 1 || nb0.b1_iseq.is_some() {
            continue;
        }
        let parent = nb0.a0.unwrap().0;
        let raw = third_raw.get(&parent).cloned().unwrap_or_default();
        let mut reduced: Vec<u32> = Vec::new();
        let mut alts: Vec<u32> = Vec::new();
        for t in raw {
            if at.altloc[t as usize].is_empty() {
                reduced.push(t);
            } else {
                alts.push(t);
            }
        }
        if !alts.is_empty() {
            reduced.extend(reduce_alternates(at, &alts));
        }
        if reduced.is_empty() {
            conn[ih] = Some(Nb { ih: ih as u32, valid: false, ..Default::default() });
            continue;
        }
        if parent_lost_heavy_neighbor(at, &fsc0, parent) {
            conn[ih] = Some(Nb { ih: ih as u32, valid: false, ..Default::default() });
            continue;
        }
        let third = reduced[0];
        let nb = conn[ih].as_mut().unwrap();
        if nb0.nh == 2 {
            nb.b1_iseq = Some(third);
            let sibling = nb0.h[0].map(|x| x.iseq);
            nb.b1_ideal = Some(if sibling.map(|s| (ih as u32) < s).unwrap_or(false) { 60.0 } else { -60.0 });
        } else {
            nb.b1_iseq = Some(third);
            nb.b1_ideal = Some(180.0);
        }
    }
    for &s in &slipped {
        conn[s as usize] = Some(Nb { ih: s, valid: false, ..Default::default() });
    }
    for ih in 0..n {
        if let Some(nb) = conn[ih].as_mut() {
            nb.is_in_plane = plane_h.contains(&(ih as u32));
        }
    }

    if at.reroute_axial_reference {
        reroute_axial_references(at, &fsc0, sites, &mut conn);
    }

    // ---------------- parameterization
    let mut axial_reference: Vec<u32> = Vec::new();
    let mut coef: Vec<Option<RidingCoef>> = vec![None; n];
    let mut unk: Vec<u32> = Vec::new();
    for ih in 0..n {
        let Some(nb) = conn[ih].clone() else { continue };
        if coef[ih].is_some() {
            continue;
        }
        if !nb.valid || nb.a0.is_none() {
            unk.push(ih as u32);
            continue;
        }
        let r = if nb.nnonh == 2 {
            process_2(&nb, &conn, sites, &mut coef)
        } else if nb.nnonh == 3 && nb.nh == 0 {
            process_3(&nb, sites, &mut coef)
        } else if nb.nnonh == 1 && (nb.nh == 0 || nb.nh == 2) {
            process_1(&nb, &conn, sites, &mut coef, use_ideal_dihedral, &mut axial_reference)
        } else if nb.nnonh == 1 && nb.nh == 1 {
            process_1_arg(&nb, &conn, sites, &mut coef, &mut axial_reference)
        } else {
            Err(())
        };
        if r.is_err() {
            unk.push(ih as u32);
        }
    }
    // ---------------- positions (i_seq order, written immediately)
    for ih in 0..n {
        if let Some(c) = coef[ih] {
            if let Some(p) = compute_h_position(&c, sites) {
                sites[ih] = p;
            }
        }
    }
    let unparameterized: Vec<u32> = (0..n as u32).filter(|&k| at.is_h[k as usize] && coef[k as usize].is_none()).collect();
    let _ = unk;
    RidingResult { coef, unparameterized, axial_reference, warnings }
}

fn parent_lost_heavy_neighbor(at: &RidingAtoms, fsc0: &[Vec<u32>], parent: u32) -> bool {
    let Some(expected) = (at.expected_heavy)(parent) else { return false };
    let palt = &at.altloc[parent as usize];
    let mut seen: FxHashSet<&str> = FxHashSet::default();
    for &j in &fsc0[parent as usize] {
        if at.is_h[j as usize] {
            continue;
        }
        let jalt = &at.altloc[j as usize];
        if !palt.is_empty() && !jalt.is_empty() && palt != jalt {
            continue;
        }
        seen.insert(at.name[j as usize].trim());
    }
    seen.len() < expected
}

fn superposed(a: Vec3, b: Vec3) -> bool {
    (a - b).length() < 0.001
}

/// Distance of `rb1` from the line through `r0` and `r1`, in Python's
/// arithmetic order (`process_1_neighbor`); Reduce2 normalizes this vector.
fn off_axis(r0: Vec3, r1: Vec3, rb1: Vec3) -> f64 {
    let u1 = (r0 - r1).normalize();
    let rb10 = rb1 - r1;
    (rb10 - u1 * rb10.dot(u1)).length()
}

/// A third neighbor closer than this to the parent bond axis does not orient
/// the group: below the coordinate precision of a model file, its direction is
/// rounding noise.
const AXIS_TOLERANCE: f64 = 0.001;

/// Give each one-neighbor group whose third neighbor lies on the parent bond
/// axis another reference off the axis: a further heavy neighbor of the first
/// neighbor when there is one (so the movers can name it too), else the
/// nearest heavy atom off the axis. The ideal torsion is kept; it is undefined
/// about the axial atom, and the groups this touches are rotatable.
fn reroute_axial_references(at: &RidingAtoms, fsc0: &[Vec<u32>], sites: &[Vec3], conn: &mut [Option<Nb>]) {
    for ih in 0..conn.len() {
        let Some(nb) = conn[ih].as_ref() else { continue };
        if !nb.valid || nb.nnonh != 1 {
            continue;
        }
        let (Some((a0, _)), Some(a1), Some(b1)) = (nb.a0, nb.a[0], nb.b1_iseq) else { continue };
        let a1 = a1.iseq;
        let (r0, r1) = (sites[a0 as usize], sites[a1 as usize]);
        if superposed(r0, r1) || off_axis(r0, r1, sites[b1 as usize]) >= AXIS_TOLERANCE {
            continue;
        }
        let altloc = &at.altloc[a0 as usize];
        let usable = |c: u32| {
            let calt = &at.altloc[c as usize];
            c != a0
                && c != a1
                && !at.is_h[c as usize]
                && (altloc.is_empty() || calt.is_empty() || calt == altloc)
                && off_axis(r0, r1, sites[c as usize]) >= AXIS_TOLERANCE
        };
        let bonded = fsc0[a1 as usize].iter().copied().find(|&c| usable(c));
        let nearest = || {
            (0..sites.len() as u32)
                .filter(|&c| usable(c))
                .min_by(|&x, &y| (sites[x as usize] - r1).length().total_cmp(&(sites[y as usize] - r1).length()))
        };
        conn[ih].as_mut().unwrap().b1_iseq = bonded.or_else(nearest);
    }
}

fn process_1(
    nb: &Nb,
    conn: &[Option<Nb>],
    sites: &[Vec3],
    coef: &mut [Option<RidingCoef>],
    use_ideal_dihedral: bool,
    axial: &mut Vec<u32>,
) -> Result<(), ()> {
    let mut nbs = nb.clone();
    if nb.nh == 2 {
        let (h1, h2) = (nb.h[0].unwrap().iseq, nb.h[1].unwrap().iseq);
        if nb.b1_ideal.is_some() {
            // use self
        } else if conn[h1 as usize].as_ref().map(|x| x.b1_ideal.is_some()).unwrap_or(false) {
            if coef[h1 as usize].is_none() {
                nbs = conn[h1 as usize].clone().unwrap();
            }
        } else if conn[h2 as usize].as_ref().map(|x| x.b1_ideal.is_some()).unwrap_or(false) && coef[h2 as usize].is_none() {
            nbs = conn[h2 as usize].clone().unwrap();
        }
    }
    let ih = nbs.ih;
    let (i_a0, disth) = nbs.a0.ok_or(())?;
    let a1 = nbs.a[0].ok_or(())?;
    let i_b1 = nbs.b1_iseq.ok_or(())?;
    let (rh, r0, r1) = (sites[ih as usize], sites[i_a0 as usize], sites[a1.iseq as usize]);
    if superposed(rh, r0) || superposed(r1, r0) {
        return Err(());
    }
    if off_axis(r0, r1, sites[i_b1 as usize]) == 0.0 {
        axial.push(ih);
        return Err(());
    }
    let dihedral = dihedral_rad_cpp(sites[ih as usize], sites[i_a0 as usize], sites[a1.iseq as usize], sites[i_b1 as usize]);
    let alpha = a1.angle.to_radians();
    let mut phi = dihedral;
    if use_ideal_dihedral {
        if let Some(d) = nbs.b1_ideal {
            phi = Some(d.to_radians());
        }
    }
    let phi = phi.ok_or(())?;
    if nbs.nh == 0 {
        coef[ih as usize] = Some(RidingCoef {
            htype: HType::Alg1b,
            ih,
            a0: i_a0,
            a1: a1.iseq,
            a2: i_b1 as i64,
            a3: -1,
            a: alpha,
            b: phi,
            h: 0.0,
            n: 0,
            disth,
        });
    } else {
        let (mut h1, mut h2) = (nbs.h[0].unwrap().iseq, nbs.h[1].unwrap().iseq);
        // check_propeller_order
        let rh2 = sites[h2 as usize];
        if !(((rh - r0).cross(rh2 - r0)).dot(r1 - r0) >= 0.0) {
            std::mem::swap(&mut h1, &mut h2);
        }
        for (n, hp) in [(0, ih), (1, h1), (2, h2)] {
            coef[hp as usize] = Some(RidingCoef {
                htype: HType::Prop,
                ih: hp,
                a0: i_a0,
                a1: a1.iseq,
                a2: i_b1 as i64,
                a3: -1,
                a: alpha,
                b: phi,
                h: 0.0,
                n,
                disth,
            });
        }
    }
    Ok(())
}

fn process_1_arg(
    nb: &Nb,
    conn: &[Option<Nb>],
    sites: &[Vec3],
    coef: &mut [Option<RidingCoef>],
    axial: &mut Vec<u32>,
) -> Result<(), ()> {
    let ih = nb.ih;
    let h1 = nb.h[0].unwrap();
    let (i_a0, disth) = nb.a0.ok_or(())?;
    let a1 = nb.a[0].ok_or(())?;
    let (ih_d, ih_nd) = if nb.b1_ideal.is_some() {
        (ih, h1.iseq)
    } else if conn[h1.iseq as usize].as_ref().map(|x| x.b1_ideal.is_some()).unwrap_or(false) {
        (h1.iseq, ih)
    } else {
        return Err(());
    };
    let nbd = conn[ih_d as usize].as_ref().unwrap();
    let i_b1 = nbd.b1_iseq.ok_or(())?;
    if h1.angle > 107.0 && h1.angle < 111.0 {
        return Err(());
    }
    let (rh, r0, r1) = (sites[ih as usize], sites[i_a0 as usize], sites[a1.iseq as usize]);
    if superposed(rh, r0) || superposed(r1, r0) {
        return Err(());
    }
    if off_axis(r0, r1, sites[i_b1 as usize]) == 0.0 {
        axial.push(ih);
        return Err(());
    }
    let alpha = a1.angle.to_radians();
    let phi = nbd.b1_ideal.unwrap().to_radians();
    for (h, ph) in [(ih_d, phi), (ih_nd, phi + PI)] {
        if coef[h as usize].is_none() {
            coef[h as usize] = Some(RidingCoef {
                htype: HType::Alg1a,
                ih: h,
                a0: i_a0,
                a1: a1.iseq,
                a2: i_b1 as i64,
                a3: -1,
                a: alpha,
                b: ph,
                h: 0.0,
                n: 0,
                disth,
            });
        }
    }
    Ok(())
}

fn process_2(nb: &Nb, conn: &[Option<Nb>], sites: &[Vec3], coef: &mut [Option<RidingCoef>]) -> Result<(), ()> {
    let _ = conn;
    let ih = nb.ih;
    let (i_a0, disth) = nb.a0.ok_or(())?;
    let a1 = nb.a[0].ok_or(())?;
    let a2 = nb.a[1].ok_or(())?;
    let (rh, r0, r1, r2) = (sites[ih as usize], sites[i_a0 as usize], sites[a1.iseq as usize], sites[a2.iseq as usize]);
    if superposed(rh, r0) || superposed(r1, r0) || superposed(r2, r0) {
        return Err(());
    }
    let uh0 = (rh - r0).normalize();
    let u10 = (r1 - r0).normalize();
    let u20 = (r2 - r0).normalize();
    let alpha0 = nb.angle_a1a0a2.ok_or(())?.to_radians();
    let alpha1 = a1.angle.to_radians();
    let alpha2 = a2.angle.to_radians();
    let (c0, c1, c2) = (alpha0.cos(), alpha1.cos(), alpha2.cos());
    let sumang = alpha0 + alpha1 + alpha2;
    let denom = 1.0 - c0 * c0;
    if denom == 0.0 {
        return Err(());
    }
    let a = (c1 - c0 * c2) / (1.0 - c0 * c0);
    let b = (c2 - c0 * c1) / (1.0 - c0 * c0);
    let flat = sumang < 2.0 * PI + 0.05 && sumang > 2.0 * PI - 0.05;
    let root = 1.0 - c1 * c1 - c2 * c2 - c0 * c0 + 2.0 * c0 * c1 * c2;
    let mut h = 0.0;
    let htype;
    if sumang > 2.0 * PI + 0.05 && root < 0.0 {
        return Err(());
    } else if flat {
        htype = HType::Flat2;
    } else if nb.nh == 1 {
        let h1 = nb.h[0].unwrap();
        let rh2 = sites[h1.iseq as usize];
        if superposed(rh2, r0) {
            return Err(());
        }
        h = h1.angle.to_radians() * 0.5;
        if u10.cross(u20).dot(uh0) < 0.0 {
            h = -h;
        }
        htype = HType::TwoTetra;
        coef[h1.iseq as usize] = Some(RidingCoef {
            htype: HType::TwoTetra,
            ih: h1.iseq,
            a0: i_a0,
            a1: a1.iseq,
            a2: a2.iseq as i64,
            a3: -1,
            a,
            b,
            h: -h,
            n: 0,
            disth,
        });
    } else {
        if root < 0.0 {
            return Err(());
        }
        let sa0 = alpha0.sin();
        if sa0 == 0.0 {
            return Err(());
        }
        h = root.sqrt() / sa0;
        if u10.cross(u20).dot(uh0) < 0.0 {
            h = -h;
        }
        htype = if nb.is_in_plane { HType::Flat2 } else { HType::TwoNeigbs };
    }
    coef[ih as usize] = Some(RidingCoef { htype, ih, a0: i_a0, a1: a1.iseq, a2: a2.iseq as i64, a3: -1, a, b, h, n: 0, disth });
    Ok(())
}

fn process_3(nb: &Nb, sites: &[Vec3], coef: &mut [Option<RidingCoef>]) -> Result<(), ()> {
    let ih = nb.ih;
    let (i_a0, disth) = nb.a0.ok_or(())?;
    let (a1, a2, a3) = (nb.a[0].ok_or(())?, nb.a[1].ok_or(())?, nb.a[2].ok_or(())?);
    let r0 = sites[i_a0 as usize];
    for x in [ih, a1.iseq, a2.iseq, a3.iseq] {
        if superposed(sites[x as usize], r0) {
            return Err(());
        }
    }
    let c1 = a1.angle.to_radians().cos();
    let c2 = a2.angle.to_radians().cos();
    let c3 = a3.angle.to_radians().cos();
    let w12 = nb.angle_a1a0a2.ok_or(())?.to_radians().cos();
    let w23 = nb.angle_a2a0a3.ok_or(())?.to_radians().cos();
    let w13 = nb.angle_a3a0a1.ok_or(())?.to_radians().cos();
    let d = det3([[1.0, w12, w13], [w12, 1.0, w23], [w13, w23, 1.0]]);
    if d == 0.0 {
        return Err(());
    }
    let dx = det3([[c1, w12, w13], [c2, 1.0, w23], [c3, w23, 1.0]]);
    let dy = det3([[1.0, c1, w13], [w12, c2, w23], [w13, c3, 1.0]]);
    let dz = det3([[1.0, w12, c1], [w12, 1.0, c2], [w13, w23, c3]]);
    coef[ih as usize] = Some(RidingCoef {
        htype: HType::ThreeNeigbs,
        ih,
        a0: i_a0,
        a1: a1.iseq,
        a2: a2.iseq as i64,
        a3: a3.iseq as i64,
        a: dx / d,
        b: dy / d,
        h: dz / d,
        n: 0,
        disth,
    });
    Ok(())
}

/// scitbx vec3 arithmetic as the arm64 cctbx build evaluates it: clang fuses
/// `x += a*b` and `a*b - c*d` written in one expression (dot products,
/// lengths, cross products); element-wise vector operators stay unfused.
pub mod cx {
    use crate::geom::{v3, Vec3};
    #[inline]
    pub fn fmadd(a: f64, b: f64, c: f64) -> f64 {
        if cfg!(target_arch = "aarch64") { a.mul_add(b, c) } else { a * b + c }
    }
    #[inline]
    pub fn dot(a: Vec3, b: Vec3) -> f64 {
        fmadd(a.z, b.z, fmadd(a.y, b.y, a.x * b.x))
    }
    #[inline]
    pub fn length(a: Vec3) -> f64 {
        dot(a, a).sqrt()
    }
    #[inline]
    pub fn normalize(a: Vec3) -> Vec3 {
        a / length(a)
    }
    #[inline]
    pub fn cross(a: Vec3, b: Vec3) -> Vec3 {
        v3(fmadd(a.y, b.z, -(b.y * a.z)), fmadd(a.z, b.x, -(b.z * a.x)), fmadd(a.x, b.y, -(b.x * a.y)))
    }
}

/// `compute_h_position` (mmtbx/hydrogens/hydrogens.h).
pub fn compute_h_position(c: &RidingCoef, sites: &[Vec3]) -> Option<Vec3> {
    let r0 = sites[c.a0 as usize];
    let r1 = sites[c.a1 as usize];
    let unit = cx::normalize;
    let (a, b, h, dh) = (c.a, c.b, c.h, c.disth);
    Some(match c.htype {
        HType::Flat2 => {
            let r2 = sites[c.a2 as usize];
            let rh0 = unit(r1 - r0) * a + unit(r2 - r0) * b;
            let l = cx::length(rh0);
            if l <= 0.0 {
                return None;
            }
            r0 + (rh0 / l) * dh
        }
        HType::TwoNeigbs => {
            let r2 = sites[c.a2 as usize];
            let u10 = unit(r1 - r0);
            let u20 = unit(r2 - r0);
            let v0 = unit(cx::cross(u10, u20));
            let rh0 = u10 * a + u20 * b + v0 * h;
            let l = cx::length(rh0);
            if l <= 0.0 {
                return None;
            }
            r0 + (rh0 / l) * dh
        }
        HType::TwoTetra => {
            let r2 = sites[c.a2 as usize];
            let u10 = unit(r1 - r0);
            let u20 = unit(r2 - r0);
            let v0 = unit(cx::cross(u10, u20));
            let d0 = unit(u10 * a + u20 * b);
            r0 + (d0 * h.cos() + v0 * h.sin()) * dh
        }
        HType::ThreeNeigbs => {
            let r2 = sites[c.a2 as usize];
            let r3 = sites[c.a3 as usize];
            let rh0 = unit(r1 - r0) * a + unit(r2 - r0) * b + unit(r3 - r0) * h;
            let l = cx::length(rh0);
            if l <= 0.0 {
                return None;
            }
            r0 + (rh0 / l) * dh
        }
        HType::Alg1a | HType::Alg1b | HType::Prop => {
            let rb1 = sites[c.a2 as usize];
            let phi = b + c.n as f64 * 2.0 * PI / 3.0;
            let (salpha, calpha) = (a.sin(), a.cos());
            let (sphi, cphi) = (phi.sin(), phi.cos());
            let u1 = unit(r0 - r1);
            let rb10 = rb1 - r1;
            let u2 = unit(rb10 - u1 * cx::dot(rb10, u1));
            let u3 = cx::cross(u1, u2);
            r0 + ((u2 * cphi + u3 * sphi) * salpha - u1 * calpha) * dh
        }
    })
}
