//! Riding hydrogen placement (port of mmtbx.hydrogens connectivity,
//! parameterization and `compute_h_position`).

use crate::geom::*;
use crate::interp::{angle_delta_deg, Adjacency, Interp};
use crate::model::ShortStr;
use rustc_hash::{FxHashMap, FxHashSet};
use smallvec::SmallVec;
use std::f64::consts::PI;

/// Neighbor list of one atom, inline for up to four.
type Neighbors = SmallVec<[u32; 4]>;

#[cfg(test)]
mod plane_tests;
mod planar;

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
    explicit_plane: bool,
    valid: bool,
}

/// The neighbor records of the atoms, found by atom index. Only hydrogens
/// get records, so each has a slot, in atom order, made when it is created.
struct Conn {
    slot: Vec<u32>,
    /// Whether each slot holds a record, so `index` need not read the
    /// large records.
    live: Vec<bool>,
    nbs: Vec<Option<Nb>>,
}

impl Conn {
    const NONE: u32 = u32::MAX;

    fn new(is_h: &[bool]) -> Conn {
        let mut next = 0u32;
        let slot = is_h.iter().map(|&h| if h { next += 1; next - 1 } else { Conn::NONE }).collect();
        let nbs = std::iter::repeat_with(|| None).take(next as usize).collect();
        Conn { slot, live: vec![false; next as usize], nbs }
    }
    /// The slot of atom `a`, if it has a record.
    #[inline]
    fn index(&self, a: u32) -> Option<usize> {
        let k = self.slot[a as usize];
        (k != Conn::NONE && self.live[k as usize]).then_some(k as usize)
    }
    #[inline]
    fn get(&self, a: u32) -> Option<&Nb> {
        let k = self.slot[a as usize];
        if k == Conn::NONE { None } else { self.nbs[k as usize].as_ref() }
    }
    #[inline]
    fn get_mut(&mut self, a: u32) -> Option<&mut Nb> {
        let k = self.slot[a as usize];
        if k == Conn::NONE { None } else { self.nbs[k as usize].as_mut() }
    }
    /// Create or replace the record of atom `a`; returns its slot.
    fn set(&mut self, a: u32, nb: Nb) -> usize {
        let mut k = self.slot[a as usize];
        if k == Conn::NONE {
            k = self.nbs.len() as u32;
            self.slot[a as usize] = k;
            self.nbs.push(None);
            self.live.push(false);
        }
        self.nbs[k as usize] = Some(nb);
        self.live[k as usize] = true;
        k as usize
    }
}

/// Inputs describing the atoms.
pub struct RidingAtoms<'a> {
    pub is_h: &'a [bool],
    pub altloc: &'a [ShortStr],
    pub name: &'a [ShortStr],
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
    /// The name of the atom a hydrogen must be cis to, when a naming
    /// convention fixes it ([`crate::names::amino_h_cis_partner`]). It
    /// overrides the dictionary torsion of that hydrogen: a period-2 torsion
    /// does not tell the two H of an NH2 apart, and the constant torsions of
    /// GeoStd's DG name its amino group the other way round. Reduce2 has none.
    pub cis_partner: &'a dyn Fn(u32) -> Option<&'static str>,
    /// Orient a one-neighbor group by another atom when its third neighbor lies
    /// on the parent bond axis (a ligand on a crystallographic two-fold, say).
    /// Reduce2 divides by zero there.
    pub reroute_axial_reference: bool,
    /// Use explicit planes when a single H has inconsistent angle restraints.
    pub use_plane_restraints: bool,
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
    // The neighbors in bond order, a stable sort of the bonds by (i, j),
    // which lists each atom's neighbors in ascending order.
    let (fsc0, fsc0_bond) = Adjacency::sorted(n, &it.bonds);
    let mut conn = Conn::new(at.is_h);
    // the parent of each record of `conn` from step 1, by slot
    let mut parent_of: Vec<u32> = vec![u32::MAX; conn.nbs.len()];
    let mut double_h = vec![false; n];
    let mut is_parent = vec![false; n];
    // 1. first neighbors: the parent of an H is the heavy atom of its first
    // bond in bond order, which is its lowest heavy neighbor
    for ih in 0..n as u32 {
        if !at.is_h[ih as usize] {
            continue;
        }
        let entries = fsc0.range(ih);
        let Some(e) = entries.clone().find(|&e| !at.is_h[fsc0.of_entry(e) as usize]) else { continue };
        let parent = fsc0.of_entry(e);
        let ideal = it.bonds[fsc0_bond[e] as usize].ideal;
        if entries.len() > 1 {
            double_h[ih as usize] = true;
        }
        let k = conn.set(ih, Nb { ih, a0: Some((parent, ideal)), valid: true, ..Default::default() });
        parent_of.resize(conn.nbs.len(), u32::MAX);
        parent_of[k] = parent;
        is_parent[parent as usize] = true;
    }
    // 2. slipped H (no bond)
    let slipped: Vec<u32> = (0..n as u32).filter(|&k| at.is_h[k as usize] && conn.index(k).is_none()).collect();
    // 3. second neighbors (angle proxy order), by record of `conn`
    let mut second_raw: Vec<Neighbors> = vec![Neighbors::new(); conn.nbs.len()];
    // the ideal angle of each entry of second_raw; an H has one parent here,
    // and a later entry for the same second neighbor replaces an earlier one
    let mut second_angle: Vec<SmallVec<[f64; 4]>> = vec![SmallVec::new(); conn.nbs.len()];
    // the angles with two heavy end atoms, in angle order, for step 7
    let mut heavy_angles: Vec<u32> = Vec::new();
    for (ai, ap) in it.angles.iter().enumerate() {
        if !at.is_h[ap.i[0] as usize] && !at.is_h[ap.i[2] as usize] {
            heavy_angles.push(ai as u32);
        }
        for &i_test in &ap.i {
            let Some(k) = conn.index(i_test) else { continue };
            let ih = i_test;
            let parent = ap.i[1];
            if !is_parent[parent as usize] {
                continue;
            }
            let second = match ap.i.iter().find(|&&x| x != ih && x != parent) {
                Some(&s) => s,
                None => continue,
            };
            if parent_of[k] != parent {
                if double_h[ih as usize] {
                    continue;
                }
                warnings.push(format!("angle and bond restraints conflict for H {}", ih));
                continue;
            }
            second_raw[k].push(second);
            second_angle[k].push(ap.ideal);
        }
    }
    // 4. Retain Reduce2's first-H placement for restraints that already work.
    let mut plane_h = vec![false; n];
    let mut explicit_plane_h = vec![false; n];
    for pl in &it.planes {
        if let Some(&h) = pl.iter().find(|&&a| at.is_h[a as usize]) {
            plane_h[h as usize] = true;
        }
        if at.use_plane_restraints {
            for &a in pl.iter().filter(|&&a| at.is_h[a as usize]) {
                explicit_plane_h[a as usize] = true;
            }
        }
    }
    // 5. process second neighbors
    let mut a1_atoms = vec![false; n];
    let mut a0a1: FxHashMap<u32, Neighbors> = FxHashMap::default();
    let (mut reduced, mut alts, mut hs, mut non_h): (Vec<u32>, Vec<u32>, Vec<u32>, Vec<u32>) = Default::default();
    for ih in 0..n as u32 {
        let Some(slot) = conn.index(ih) else { continue };
        let Some(nb) = conn.nbs[slot].as_mut() else { continue };
        let parent = parent_of[slot];
        reduced.clear();
        alts.clear();
        for &s in &second_raw[slot] {
            if at.altloc[s as usize].is_empty() {
                reduced.push(s);
            } else {
                alts.push(s);
            }
        }
        if !alts.is_empty() {
            reduced.extend(reduce_alternates(at, &alts));
        }
        hs.clear();
        hs.extend(reduced.iter().copied().filter(|&x| at.is_h[x as usize]));
        // non-H: Python set order (iteration by hash); positions do not depend on it
        non_h.clear();
        for &x in &reduced {
            if !at.is_h[x as usize] && !non_h.contains(&x) {
                non_h.push(x);
            }
        }
        let angle_to = |x: u32| {
            let k = second_raw[slot].iter().rposition(|&s| s == x).unwrap_or(0);
            second_angle[slot][k]
        };
        let mut cnt = 0;
        for (k, &x) in non_h.iter().take(3).enumerate() {
            nb.a[k] = Some(NbAngle { iseq: x, angle: angle_to(x) });
            cnt += 1;
        }
        let mut hc = 0;
        for (k, &x) in hs.iter().take(2).enumerate() {
            nb.h[k] = Some(NbAngle { iseq: x, angle: angle_to(x) });
            hc += 1;
        }
        nb.nh = hc;
        nb.nnonh = cnt;
        if cnt == 1 {
            let a1 = nb.a[0].unwrap().iseq;
            a1_atoms[a1 as usize] = true;
            a0a1.entry(parent).or_default().push(a1);
        }
    }
    // 6. third neighbors from dihedral proxies (const first, then var; last write wins)
    let assign_group = |conn: &mut Conn, ih: u32, third: u32| {
        let Some(nb) = conn.get(ih) else { return };
        let (nh, h1, h2) = (nb.nh, nb.h[0].map(|x| x.iseq), nb.h[1].map(|x| x.iseq));
        if nh > 0 {
            if let Some(h1) = h1 {
                if let Some(o) = conn.get_mut(h1) {
                    o.b1_iseq = Some(third);
                    o.b1_ideal = None;
                }
            }
            if nh == 2 {
                if let Some(h2) = h2 {
                    if let Some(o) = conn.get_mut(h2) {
                        o.b1_iseq = Some(third);
                        o.b1_ideal = None;
                    }
                }
            }
        }
    };
    // the torsion that puts H `t` cis (0) or trans (180) to the third atom of
    // its dihedral by the naming convention, if one applies
    let named = |t: u32, third: u32| {
        (at.cis_partner)(t).map(|partner| if at.name[third as usize].trim() == partner { 0.0 } else { 180.0 })
    };
    for dp in &it.const_dihedrals {
        for &t in &dp.i {
            if conn.index(t).is_none() {
                continue;
            }
            let (third, ideal) = if t == dp.i[0] {
                (dp.i[3], dp.ideal)
            } else if t == dp.i[3] {
                (dp.i[0], -dp.ideal)
            } else {
                continue;
            };
            if let Some(nb) = conn.get_mut(t) {
                nb.b1_iseq = Some(third);
                nb.b1_ideal = Some(named(t, third).unwrap_or(ideal));
            }
            assign_group(&mut conn, t, third);
        }
    }
    'dp: for dp in &it.dihedrals {
        for &t in &dp.i {
            if conn.index(t).is_none() {
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
            if let Some(nb) = conn.get_mut(t) {
                nb.b1_iseq = Some(third);
                let nh2 = at.dictionary_nh2_torsion && dp.period == 2 && nb.nh == 1;
                // Both images of a period-2 torsion fit the dictionary. Groups
                // without a naming convention keep the image opposite its
                // value, which names the amino groups of GeoStd's RNA A and C
                // as IUPAC does.
                nb.b1_ideal = Some(match named(t, third) {
                    Some(torsion) => torsion,
                    None if nh2 => dp.ideal - 180.0,
                    None => model_deg + delta,
                });
            }
            assign_group(&mut conn, t, third);
        }
    }
    // 7. angles at the parent and fallback third neighbors
    // heavy-atom angles at each parent, in angle order: (first, third atom,
    // ideal), grouped by the parent
    let mut pa_start = vec![0u32; n + 1];
    for &ai in &heavy_angles {
        let iy = it.angles[ai as usize].i[1];
        if is_parent[iy as usize] {
            pa_start[iy as usize + 1] += 1;
        }
    }
    for k in 0..n {
        pa_start[k + 1] += pa_start[k];
    }
    let mut pa_fill = pa_start.clone();
    let mut pa_list: Vec<(u32, u32, f64)> = vec![(0, 0, 0.0); pa_start[n] as usize];
    let mut third_raw: FxHashMap<u32, Neighbors> = FxHashMap::default();
    for &ai in &heavy_angles {
        let ap = &it.angles[ai as usize];
        let [ix, iy, iz] = ap.i;
        if is_parent[iy as usize] {
            pa_list[pa_fill[iy as usize] as usize] = (ix, iz, ap.ideal);
            pa_fill[iy as usize] += 1;
        }
        if !a1_atoms[iy as usize] {
            continue;
        }
        let (parent, third) = if is_parent[ix as usize] && a0a1.contains_key(&ix) {
            if a0a1[&ix].contains(&iy) {
                (ix, iz)
            } else {
                continue;
            }
        } else if is_parent[iz as usize] && a0a1.contains_key(&iz) {
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
            if is_parent[third as usize] {
                third_raw.entry(third).or_default().push(parent);
            }
        }
    }
    // the last angle (a0; x, z) recorded, else the last (a0; z, x)
    let get_pa = |a0: u32, x: u32, z: u32| -> Option<f64> {
        let list = &pa_list[pa_start[a0 as usize] as usize..pa_start[a0 as usize + 1] as usize];
        let last = |p: u32, q: u32| list.iter().rev().find(|e| e.0 == p && e.1 == q).map(|e| e.2);
        last(x, z).or_else(|| last(z, x))
    };
    for ih in 0..n {
        let Some(slot) = conn.index(ih as u32) else { continue };
        let Some(nb) = conn.nbs[slot].as_mut() else { continue };
        // the record as it was before this step
        let (nnonh, nh, b1_iseq, h0) = (nb.nnonh, nb.nh, nb.b1_iseq, nb.h[0]);
        let parent = parent_of[slot];
        // assign_a0_angles
        if let (true, Some(a1), Some(a2)) = (nnonh > 1, nb.a[0], nb.a[1]) {
            let (a1, a2) = (a1.iseq, a2.iseq);
            let ang = get_pa(parent, a1, a2);
            nb.angle_a1a0a2 = ang;
            if ang.is_none() {
                *nb = Nb { ih: ih as u32, valid: false, ..Default::default() };
                continue;
            }
            if let (3, Some(a3)) = (nnonh, nb.a[2]) {
                let a3 = a3.iseq;
                nb.angle_a2a0a3 = get_pa(parent, a2, a3);
                nb.angle_a3a0a1 = get_pa(parent, a3, a1);
                if nb.angle_a2a0a3.is_none() || nb.angle_a3a0a1.is_none() {
                    *nb = Nb { ih: ih as u32, valid: false, ..Default::default() };
                    continue;
                }
            }
        }
        if nnonh != 1 || b1_iseq.is_some() {
            continue;
        }
        let raw = third_raw.get(&parent).map_or(&[][..], |v| v.as_slice());
        reduced.clear();
        alts.clear();
        for &t in raw {
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
            conn.nbs[slot] = Some(Nb { ih: ih as u32, valid: false, ..Default::default() });
            continue;
        }
        if parent_lost_heavy_neighbor(at, &fsc0, parent) {
            conn.nbs[slot] = Some(Nb { ih: ih as u32, valid: false, ..Default::default() });
            continue;
        }
        let third = reduced[0];
        let Some(nb) = conn.nbs[slot].as_mut() else { continue };
        if nh == 2 {
            nb.b1_iseq = Some(third);
            let sibling = h0.map(|x| x.iseq);
            nb.b1_ideal = Some(if sibling.map(|s| (ih as u32) < s).unwrap_or(false) { 60.0 } else { -60.0 });
        } else {
            nb.b1_iseq = Some(third);
            nb.b1_ideal = Some(180.0);
        }
    }
    for &s in &slipped {
        conn.set(s, Nb { ih: s, valid: false, ..Default::default() });
    }
    for ih in 0..n {
        if let Some(nb) = conn.get_mut(ih as u32) {
            nb.is_in_plane = plane_h[ih];
            nb.explicit_plane = explicit_plane_h[ih];
        }
    }

    if at.reroute_axial_reference {
        reroute_axial_references(at, &fsc0, sites, &mut conn);
    }

    // ---------------- parameterization
    let mut axial_reference: Vec<u32> = Vec::new();
    let mut coef: Vec<Option<RidingCoef>> = std::iter::repeat_with(|| None).take(n).collect();
    let mut unk: Vec<u32> = Vec::new();
    for ih in 0..n {
        let Some(nb) = conn.get(ih as u32) else { continue };
        if coef[ih].is_some() {
            continue;
        }
        if !nb.valid || nb.a0.is_none() {
            unk.push(ih as u32);
            continue;
        }
        let r = if nb.nnonh == 2 {
            process_2(nb, sites, &mut coef, at.use_plane_restraints,
                &|| planar::dictionary_carbon(it, at, &fsc0, nb))
        } else if nb.nnonh == 3 && nb.nh == 0 {
            process_3(nb, sites, &mut coef)
        } else if nb.nnonh == 1 && (nb.nh == 0 || nb.nh == 2) {
            process_1(nb, &conn, sites, &mut coef, use_ideal_dihedral, &mut axial_reference)
        } else if nb.nnonh == 1 && nb.nh == 1 {
            process_1_arg(nb, &conn, sites, &mut coef, &mut axial_reference)
        } else {
            Err(())
        };
        if r.is_err() {
            unk.push(ih as u32);
        }
    }
    // ---------------- positions (i_seq order, written immediately)
    // only hydrogens get coefficients
    let mut unparameterized: Vec<u32> = Vec::new();
    for ih in 0..n {
        if !at.is_h[ih] {
            continue;
        }
        match coef[ih] {
            Some(c) => {
                if let Some(p) = compute_h_position(&c, sites) {
                    sites[ih] = p;
                }
            }
            None => unparameterized.push(ih as u32),
        }
    }
    let _ = unk;
    RidingResult { coef, unparameterized, axial_reference, warnings }
}

fn parent_lost_heavy_neighbor(at: &RidingAtoms, fsc0: &Adjacency, parent: u32) -> bool {
    let Some(expected) = (at.expected_heavy)(parent) else { return false };
    let palt = &at.altloc[parent as usize];
    let mut seen: FxHashSet<&str> = FxHashSet::default();
    for &j in fsc0.of(parent) {
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
fn reroute_axial_references(at: &RidingAtoms, fsc0: &Adjacency, sites: &[Vec3], conn: &mut Conn) {
    for ih in 0..conn.slot.len() as u32 {
        let Some(nb) = conn.get(ih) else { continue };
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
        let bonded = fsc0.of(a1).iter().copied().find(|&c| usable(c));
        let nearest = || {
            (0..sites.len() as u32)
                .filter(|&c| usable(c))
                .min_by(|&x, &y| (sites[x as usize] - r1).length().total_cmp(&(sites[y as usize] - r1).length()))
        };
        let b1 = bonded.or_else(nearest);
        if let Some(nb) = conn.get_mut(ih) {
            nb.b1_iseq = b1;
        }
    }
}

fn process_1(
    nb: &Nb,
    conn: &Conn,
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
        } else if let Some(o) = conn.get(h1).filter(|x| x.b1_ideal.is_some()) {
            if coef[h1 as usize].is_none() {
                nbs = o.clone();
            }
        } else if let Some(o) = conn.get(h2).filter(|x| x.b1_ideal.is_some() && coef[h2 as usize].is_none()) {
            nbs = o.clone();
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
    conn: &Conn,
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
    } else if conn.get(h1.iseq).is_some_and(|x| x.b1_ideal.is_some()) {
        (h1.iseq, ih)
    } else {
        return Err(());
    };
    let nbd = conn.get(ih_d).ok_or(())?;
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

fn process_2(nb: &Nb, sites: &[Vec3], coef: &mut [Option<RidingCoef>], use_plane_restraints: bool, dictionary_plane: &dyn Fn() -> bool) -> Result<(), ()> {
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
    if use_plane_restraints && nb.nh == 0 && sumang > 2.0 * PI + 0.05 && root < 0.0
        && (nb.explicit_plane || dictionary_plane()) {
        htype = HType::Flat2;
    } else if sumang > 2.0 * PI + 0.05 && root < 0.0 {
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
