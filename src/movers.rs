//! Movers: descriptions of the alternative positions of groups of atoms
//! (port of mmtbx/reduce/Movers.py).
//!
//! A Mover never moves atoms by itself during optimization; it describes
//! coarse and fine states (positions, optional per-state atom info and
//! deletions, preference energies) and a final fix-up. Like the original,
//! constructors do move atoms into a canonical starting orientation.

use crate::geom::*;
use crate::probe::AtomInfo;
use crate::world::World;

#[derive(Clone, Debug, PartialEq)]
pub enum MoverKind {
    SingleHydrogenRotator,
    NH3Rotator,
    AromaticMethylRotator,
    AmideFlip,
    HisFlip { enabled: u8, enable_fixup: bool },
}

impl MoverKind {
    pub fn class_name(&self) -> &'static str {
        match self {
            MoverKind::SingleHydrogenRotator => "MoverSingleHydrogenRotator",
            MoverKind::NH3Rotator => "MoverNH3Rotator",
            MoverKind::AromaticMethylRotator => "MoverAromaticMethylRotator",
            MoverKind::AmideFlip => "MoverAmideFlip",
            MoverKind::HisFlip { .. } => "MoverHisFlip",
        }
    }
    pub fn is_flip(&self) -> bool {
        matches!(self, MoverKind::AmideFlip | MoverKind::HisFlip { .. })
    }
}

#[derive(Clone, Debug, Default)]
pub struct RotatorMeta {
    pub offset: f64,
    pub coarse_angles: Vec<f64>,
    pub fine_angles: Vec<f64>,
    pub axis_origin: Vec3,
    pub axis_dir: Vec3,
}

#[derive(Clone, Debug)]
pub struct Mover {
    pub kind: MoverKind,
    /// `CoarsePositions().atoms`.
    pub atoms: Vec<u32>,
    /// Number of leading atoms that have positions in the coarse/fine states.
    pub n_moved: usize,
    pub coarse_pos: Vec<Vec<Vec3>>,
    /// Per coarse state: atom info for the leading atoms (empty when unchanged).
    pub coarse_info: Vec<Vec<AtomInfo>>,
    /// Per coarse state: deletion flags for the leading atoms (empty = none).
    pub coarse_del: Vec<Vec<bool>>,
    pub coarse_pref: Vec<f64>,
    pub fine_pos: Vec<Vec<Vec<Vec3>>>,
    pub fine_pref: Vec<Vec<f64>>,
    pub fixup_pos: Vec<Vec<Vec3>>,
    pub fixup_info: Vec<Vec<AtomInfo>>,
    pub fixup_del: Vec<Vec<bool>>,
    pub rot: Option<RotatorMeta>,
}

impl Mover {
    pub fn num_coarse(&self) -> usize {
        self.coarse_pos.len()
    }

    /// `PoseDescription(coarseIndex, fineIndex, fixedUp)`.
    pub fn pose_description(&self, coarse: usize, fine: Option<usize>, fixed_up: bool) -> String {
        let nfine0 = self.fine_pos.first().map(|f| f.len()).unwrap_or(0);
        let bad = coarse >= self.num_coarse() || matches!(fine, Some(f) if f > 0 && f >= nfine0);
        match &self.kind {
            MoverKind::SingleHydrogenRotator | MoverKind::NH3Rotator | MoverKind::AromaticMethylRotator => {
                if bad {
                    return "Unrecognized state . .".into();
                }
                let r = self.rot.as_ref().unwrap();
                let fine_off = fine.map(|f| r.fine_angles[f]).unwrap_or(0.0);
                let mut angle = r.offset + r.coarse_angles[coarse] + fine_off;
                while angle > 180.0 {
                    angle -= 360.0;
                }
                while angle < -180.0 {
                    angle += 360.0;
                }
                format!("Angle {} deg .", fmt_f(angle, 1))
            }
            MoverKind::AmideFlip => {
                let fstr = if coarse == 1 {
                    if fixed_up { "AnglesAdjusted" } else { "AnglesNotAdjusted" }
                } else {
                    "."
                };
                if bad {
                    "Unrecognized state .".into()
                } else if coarse == 0 {
                    format!("Unflipped . . {}", fstr)
                } else {
                    format!("Flipped . . {}", fstr)
                }
            }
            MoverKind::HisFlip { enabled, enable_fixup } => {
                let mut ret = match enabled {
                    3 => {
                        if bad {
                            return "Unrecognized state . .".into();
                        }
                        if coarse < 4 { "Unflipped".to_string() } else { "Flipped".to_string() }
                    }
                    1 => "Unflipped".to_string(),
                    2 => "Flipped".to_string(),
                    _ => return "Unrecognized flip states .".into(),
                };
                ret += if coarse % 4 == 0 || coarse % 4 == 1 { " HD1Placed" } else { " HD1NotPlaced" };
                ret += if coarse % 4 == 0 || coarse % 4 == 2 { " HE2Placed" } else { " HE2NotPlaced" };
                let half = self.num_coarse() as f64 / 2.0;
                if *enabled == 2 || (*enabled == 3 && (coarse as f64) >= half) {
                    if *enable_fixup && fixed_up {
                        ret += " AnglesAdjusted";
                    } else {
                        ret += " AnglesNotAdjusted";
                    }
                } else {
                    ret += " .";
                }
                ret
            }
        }
    }
}

/// Python-compatible fixed-point formatting ("{:.Nf}"), which rounds half to
/// even on the exact binary value just like Rust's formatter.
pub fn fmt_f(v: f64, prec: usize) -> String {
    let s = format!("{:.*}", prec, v);
    // Python prints "-0.0" for negative values that round to zero; Rust does too.
    s
}

pub type MResult<T> = Result<T, String>;

// ----------------------------------------------------------------------------
// Rotators

pub struct RotatorSpec {
    pub atoms: Vec<u32>,
    pub axis_origin: Vec3,
    pub axis_dir: Vec3,
    pub dihedral: f64,
    pub offset: f64,
    pub coarse_range: f64,
    pub coarse_step: f64,
    pub do_fine: bool,
    pub fine_step: f64,
    pub preference: Option<fn(f64) -> f64>,
    pub pref_scale: f64,
}

fn rotator_angles(range: f64, step: f64) -> Vec<f64> {
    let mut v = vec![0.0];
    let mut cur = step;
    while cur <= range {
        v.push(-cur);
        if cur < range {
            v.push(cur);
        }
        cur += step;
    }
    v
}

fn fine_angles(coarse_step: f64, fine_step: f64) -> Vec<f64> {
    let mut v = Vec::new();
    let mut cur = fine_step;
    let range = coarse_step / 2.0;
    while cur <= range {
        v.push(-cur);
        if cur < range {
            v.push(cur);
        }
        cur += fine_step;
    }
    v
}

fn poses_for(w: &World, atoms: &[u32], origin: Vec3, dir: Vec3, angles: &[f64]) -> Vec<Vec<Vec3>> {
    angles
        .iter()
        .map(|&ang| atoms.iter().map(|&a| rotate_deg_axis_dir(origin, dir, w.pos[a as usize], ang)).collect())
        .collect()
}

/// `_MoverRotator.__init__`: rotate the atoms into the canonical start and
/// build coarse and fine poses.
fn build_rotator(w: &mut World, kind: MoverKind, spec: RotatorSpec) -> Mover {
    for &a in &spec.atoms {
        let p = w.pos[a as usize];
        w.set_pos(a, rotate_deg_axis_dir(spec.axis_origin, spec.axis_dir, p, spec.offset + spec.dihedral));
    }
    let coarse_angles = rotator_angles(spec.coarse_range, spec.coarse_step);
    let fine = fine_angles(spec.coarse_step, spec.fine_step);
    let mut m = Mover {
        kind,
        atoms: spec.atoms.clone(),
        n_moved: spec.atoms.len(),
        coarse_pos: vec![],
        coarse_info: vec![],
        coarse_del: vec![],
        coarse_pref: vec![],
        fine_pos: vec![],
        fine_pref: vec![],
        fixup_pos: vec![],
        fixup_info: vec![],
        fixup_del: vec![],
        rot: Some(RotatorMeta {
            offset: spec.offset,
            coarse_angles,
            fine_angles: if spec.do_fine { fine } else { vec![] },
            axis_origin: spec.axis_origin,
            axis_dir: spec.axis_dir,
        }),
    };
    recompute_rotator_positions(w, &mut m, spec.preference, spec.pref_scale, spec.do_fine);
    m
}

fn recompute_rotator_positions(w: &World, m: &mut Mover, pref: Option<fn(f64) -> f64>, scale: f64, do_fine: bool) {
    let r = m.rot.as_ref().unwrap();
    let pf = |angles: &[f64]| -> Vec<f64> {
        match pref {
            Some(f) => angles.iter().map(|&a| f(a) * scale).collect(),
            None => vec![0.0; angles.len()],
        }
    };
    m.coarse_pos = poses_for(w, &m.atoms, r.axis_origin, r.axis_dir, &r.coarse_angles);
    m.coarse_pref = pf(&r.coarse_angles);
    let n = r.coarse_angles.len();
    m.coarse_info = vec![vec![]; n];
    m.coarse_del = vec![vec![]; n];
    m.fine_pos = Vec::with_capacity(n);
    m.fine_pref = Vec::with_capacity(n);
    for &ca in &r.coarse_angles {
        if !do_fine {
            m.fine_pos.push(vec![]);
            m.fine_pref.push(vec![]);
            continue;
        }
        let angles: Vec<f64> = r.fine_angles.iter().map(|fa| fa + ca).collect();
        m.fine_pos.push(poses_for(w, &m.atoms, r.axis_origin, r.axis_dir, &angles));
        m.fine_pref.push(pf(&angles));
    }
    m.fixup_pos = vec![vec![]; n];
    m.fixup_info = vec![vec![]; n];
    m.fixup_del = vec![vec![]; n];
}

fn pref_120(deg: f64) -> f64 {
    0.1 + 0.1 * (deg * (std::f64::consts::PI / 180.0) * (360.0 / 120.0)).cos()
}

/// Riding-parameter lookup needed to choose the conventional dihedral.
#[derive(Clone, Copy, Debug)]
pub struct RidingRef {
    pub n: i32,
    pub a2: i64,
}

/// `dihedralChoicesForRotatableHydrogens`: the hydrogen whose riding `n` is 0
/// and the potential third-neighbor named by its `a2` (last match wins).
/// With `fallback`, a group whose hydrogens name none of the potentials (a
/// methanol methyl, whose dictionary has no torsion and no heavy third
/// neighbor) is measured from the first pair that defines a dihedral.
fn dihedral_choice_or_any(w: &World, hydrogens: &[u32], potentials: &[u32], partner: u32, atom: u32, fallback: bool) -> MResult<(u32, u32)> {
    let named = dihedral_choice(w, hydrogens, potentials);
    if named.is_ok() || !fallback {
        return named;
    }
    hydrogens
        .iter()
        .flat_map(|&h| potentials.iter().map(move |&p| (h, p)))
        .find(|&(h, p)| dihedral_of(w, h, partner, atom, p).is_ok())
        .ok_or_else(|| named.unwrap_err())
}

fn dihedral_choice(w: &World, hydrogens: &[u32], potentials: &[u32]) -> MResult<(u32, u32)> {
    let mut res = None;
    for &h in hydrogens {
        if let Some(Some(item)) = w.riding.get(h as usize) {
            if item.n == 0 {
                for &p in potentials {
                    if p as i64 == item.a2 {
                        res = Some((h, p));
                    }
                }
            }
        }
    }
    res.ok_or_else(|| "mmtbx.probe.Helpers.dihedralChoicesForRotatableHydrogens(): Could not determine atoms to use".to_string())
}

fn dihedral_of(w: &World, a: u32, b: u32, c: u32, d: u32) -> MResult<f64> {
    dihedral_deg(w.pos[a as usize], w.pos[b as usize], w.pos[c as usize], w.pos[d as usize])
        .ok_or_else(|| "unsupported operand type(s) for +: 'int' and 'NoneType'".to_string())
}

pub struct SingleHOptions {
    pub circular_angle_spacing: bool,
    /// Rotate a hydrogen whose partner has any number of other bonds
    /// (S-OH, P-OH, a metal oxo-hydroxide); Reduce2 needs two or three.
    pub any_partner_valence: bool,
}

/// `MoverSingleHydrogenRotator`.
pub fn single_hydrogen_rotator(
    w: &mut World,
    atom: u32,
    potential_acceptors: &[u32],
    potential_touches: &[u32],
    opts: &SingleHOptions,
) -> MResult<Mover> {
    if w.elem(atom) != "H" {
        return Err("MoverSingleHydrogenRotator(): Atom is not a hydrogen".into());
    }
    let neighbors = w.bonded[atom as usize].clone();
    if neighbors.len() != 1 {
        return Err("MoverSingleHydrogenRotator(): Atom does not have a single bonded neighbor".into());
    }
    let neighbor = neighbors[0];
    let partners = w.bonded[neighbor as usize].clone();
    if partners.len() != 2 {
        return Err("MoverSingleHydrogenRotator(): Atom does not have a single bonded neighbor's neighbor".into());
    }
    let mut partner = partners[0];
    if partner == atom {
        partner = partners[1];
    }
    let friends: Vec<u32> = w.bonded[partner as usize].iter().copied().filter(|&b| b != neighbor).collect();
    let valence_ok = if opts.any_partner_valence { !friends.is_empty() } else { friends.len() == 2 || friends.len() == 3 };
    if !valence_ok {
        return Err(format!(
            "MoverSingleHydrogenRotator(): Atom's bonded neighbor's neighbor does not have 2-3 other bonds it has {}",
            friends.len()
        ));
    }
    let normal = (w.pos[neighbor as usize] - w.pos[partner as usize]).normalize();
    let origin = w.pos[partner as usize];
    let atoms = vec![atom, neighbor];
    let (conv_h, conv_friend) = dihedral_choice_or_any(w, &atoms, &friends, partner, neighbor, opts.any_partner_valence)?;
    let dihedral = dihedral_of(w, conv_h, partner, neighbor, conv_friend)?;
    let mut m = build_rotator(
        w,
        MoverKind::SingleHydrogenRotator,
        RotatorSpec {
            atoms,
            axis_origin: origin,
            axis_dir: normal,
            dihedral,
            offset: 180.0,
            coarse_range: 180.0,
            coarse_step: 10.0,
            do_fine: true,
            fine_step: 1.0,
            preference: None,
            pref_scale: 1.0,
        },
    );

    // Orientations toward potential acceptors, sorted by an atom-ID string.
    let mut acc: Vec<u32> = potential_acceptors.to_vec();
    acc.sort_by_cached_key(|&a| w.atom_id_string(a));
    let mut acceptor_angles = Vec::new();
    for a in acc {
        if let Some(deg) = dihedral_deg(w.pos[atom as usize], w.pos[partner as usize], w.pos[neighbor as usize], w.pos[a as usize]) {
            acceptor_angles.push(deg);
        }
    }
    // Best "just touching" angle among the coarse angles.
    let mut best_touch_angle = 0.0;
    let mut best_touch_gap = 1e100;
    let ra = w.info[atom as usize].vdw_radius;
    let coarse_angles = m.rot.as_ref().unwrap().coarse_angles.clone();
    for (i, &ang) in coarse_angles.iter().enumerate() {
        let mut min_gap = 1e100;
        for &pt in potential_touches {
            let rt = w.info[pt as usize].vdw_radius;
            let dist = m.coarse_pos[i][0].dist(w.pos[pt as usize]);
            let gap = dist - (ra + rt);
            if gap < min_gap {
                min_gap = gap;
            }
        }
        if min_gap.abs() < best_touch_gap {
            best_touch_gap = min_gap.abs();
            best_touch_angle = ang;
        }
    }
    let mut sofar = vec![best_touch_angle];
    sofar.extend(acceptor_angles);
    for &ang in &coarse_angles {
        let mut min_ang = 360.0f64;
        for &a in &sofar {
            let mut diff = (a - ang).abs();
            if opts.circular_angle_spacing {
                diff = diff % 360.0;
                diff = diff.min(360.0 - diff);
            }
            if diff < min_ang {
                min_ang = diff;
            }
        }
        if min_ang >= 45.0 {
            sofar.push(ang);
        }
    }
    m.rot.as_mut().unwrap().coarse_angles = sofar;
    recompute_rotator_positions(w, &mut m, None, 1.0, true);
    Ok(m)
}

/// Shared checks for the three-hydrogen rotators. Returns (hydrogens, partner, friends).
fn three_h_group(w: &World, atom: u32, elem: &str, who: &str) -> MResult<(Vec<u32>, u32, Vec<u32>)> {
    let atomword = if elem == "N" { "Nitrogen" } else { "Carbon" };
    if w.elem(atom) != elem {
        return Err(format!("{}(): atom is not a {}", who, atomword));
    }
    let partners = &w.bonded[atom as usize];
    if partners.len() != 4 {
        return Err(format!("{}(): atom does not have four bonded neighbors", who));
    }
    let mut hydrogens = Vec::new();
    let mut partner = None;
    for &a in partners {
        if w.is_h(a) {
            hydrogens.push(a);
        } else {
            partner = Some(a);
        }
    }
    if hydrogens.len() != 3 {
        return Err(format!("{}(): atom does not have three bonded hydrogens", who));
    }
    let partner = partner.unwrap();
    let friends: Vec<u32> = w.bonded[partner as usize].iter().copied().filter(|&b| b != atom).collect();
    Ok((hydrogens, partner, friends))
}

/// `MoverNH3Rotator`.
/// With `any_partner_valence`, the partner may have any number of other bonds
/// (an ammine on a two-coordinate metal, say); Reduce2 needs at least three.
pub fn nh3_rotator(w: &mut World, atom: u32, any_partner_valence: bool) -> MResult<Mover> {
    let (hydrogens, partner, friends) = three_h_group(w, atom, "N", "MoverNH3Rotator")?;
    if friends.is_empty() || (friends.len() < 3 && !any_partner_valence) {
        return Err("MoverNH3Rotator(): Partner does not have at least three bonded friends".into());
    }
    let preference: Option<fn(f64) -> f64> = if friends.len() == 3 { Some(pref_120) } else { None };
    let normal = (w.pos[atom as usize] - w.pos[partner as usize]).normalize();
    let origin = w.pos[partner as usize];
    let (conv_h, conv_friend) = dihedral_choice_or_any(w, &hydrogens, &friends, partner, atom, any_partner_valence)?;
    let dihedral = dihedral_of(w, conv_h, partner, atom, conv_friend)?;
    let mut atoms = vec![atom];
    atoms.extend(&hydrogens);
    Ok(build_rotator(
        w,
        MoverKind::NH3Rotator,
        RotatorSpec {
            atoms,
            axis_origin: origin,
            axis_dir: normal,
            dihedral,
            offset: 180.0,
            coarse_range: 60.0,
            coarse_step: 15.0,
            do_fine: true,
            fine_step: 1.0,
            preference,
            pref_scale: 1.0,
        },
    ))
}

/// `MoverAromaticMethylRotator`.
pub fn aromatic_methyl_rotator(w: &mut World, atom: u32, any_reference: bool) -> MResult<Mover> {
    let (hydrogens, partner, friends) = three_h_group(w, atom, "C", "MoverAromaticMethylRotator")?;
    if friends.len() != 2 {
        return Err("MoverAromaticMethylRotator(): Partner does not have two bonded friends".into());
    }
    let normal = (w.pos[atom as usize] - w.pos[partner as usize]).normalize();
    let origin = w.pos[partner as usize];
    let (conv_h, conv_friend) = dihedral_choice_or_any(w, &hydrogens, &friends, partner, atom, any_reference)?;
    let dihedral = dihedral_of(w, conv_h, partner, atom, conv_friend)?;
    let mut atoms = vec![atom];
    atoms.extend(&hydrogens);
    Ok(build_rotator(
        w,
        MoverKind::AromaticMethylRotator,
        RotatorSpec {
            atoms,
            axis_origin: origin,
            axis_dir: normal,
            dihedral,
            offset: 180.0 + 90.0,
            coarse_range: 180.0,
            coarse_step: 180.0,
            do_fine: false,
            fine_step: 1.0,
            preference: None,
            pref_scale: 1.0,
        },
    ))
}

/// `MoverTetrahedralMethylRotator`: only its side effect (staggering the
/// hydrogens) is used by Reduce2.
pub fn stagger_tetrahedral_methyl(w: &mut World, atom: u32, any_reference: bool) -> MResult<()> {
    let (hydrogens, partner, friends) = three_h_group(w, atom, "C", "MoverTetrahedralMethylRotator")?;
    if friends.len() != 1 && friends.len() != 3 {
        return Err("MoverTetrahedralMethylRotator(): Partner does not have one or three bonded friends".into());
    }
    let normal = (w.pos[atom as usize] - w.pos[partner as usize]).normalize();
    let origin = w.pos[partner as usize];
    let (conv_h, conv_friend) = dihedral_choice_or_any(w, &hydrogens, &friends, partner, atom, any_reference)?;
    let dihedral = dihedral_of(w, conv_h, partner, atom, conv_friend)?;
    let mut atoms = vec![atom];
    atoms.extend(&hydrogens);
    for &a in &atoms {
        let p = w.pos[a as usize];
        w.set_pos(a, rotate_deg_axis_dir(origin, normal, p, 180.0 + dihedral));
    }
    Ok(())
}

// ----------------------------------------------------------------------------
// Flips

/// `_rotateHingeDock` (three-point dock after a flip). `clamp_acos` guards
/// against rounding pushing a cosine just outside [-1, 1] (the original raised
/// a math domain error and silently dropped the Mover).
fn rotate_hinge_dock(
    w: &World,
    movable_atoms: &[u32],
    hinge_index: usize,
    first_dock: usize,
    second_dock: usize,
    alpha_carbon: u32,
    clamp_acos: bool,
) -> MResult<Vec<Vec3>> {
    let acos_deg = |x: f64| -> MResult<f64> {
        if clamp_acos {
            Ok(acos_deg_clamped(x))
        } else if (-1.0..=1.0).contains(&x) {
            Ok(x.acos() * (180.0 / std::f64::consts::PI))
        } else {
            Err("math domain error".into())
        }
    };
    let p = |a: u32| w.pos[a as usize];
    let first = p(movable_atoms[first_dock]);
    let second = p(movable_atoms[second_dock]);
    let hinge_atom = p(movable_atoms[hinge_index]);
    let pivot_atom = p(movable_atoms[hinge_index + 1]);
    let ca = p(alpha_carbon);
    let mut movable: Vec<Vec3> = movable_atoms.iter().map(|&a| p(a)).collect();

    // A) rotate the atoms before the hinge 180 degrees about hinge->pivot
    let normal = (pivot_atom - hinge_atom).normalize();
    for m in movable.iter_mut().take(hinge_index) {
        *m = rotate_deg_axis_dir(hinge_atom, normal, *m, 180.0);
    }
    // B) hinge back into the original plane
    let c_to_o = first - hinge_atom;
    let n_to_o = second - hinge_atom;
    let old_normal = c_to_o.cross(n_to_o).normalize();
    let new_n_to_o = movable[second_dock] - movable[hinge_index];
    let new_c_to_o = movable[first_dock] - movable[hinge_index];
    let new_normal = new_n_to_o.cross(new_c_to_o).normalize();
    let hinge = old_normal.cross(new_normal);
    if hinge.length() > 0.0 {
        let hinge = hinge / hinge.length();
        let degrees = acos_deg(old_normal.dot(new_normal))?;
        for m in movable.iter_mut().take(hinge_index) {
            *m = rotate_deg_axis_dir(hinge_atom, hinge, *m, -degrees);
        }
    }
    // C1) rotate everything about the alpha carbon to put the first docked atom on CA->old second
    let a_to_new_o = movable[first_dock] - ca;
    let a_to_old_n = second - ca;
    let normal = a_to_new_o.cross(a_to_old_n);
    if normal.length() > 0.0 {
        let normal = normal / normal.length();
        let degrees = acos_deg((a_to_old_n / a_to_old_n.length()).dot(a_to_new_o / a_to_new_o.length()))?;
        for m in movable.iter_mut() {
            *m = rotate_deg_axis_dir(ca, normal, *m, degrees);
        }
    }
    // C2) twist about CA->old second to restore the plane
    let degrees = match dihedral_deg(movable[second_dock], ca, second, first) {
        Some(d) => d,
        None => {
            if clamp_acos {
                0.0
            } else {
                return Err("bad operand type for unary -: 'NoneType'".into());
            }
        }
    };
    let hinge = ca - second;
    if hinge.length() > 0.0 {
        let hinge = hinge / hinge.length();
        for m in movable.iter_mut() {
            *m = rotate_deg_axis_dir(ca, hinge, *m, -degrees);
        }
    }
    Ok(movable)
}

/// `MoverAmideFlip`.
pub fn amide_flip(w: &mut World, nh2: u32, ca_name: &str, non_flip_preference: f64, clamp: bool) -> MResult<Mover> {
    if w.elem(nh2) != "N" {
        return Err("MoverAmideFlip(): nh2Atom is not a Nitrogen".into());
    }
    let partners = w.bonded[nh2 as usize].clone();
    if partners.len() != 3 {
        return Err("MoverAmideFlip(): nh2Atom does not have three bonded neighbors".into());
    }
    let mut hs = Vec::new();
    let mut hinge = None;
    for &a in &partners {
        if w.is_h(a) {
            hs.push(a);
        } else {
            hinge = Some(a);
        }
    }
    if hs.len() != 2 {
        return Err("MoverAmideFlip(): nh2Atom does not have two bonded hydrogens".into());
    }
    let hinge = hinge.ok_or("MoverAmideFlip(): nh2Atom does not have bonded (hinge) Carbon friend")?;
    let mut oxygen = None;
    let mut pivot = None;
    for &b in &w.bonded[hinge as usize] {
        if w.elem(b) == "O" {
            oxygen = Some(b);
        } else if w.elem(b) == "C" {
            pivot = Some(b);
        }
    }
    let pivot = pivot.ok_or("MoverAmideFlip(): Hinge does not have bonded (pivot) Carbon friend")?;
    let oxygen = oxygen.ok_or("MoverAmideFlip(): Hinge does not have bonded oxygen friend")?;
    if w.bonded[oxygen as usize].len() != 1 {
        return Err("MoverAmideFlip(): Oxygen has more than one bonded neighbor".into());
    }
    let mut atoms = vec![hs[0], hs[1], nh2, oxygen, hinge, pivot];

    // linkers from the pivot to the alpha carbon
    let mut linkers = Vec::new();
    let mut linker_h = Vec::new();
    let mut prev = hinge;
    let mut cur = pivot;
    let ca_atom;
    loop {
        let bonded = w.bonded[cur as usize].clone();
        if bonded.len() != 4 {
            return Err("MoverAmideFlip(): Linker chain has an element with other than four bonds".into());
        }
        let mut link = None;
        for &b in &bonded {
            if b == prev {
                continue;
            } else if w.elem(b) == "C" {
                link = Some(b);
            } else if w.is_h(b) {
                linker_h.push(b);
            }
        }
        let link = link.ok_or("MoverAmideFlip(): Did not find Carbon in linker chain step")?;
        if w.name(link) == ca_name.trim().to_ascii_uppercase() {
            ca_atom = link;
            break;
        }
        linkers.push(link);
        prev = cur;
        cur = link;
    }
    if linker_h.len() != 2 * (linkers.len() + 1) {
        return Err(format!(
            "MoverAmideFlip(): Linker carbons do not have two Hydrogens each: {},{}",
            linkers.len(),
            linker_h.len()
        ));
    }
    atoms.extend(&linkers);
    atoms.extend(&linker_h);

    let pn = w.pos[nh2 as usize];
    let ph = w.pos[hinge as usize];
    let po = w.pos[oxygen as usize];
    // Place the H in the N-C-O plane at +/-120 degrees from C->N.
    let c_to_n = pn - ph;
    let c_to_o = po - ph;
    let normal = c_to_n.cross(c_to_o);
    let normal = normal / normal.length();
    let h0len = w.pos[hs[0] as usize].dist(pn);
    let h1len = w.pos[hs[1] as usize].dist(pn);
    let u = -(c_to_n / c_to_n.length());
    let h0 = pn + rotate_around_origin_scitbx(u * h0len, normal, 120.0 * std::f64::consts::PI / 180.0);
    let h1 = pn + rotate_around_origin_scitbx(u * h1len, normal, -120.0 * std::f64::consts::PI / 180.0);
    w.set_pos(hs[0], h0);
    w.set_pos(hs[1], h1);
    // flipped H positions around the oxygen
    let c_to_o = po - ph;
    let c_to_n = pn - ph;
    let normal = c_to_o.cross(c_to_n);
    let normal = normal / normal.length();
    let h0len = w.pos[hs[0] as usize].dist(pn);
    let h1len = w.pos[hs[1] as usize].dist(pn);
    let u = -(c_to_o / c_to_o.length());
    let new_h0 = po + rotate_around_origin_scitbx(u * h0len, normal, 120.0 * std::f64::consts::PI / 180.0);
    let new_h1 = po + rotate_around_origin_scitbx(u * h1len, normal, -120.0 * std::f64::consts::PI / 180.0);

    let start: Vec<Vec3> = atoms.iter().map(|&a| w.pos[a as usize]).collect();
    let mut newp = start.clone();
    newp[0] = new_h0;
    newp[1] = new_h1;
    newp[2] = po;
    newp[3] = pn;
    let coarse_pos = vec![start[..5].to_vec(), newp[..5].to_vec()];
    let movable = rotate_hinge_dock(w, &atoms, 4, 3, 2, ca_atom, clamp)?;
    Ok(Mover {
        kind: MoverKind::AmideFlip,
        atoms,
        n_moved: 5,
        coarse_pos,
        coarse_info: vec![vec![], vec![]],
        coarse_del: vec![vec![], vec![]],
        coarse_pref: vec![0.0, -non_flip_preference],
        fine_pos: vec![vec![], vec![]],
        fine_pref: vec![vec![], vec![]],
        fixup_pos: vec![vec![], movable],
        fixup_info: vec![vec![], vec![]],
        fixup_del: vec![vec![], vec![]],
        rot: None,
    })
}

/// `MoverHisFlip`.
pub fn his_flip(
    w: &mut World,
    ne2: u32,
    non_flip_preference: f64,
    enabled: u8,
    enable_fixup: bool,
    clamp: bool,
) -> MResult<Mover> {
    if w.elem(ne2) != "N" {
        return Err("MoverHisFlip(): ne2Atom is not a Nitrogen".into());
    }
    let partners = w.bonded[ne2 as usize].clone();
    if partners.len() < 3 {
        return Err("MoverHisFlip(): NE2 does not have three bonded neighbors".into());
    }
    let mut hyd = Vec::new();
    let mut carb = Vec::new();
    for &a in &partners {
        if w.is_h(a) {
            hyd.push(a);
        } else if w.elem(a) == "C" {
            carb.push(a);
        }
    }
    if hyd.len() != 1 {
        return Err("MoverHisFlip(): NE2 does not have one bonded hydrogen (probably ionically bound)".into());
    }
    if carb.len() != 2 {
        return Err("MoverHisFlip(): NE2 does not have two bonded carbons".into());
    }
    let ne2h = hyd[0];
    let ctest = carb[0];
    if w.bonded[ctest as usize].len() != 3 {
        return Err("MoverHisFlip(): NE2 neighbor does not have three bonded neighbors".into());
    }
    let cs = w.bonded[ctest as usize].iter().filter(|&&b| w.elem(b) == "C").count();
    let (ce1, cd2) = if cs == 0 { (carb[0], carb[1]) } else { (carb[1], carb[0]) };
    if w.bonded[ce1 as usize].len() != 3 {
        return Err("MoverHisFlip(): CE1 does not have three bonded neighbors".into());
    }
    let ce1h = w.bonded[ce1 as usize].iter().copied().filter(|&b| w.is_h(b)).last();
    let ce1h = ce1h.ok_or("MoverHisFlip(): Could not find Hydrogen attached to CE1")?;
    if w.bonded[cd2 as usize].len() != 3 {
        return Err("MoverHisFlip(): CD2 does not have three bonded neighbors".into());
    }
    let cd2h = w.bonded[cd2 as usize].iter().copied().filter(|&b| w.is_h(b)).last();
    let cd2h = cd2h.ok_or("MoverHisFlip(): CD2 does not have a bonded hydrogen")?;
    let cg = w.bonded[cd2 as usize].iter().copied().filter(|&b| b != cd2 && w.elem(b) == "C").last();
    let cg = cg.ok_or("MoverHisFlip(): Could not find CG")?;
    if w.bonded[cg as usize].len() != 3 {
        return Err("MoverHisFlip(): CG does not have three bonded neighbors".into());
    }
    let mut cb = None;
    let mut nd1 = None;
    for &b in &w.bonded[cg as usize] {
        if b == cd2 {
            continue;
        } else if w.elem(b) == "N" {
            nd1 = Some(b);
        } else if w.elem(b) == "C" {
            cb = Some(b);
        }
    }
    let nd1 = nd1.ok_or("MoverHisFlip(): Could not find ND1")?;
    let cb = cb.ok_or("MoverHisFlip(): Could not find CB")?;
    let p2 = &w.bonded[nd1 as usize];
    if p2.len() < 3 {
        return Err("MoverHisFlip(): ND1 does not have three bonded neighbors".into());
    }
    let hyd2: Vec<u32> = p2.iter().copied().filter(|&a| w.is_h(a)).collect();
    let carb2 = p2.iter().filter(|&&a| !w.is_h(a) && w.elem(a) == "C").count();
    if hyd2.len() != 1 {
        return Err("MoverHisFlip(): ND1 does not have one bonded hydrogen (probably ionically bound)".into());
    }
    if carb2 != 2 {
        return Err("MoverHisFlip(): ND1 does not have two bonded carbons".into());
    }
    let nd1h = hyd2[0];
    let cbb = &w.bonded[cb as usize];
    if cbb.len() != 4 {
        return Err(format!("MoverHisFlip(): CB does not have four bonded neighbors, has {}", cbb.len()));
    }
    let mut ca = None;
    let mut cbh = Vec::new();
    for &b in cbb {
        if b == cg {
            continue;
        } else if w.elem(b) == "C" {
            ca = Some(b);
        } else if w.is_h(b) {
            cbh.push(b);
        }
    }
    let ca = ca.ok_or("MoverHisFlip(): Could not find CA")?;
    if cbh.len() != 2 {
        return Err("MoverHisFlip(): Could not find Hydrogens on CB".into());
    }
    let atoms = vec![ne2, ne2h, ce1, ce1h, nd1, nd1h, cd2, cd2h, cg, cb, cbh[0], cbh[1]];
    let p = |a: u32| w.pos[a as usize];
    let nd1hv = p(nd1h) - p(nd1);
    let ne2hv = p(ne2h) - p(ne2);
    let ce1hv = p(ce1h) - p(ce1);
    let cd2hv = p(cd2h) - p(cd2);
    let unit = |v: Vec3| v / v.length();
    let nd1h_new = p(cd2) + unit(cd2hv) * nd1hv.length();
    let cd2h_new = p(nd1) + unit(nd1hv) * cd2hv.length();
    let ce1h_new = p(ne2) + unit(ne2hv) * ce1hv.length();
    let ne2h_new = p(ce1) + unit(ce1hv) * ne2hv.length();
    let start: Vec<Vec3> = atoms.iter().map(|&a| p(a)).collect();
    let mut newp = start.clone();
    newp[0] = p(ce1);
    newp[1] = ne2h_new;
    newp[2] = p(ne2);
    newp[3] = ce1h_new;
    newp[4] = p(cd2);
    newp[5] = nd1h_new;
    newp[6] = p(nd1);
    newp[7] = cd2h_new;
    let mut coarse_pos = Vec::new();
    if enabled & 1 != 0 {
        for _ in 0..4 {
            coarse_pos.push(start[..9].to_vec());
        }
    }
    if enabled & 2 != 0 {
        for _ in 0..4 {
            coarse_pos.push(newp[..9].to_vec());
        }
    }
    let fixed = rotate_hinge_dock(w, &atoms, 8, 0, 2, ca, clamp)?;
    let mut fixup_pos = Vec::new();
    if enabled & 1 != 0 {
        for _ in 0..4 {
            fixup_pos.push(vec![]);
        }
    }
    if enabled & 2 != 0 {
        for _ in 0..4 {
            fixup_pos.push(if enable_fixup { fixed.clone() } else { vec![] });
        }
    }
    let n = coarse_pos.len();
    let mut infos = Vec::with_capacity(n);
    let mut dels = Vec::with_capacity(n);
    for i in 0..n {
        let mut ex: Vec<AtomInfo> = atoms.iter().map(|&a| w.info[a as usize]).collect();
        let mut de = vec![false; atoms.len()];
        if i % 4 == 1 || i % 4 == 3 {
            ex[0].is_acceptor = true;
            de[1] = true;
        }
        if i % 4 == 2 || i % 4 == 3 {
            ex[4].is_acceptor = true;
            de[5] = true;
        }
        infos.push(ex);
        dels.push(de);
    }
    let mut pref = Vec::new();
    if enabled & 1 != 0 {
        pref.extend([0.0 - 0.05, 0.0, 0.0, 0.0 - 1.0]);
    }
    if enabled & 2 != 0 {
        pref.extend([-non_flip_preference - 0.05, -non_flip_preference, -non_flip_preference, -non_flip_preference - 1.0]);
    }
    Ok(Mover {
        kind: MoverKind::HisFlip { enabled, enable_fixup },
        atoms,
        n_moved: 9,
        coarse_pos,
        coarse_info: infos.clone(),
        coarse_del: dels.clone(),
        coarse_pref: pref,
        fine_pos: vec![vec![]; n],
        fine_pref: vec![vec![]; n],
        fixup_pos,
        fixup_info: infos,
        fixup_del: dels,
        rot: None,
    })
}
