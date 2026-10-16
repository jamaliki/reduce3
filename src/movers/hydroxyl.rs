//! Local geometry shared by conjugated O-H search and output validation.

use crate::geom::{Vec3, dihedral_deg, rotate_deg_axis_dir};
use crate::world::World;
use super::RotatorMeta;
use smallvec::SmallVec;

/// Fixed-mode geometry guard in A, not a fitted energy or an experimental limit.
/// It admits the checked 8E1W Tyr164 D position (0.089-0.100 A).
pub(crate) const MAX_HEIGHT: f64 = 0.12;

/// Distinguish a planar three-coordinate partner from a pyramidal atom.
const PLANAR_PARTNER_TOLERANCE: f64 = 0.2;

pub(super) struct PlanarHydroxyl {
    pub planes: [u32; 2],
    pub carbonyl: SmallVec<[u32; 2]>,
    oxygen: u32,
    partner: u32,
}

pub(super) fn planar_hydroxyl(w: &World, oxygen: u32, partner: u32, friends: &[u32]) -> Option<PlanarHydroxyl> {
    if w.elem(oxygen) != "O" || friends.len() != 2 {
        return None;
    }
    let o = w.pos[oxygen as usize];
    let normal = (w.pos[friends[0] as usize] - o).cross(w.pos[friends[1] as usize] - o);
    if normal.length() < 1e-6
        || ((w.pos[partner as usize] - o).dot(normal) / normal.length()).abs() > PLANAR_PARTNER_TOLERANCE {
        return None;
    }
    let carbonyl = friends.iter().copied()
        .filter(|&f| w.elem(f) == "O" && w.bonded[f as usize].len() == 1).collect();
    Some(PlanarHydroxyl { planes: [friends[0], friends[1]], carbonyl, oxygen, partner })
}

impl PlanarHydroxyl {
    fn height(&self, w: &World, hydrogen: Vec3) -> f64 {
        if !hydrogen.as_array().iter().all(|x| x.is_finite()) { return f64::INFINITY; }
        let o = w.pos[self.oxygen as usize];
        let axis = w.pos[self.partner as usize] - o;
        let mut height: f64 = 0.0;
        for &f in &self.planes {
            let normal = axis.cross(w.pos[f as usize] - o);
            let length = normal.length();
            if !length.is_finite() || length < 1e-6 { return f64::INFINITY; }
            height = height.max((hydrogen - o).dot(normal / length).abs());
        }
        height
    }
}

/// Keep the normal fine-angle grid. Add legal band endpoints so contact scoring
/// can use the allowed tilt without adding a dense coarse search.
pub(super) fn constrain_rotator(w: &World, h: u32, plane: &PlanarHydroxyl, rot: &mut RotatorMeta) -> Result<(), String> {
    let position = |angle| rotate_deg_axis_dir(rot.axis_origin, rot.axis_dir, w.pos[h as usize], angle);
    let legal = |angle: f64, fine: &[f64]| {
        std::iter::once(0.0).chain(fine.iter().copied())
            .all(|offset| plane.height(w, position(angle + offset)) <= MAX_HEIGHT + 1e-9)
    };
    let reference = plane.carbonyl.first().copied().unwrap_or(plane.planes[0]);
    let center = dihedral_deg(w.pos[h as usize], w.pos[plane.partner as usize],
        w.pos[plane.oxygen as usize], w.pos[reference as usize])
        .ok_or("Cannot define the conjugated O-H plane")?;
    let mut angles = rot.coarse_angles.clone();
    for base in [center, center + 180.0] {
        let band: Vec<_> = (-90..=90).map(|offset| base + offset as f64)
            .filter(|&angle| legal(angle, &rot.fine_angles)).collect();
        angles.extend(band.first().copied());
        angles.extend(band.last().copied());
    }
    angles.retain(|&angle| legal(angle, &rot.fine_angles));
    if angles.is_empty() {
        // A distorted heavy frame can leave less room than the fine grid spans.
        // Search its legal one-degree states without a fine step.
        rot.fine_angles.clear();
        angles = (-180..180).map(|offset| center + offset as f64)
            .filter(|&angle| legal(angle, &[])).collect();
    }
    if angles.is_empty() {
        return Err(format!("No conjugated O-H orientation fits the {MAX_HEIGHT:.2} A local-plane limit"));
    }
    for angle in &mut angles {
        *angle = (*angle + 180.0).rem_euclid(360.0) - 180.0;
    }
    angles.sort_by(f64::total_cmp);
    angles.dedup_by(|a, b| (*a - *b).abs() < 1e-8);
    rot.coarse_angles = angles;
    Ok(())
}

/// Check the final coordinates without relying on a Mover's scores or states.
pub(crate) fn validate(w: &World, hydrogens: impl IntoIterator<Item = u32>) -> Result<usize, String> {
    let mut checked = 0;
    for h in hydrogens {
        let [oxygen] = w.bonded[h as usize].as_slice() else { continue };
        let neighbors = &w.bonded[*oxygen as usize];
        if neighbors.len() != 2 { continue; }
        let Some(&partner) = neighbors.iter().find(|&&a| a != h) else { continue };
        let friends: Vec<_> = w.bonded[partner as usize].iter().copied().filter(|&a| a != *oxygen).collect();
        let Some(plane) = planar_hydroxyl(w, *oxygen, partner, &friends) else { continue };
        let height = plane.height(w, w.pos[h as usize]);
        if !height.is_finite() || height > MAX_HEIGHT + 1e-8 {
            return Err(format!("Conjugated O-H geometry failed for {} {}: local-plane distance {height:.4} A exceeds {MAX_HEIGHT:.2} A",
                w.res_name_and_id(h), w.name(h)));
        }
        checked += 1;
    }
    Ok(checked)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::v3;
    use crate::movers::{RidingRef, SingleHOptions, single_hydrogen_rotator};
    use crate::probe::AtomInfo;
    use crate::world::AtomLabels;

    fn world() -> World {
        let tilt = 42_f64.to_radians();
        let pos = vec![v3(-0.287, 0.789 * tilt.cos(), 0.789 * tilt.sin()),
            v3(0., 0., 0.), v3(1.36, 0., 0.), v3(2.05, 1.19, 0.),
            v3(2.05, -1.19, 0.), v3(-1., 2., 2.)];
        let labels = [("HH", "H"), ("OH", "O"), ("C", "C"), ("C1", "C"), ("C2", "C"), ("OA", "O")]
            .map(|(name, element)| AtomLabels { name: name.into(), element: element.into(),
                resname: "PHN".into(), chain: "A".into(), ..Default::default() }).to_vec();
        World { labels, pos, occ: vec![1.; 6], b: vec![10.; 6],
            info: vec![AtomInfo::default(); 6],
            bonded: vec![vec![1], vec![0, 2], vec![1, 3, 4], vec![2], vec![2], vec![]],
            riding: vec![Some(RidingRef { n: 0, a2: 3 }), None, None, None, None, None], n_real: 6 }
    }

    fn options(constrain: bool) -> SingleHOptions {
        SingleHOptions { circular_angle_spacing: true, any_partner_valence: true,
            constrain_planar_hydroxyls: constrain, planar_preference: 0., syn_preference: 0. }
    }

    #[test]
    fn every_coarse_and_fine_pose_respects_both_local_planes_without_preferences() {
        let mut w = world();
        let heavy = w.pos[1..].to_vec();
        let mover = single_hydrogen_rotator(&mut w, 0, &[5], &[], &options(true)).unwrap();
        let plane = planar_hydroxyl(&w, 1, 2, &[3, 4]).unwrap();
        assert!(mover.num_coarse() <= 12);
        assert!(!mover.fine_pos[0].is_empty());
        for pose in mover.coarse_pos.iter().chain(mover.fine_pos.iter().flatten()) {
            assert!(plane.height(&w, pose[0]) <= MAX_HEIGHT + 1e-9);
            assert!(pose[1].dist(heavy[0]) < 1e-12);
        }
        assert_eq!(w.pos[1..], heavy);
        assert!(mover.coarse_pref.iter().all(|&p| p == 0.));
    }

    #[test]
    fn legacy_options_leave_nonplanar_states_available() {
        let mut w = world();
        let mover = single_hydrogen_rotator(&mut w, 0, &[5], &[], &options(false)).unwrap();
        let plane = planar_hydroxyl(&w, 1, 2, &[3, 4]).unwrap();
        assert!(mover.coarse_pos.iter().any(|pose| plane.height(&w, pose[0]) > 0.4));
    }

    #[test]
    fn tetrahedral_alcohols_and_thiols_are_not_constrained() {
        let mut w = world();
        w.bonded[2].push(5);
        assert_eq!(validate(&w, [0]).unwrap(), 0);
        let mover = single_hydrogen_rotator(&mut w, 0, &[], &[], &options(true)).unwrap();
        assert!(mover.coarse_pos.iter().any(|pose| pose[0].z.abs() > 0.4));
        let mut thiol = world();
        thiol.labels[1].element = "S".into();
        assert_eq!(validate(&thiol, [0]).unwrap(), 0);
    }

    #[test]
    fn final_validation_rejects_bad_coordinates_without_a_mover() {
        let mut w = world();
        let error = validate(&w, [0]).unwrap_err();
        assert!(error.contains("PHN") && error.contains("HH") && error.contains("0.12 A"));
        w.pos[0] = v3(-0.287, 0.789, 0.);
        assert_eq!(validate(&w, [0]).unwrap(), 1);
        w.pos[0].z = f64::NAN;
        assert!(validate(&w, [0]).is_err());
    }

    #[test]
    fn both_neighbor_planes_are_checked() {
        let w = world();
        let forward = planar_hydroxyl(&w, 1, 2, &[3, 4]).unwrap();
        let reversed = planar_hydroxyl(&w, 1, 2, &[4, 3]).unwrap();
        assert_eq!(forward.height(&w, w.pos[0]), reversed.height(&w, w.pos[0]));
    }

    #[test]
    fn narrow_legal_bands_do_not_use_invalid_fine_steps() {
        let mut w = world();
        w.pos[3].z = 0.15;
        w.pos[4].z = 0.15;
        let mover = single_hydrogen_rotator(&mut w, 0, &[5], &[], &options(true)).unwrap();
        assert!(mover.fine_pos.iter().all(Vec::is_empty));
        let plane = planar_hydroxyl(&w, 1, 2, &[3, 4]).unwrap();
        assert!(mover.coarse_pos.iter().all(|pose| plane.height(&w, pose[0]) <= MAX_HEIGHT + 1e-9));
    }

    #[test]
    fn incompatible_heavy_planes_fail_instead_of_waiving_the_limit() {
        let mut w = world();
        w.pos[3].z = 0.3;
        w.pos[4].z = 0.3;
        let error = single_hydrogen_rotator(&mut w, 0, &[5], &[], &options(true)).unwrap_err();
        assert!(error.contains("No conjugated O-H orientation"));
        assert!(validate(&w, [0]).is_err());
    }

    #[test]
    fn experimental_8e1w_tyr164_deuterium_remains_within_the_guard() {
        let mut w = world();
        // Deposited DH, OH, CZ, CE1 and CE2 coordinates from PDB 8E1W.
        w.pos[..5].copy_from_slice(&[v3(-9.078, -5.290, 19.196), v3(-8.828, -5.110, 18.396),
            v3(-8.638, -6.289, 17.697), v3(-8.363, -6.241, 16.363), v3(-8.708, -7.473, 18.328)]);
        w.labels[0].element = "D".into();
        let plane = planar_hydroxyl(&w, 1, 2, &[3, 4]).unwrap();
        let height = plane.height(&w, w.pos[0]);
        assert!((height - 0.0996122786).abs() < 1e-9);
        assert_eq!(validate(&w, [0]).unwrap(), 1);
    }
}
