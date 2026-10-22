use super::*;

fn fixture() -> (Interp, Vec<Vec3>) {
    let mut it = Interp::empty(5);
    it.add_bond(0, 1, 1.4, 0);
    it.add_bond(0, 2, 1.4, 0);
    it.add_bond(0, 4, 0.86, 0);
    for (i, value) in [(1, 120.13), (2, 120.09)] {
        it.push_angle_unchecked(4, 0, i, value, 3.0);
    }
    it.push_angle_unchecked(1, 0, 2, 126.93, 3.0);
    it.planes.push(vec![3, 4, 0, 1, 2]);
    let angle = 126.93_f64.to_radians();
    (it, vec![Vec3::ZERO, v3(1.4, 0., 0.), v3(1.4 * angle.cos(), 1.4 * angle.sin(), 0.),
        v3(0., 0., 1.), v3(0.5, 0.5, 0.5)])
}

fn parameterize(it: &Interp, sites: &mut [Vec3], fixed: bool) -> RidingResult {
    let names = ["N3", "C4", "C2", "H5", "HN3"].map(ShortStr::from);
    let alts = vec![ShortStr::new(); 5];
    riding(it, sites, &RidingAtoms {
        is_h: &[false, false, false, true, true], altloc: &alts,
        name: &names, occ: &[1.; 5], rg_of: &[0; 5], expected_heavy: &|_| None,
        dictionary_nh2_torsion: fixed, cis_partner: &|_| None, reroute_axial_reference: fixed, use_plane_restraints: fixed,
    }, true)
}

#[test]
fn a_later_h_in_an_explicit_plane_keeps_its_plane_and_bond_length() {
    let (it, mut sites) = fixture();
    let before = sites[..3].to_vec();
    let result = parameterize(&it, &mut sites, true);
    assert_eq!(result.coef[4].unwrap().htype, HType::Flat2);
    assert!((sites[4].dist(sites[0]) - 0.86).abs() < 1e-12);
    assert!(sites[4].z.abs() < 1e-12);
    assert_eq!(sites[..3], before);
}

#[test]
fn invalid_angles_without_explicit_plane_and_compat_still_fail() {
    let (mut it, mut sites) = fixture();
    assert!(parameterize(&it, &mut sites, false).coef[4].is_none());
    it.planes.clear();
    assert!(parameterize(&it, &mut sites, true).coef[4].is_none());
}

#[test]
fn explicit_plane_does_not_allow_an_extra_hydrogen() {
    let (mut it, mut sites) = fixture();
    it.add_bond(0, 3, 0.86, 0);
    it.push_angle_unchecked(4, 0, 3, 120., 3.);
    assert!(parameterize(&it, &mut sites, true).coef[4].is_none());
}
