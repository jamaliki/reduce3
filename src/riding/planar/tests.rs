use super::*;
use crate::monlib::{Comp, CompAtom, CompBond};

fn fixture(orders: [&str; 2]) -> (Interp, Vec<Vec3>) {
    let mut it = Interp::empty(4);
    let names = ["CHB", "C4A", "C1B", "HHB"];
    let mut comp = Comp { id: "BCL".into(), ..Default::default() };
    for (i, name) in names.iter().enumerate() {
        comp.atoms.push(CompAtom { id: (*name).into(), type_symbol: if i == 3 { "H" } else { "C" }.into(), type_energy: None });
    }
    for (i, order) in [orders[0], orders[1], "single"].iter().enumerate() {
        comp.bonds.push(CompBond { a1: "CHB".into(), a2: names[i + 1].into(), type_: (*order).into(), ..Default::default() });
        it.add_bond(0, (i + 1) as u32, if i == 2 { 1.09 } else { 1.4 }, 0);
    }
    let comp = Arc::new(comp);
    for (i, name) in names.iter().enumerate() {
        it.dict[i] = crate::interp::AtomDictInfo { comp: Some(comp.clone()), dict_id: Some((*name).into()), instance: 0 };
    }
    it.push_angle_unchecked(3, 0, 1, 120., 3.);
    it.push_angle_unchecked(3, 0, 2, 120., 3.);
    it.push_angle_unchecked(1, 0, 2, 131.496, 2.737);
    let angle = 131.496_f64.to_radians();
    (it, vec![Vec3::ZERO, v3(1.4, 0., 0.), v3(1.4 * angle.cos(), 1.4 * angle.sin(), 0.), v3(0.5, 0.5, 0.5)])
}

fn parameterize(it: &Interp, sites: &mut [Vec3], fixed: bool) -> RidingResult {
    let count = sites.len();
    let names: Vec<_> = ["CHB", "C4A", "C1B", "HHB", "EXTRA"].iter().take(count).map(|name| (*name).into()).collect();
    let mut is_h = vec![false; count];
    is_h[3] = true;
    riding(it, sites, &RidingAtoms {
        is_h: &is_h, altloc: &vec![ShortStr::new(); count],
        name: &names, occ: &vec![1.; count], rg_of: &vec![0; count], expected_heavy: &|_| None,
        dictionary_nh2_torsion: fixed, reroute_axial_reference: fixed, use_plane_restraints: fixed,
    }, true)
}

#[test]
fn complete_source_sp2_carbon_recovers_only_in_fixed_mode() {
    for orders in [["single", "double"], ["deloc", "deloc"], ["aromatic", "aromatic"]] {
        let (it, mut sites) = fixture(orders);
        let original = sites[..3].to_vec();
        assert!(parameterize(&it, &mut sites, false).coef[3].is_none());
        assert_eq!(parameterize(&it, &mut sites, true).coef[3].unwrap().htype, HType::Flat2);
        assert!(sites[3].z.abs() < 1e-12);
        assert!((sites[3].dist(sites[0]) - 1.09).abs() < 1e-12);
        assert_eq!(sites[..3], original);
    }
}

#[test]
fn tetrahedral_unknown_order_overvalent_or_missing_source_still_fails() {
    for orders in [["single", "single"], ["triple", "single"], ["deloc", "single"], ["unknown", "double"]] {
        let (it, mut sites) = fixture(orders);
        assert!(parameterize(&it, &mut sites, true).coef[3].is_none());
    }
    for invalid in 0..7 {
        let (mut it, mut sites) = fixture(["deloc", "deloc"]);
        match invalid {
            0 => it.dict[0].comp = None,
            1 => { Arc::make_mut(it.dict[0].comp.as_mut().unwrap()).atoms[0].type_symbol = "N".into(); }
            2 => it.dict[2].dict_id = None,
            3 => { Arc::make_mut(it.dict[0].comp.as_mut().unwrap()).bonds.push(CompBond {
                a1: "CHB".into(), a2: "EXTRA".into(), type_: "single".into(), ..Default::default() }); }
            4 => it.dict[3].comp = None,
            5 => it.bonds.retain(|b| b.i != 2 && b.j != 2),
            6 => {
                sites.push(v3(0., 0., 1.4));
                it.dict.push(crate::interp::AtomDictInfo::default());
                it.add_bond(0, 4, 1.4, 0);
            }
            _ => unreachable!(),
        }
        assert!(parameterize(&it, &mut sites, true).coef[3].is_none());
    }
}
