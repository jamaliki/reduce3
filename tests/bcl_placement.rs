use reduce3::{hplace, mmcif, MonLib};
use reduce3::model::Structure;

fn library() -> Option<MonLib> {
    let Some(root) = MonLib::locate(None) else {
        eprintln!("skipped: chem_data not found (set REDUCE3_CHEM_DATA)");
        return None;
    };
    Some(MonLib::load(&root).unwrap())
}

fn position(st: &Structure, name: &str) -> reduce3::geom::Vec3 {
    st.atom_paths().into_iter().find(|&p| st.atom(p).name.trim() == name).map(|p| st.atom(p).xyz).unwrap()
}

#[test]
fn complete_indexed_bcl_keeps_heavy_coordinates_and_places_source_planar_h() {
    let Some(library) = library() else { return };
    let mut st = mmcif::read_mmcif(include_str!("data/bcl_7vzg_902.cif")).unwrap();
    let original: Vec<_> = st.atom_paths().into_iter().map(|p| (st.atom(p).name.clone(), st.atom(p).xyz)).collect();
    assert_eq!(original.len(), 66);
    let placed = hplace::place_hydrogens(&mut st, &library, &hplace::HPlaceParams {
        neutron: false, n_terminal_charge: hplace::NTermCharge::ResidueOne,
        exclude_water: true, keep_existing_h: false, adp_scale: 1.0, compat: false, cell: None,
    });
    assert!(placed.site_labels_no_para.is_empty(), "{:?}", placed.site_labels_no_para);
    for (name, xyz) in original { assert_eq!(position(&st, name.trim()), xyz); }
    for (h, parent, first, second) in [
        ("H2", "C2", "C1", "C3"), ("HHB", "CHB", "C4A", "C1B"),
        ("HHC", "CHC", "C4B", "C1C"), ("HHD", "CHD", "C4C", "C1D"),
    ] {
        let (hydrogen, center, a, b) = (position(&st, h), position(&st, parent), position(&st, first), position(&st, second));
        assert!((hydrogen.dist(center) - 1.09).abs() < 1e-12);
        let normal = (a - center).cross(b - center).normalize();
        assert!((hydrogen - center).dot(normal).abs() < 1e-12);
    }
}

#[test]
fn compat_keeps_the_existing_bcl_restraint_failure() {
    let Some(library) = library() else { return };
    let mut st = mmcif::read_mmcif(include_str!("data/bcl_7vzg_902.cif")).unwrap();
    let placed = hplace::place_hydrogens(&mut st, &library, &hplace::HPlaceParams {
        neutron: false, n_terminal_charge: hplace::NTermCharge::ResidueOne,
        exclude_water: true, keep_existing_h: false, adp_scale: 1.0, compat: true, cell: None,
    });
    for name in ["H2", "HHB", "HHC", "HHD"] {
        assert!(placed.site_labels_no_para.iter().any(|label| label.trim_start().starts_with(name)));
    }
}
