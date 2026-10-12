#[allow(dead_code)]
mod common {
    pub mod validation_fixture;
}

use reduce3::atominfo::{extra_atom_info, extra_atom_info_validation};
use reduce3::interp::{EType, FlatAtoms};
use reduce3::monlib::{EnerAtom, MonLib};

#[test]
fn validation_neutron_radii_keep_polar_override_and_legacy_typing() {
    let fixture = common::validation_fixture::library();
    let mut library = MonLib::load_validation(&fixture.0).unwrap();
    for (name, electron, neutron, element) in [
        ("HNH1", 1.05, Some(1.0), "H"),
        ("HCH", 1.22, Some(1.17), "H"),
        ("HAR", 1.05, Some(1.0), "H"),
        ("NH1", 1.55, None, "N"),
    ] {
        library.ener.insert(
            name.into(),
            EnerAtom {
                type_: name.into(),
                hb_type: "N".into(),
                vdw_radius: Some(electron),
                vdw_radius_neutron: neutron,
                element: element.into(),
                ..Default::default()
            },
        );
    }
    let structure = reduce3::pdbio::read_pdb(concat!(
        "ATOM      1  H   ALA A   1       0.000   0.000   0.000  1.00  0.00           H\n",
        "ATOM      2  HB1 ALA A   1       1.000   0.000   0.000  1.00  0.00           H\n",
        "ATOM      3  HE1 TYR A   2       2.000   0.000   0.000  1.00  0.00           H\n",
        "ATOM      4  N   ALA A   3       3.000   0.000   0.000  1.00  0.00           N\n",
    ));
    let flat = FlatAtoms::from_structure(&structure);
    let types: Vec<_> = ["HNH1", "HCH", "HAR", "NH1"]
        .into_iter()
        .map(|name| EType::Typed(name.into()))
        .collect();
    let bonded = vec![vec![3], vec![], vec![], vec![0]];
    let radii = |info: reduce3::atominfo::ExtraInfo| {
        assert!(!info.warnings.contains("Warning:"));
        info.info
            .into_iter()
            .map(|atom| atom.vdw_radius)
            .collect::<Vec<_>>()
    };
    assert_eq!(
        radii(extra_atom_info_validation(
            &structure, &flat, &types, &library, &bonded, true
        )),
        [1.05, 1.17, 1.0, 1.55]
    );
    assert_eq!(
        radii(extra_atom_info_validation(
            &structure, &flat, &types, &library, &bonded, false
        )),
        [1.05, 1.22, 1.05, 1.55]
    );
    assert_eq!(
        radii(extra_atom_info(
            &structure, &flat, &types, &library, &bonded, true
        )),
        [1.05, 1.22, 1.05, 1.55]
    );
}
