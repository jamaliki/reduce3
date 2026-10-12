use reduce3::MonLib;
use reduce3::interp::{
    FlatAtoms,
    validation::{ValidationOptions, build},
};
#[path = "common/validation_fixture.rs"]
mod fixture;
use fixture::{library, one_residue, optional_library, structure};

#[test]
fn oxygen_only_water_checks_supplied_atoms_without_requiring_hydrogens() {
    let fixture = library();
    optional_library(
        &fixture,
        "HOH",
        &[("O", "O"), ("H1", "H"), ("H2", "H")],
        "O H1 single .96 .02\nO H2 single .96 .02\n",
        "H1 O H2 104.5 3.0\n",
    );
    let ml = MonLib::load_validation(&fixture.0).unwrap();
    let st = one_residue("HOH", &[("O", "O")]);
    let graph = build(
        &st,
        &FlatAtoms::from_structure(&st),
        &ml,
        &ValidationOptions::default(),
    )
    .unwrap();
    assert!(graph.complete());
    assert!(graph.bonds.is_empty());
    assert!(graph.angles.is_empty());
    assert_eq!(graph.sites.len(), 1);
}

#[test]
fn neutral_histidine_may_omit_only_the_other_ring_proton() {
    let fixture = library();
    optional_library(
        &fixture,
        "HIS",
        &[("ND1", "N"), ("NE2", "N"), ("HD1", "H"), ("HE2", "H")],
        "ND1 HD1 single 1.0 .02\nNE2 HE2 single 1.0 .02\n",
        "HD1 ND1 NE2 120.0 3.0\nHE2 NE2 ND1 120.0 3.0\n",
    );
    let ml = MonLib::load_validation(&fixture.0).unwrap();
    for proton in ["HD1", "HE2"] {
        let st = one_residue("HIS", &[("ND1", "N"), ("NE2", "N"), (proton, "H")]);
        let graph = build(
            &st,
            &FlatAtoms::from_structure(&st),
            &ml,
            &ValidationOptions::default(),
        )
        .unwrap();
        assert!(graph.complete(), "{:?}", graph.coverage);
        assert_eq!(graph.bonds.len(), 1);
        assert_eq!(graph.angles.len(), 1);
    }
    let st = one_residue("HIS", &[("ND1", "N"), ("NE2", "N")]);
    let graph = build(
        &st,
        &FlatAtoms::from_structure(&st),
        &ml,
        &ValidationOptions::default(),
    )
    .unwrap();
    assert!(!graph.complete());
    assert_eq!(
        graph
            .coverage
            .iter()
            .find(|c| c.check == "hydrogens")
            .unwrap()
            .unresolved,
        2
    );
}

#[test]
fn validation_uses_effective_library_and_sequential_link_modifications() {
    let fixture = library();
    let native = MonLib::load_validation(&fixture.0).unwrap();
    let legacy = MonLib::load(&fixture.0).unwrap();
    let st = structure();
    let flat = FlatAtoms::from_structure(&st);
    let options = ValidationOptions {
        use_cdl: false,
        ..Default::default()
    };
    let current = build(&st, &flat, &native, &options).unwrap();
    let old = build(&st, &flat, &legacy, &options).unwrap();
    let first_co = current.bonds.iter().find(|b| b.indices == [2, 3]).unwrap();
    let second_nh = current.bonds.iter().find(|b| b.indices == [5, 9]).unwrap();
    let old_nh = old.bonds.iter().find(|b| b.indices == [5, 9]).unwrap();
    assert_eq!(first_co.ideal, 1.23);
    assert_eq!(first_co.weight, 2500.);
    assert_eq!(second_nh.ideal, 0.914);
    assert_eq!(second_nh.weight, 10000.);
    assert_eq!(old_nh.ideal, 0.860);
    assert_eq!(old_nh.weight, 2500.);
    assert_eq!(current.sites, flat.pos);
}

#[test]
fn nuclear_link_targets_keep_the_ordinary_esd() {
    let fixture = library();
    let file = fixture.0.join("mon_lib/list/mon_lib_list.cif");
    let text = std::fs::read_to_string(&file).unwrap().replace(
        "_chem_mod_bond.new_value_dist_esd\nchange N H .914 .010",
        "_chem_mod_bond.new_value_dist_esd\n_chem_mod_bond.new_value_dist_nucleus\n_chem_mod_bond.new_value_dist_nucleus_esd\nchange N H .914 .010 1.023 .123",
    );
    std::fs::write(file, text).unwrap();
    let ml = MonLib::load_validation(&fixture.0).unwrap();
    let st = structure();
    let graph = build(
        &st,
        &FlatAtoms::from_structure(&st),
        &ml,
        &ValidationOptions {
            nuclear: true,
            use_cdl: false,
            ..Default::default()
        },
    )
    .unwrap();
    let second_nh = graph.bonds.iter().find(|b| b.indices == [5, 9]).unwrap();
    assert_eq!(second_nh.ideal, 1.023);
    assert_eq!(second_nh.weight, 10000.);
}

#[test]
fn disulfide_targets_use_the_reference_ss_link() {
    let fixture = library();
    optional_library(
        &fixture,
        "CYS",
        &[("CB", "C"), ("SG", "S"), ("HG", "H")],
        "CB SG single 1.81 .02\nSG HG single 1.20 .02\n",
        "CB SG HG 109.0 3.0\n",
    );
    let file = fixture.0.join("mon_lib/list/mon_lib_list.cif");
    let mut text = std::fs::read_to_string(&file).unwrap();
    text.push_str(
        r#"
data_link_SS
loop_
_chem_link_bond.atom_1_comp_id
_chem_link_bond.atom_id_1
_chem_link_bond.atom_2_comp_id
_chem_link_bond.atom_id_2
_chem_link_bond.value_dist
_chem_link_bond.value_dist_esd
1 SG 2 SG 2.031 .020
loop_
_chem_link_angle.atom_1_comp_id
_chem_link_angle.atom_id_1
_chem_link_angle.atom_2_comp_id
_chem_link_angle.atom_id_2
_chem_link_angle.atom_3_comp_id
_chem_link_angle.atom_id_3
_chem_link_angle.value_angle
_chem_link_angle.value_angle_esd
1 CB 1 SG 2 SG 104.2 2.1
1 SG 2 SG 2 CB 104.2 2.1
"#,
    );
    std::fs::write(file, text).unwrap();
    let ml = MonLib::load_validation(&fixture.0).unwrap();
    let mut st = one_residue("CYS", &[("CB", "C"), ("SG", "S")]);
    let mut second = st.models[0].chains[0].residue_groups[0].clone();
    second.resseq = "2".into();
    for atom in &mut second.atom_groups[0].atoms {
        atom.xyz.x += 2.031;
    }
    st.models[0].chains[0].residue_groups.push(second);
    let graph = build(
        &st,
        &FlatAtoms::from_structure(&st),
        &ml,
        &ValidationOptions {
            use_cdl: false,
            ..Default::default()
        },
    )
    .unwrap();
    let bond = graph.bonds.iter().find(|b| b.indices == [1, 3]).unwrap();
    assert_eq!(bond.origin, reduce3::interp::ORIGIN_SS);
    assert_eq!(bond.ideal, 2.031);
    assert_eq!(bond.weight, 2500.);
    assert!(graph.complete(), "{:?}", graph.coverage);
    assert_eq!(graph.angles.len(), 2);
    for angle in graph.angles {
        assert_eq!(angle.ideal, 104.2);
        assert_eq!(angle.weight, 1. / 2.1_f64.powi(2));
    }
}
