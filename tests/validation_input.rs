use reduce3::MonLib;
use std::path::Path;
#[path = "common/validation_fixture.rs"]
#[allow(dead_code)]
mod fixture;

#[test]
fn energy_overlays_copy_only_changed_tables_and_leave_other_overlays_unchanged() {
    let root = fixture::library();
    let base = MonLib::load_validation(&root.0).unwrap();
    let dictionary =
        "data_energy\nloop_\n_lib_atom.type\n_lib_atom.element\n_lib_atom.vdw_radius\nTEST C 1.7\n";
    let first = base
        .with_validation_cif(&[(dictionary, Path::new("energy"))])
        .unwrap();
    let changed = dictionary.replace("1.7", "2.0");
    let second = base
        .with_validation_cif(&[(&changed, Path::new("changed-energy"))])
        .unwrap();
    assert!(!std::ptr::eq(&*base.ener, &*first.ener));
    assert!(std::ptr::eq(&*base.mods, &*first.mods));
    assert!(base.ener.get("TEST").is_none());
    assert_eq!(first.ener["TEST"].vdw_radius, Some(1.7));
    assert_eq!(second.ener["TEST"].vdw_radius, Some(2.0));
}

#[test]
fn malformed_link_sides_and_unknown_energy_tables_fail() {
    let root = fixture::library();
    let base = MonLib::load_validation(&root.0).unwrap();
    for text in [
        "data_link_BAD\nloop_\n_chem_link_bond.atom_1_comp_id\n_chem_link_bond.atom_id_1\n_chem_link_bond.atom_2_comp_id\n_chem_link_bond.atom_id_2\n3 C 2 N\n",
        "data_energy\nloop_\n_lib_bond.atom_type_1\n_lib_bond.atom_type_2\nC N\n",
        "data_mod_BAD\nloop_\n_chem_mod_atom.function\n_chem_mod_atom.atom_id\nignore C\n",
    ] {
        assert!(
            base.with_validation_cif(&[(text, Path::new("bad"))])
                .is_err()
        );
    }
}

#[test]
fn validation_aliases_keep_reference_case_without_changing_preparation() {
    let root = fixture::library();
    let file = root.0.join("mon_lib/list/mon_lib_list.cif");
    let mut text = std::fs::read_to_string(&file).unwrap();
    text.push_str("\ndata_comp_synonym_list\nloop_\n_chem_comp_synonym.comp_id\n_chem_comp_synonym.comp_alternative_id\nDA Ad\n");
    std::fs::write(file, &text).unwrap();
    std::fs::write(root.0.join("geostd/list/mon_lib_list.cif"), text).unwrap();
    let validation = MonLib::load_validation(&root.0).unwrap();
    let preparation = MonLib::load(&root.0).unwrap();
    assert_eq!(
        validation.comp_synonyms.get("Ad").map(String::as_str),
        Some("DA")
    );
    assert!(validation.comp_synonyms.get("AD").is_none());
    assert!(preparation.comp_synonyms.get("Ad").is_none());
    assert_eq!(
        preparation.comp_synonyms.get("AD").map(String::as_str),
        Some("DA")
    );
}

#[test]
fn private_dictionary_overlay_cannot_change_the_base_component() {
    let root = fixture::library();
    let base = MonLib::load_validation(&root.0).unwrap();
    let text = "data_comp_XYZ\nloop_\n_chem_comp.id\n_chem_comp.group\nXYZ non-polymer\nloop_\n_chem_comp_atom.atom_id\n_chem_comp_atom.type_symbol\n_chem_comp_atom.type_energy\nC C CH1\nH H HCH1\nloop_\n_chem_comp_bond.atom_id_1\n_chem_comp_bond.atom_id_2\n_chem_comp_bond.value_dist\n_chem_comp_bond.value_dist_esd\nC H .97 .02\n";
    let overlay = base
        .with_validation_cif(&[(text, Path::new("test-dictionary"))])
        .unwrap();
    assert!(std::ptr::eq(&*base.ener, &*overlay.ener));
    assert!(std::ptr::eq(&*base.mods, &*overlay.mods));
    assert!(std::ptr::eq(&*base.comp_synonyms, &*overlay.comp_synonyms));
    assert!(base.comp("XYZ").is_none());
    let comp = overlay.comp("XYZ").unwrap();
    assert_eq!(comp.bonds[0].value_dist, Some(0.97));
    assert_eq!(comp.bonds[0].esd, Some(0.02));
    let duplicate = text.replace("C C CH1\n", "C C CH1\nC C CH1\n");
    assert!(
        base.with_validation_cif(&[(&duplicate, Path::new("bad-dictionary"))])
            .is_err()
    );
}
