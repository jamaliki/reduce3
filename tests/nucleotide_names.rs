use reduce3::{names, MonLib};

#[test]
fn thymine_methyl_names_share_the_v3_reference() {
    for (given, reference) in [("C5M", "C7"), ("H5M1", "H71"), ("H5M2", "H72"), ("H5M3", "H73")] {
        assert_eq!(names::na_reference_name(given), reference);
        assert_eq!(names::na_reference_name(reference), reference);
    }
    assert_eq!(names::na_reference_name("C1*"), "C1'");
    assert_eq!(names::na_reference_name("O1P"), "OP1");
}

#[test]
fn dictionary_mapping_accepts_both_thymine_methyl_spellings() {
    let Some(root) = MonLib::locate(None) else { return };
    let library = MonLib::load(&root).unwrap();
    let dictionary = library.comp("TD").unwrap();
    let v2 = ["C5M", "H5M1", "H5M2", "H5M3"].map(String::from);
    let v3 = ["C7", "H71", "H72", "H73"].map(String::from);
    let old = names::map_atoms(&dictionary, &v2, None, true, None);
    let current = names::map_atoms(&dictionary, &v3, None, true, None);
    assert!(old.unexpected.is_empty() && current.unexpected.is_empty());
    assert_eq!(old.ids, current.ids);
    assert!(old.ids.iter().all(Option::is_some));
}
