use super::*;

fn fixture() -> (MonLib, Comp, crate::monlib::CompBond) {
    let library = MonLib::load(&MonLib::locate(None).expect("set REDUCE3_CHEM_DATA for source-order tests")).unwrap();
    let comp = (*library.comp("TD").unwrap()).clone();
    let bond = comp.bonds.iter().find(|b| b.a1 == "C5'" && b.a2 == "O5'").unwrap().clone();
    assert_eq!(bond.type_, "coval");
    (library, comp, bond)
}

#[test]
fn actual_dt_proves_aliased_coval_bonds_without_changing_the_dictionary() {
    let (library, comp, bond) = fixture();
    assert_eq!(source_bond_order(&library, &comp, &bond, "DT"), Some(2));
    let terminal = library.apply_mod(&comp, "5*END").unwrap();
    assert!(terminal.test_for_rna_dna().is_none());
    assert_eq!(source_bond_order(&library, &terminal, &bond, " DT"), Some(2));
    assert_eq!(comp.id, "Td");
    assert_eq!(bond.type_, "coval");
    assert_eq!(bond.value_dist, Some(1.418));
    let double = comp.bonds.iter().find(|b| b.a1 == "P" && b.a2 == "OP1").unwrap();
    assert_eq!(source_bond_order(&library, &comp, double, "DT"), Some(4));
}

#[test]
fn unknown_components_atoms_elements_or_orders_cannot_prove_a_single_bond() {
    for invalid in 0..9 {
        let (library, mut comp, mut bond) = fixture();
        let actual = match invalid {
            0 => "TD",
            1 => "UNKNOWN",
            2 => "DA",
            3 => { comp.atoms.iter_mut().find(|a| a.id == "O5'").unwrap().type_symbol = "N".into(); "DT" }
            4 => { bond.a2 = "ABSENT".into(); "DT" }
            5 => { bond.type_ = "unknown".into(); "DT" }
            6 => { comp.source = "unknown".into(); "DT" }
            7 => { comp.from_ccd = true; "DT" }
            8 => { comp.atoms.push(comp.atom("O5'").unwrap().clone()); "DT" }
            _ => unreachable!(),
        };
        assert_eq!(source_bond_order(&library, &comp, &bond, actual), None, "case {invalid}");
    }
}
