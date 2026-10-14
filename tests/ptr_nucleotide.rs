use reduce3::{hplace, mmcif, pipeline, MonLib, Params};
use reduce3::model::Structure;

fn library() -> MonLib {
    MonLib::load(&MonLib::locate(None).expect("set REDUCE3_CHEM_DATA for the PTR nucleotide test")).unwrap()
}

fn has(st: &Structure, name: &str) -> bool {
    st.atom_paths().into_iter().any(|p| st.chain(p).id == "B" && st.atom_group(p).resname.trim() == "DT"
        && st.atom(p).name.trim() == name)
}

#[test]
fn declared_ptr_nucleotide_product_keeps_both_named_c5_hydrogens() {
    let input = mmcif::read_mmcif(include_str!("data/ptr_dt_3px7.cif")).unwrap();
    assert_eq!(input.links.len(), 1);
    let required = input.atom_paths().into_iter().filter(|&p| input.atom_group(p).resname.trim() == "DT")
        .map(|p| input.atom(p).src).collect();
    let original: Vec<_> = input.atom_paths().into_iter().map(|p| (input.chain(p).id.clone(), input.atom(p).name.clone(), input.atom(p).xyz, input.atom(p).charge.clone())).collect();
    let output = pipeline::run(input, &library(), &Params {
        stop_on_any_missing_hydrogen: true, required_hydrogen_parent_sources: Some(required), ..Params::default()
    }).unwrap();
    assert!(has(&output.structure, "H5'") && has(&output.structure, "H5''"));
    assert!(!has(&output.structure, "O5'"));
    for (chain, name, xyz, charge) in original {
        let p = output.structure.atom_paths().into_iter().find(|&p| output.structure.chain(p).id == chain
            && output.structure.atom(p).name == name).unwrap();
        assert!(output.structure.atom(p).xyz.dist(xyz) < 1e-12);
        assert_eq!(output.structure.atom(p).charge, charge);
    }
}

#[test]
fn missing_ambiguous_or_invalid_declaration_does_not_supply_c5_angles() {
    for invalid in 0..8 {
        let mut input = mmcif::read_mmcif(include_str!("data/ptr_dt_3px7.cif")).unwrap();
        match invalid {
            0 => input.links.clear(),
            1 => input.links.push(input.links[0].clone()),
            2 => input.links[0].order = Some("unknown".into()),
            3 => input.links[0].atom1.chain = "X".into(),
            4 => input.links[0].atom2.resseq = "2".into(),
            5 => input.links[0].order = Some("doub".into()),
            6 => {
                for chain in &mut input.models[0].chains {
                    for residue in &mut chain.residue_groups {
                        for group in &mut residue.atom_groups {
                            group.altloc = if chain.id == "A" { "A" } else { "B" }.into();
                        }
                    }
                }
            }
            7 => {
                let other = input.models[0].chains.pop().unwrap();
                let mut model = input.models[0].clone();
                model.id = "2".into();
                model.chains = vec![other];
                input.models.push(model);
            }
            _ => unreachable!(),
        }
        let placed = hplace::place_hydrogens(&mut input, &library(), &hplace::HPlaceParams {
            neutron: false, n_terminal_charge: hplace::NTermCharge::ResidueOne,
            exclude_water: true, keep_existing_h: false, adp_scale: 1.0, compat: false, cell: None,
        });
        if invalid < 6 {
            assert!(placed.site_labels_no_para.iter().any(|name| name.trim().starts_with("H5'")), "case {invalid}");
        }
        assert!(!has(&input, "H5'") && !has(&input, "H5''"), "case {invalid}");
    }
}

#[test]
fn compatibility_mode_keeps_its_existing_missing_parameter_result() {
    let mut input = mmcif::read_mmcif(include_str!("data/ptr_dt_3px7.cif")).unwrap();
    let placed = hplace::place_hydrogens(&mut input, &library(), &hplace::HPlaceParams {
        neutron: false, n_terminal_charge: hplace::NTermCharge::ResidueOne,
        exclude_water: true, keep_existing_h: false, adp_scale: 1.0, compat: true, cell: None,
    });
    assert!(placed.site_labels_no_para.iter().any(|name| name.trim().starts_with("H5'")));
    assert!(!has(&input, "H5'") && !has(&input, "H5''"));
}
