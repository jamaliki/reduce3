//! Explicit metal hydrides keep both original parent rows.

use reduce3::geom::v3;
use reduce3::model::{Atom, AtomPath, Structure};
use reduce3::{mmcif, pipeline, MonLib, Params};

fn library() -> MonLib {
    MonLib::load(&MonLib::locate(None).expect("Bridge tests require REDUCE3_CHEM_DATA")).unwrap()
}

fn fixture(component: &str) -> Structure {
    let rows = match component {
        "NWN" => concat!(
            "FE FE 19.249 5.075 -13.055\nNI NI 17.788 5.315 -15.164\n",
            "C C1 21.026 5.409 -13.514\nN N1 22.114 5.592 -13.838\n",
            "C C2 19.620 3.270 -12.628\nN N2 19.817 2.145 -12.312\n",
            "C C3 19.499 5.427 -11.383\nO O3 19.663 5.636 -10.236\n",
        ),
        "EJ2" => concat!(
            "FE FE -4.600 -5.078 26.969\nNI NI -5.053 -2.663 26.213\n",
            "C C1 -3.495 -5.027 28.482\nN N1 -2.796 -4.998 29.405\n",
            "C C2 -3.227 -5.893 25.921\nN N2 -2.353 -6.357 25.319\n",
            "C C3 -5.225 -6.686 27.290\nO O3 -5.613 -7.733 27.462\n",
        ),
        _ => panic!("Unknown bridge source component"),
    }.lines().enumerate().map(|(i, row)| format!("HETATM {} {row} {component} L 601 . ? 1.0 20.0 1\n", i + 1)).collect::<String>();
    let header = concat!(
        "data_bridge\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n",
        "_atom_site.label_atom_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n",
        "_atom_site.label_comp_id\n_atom_site.auth_asym_id\n_atom_site.auth_seq_id\n",
        "_atom_site.label_alt_id\n_atom_site.pdbx_PDB_ins_code\n_atom_site.occupancy\n",
        "_atom_site.B_iso_or_equiv\n_atom_site.pdbx_PDB_model_num\n",
    );
    mmcif::read_mmcif(&format!("{header}{rows}")).unwrap()
}

fn find(st: &Structure, name: &str) -> AtomPath {
    st.atom_paths().into_iter().find(|p| st.atom(*p).name.trim() == name).unwrap()
}

fn params(st: &Structure) -> Params {
    let rows = st.atom_paths().into_iter().filter(|p| !st.atom(*p).is_hydrogen())
        .map(|p| st.atom(p).src).collect::<Vec<_>>();
    Params {
        stop_on_any_missing_hydrogen: true,
        required_hydrogen_parent_sources: Some(rows.clone()),
        allowed_unplaced_hydrogen_parent_sources: rows,
        allowed_unplaced_bridge_hydrogens: vec![("H".into(), [st.atom(find(st, "FE")).src, st.atom(find(st, "NI")).src])],
        ..Params::default()
    }
}

#[test]
fn exact_source_hydrides_keep_both_rows_without_inventing_coordinates() {
    let ml = library();
    for component in ["NWN", "EJ2"] {
        let input = fixture(component);
        let output = pipeline::run(input.clone(), &ml, &params(&input)).unwrap();
        assert_eq!(output.hydrogen_uncertainties.len(), 1);
        let site = &output.hydrogen_uncertainties[0];
        assert_eq!(site.hydrogen_name, "H");
        assert_eq!(site.parent_sources, [0, 1]);
        assert_eq!(site.reason, pipeline::HydrogenUncertaintyReason::MissingObservedBridge);
        assert!(output.structure.atom_paths().into_iter().all(|p| output.structure.atom(p).name.trim() != "H"));
        for p in input.atom_paths() {
            let atom = input.atom(p);
            let after = output.structure.atom_paths().into_iter().find(|p| output.structure.atom(*p).src == atom.src).unwrap();
            let xyz = |atom: &Atom| [atom.xyz.x as f32, atom.xyz.y as f32, atom.xyz.z as f32];
            assert_eq!(xyz(output.structure.atom(after)), xyz(atom));
        }
        let mut undeclared = params(&input);
        undeclared.allowed_unplaced_bridge_hydrogens.clear();
        let strict = pipeline::run(input.clone(), &ml, &undeclared);
        match strict {
            Ok(strict) => {
                assert!(!strict.hydrogen_uncertainties.iter().any(|site| site.parent_sources.len() == 2));
                assert!(strict.structure.atom_paths().into_iter().any(|p| strict.structure.atom(p).name.trim() == "H"));
            }
            Err(error) => assert!(error.contains("Insufficient restraints"), "{error}"),
        }
    }
}

#[test]
fn observed_hydride_keeps_its_original_coordinate_charge_and_row() {
    let ml = library();
    for (name, element) in [("H", "H"), ("H", "D"), ("D", "D")] {
    let mut input = fixture("EJ2");
    let mut hydrogen = Atom::new(name, element, v3(-3.806, -3.506, 26.700));
    hydrogen.src = 8;
    hydrogen.charge = "1-".into();
    hydrogen.occ = 0.61;
    hydrogen.b = 12.37;
    input.models[0].chains[0].residue_groups[0].atom_groups[0].atoms.push(hydrogen.clone());
    for keep_existing_h in [false, true] {
        let mut policy = params(&input);
        policy.keep_existing_h = keep_existing_h;
        let output = pipeline::run(input.clone(), &ml, &policy).unwrap();
        assert!(output.hydrogen_uncertainties.is_empty());
        let atom = output.structure.atom(find(&output.structure, name));
        assert_eq!(atom.src, hydrogen.src);
        assert_eq!(atom.xyz, hydrogen.xyz);
        assert_eq!(atom.charge, hydrogen.charge);
        assert_eq!(atom.occ, hydrogen.occ);
        assert_eq!(atom.b, hydrogen.b);
        assert_eq!(atom.name, hydrogen.name);
        assert_eq!(atom.element, hydrogen.element);
    }
    }
}

#[test]
fn unsupported_or_ambiguous_bridge_isotopes_fail_before_placement() {
    let ml = library();
    for isotopes in [vec![("T", "T")], vec![("H", "T")], vec![("D", "H")],
        vec![("H", "H"), ("D", "D")]] {
        let mut input = fixture("EJ2");
        for (offset, (name, element)) in isotopes.iter().enumerate() {
            let mut atom = Atom::new(name, element, v3(-3.806, -3.506, 26.700));
            atom.src = 8 + offset as u32;
            input.models[0].chains[0].residue_groups[0].atom_groups[0].atoms.push(atom);
        }
        let error = pipeline::run(input.clone(), &ml, &params(&input)).err().unwrap();
        assert!(error.contains("tritium") || error.contains("isotope identity"), "{error}");
    }
}

#[test]
fn invalid_hydride_source_declarations_stay_strict() {
    let ml = library();
    for case in 0..12 {
        let mut input = fixture("EJ2");
        let mut policy = params(&input);
        let ni = find(&input, "NI");
        match case {
            0 => policy.allowed_unplaced_bridge_hydrogens[0].1 = [0, 0],
            1 => policy.allowed_unplaced_bridge_hydrogens[0].1 = [0, 99],
            2 => policy.allowed_unplaced_bridge_hydrogens[0].1 = [0, 2],
            3 => policy.allowed_unplaced_bridge_hydrogens[0].0 = "H2".into(),
            4 => input.atom_mut(ni).element = "C".into(),
            5 => input.atom_mut(ni).xyz.x = f64::NAN,
            6 => policy.allowed_unplaced_bridge_hydrogens.push(policy.allowed_unplaced_bridge_hydrogens[0].clone()),
            7 => {
                let mut other = input.models[0].chains[0].residue_groups[0].atom_groups[0].clone();
                other.resname = "NWN".into();
                input.models[0].chains[0].residue_groups[0].atom_groups.push(other);
            }
            8 => {
                let mut other = input.atom_group(ni).clone();
                other.atoms.retain(|atom| atom.name.trim() == "NI");
                other.atoms[0].src = 100;
                other.altloc = "A".into();
                input.models[0].chains[0].residue_groups[0].atom_groups.push(other);
                policy.allowed_unplaced_bridge_hydrogens[0].1 = [0, 100];
            }
            9 => input.atom_mut(ni).element = "FE".into(),
            10 => {
                let mut other = input.atom(ni).clone();
                other.src = 100;
                input.models[0].chains[0].residue_groups[0].atom_groups[0].atoms.push(other);
            }
            _ => {
                let mut other = input.models[0].clone();
                other.id = "2".into();
                for group in &mut other.chains[0].residue_groups[0].atom_groups {
                    for atom in &mut group.atoms { atom.src += 100; }
                }
                input.models.push(other);
                policy.model_id = Some(1);
                policy.allowed_unplaced_bridge_hydrogens[0].1 = [0, 101];
            }
        }
        assert!(pipeline::run(input, &ml, &policy).is_err(), "invalid declaration {case}");
    }
}
