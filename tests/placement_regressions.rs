use reduce3::{hplace, mmcif, pipeline, MonLib, Params};
use reduce3::model::{AtomPath, Structure};

fn library() -> Option<MonLib> {
    let Some(root) = MonLib::locate(None) else {
        eprintln!("skipped: chem_data not found (set REDUCE3_CHEM_DATA)");
        return None;
    };
    Some(MonLib::load(&root).unwrap())
}

fn find(st: &Structure, component: &str, sequence: &str, name: &str) -> AtomPath {
    st.atom_paths().into_iter().find(|&p| st.atom_group(p).resname.trim() == component
        && st.residue_group(p).resseq.trim() == sequence && st.atom(p).name.trim() == name).unwrap()
}

#[test]
fn partial_esterified_sugar_removes_its_h_before_strict_parameter_check() {
    let Some(library) = library() else { return };
    let input = mmcif::read_mmcif(include_str!("data/linked_xyp_1ux7.cif")).unwrap();
    let before = input.atom(find(&input, "XYP", "1", "O1")).xyz;
    let output = pipeline::run(input.clone(), &library, &Params {
        stop_on_any_missing_hydrogen: true, ..Params::default()
    }).unwrap();
    assert!(output.structure.atom_paths().into_iter().all(|p| st_sequence(&output.structure, p) != "1"
        || output.structure.atom(p).name.trim() != "HO1"));
    assert!(output.structure.atom(find(&output.structure, "XYP", "1", "O1")).xyz.dist(before) < 1e-10);
    let lenient = pipeline::run(input, &library, &Params::default()).unwrap();
    assert_eq!(mmcif::write_mmcif(&output.structure), mmcif::write_mmcif(&lenient.structure));
}

fn st_sequence(st: &Structure, p: AtomPath) -> &str { st.residue_group(p).resseq.trim() }

#[test]
fn alternate_conformer_cannot_link_across_a_skipped_residue_group() {
    let Some(library) = library() else { return };
    let mut input = mmcif::read_mmcif(include_str!("data/altloc_peptide_5z3k.cif")).unwrap();
    let placed = hplace::place_hydrogens(&mut input, &library, &hplace::HPlaceParams {
        neutron: false, n_terminal_charge: hplace::NTermCharge::ResidueOne,
        exclude_water: true, keep_existing_h: false, adp_scale: 1.0, compat: false, cell: None,
    });
    assert!(placed.site_labels_no_para.is_empty());
    let nitrogen = input.atom(find(&input, "ASP", "235", "N")).i_seq as u32;
    let carbon = input.atom(find(&input, "ARG", "234", "C")).i_seq as u32;
    let skipped = input.atom(find(&input, "SER", "233", "C")).i_seq as u32;
    let bonded = |i, j| placed.bonds.iter().any(|&(a, b, _)| (a == i && b == j) || (b == i && a == j));
    assert!(bonded(nitrogen, carbon));
    assert!(!bonded(nitrogen, skipped));
    let hydrogen = input.atom(find(&input, "ASP", "235", "H")).xyz;
    assert!((hydrogen.dist(input.atom(find(&input, "ASP", "235", "N")).xyz) - 0.86).abs() < 1e-10);
}
