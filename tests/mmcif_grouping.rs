use reduce3::interp::target_uncertainty::BondUncertaintyReason;
use reduce3::model::{Atom, AtomPath, Structure};
use reduce3::{cif, mmcif, pipeline, MonLib, Params};

// Checked binding input SHA256: ddd3b8d5d99dbb213b8698da2d764197202b12780a9ad003cfd379b6c7fb08f8.
// The added observed HG follows the original 52 heavy rows from the 9HJ3 fixture.
const OBSERVED: &str = include_str!("data/interleaved_observed_hg.cif");

fn reordered(order: &[usize]) -> Structure {
    let mut document = cif::parse(OBSERVED);
    let table = document.blocks[0].categories.iter_mut().find(|table| table.name == "_atom_site").unwrap();
    let old = table.values.clone();
    let width = table.ncols();
    table.values = order.iter().flat_map(|&row| old[row * width..(row + 1) * width].iter().copied()).collect();
    mmcif::structure_from_cif(&document.blocks[0]).unwrap()
}

fn source(st: &Structure, row: u32) -> AtomPath {
    let found = st.atom_paths().into_iter().filter(|p| st.atom(*p).src == row).collect::<Vec<_>>();
    assert_eq!(found.len(), 1, "source row {row}");
    found[0]
}

fn metadata(actual: &Atom, original: &Atom) {
    assert_eq!(actual.src, original.src);
    assert_eq!(actual.name, original.name);
    assert_eq!(actual.element, original.element);
    assert_eq!(actual.charge, original.charge);
    assert_eq!(actual.xyz.as_array().map(|value| (value as f32).to_bits()),
        original.xyz.as_array().map(|value| (value as f32).to_bits()),
        "source row {} ({})", original.src, original.name.trim());
    assert_eq!(actual.occ.to_bits(), original.occ.to_bits());
    assert_eq!(actual.b.to_bits(), original.b.to_bits());
    assert_eq!(actual.segid, original.segid);
    assert_eq!(actual.uij, original.uij);
    assert_eq!(actual.hetero, original.hetero);
}

#[test]
fn appended_observed_hg_joins_its_original_heavy_group_and_keeps_row_metadata() {
    let st = mmcif::read_mmcif(OBSERVED).unwrap();
    let hg = source(&st, 52);
    let sg = source(&st, 5);
    assert_eq!((hg.model, hg.chain, hg.rg, hg.ag), (sg.model, sg.chain, sg.rg, sg.ag));
    assert_eq!(st.models.len(), 1);
    assert!(st.models[0].id.is_empty());
    assert_eq!(st.models[0].chains.len(), 2);
    assert_eq!(st.models[0].chains[0].residue_groups.len(), 1);
    assert_eq!(st.models[0].chains[1].residue_groups.iter().map(|rg| rg.resseq.trim()).collect::<Vec<_>>(), ["601", "602"]);
    assert_eq!(st.atom_group(hg).atoms.iter().map(|atom| atom.src).collect::<Vec<_>>(), [0, 1, 2, 3, 4, 5, 52]);
    let atom = st.atom(hg);
    assert_eq!(atom.xyz.as_array(), [151.99800000000002, 173.434, 173.09900000000002]);
    assert_eq!(atom.name.trim(), "HG");
    assert_eq!(atom.charge, "1-");
    assert_eq!(atom.occ, 0.43);
    assert_eq!(atom.b, 27.1);
    assert!(!atom.hetero);
    assert_eq!(st.atoms_size(), 53);
}

#[test]
fn interleaved_heavy_and_h_rows_match_the_contiguous_layout_with_exact_source_rows() {
    let original = mmcif::read_mmcif(OBSERVED).unwrap();
    let order = [0, 6, 1, 7, 2, 8, 3, 9, 4, 10, 5, 11].into_iter().chain(12..53).collect::<Vec<_>>();
    let st = reordered(&order);
    assert_eq!(mmcif::write_mmcif(&st), mmcif::write_mmcif(&original));
    for (new_row, &old_row) in order.iter().enumerate() {
        let actual = st.atom(source(&st, new_row as u32));
        let mut expected = original.atom(source(&original, old_row as u32)).clone();
        expected.src = new_row as u32;
        metadata(actual, &expected);
        assert_eq!(actual.xyz.as_array().map(f64::to_bits), expected.xyz.as_array().map(f64::to_bits));
    }
    let document = cif::parse(OBSERVED);
    assert_eq!(mmcif::write_mmcif(&mmcif::structure_from_cif(&document.blocks[0]).unwrap()), mmcif::write_mmcif(&original));
}

#[test]
fn both_observed_h_orders_preserve_targets_without_relaxing_strict_defaults() {
    let ml = MonLib::load(&MonLib::locate(None).expect("Set REDUCE3_CHEM_DATA")).unwrap();
    let contiguous = (0..6).chain([52]).chain(6..52).collect::<Vec<_>>();
    for input in [mmcif::read_mmcif(OBSERVED).unwrap(), reordered(&contiguous)] {
        let sg = input.atom_paths().into_iter().find(|p| input.atom(*p).name.trim() == "SG").unwrap();
        let hg = input.atom_paths().into_iter().find(|p| input.atom(*p).name.trim() == "HG").unwrap();
        let observed = input.atom(hg).clone();
        let required: Vec<u32> = input.atom_group(sg).atoms.iter().filter(|atom| !atom.is_hydrogen()).map(|atom| atom.src).collect();
        for keep_existing_h in [false, true] {
            let params = Params { stop_on_any_missing_hydrogen: true, keep_existing_h,
                required_hydrogen_parent_sources: Some(required.clone()), allow_unresolved_bond_targets: true,
                ..Params::default() };
            let output = pipeline::run(input.clone(), &ml, &params).unwrap();
            let actual = output.structure.atom(source(&output.structure, observed.src));
            metadata(actual, &observed);
            assert_eq!(actual.xyz.as_array().map(f64::to_bits), observed.xyz.as_array().map(f64::to_bits));
            assert!(output.hydrogen_uncertainties.is_empty());
            assert_eq!(output.bond_uncertainties.len(), 1);
            let report = &output.bond_uncertainties[0];
            assert_eq!(report.parent_source, input.atom(sg).src);
            assert!(!report.hydrogens_unknown);
            assert_eq!(report.reason, BondUncertaintyReason::UnsupportedResidueEnvironment);
            for p in input.atom_paths().into_iter().filter(|p| !input.atom(*p).is_hydrogen()) {
                let before = input.atom(p);
                metadata(output.structure.atom(source(&output.structure, before.src)), before);
            }
            let mut strict = params;
            strict.allow_unresolved_bond_targets = false;
            let output = pipeline::run(input.clone(), &ml, &strict).unwrap();
            assert!(output.structure.atom_paths().into_iter().all(|p| output.structure.atom(p).src != observed.src));
        }
    }
}

fn small(rows: &str) -> Structure {
    let header = concat!("data_groups\nloop_\n_atom_site.pdbx_PDB_model_num\n_atom_site.label_asym_id\n",
        "_atom_site.auth_asym_id\n_atom_site.auth_seq_id\n_atom_site.pdbx_PDB_ins_code\n",
        "_atom_site.label_comp_id\n_atom_site.label_alt_id\n_atom_site.label_atom_id\n",
        "_atom_site.type_symbol\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n",
        "_atom_site.occupancy\n_atom_site.B_iso_or_equiv\n_atom_site.auth_break\n");
    mmcif::read_mmcif(&format!("{header}{rows}")).unwrap()
}

#[test]
fn first_seen_residue_order_altlocs_components_insertions_and_breaks_stay_distinct() {
    let st = small(concat!(
        "18 A G 18 ? CYS A SG S 0 0 0 1 10 0\n",
        "18 A G 19 ? ALA . N N 1 0 0 1 10 1\n",
        "18 A G 18 ? CYS B SG S 2 0 0 1 10 0\n",
        "18 A G 18 ? CYX A SG S 3 0 0 1 10 0\n",
        "18 A G 18 ? CYS . N N 4 0 0 1 10 0\n",
        "18 A G 18 A CYS . N N 5 0 0 1 10 0\n",
        "18 A G 19 ? ALA . CA C 6 0 0 1 10 0\n",
    ));
    let residues = &st.models[0].chains[0].residue_groups;
    assert_eq!(residues.iter().map(|rg| (rg.resseq.trim(), rg.icode.as_str())).collect::<Vec<_>>(),
        [("18", " "), ("19", " "), ("18", "A")]);
    assert!(!residues[1].link_to_previous);
    assert_eq!(residues[1].atom_groups[0].atoms.iter().map(|atom| atom.src).collect::<Vec<_>>(), [1, 6]);
    assert_eq!(residues[0].atom_groups.iter().map(|group| (group.altloc.as_str(), group.resname.as_str())).collect::<Vec<_>>(),
        [("", "CYS"), ("A", "CYS"), ("B", "CYS"), ("A", "CYX")]);
}

#[test]
fn auth_chain_returns_keep_first_seen_chain_and_residue_order() {
    let st = small(concat!(
        "1 A X 5 ? ALA . N N 0 0 0 1 10 0\n",
        "1 B Y 9 ? ALA . N N 1 0 0 1 10 0\n",
        "1 A X 3 ? ALA . N N 2 0 0 1 10 1\n",
        "1 B Y 9 ? ALA . CA C 3 0 0 1 10 0\n",
        "1 A X 5 ? ALA . CA C 4 0 0 1 10 0\n",
    ));
    let chains = &st.models[0].chains;
    assert_eq!(chains.iter().map(|chain| chain.id.as_str()).collect::<Vec<_>>(), ["X", "Y"]);
    assert_eq!(chains[0].residue_groups.iter().map(|rg| rg.resseq.trim()).collect::<Vec<_>>(), ["5", "3"]);
    assert!(!chains[0].residue_groups[1].link_to_previous);
    assert_eq!(chains[0].residue_groups[0].atom_groups[0].atoms.iter().map(|atom| atom.src).collect::<Vec<_>>(), [0, 4]);
    assert_eq!(chains[1].residue_groups[0].atom_groups[0].atoms.iter().map(|atom| atom.src).collect::<Vec<_>>(), [1, 3]);
}

#[test]
fn repeated_and_numeric_alias_model_blocks_keep_the_existing_identity_contract() {
    for ids in [vec!["35", "18"], vec!["18", "35", "18"], vec!["18", "018"]] {
        let rows = ids.iter().map(|id| format!("{id} A G 18 ? CYS . SG S 0 0 0 1 10 0\n")).collect::<String>();
        let st = small(&rows);
        assert_eq!(st.models.iter().map(|model| model.id.as_str()).collect::<Vec<_>>(), ids);
        assert_eq!(st.models.iter().map(|model| model.chains[0].residue_groups[0].atom_groups[0].atoms.len()).collect::<Vec<_>>(), vec![1; ids.len()]);
        for (row, p) in st.atom_paths().into_iter().enumerate() { assert_eq!(st.atom(p).src, row as u32); }
    }
}

#[test]
fn repeated_atom_identity_remains_visible_and_cannot_authorize_observed_h() {
    let order = (0..53).chain([52]).collect::<Vec<_>>();
    let input = reordered(&order);
    let group = input.atom_group(source(&input, 5));
    assert_eq!(group.atoms.iter().filter(|atom| atom.name.trim() == "HG").map(|atom| atom.src).collect::<Vec<_>>(), [52, 53]);
    let ml = MonLib::load(&MonLib::locate(None).expect("Set REDUCE3_CHEM_DATA")).unwrap();
    let error = pipeline::run(input, &ml, &Params { stop_on_any_missing_hydrogen: true,
        allow_unresolved_bond_targets: true, ..Params::default() }).err().unwrap();
    assert!(error.contains("no unique retained dictionary parent"), "{error}");
}
