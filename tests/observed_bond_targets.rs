//! Exact deposited bond failures retain heavy atoms and observed H.

use reduce3::geom::v3;
use reduce3::interp::target_uncertainty::{self, BondUncertaintyReason as Reason};
use reduce3::model::{Atom, AtomPath, Structure};
use reduce3::{mmcif, pipeline, MonLib, Params};

fn library() -> MonLib {
    MonLib::load(&MonLib::locate(None).expect("Set REDUCE3_CHEM_DATA")).unwrap()
}

fn find(st: &Structure, component: &str, name: &str) -> AtomPath {
    st.atom_paths().into_iter().find(|p| st.atom_group(*p).resname.trim() == component
        && st.atom(*p).name.trim() == name).unwrap()
}

fn required(st: &Structure) -> Vec<u32> {
    st.atom_paths().into_iter().filter(|p| !st.atom(*p).is_hydrogen()).map(|p| st.atom(p).src).collect()
}

fn policy(st: &Structure) -> Params {
    Params { stop_on_any_missing_hydrogen: true, required_hydrogen_parent_sources: Some(required(st)),
        allow_unresolved_bond_targets: true, ..Params::default() }
}

fn observed(st: &mut Structure, component: &str, parent: &str, name: &str) -> Atom {
    let p = find(st, component, parent);
    let mut atom = Atom::new(name, "H", st.atom(p).xyz + v3(0.27, 0.31, 0.87));
    atom.src = 999_001;
    atom.charge = "1-".into();
    atom.occ = 0.43;
    atom.b = 27.1;
    atom.segid = "OBS".into();
    atom.uij = Some([0.1, 0.2, 0.3, 0.01, 0.02, 0.03]);
    atom.hetero = true;
    st.models[p.model as usize].chains[p.chain as usize].residue_groups[p.rg as usize]
        .atom_groups[p.ag as usize].atoms.push(atom.clone());
    atom
}

fn assert_observed(st: &Structure, original: &Atom) {
    let found = st.atom_paths().into_iter().filter(|p| st.atom(*p).src == original.src).collect::<Vec<_>>();
    assert_eq!(found.len(), 1);
    let atom = st.atom(found[0]);
    assert_eq!(atom.xyz, original.xyz);
    assert_eq!(atom.name, original.name);
    assert_eq!(atom.element, original.element);
    assert_eq!(atom.charge, original.charge);
    assert_eq!(atom.occ, original.occ);
    assert_eq!(atom.b, original.b);
    assert_eq!(atom.segid, original.segid);
    assert_eq!(atom.uij, original.uij);
    assert_eq!(atom.hetero, original.hetero);
}

#[test]
fn accepted_deposited_sg_link_retains_observed_hg_as_chemical_only_uncertainty() {
    let ml = library();
    let mut input = mmcif::read_mmcif(include_str!("data/linked_9hj3.cif")).unwrap();
    let hydrogen = observed(&mut input, "CYS", "SG", "HG");
    let sg = input.atom(find(&input, "CYS", "SG")).src;
    for keep_existing_h in [false, true] {
        let mut params = policy(&input);
        params.required_hydrogen_parent_sources = Some(input.atom_group(find(&input, "CYS", "SG"))
            .atoms.iter().filter(|atom| !atom.is_hydrogen()).map(|atom| atom.src).collect());
        params.keep_existing_h = keep_existing_h;
        let output = pipeline::run(input.clone(), &ml, &params).unwrap();
        assert!(output.bond_uncertainties.iter().any(|report| report.parent_source == sg
            && !report.hydrogens_unknown && report.reason == Reason::UnsupportedResidueEnvironment));
        assert_observed(&output.structure, &hydrogen);
        let mut strict = params;
        strict.allow_unresolved_bond_targets = false;
        let control = pipeline::run(input.clone(), &ml, &strict).unwrap();
        assert!(control.structure.atom_paths().into_iter().all(|p| control.structure.atom(p).src != hydrogen.src));
    }
}

fn degenerate_fixture() -> Structure {
    let mut st = mmcif::read_mmcif(concat!(
        "data_1nlm\nloop_\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_alt_id\n",
        "_atom_site.auth_comp_id\n_atom_site.auth_asym_id\n_atom_site.auth_seq_id\n",
        "_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n",
        "_atom_site.occupancy\n_atom_site.B_iso_or_equiv\n",
        "N N . LEU B 198 83.244 -59.862 7.138 1 10\nC CA . LEU B 198 84.592 -59.773 7.696 1 10\n",
        "C C . LEU B 198 85.145 -61.191 7.934 1 10\nO O . LEU B 198 86.363 -61.421 7.870 1 10\n",
        "C CB . LEU B 198 84.141 -57.491 8.888 1 10\nC CG . LEU B 198 84.141 -57.491 8.888 1 10\n",
        "C CD1 . LEU B 198 84.309 -56.805 10.239 1 10\nC CD2 . LEU B 198 84.958 -56.771 7.839 1 10\n",
        "N N . ASP B 211 103.295 -62.225 7.998 1 10\nC CA . ASP B 211 103.909 -61.927 9.288 1 10\n",
        "C C . ASP B 211 103.974 -60.462 9.698 1 10\nO O . ASP B 211 104.476 -60.149 10.779 1 10\n",
        "C CB . ASP B 211 105.347 -63.948 9.842 1 10\nC CG . ASP B 211 105.347 -63.948 9.842 1 10\n",
        "O OD1 . ASP B 211 106.450 -64.534 9.918 1 10\nO OD2 . ASP B 211 104.264 -64.487 10.160 1 10\n",
    )).unwrap();
    for p in st.atom_paths() {
        let offset = if st.atom_group(p).resname == "LEU" { 4059 } else { 4153 - 8 };
        st.atom_mut(p).src += offset;
    }
    st
}

#[test]
fn exact_1nlm_coincident_intrinsic_bonds_mask_only_four_original_rows() {
    let ml = library();
    let input = degenerate_fixture();
    let params = policy(&input);
    let output = pipeline::run(input.clone(), &ml, &params).unwrap();
    assert_eq!(output.bond_uncertainties.iter().map(|report| report.parent_source).collect::<Vec<_>>(),
        [4063, 4064, 4157, 4158]);
    assert!(output.bond_uncertainties.iter().all(|report| report.hydrogens_unknown
        && report.reason == Reason::DegenerateHeavyBond));
    assert_eq!(output.hydrogen_uncertainties.len(), 5);
    let mut before_positions = Vec::new();
    let mut after_positions = Vec::new();
    for p in input.atom_paths() {
        let atom = input.atom(p);
        let actual = output.structure.atom_paths().into_iter().find(|p| output.structure.atom(*p).src == atom.src).unwrap();
        let xyz = |atom: &Atom| [atom.xyz.x as f32, atom.xyz.y as f32, atom.xyz.z as f32];
        before_positions.push(xyz(atom));
        after_positions.push(xyz(output.structure.atom(actual)));
        if [4063, 4064, 4157, 4158].contains(&atom.src) { assert_eq!(xyz(atom), xyz(output.structure.atom(actual))); }
    }
    let order = |a: &[f32; 3], b: &[f32; 3]| a[0].total_cmp(&b[0]).then(a[1].total_cmp(&b[1])).then(a[2].total_cmp(&b[2]));
    before_positions.sort_by(order);
    after_positions.sort_by(order);
    assert_eq!(before_positions, after_positions);
    let mut strict = params;
    strict.allow_unresolved_bond_targets = false;
    assert!(pipeline::run(input.clone(), &ml, &strict).err().unwrap().contains("Insufficient restraints"));
    let mut observed_input = input;
    let hydrogen = observed(&mut observed_input, "LEU", "CB", "HB2");
    let output = pipeline::run(observed_input.clone(), &ml, &policy(&observed_input)).unwrap();
    assert_observed(&output.structure, &hydrogen);
}

#[test]
fn coincident_helper_requires_an_accepted_finite_original_heavy_pair() {
    use reduce3::interp::ORIGIN_COVALENT;
    let mut st = degenerate_fixture();
    let (cb, cg) = (find(&st, "LEU", "CB"), find(&st, "LEU", "CG"));
    let paths = st.atom_paths();
    let edge = (paths.iter().position(|p| *p == cb).unwrap() as u32,
        paths.iter().position(|p| *p == cg).unwrap() as u32, ORIGIN_COVALENT);
    for case in 0..5 {
        let mut input = st.clone();
        let mut bonds = vec![edge];
        match case {
            0 => bonds.clear(),
            1 => input.atom_mut(cg).xyz.x += 1e-9,
            2 => input.atom_mut(cg).src = Atom::NEW,
            3 => input.atom_mut(cg).element = "H".into(),
            _ => { input.atom_mut(cb).xyz.x = f64::INFINITY; input.atom_mut(cg).xyz.x = f64::INFINITY; }
        }
        let mut reports = Vec::new();
        target_uncertainty::add_degenerate_heavy_bonds(&input, &bonds, &mut reports);
        assert!(reports.is_empty(), "case {case}");
    }
    st.atom_mut(cg).xyz.x += 1.0;
    let ml = library();
    st.atom_mut(find(&st.clone(), "ASP", "CG")).xyz.x += 1.0;
    let output = pipeline::run(st.clone(), &ml, &policy(&st)).unwrap();
    assert!(output.bond_uncertainties.is_empty());
}
