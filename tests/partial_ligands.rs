//! Undefined ligand H targets keep their parent rows and observed H.

#[path = "common/chem_data.rs"]
mod chem_data;

use reduce3::geom::v3;
use reduce3::model::{Atom, AtomPath, Structure};
use reduce3::pipeline::{HydrogenUncertaintyReason, Output};
use reduce3::{hplace, mmcif, pipeline, MonLib, Params};

fn library() -> Option<MonLib> {
    let root = chem_data::locate()?;
    Some(MonLib::load(&root).unwrap())
}

fn find(st: &Structure, component: &str, name: &str) -> AtomPath {
    st.atom_paths().into_iter().find(|p| st.atom_group(*p).resname.trim() == component
        && st.atom(*p).name.trim() == name).unwrap()
}

fn rows(st: &Structure, component: &str) -> Vec<u32> {
    st.atom_paths().into_iter().filter_map(|p| {
        let atom = st.atom(p);
        (st.atom_group(p).resname.trim() == component && !atom.is_hydrogen()).then_some(atom.src)
    }).collect()
}

fn params(st: &Structure, component: &str) -> Params {
    let ligand = rows(st, component);
    Params {
        stop_on_any_missing_hydrogen: true,
        required_hydrogen_parent_sources: Some(ligand.clone()),
        allowed_unplaced_hydrogen_parent_sources: ligand,
        ..Params::default()
    }
}

fn heavy_unchanged(before: &Structure, output: &Output) {
    let after = &output.structure;
    for p in before.atom_paths().into_iter().filter(|p| !before.atom(*p).is_hydrogen()) {
        let atom = before.atom(p);
        let actual = after.atom_paths().into_iter().find(|p| after.atom(*p).src == atom.src).unwrap();
        let xyz = |atom: &Atom| [atom.xyz.x as f32, atom.xyz.y as f32, atom.xyz.z as f32];
        assert_eq!(xyz(atom), xyz(after.atom(actual)));
        assert_eq!(atom.charge, after.atom(actual).charge);
        assert_eq!(atom.occ, after.atom(actual).occ);
        assert_eq!(atom.b, after.atom(actual).b);
    }
}

#[test]
fn exact_9l0k_tail_reports_only_retained_parent_c6() {
    let Some(ml) = library() else { return };
    let input = mmcif::read_mmcif(include_str!("data/partial_cla_9l0k.cif")).unwrap();
    let c6 = input.atom(find(&input, "CLA", "C6")).src;
    let output = pipeline::run(input.clone(), &ml, &params(&input, "CLA")).unwrap();
    heavy_unchanged(&input, &output);
    let sites = &output.hydrogen_uncertainties;
    assert_eq!(sites.len(), 2);
    assert!(sites.iter().all(|site| site.parent_sources == [c6]
        && site.reason == HydrogenUncertaintyReason::MissingHeavyNeighbor));
    assert_eq!(sites.iter().map(|site| site.hydrogen_name.as_str()).collect::<Vec<_>>(), ["H61", "H62"]);
    assert!(output.structure.atom_paths().into_iter().all(|p| !["H61", "H62"].contains(&output.structure.atom(p).name.trim())));
    let strict = pipeline::run(input.clone(), &ml, &Params {
        stop_on_any_missing_hydrogen: true, ..Params::default()
    }).unwrap();
    assert!(strict.hydrogen_uncertainties.is_empty());
    assert_eq!(mmcif::write_mmcif(&strict.structure), mmcif::write_mmcif(&output.structure));
}

#[test]
fn exact_5lbm_formaldehyde_keeps_strict_default_and_reports_allowed_parent() {
    let Some(ml) = library() else { return };
    let input = mmcif::read_mmcif(include_str!("data/partial_for_5lbm.cif")).unwrap();
    let mut strict = params(&input, "FOR");
    strict.allowed_unplaced_hydrogen_parent_sources.clear();
    let error = pipeline::run(input.clone(), &ml, &strict).err().unwrap();
    assert!(error.contains("H1") && error.contains("H2") && error.contains("FOR"));
    let output = pipeline::run(input.clone(), &ml, &params(&input, "FOR")).unwrap();
    heavy_unchanged(&input, &output);
    assert_eq!(output.hydrogen_uncertainties.len(), 2);
    assert!(output.hydrogen_uncertainties.iter().all(|site| site.parent_sources == [rows(&input, "FOR")[0]]
        && site.reason == HydrogenUncertaintyReason::InsufficientRestraints));
}

#[test]
fn observed_hydrogens_on_uncertain_parents_keep_source_metadata() {
    let Some(ml) = library() else { return };
    for (component, file, parent, name, element) in [
        ("CLA", include_str!("data/partial_cla_9l0k.cif"), "C6", "H61", "H"),
        ("CLA", include_str!("data/partial_cla_9l0k.cif"), "C6", "D61", "D"),
        ("FOR", include_str!("data/partial_for_5lbm.cif"), "C", "H1", "H"),
    ] {
        let mut input = mmcif::read_mmcif(file).unwrap();
        let p = find(&input, component, parent);
        let mut hydrogen = Atom::new(name, element, input.atom(p).xyz + v3(0.27, 0.31, 0.87));
        hydrogen.src = 999_001;
        hydrogen.charge = "1-".into();
        hydrogen.occ = 0.47;
        hydrogen.b = 21.3;
        hydrogen.segid = "OBS".into();
        hydrogen.uij = Some([0.1, 0.2, 0.3, 0.01, 0.02, 0.03]);
        hydrogen.hetero = true;
        input.models[p.model as usize].chains[p.chain as usize].residue_groups[p.rg as usize]
            .atom_groups[p.ag as usize].atoms.push(hydrogen.clone());
        for keep_existing_h in [false, true] {
            let mut policy = params(&input, component);
            policy.keep_existing_h = keep_existing_h;
            let output = pipeline::run(input.clone(), &ml, &policy).unwrap();
            heavy_unchanged(&input, &output);
            let observed = output.structure.atom_paths().into_iter()
                .filter(|p| output.structure.atom(*p).src == hydrogen.src).collect::<Vec<_>>();
            assert_eq!(observed.len(), 1);
            let actual = output.structure.atom(observed[0]);
            assert_eq!(actual.name, hydrogen.name);
            assert_eq!(actual.element, hydrogen.element);
            assert_eq!(actual.xyz, hydrogen.xyz);
            assert_eq!(actual.charge, hydrogen.charge);
            assert_eq!(actual.occ, hydrogen.occ);
            assert_eq!(actual.b, hydrogen.b);
            assert_eq!(actual.segid, hydrogen.segid);
            assert_eq!(actual.uij, hydrogen.uij);
            assert_eq!(actual.hetero, hydrogen.hetero);
        }
    }
}

#[test]
fn allowed_ligands_do_not_relax_missing_neighbor_protein_targets() {
    let Some(ml) = library() else { return };
    let mut input = mmcif::read_mmcif(include_str!("data/linked_9hj3.cif")).unwrap();
    input.links.clear();
    let keep = input.atom_paths().into_iter().filter(|p| st_component(&input, *p) != "CYS"
        || input.atom(*p).name.trim() == "SG").map(|p| input.atom(p).src).collect::<Vec<_>>();
    input.retain_atoms(|atom| keep.contains(&atom.src));
    let mut policy = params(&input, "Z41");
    policy.required_hydrogen_parent_sources = Some(rows(&input, "CYS"));
    let error = pipeline::run(input, &ml, &policy).err().unwrap();
    assert!(error.contains("HG") && error.contains("CYS"), "{error}");
}

fn st_component(st: &Structure, p: AtomPath) -> &str { st.atom_group(p).resname.trim() }

#[test]
fn invalid_parent_sources_and_modes_do_not_allow_uncertainty() {
    let Some(ml) = library() else { return };
    let input = mmcif::read_mmcif(include_str!("data/partial_cla_9l0k.cif")).unwrap();
    let row = rows(&input, "CLA")[0];
    for allowed in [vec![row, row], vec![Atom::NEW], vec![999_999]] {
        let mut policy = params(&input, "CLA");
        policy.allowed_unplaced_hydrogen_parent_sources = allowed;
        assert!(pipeline::run(input.clone(), &ml, &policy).err().unwrap().contains("unique heavy source rows"));
    }
    for mode in 0..3 {
        let mut policy = params(&input, "CLA");
        match mode {
            0 => policy.compat = true,
            1 => policy.approach = pipeline::Approach::Optimize,
            _ => policy.stop_on_any_missing_hydrogen = false,
        }
        assert!(pipeline::run(input.clone(), &ml, &policy).err().unwrap().contains("strict fixed-mode"));
    }
}

#[test]
fn typed_placement_reports_missing_neighbors_that_old_logs_omit() {
    let Some(ml) = library() else { return };
    let mut input = mmcif::read_mmcif(include_str!("data/partial_cla_9l0k.cif")).unwrap();
    let c6 = input.atom(find(&input, "CLA", "C6")).src;
    let placed = hplace::place_hydrogens(&mut input, &ml, &hplace::HPlaceParams {
        neutron: false, n_terminal_charge: hplace::NTermCharge::ResidueOne,
        exclude_water: true, keep_existing_h: false, adp_scale: 1.0, compat: false, cell: None,
    });
    assert!(placed.site_labels_no_para.is_empty());
    assert_eq!(placed.unplaced_hydrogens.iter().filter(|site| site.parent_sources == [c6]).count(), 2);
    assert!(placed.unplaced_hydrogens.iter().any(|site| site.parent_sources.is_empty()));
}

#[test]
fn unresolved_source_bonds_restore_only_the_affected_observed_hydrogen() {
    let Some(ml) = library() else { return };
    let mut input = mmcif::read_mmcif(include_str!("data/linked_9hj3.cif")).unwrap();
    let required = rows(&input, "CYS");
    let sg = find(&input, "CYS", "SG");
    let source = input.atom(sg).src;
    let mut observed = Atom::new("HG", "H", input.atom(sg).xyz + v3(0.27, 0.31, 0.87));
    observed.src = 999_001;
    observed.occ = 0.63;
    observed.b = 37.2;
    input.models[0].chains[0].residue_groups[0].atom_groups[0].atoms.push(observed.clone());
    let mut other = input.models[0].chains[0].residue_groups[0].clone();
    other.resseq = "  19".into();
    other.link_to_previous = false;
    other.atom_groups[0].atoms.retain(|atom| atom.name.trim() == "SG");
    other.atom_groups[0].atoms[0].src = 100;
    other.atom_groups[0].atoms[0].xyz += v3(3.5, 0.0, 0.0);
    input.models[0].chains[0].residue_groups.push(other);
    let mut second = input.links[0].clone();
    second.atom2 = second.atom1.clone();
    second.atom2.resseq = "19".into();
    second.kind = "disulf".into();
    second.distance = Some(3.5);
    second.order = None;
    input.links.push(second);
    let policy = Params { stop_on_any_missing_hydrogen: true,
        required_hydrogen_parent_sources: Some(required), allow_unresolved_bond_targets: true,
        ..Params::default() };
    let output = pipeline::run(input.clone(), &ml, &policy).unwrap();
    assert!(output.bond_uncertainties.iter().any(|bond| bond.parent_source == source && bond.hydrogens_unknown));
    let p = output.structure.atom_paths().into_iter().find(|p| output.structure.atom(*p).src == observed.src).unwrap();
    let actual = output.structure.atom(p);
    assert_eq!(actual.xyz, observed.xyz);
    assert_eq!(actual.name, observed.name);
    assert_eq!(actual.occ, observed.occ);
    assert_eq!(actual.b, observed.b);
    let mut strict = policy;
    strict.allow_unresolved_bond_targets = false;
    let strict = pipeline::run(input, &ml, &strict).unwrap();
    assert!(strict.structure.atom_paths().into_iter().all(|p| strict.structure.atom(p).src != observed.src));
}

#[test]
fn unresolved_bond_guesses_leave_the_graph_before_optimization() {
    let Some(ml) = library() else { return };
    let mut input = mmcif::read_mmcif(include_str!("data/linked_9hj3.cif")).unwrap();
    let required = rows(&input, "CYS");
    input.retain_atoms(|atom| required.contains(&atom.src));
    let mut second = input.links.remove(0);
    second.atom2 = second.atom1.clone();
    second.atom2.resseq = "19".into();
    second.kind = "disulf".into();
    second.distance = Some(3.5);
    second.order = None;
    input.links = vec![second];
    let mut other = input.models[0].chains[0].residue_groups[0].clone();
    other.resseq = "  19".into();
    other.link_to_previous = false;
    other.atom_groups[0].atoms.retain(|atom| atom.name.trim() == "SG");
    other.atom_groups[0].atoms[0].src = 100;
    other.atom_groups[0].atoms[0].xyz += v3(3.5, 0.0, 0.0);
    input.models[0].chains[0].residue_groups.push(other);
    let policy = Params { stop_on_any_missing_hydrogen: true,
        required_hydrogen_parent_sources: Some(required), allow_unresolved_bond_targets: true,
        ..Params::default() };
    let output = pipeline::run(input.clone(), &ml, &policy).unwrap();
    assert!(output.bond_uncertainties.iter().any(|bond| bond.hydrogens_unknown));
    assert!(output.structure.atom_paths().into_iter().all(|p| output.structure.atom(p).name.trim() != "HG"));
    let mut strict = policy;
    strict.allow_unresolved_bond_targets = false;
    let strict = pipeline::run(input, &ml, &strict).unwrap();
    assert!(strict.structure.atom_paths().into_iter().any(|p| strict.structure.atom(p).name.trim() == "HG"));
}

#[test]
fn observed_terminal_aliases_keep_their_named_n_parent() {
    let Some(ml) = library() else { return };
    for (name, element) in [("H1", "H"), ("D1", "D")] {
        let mut input = mmcif::read_mmcif(include_str!("data/linked_9hj3.cif")).unwrap();
        let required = rows(&input, "CYS");
        input.retain_atoms(|atom| required.contains(&atom.src));
        let n = find(&input, "CYS", "N");
        let observed = Atom { src: 999_001, ..Atom::new(name, element, input.atom(n).xyz + v3(0.27, 0.31, 0.87)) };
        input.models[0].chains[0].residue_groups[0].atom_groups[0].atoms.push(observed.clone());
        let mut other = input.models[0].chains[0].residue_groups[0].clone();
        other.resseq = " 900".into();
        other.link_to_previous = false;
        other.atom_groups[0].resname = "HOH".into();
        other.atom_groups[0].atoms = vec![Atom { src: 100, ..Atom::new("O", "O", input.atom(n).xyz + v3(10.0, 0.0, 0.0)) }];
        input.models[0].chains[0].residue_groups.push(other);
        let mut link = input.links.remove(0);
        link.atom1.name = "N".into();
        link.atom2 = link.atom1.clone();
        link.atom2.name = "O".into();
        link.atom2.resname = "HOH".into();
        link.atom2.resseq = "900".into();
        link.distance = Some(10.0);
        input.links = vec![link];
        let policy = Params { stop_on_any_missing_hydrogen: true,
            required_hydrogen_parent_sources: Some(required), allow_unresolved_bond_targets: true, ..Params::default() };
        let output = pipeline::run(input, &ml, &policy).unwrap();
        let p = output.structure.atom_paths().into_iter().find(|p| output.structure.atom(*p).src == observed.src).unwrap();
        assert_eq!(output.structure.atom(p).name, observed.name);
        assert_eq!(output.structure.atom(p).element, observed.element);
        assert_eq!(output.structure.atom(p).xyz, observed.xyz);
    }
}
