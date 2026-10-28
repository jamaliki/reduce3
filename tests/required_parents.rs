#[path = "common/chem_data.rs"]
mod chem_data;

use reduce3::{mmcif, pipeline, MonLib, Params};
use reduce3::model::{Atom, Structure};

fn library() -> Option<MonLib> {
    let Some(root) = chem_data::locate() else { return None };
    Some(MonLib::load(&root).unwrap())
}

fn fixture() -> Structure {
    let mut st = mmcif::read_mmcif(include_str!("data/linked_9hj3.cif")).unwrap();
    st.links.clear();
    st
}

fn rows(st: &Structure, component: &str) -> Vec<u32> {
    st.atom_paths().into_iter().filter_map(|p| {
        let atom = st.atom(p);
        (st.atom_group(p).resname == component && !atom.is_hydrogen()).then_some(atom.src)
    }).collect()
}

#[test]
fn only_context_unique_parents_allow_unplaced_hydrogens() {
    let Some(library) = library() else { return };
    let input = fixture();
    let output = pipeline::run(input.clone(), &library, &Params {
        stop_on_any_missing_hydrogen: true,
        required_hydrogen_parent_sources: Some(rows(&input, "CYS")), ..Params::default()
    }).unwrap();
    assert!(output.description.contains("Removed unplaced context hydrogens")
        && output.description.contains("H351") && output.description.contains("H352"));
    for p in input.atom_paths() {
        let source = input.atom(p).src;
        let after = output.structure.atom_paths().into_iter().find(|p| output.structure.atom(*p).src == source).unwrap();
        assert!(input.atom(p).xyz.dist(output.structure.atom(after).xyz) < 1e-10);
    }
    let control = pipeline::run(input.clone(), &library, &Params::default()).unwrap();
    assert_eq!(mmcif::write_mmcif(&output.structure), mmcif::write_mmcif(&control.structure));
    let error = pipeline::run(input.clone(), &library, &Params {
        stop_on_any_missing_hydrogen: true,
        required_hydrogen_parent_sources: Some(rows(&input, "Z41")), ..Params::default()
    }).err().unwrap();
    assert!(error.contains("H351") && error.contains("H352"));

    let mut unknown = input.clone();
    let parent = unknown.atom_paths().into_iter().find(|p| unknown.atom(*p).name.trim() == "C35").unwrap();
    unknown.atom_mut(parent).src = Atom::NEW;
    let error = pipeline::run(unknown, &library, &Params {
        stop_on_any_missing_hydrogen: true,
        required_hydrogen_parent_sources: Some(rows(&input, "CYS")), ..Params::default()
    }).err().unwrap();
    assert!(error.contains("H351") && error.contains("H352"));
}

#[test]
fn stale_duplicate_hydrogen_or_unsupported_parent_rows_fail() {
    let Some(library) = library() else { return };
    let input = fixture();
    let row = rows(&input, "CYS")[0];
    for rows in [vec![], vec![row, row], vec![Atom::NEW]] {
        let error = pipeline::run(input.clone(), &library, &Params {
            stop_on_any_missing_hydrogen: true, required_hydrogen_parent_sources: Some(rows), ..Params::default()
        }).err().unwrap();
        assert!(error.contains("unique heavy source rows"));
    }
    for incompatible in 0..3 {
        let mut params = Params { stop_on_any_missing_hydrogen: true,
            required_hydrogen_parent_sources: Some(vec![row]), ..Params::default() };
        match incompatible {
            0 => params.compat = true,
            1 => params.stop_on_any_missing_hydrogen = false,
            2 => params.approach = pipeline::Approach::Optimize,
            _ => unreachable!(),
        }
        assert!(pipeline::run(input.clone(), &library, &params).err().unwrap().contains("strict fixed-mode"));
    }
}
