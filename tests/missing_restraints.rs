//! Scoped missing dictionaries keep heavy contact context and strict targets.

#[path = "common/chem_data.rs"]
mod chem_data;

use reduce3::geom::{v3, Vec3};
use reduce3::model::{Atom, AtomGroup, Chain, ResidueGroup, Structure};
use reduce3::{mmcif, pipeline, MonLib, Params};

const TYROSINE: &str = "data_tyr\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_alt_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_entity_id\n_atom_site.label_seq_id\n_atom_site.pdbx_PDB_ins_code\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n_atom_site.occupancy\n_atom_site.B_iso_or_equiv\n_atom_site.pdbx_formal_charge\n_atom_site.auth_seq_id\n_atom_site.auth_comp_id\n_atom_site.auth_asym_id\n_atom_site.auth_atom_id\n_atom_site.pdbx_PDB_model_num\n\
ATOM 45 N N . TYR A 1 7 ? 4.981 -6.717 10.785 1.00 12.74 ? 7 TYR A N 1\n\
ATOM 46 C CA . TYR A 1 7 ? 4.621 -7.832 11.650 1.00 13.75 ? 7 TYR A CA 1\n\
ATOM 47 C C . TYR A 1 7 ? 5.516 -9.043 11.418 1.00 16.26 ? 7 TYR A C 1\n\
ATOM 48 O O . TYR A 1 7 ? 5.211 -10.147 11.865 1.00 17.99 ? 7 TYR A O 1\n\
ATOM 49 C CB . TYR A 1 7 ? 3.154 -8.214 11.434 1.00 12.91 ? 7 TYR A CB 1\n\
ATOM 50 C CG . TYR A 1 7 ? 2.185 -7.141 11.878 1.00 10.04 ? 7 TYR A CG 1\n\
ATOM 51 C CD1 . TYR A 1 7 ? 1.745 -6.156 10.994 1.00 9.16 ? 7 TYR A CD1 1\n\
ATOM 52 C CD2 . TYR A 1 7 ? 1.731 -7.094 13.195 1.00 9.58 ? 7 TYR A CD2 1\n\
ATOM 53 C CE1 . TYR A 1 7 ? 0.875 -5.148 11.410 1.00 7.70 ? 7 TYR A CE1 1\n\
ATOM 54 C CE2 . TYR A 1 7 ? 0.866 -6.091 13.623 1.00 8.82 ? 7 TYR A CE2 1\n\
ATOM 55 C CZ . TYR A 1 7 ? 0.441 -5.122 12.728 1.00 7.68 ? 7 TYR A CZ 1\n\
ATOM 56 O OH . TYR A 1 7 ? -0.415 -4.135 13.154 1.00 8.87 ? 7 TYR A OH 1\n";

fn library() -> Option<MonLib> {
    let root = chem_data::locate()?;
    Some(MonLib::load(&root).unwrap())
}

fn find_atom<'a>(structure: &'a Structure, component: &str, name: &str) -> Option<&'a Atom> {
    structure.models.iter().flat_map(|m| &m.chains).flat_map(|c| &c.residue_groups)
        .flat_map(|r| &r.atom_groups).filter(|g| g.resname.trim() == component)
        .flat_map(|g| &g.atoms).find(|a| a.name.trim() == name)
}

fn atom<'a>(structure: &'a Structure, component: &str, name: &str) -> &'a Atom {
    find_atom(structure, component, name).unwrap()
}

fn add_context(structure: &mut Structure, component: &str, position: Vec3) {
    let mut carbon = Atom::new("C1", "C", position);
    carbon.hetero = true;
    structure.models[0].chains.push(Chain { id: component.into(), residue_groups: vec![ResidueGroup {
        resseq: " 900".into(), icode: " ".into(), link_to_previous: false,
        atom_groups: vec![AtomGroup { altloc: String::new(), resname: component.into(), atoms: vec![carbon] }],
    }] });
}

fn strict(allowed: &[&str]) -> Params {
    Params { stop_on_any_missing_hydrogen: true,
        allow_missing_restraints_for: allowed.iter().map(|s| s.to_string()).collect(), ..Params::default() }
}

#[test]
fn only_named_missing_dictionaries_are_allowed() {
    let Some(library) = library() else { return };
    let mut input = mmcif::read_mmcif(TYROSINE).unwrap();
    add_context(&mut input, "UNL", v3(20., 20., 20.));
    assert!(pipeline::run(input.clone(), &library, &strict(&[])).err().unwrap().contains("UNL"));
    assert!(pipeline::run(input.clone(), &library, &strict(&["ZZZ"])).err().unwrap().contains("UNL"));
    let output = pipeline::run(input.clone(), &library, &strict(&[" unl "])).unwrap();
    assert!(output.description.contains("Allowed missing restraints for components UNL"));
    assert_eq!(atom(&output.structure, "UNL", "C1").xyz, v3(20., 20., 20.));
    let control = pipeline::run(input.clone(), &library, &Params::default()).unwrap();
    assert_eq!(mmcif::write_mmcif(&output.structure), mmcif::write_mmcif(&control.structure));
    add_context(&mut input, "NO_SUCH_COMPONENT", v3(30., 30., 30.));
    let error = pipeline::run(input, &library, &strict(&["UNL"])).err().unwrap();
    assert!(error.contains("NO_SUCH_COMPONENT") && !error.contains("UNL"), "{error}");
}

#[test]
fn retained_unl_heavy_atoms_affect_neighbor_optimization() {
    let Some(library) = library() else { return };
    let mut input = mmcif::read_mmcif(TYROSINE).unwrap();
    let without = pipeline::run(input.clone(), &library, &strict(&[])).unwrap();
    let initial = atom(&without.structure, "TYR", "HH").xyz;
    let context = initial + v3(0.15, 0.1, 0.1);
    add_context(&mut input, "UNL", context);
    let with = pipeline::run(input.clone(), &library, &strict(&["UNL"])).unwrap();
    assert_eq!(atom(&with.structure, "UNL", "C1").xyz, context);
    assert!(find_atom(&with.structure, "TYR", "HH").is_none_or(|a| a.xyz.dist(initial) > 0.5));
    let reference = pipeline::run(input, &library, &Params::default()).unwrap();
    assert_eq!(mmcif::write_mmcif(&with.structure), mmcif::write_mmcif(&reference.structure));
}

#[test]
fn known_dictionary_outputs_do_not_change_when_listed() {
    let Some(library) = library() else { return };
    let input = mmcif::read_mmcif(TYROSINE).unwrap();
    let before = pipeline::run(input.clone(), &library, &strict(&[])).unwrap();
    let after = pipeline::run(input, &library, &strict(&["TYR", "UNL", "D8U", "HOH", "DOD"])).unwrap();
    assert_eq!(mmcif::write_mmcif(&before.structure), mmcif::write_mmcif(&after.structure));
    assert!(!after.description.contains("Allowed missing restraints"));
}

#[test]
fn allowed_missing_dictionary_does_not_allow_unparameterized_hydrogens() {
    let Some(library) = library() else { return };
    let mut input = mmcif::read_mmcif(include_str!("data/linked_9hj3.cif")).unwrap();
    input.links.clear();
    let error = pipeline::run(input, &library, &strict(&["Z41", "UNL"])).err().unwrap();
    assert!(error.contains("H351") && error.contains("H352"), "{error}");
}

#[test]
fn compat_and_invalid_component_ids_fail() {
    let Some(library) = library() else { return };
    for allowed in ["", "?", ".", "UNL ZZZ"] {
        assert!(pipeline::run(mmcif::read_mmcif(TYROSINE).unwrap(), &library, &strict(&[allowed]))
            .err().unwrap().contains("component identifiers"));
    }
    let mut params = strict(&["UNL"]);
    params.compat = true;
    assert!(pipeline::run(mmcif::read_mmcif(TYROSINE).unwrap(), &library, &params)
        .err().unwrap().contains("require fixed mode"));
}
