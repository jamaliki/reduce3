use reduce3::{hplace::NTermCharge, mmcif, pipeline, MonLib, Params};
use reduce3::model::{Atom, AtomPath, Structure};
use reduce3::geom::v3;

fn library() -> Option<MonLib> {
    let root = MonLib::locate(None)?;
    Some(MonLib::load(&root).unwrap())
}

fn fixture() -> Structure {
    mmcif::read_mmcif(include_str!("data/terminal_0g7_1hut.cif")).unwrap()
}

fn atom(structure: &Structure, chain: &str, name: &str) -> AtomPath {
    structure.atom_paths().into_iter().find(|&p| structure.chain(p).id == chain
        && structure.atom(p).name.trim() == name).unwrap()
}

fn sites(structure: &Structure) -> Vec<(String, String, [u32; 3], String)> {
    let mut sites: Vec<_> = structure.atom_paths().into_iter().filter(|&p| !structure.atom(p).is_hydrogen()).map(|p| (
        structure.chain(p).id.clone(), structure.atom(p).name.to_string(),
        structure.atom(p).xyz.as_array().map(|v| (v as f32).to_bits()), structure.atom(p).charge.to_string(),
    )).collect();
    sites.sort_by(|a, b| (&a.0, &a.1).cmp(&(&b.0, &b.1)));
    sites
}

fn hydrogens(structure: &Structure, chain: &str, parent: &str) -> Vec<String> {
    let center = structure.atom(atom(structure, chain, parent)).xyz;
    structure.atom_paths().into_iter().filter(|&p| structure.chain(p).id == chain
        && structure.atom(p).is_hydrogen() && structure.atom(p).xyz.dist(center) < 1.25)
        .map(|p| structure.atom(p).name.trim().to_owned()).collect()
}

#[test]
fn charged_modified_terminus_preserves_internal_amide_and_heavy_atoms() {
    let Some(library) = library() else { return };
    let input = fixture();
    let original = sites(&input);
    let output = pipeline::run(input, &library, &Params {
        stop_on_any_missing_hydrogen: true, ..Params::default()
    }).unwrap().structure;
    assert_eq!(sites(&output), original);
    assert_eq!(hydrogens(&output, "H", "N"), ["HNT1", "HNT2", "HNT3"]);
    assert_eq!(hydrogens(&output, "H", "N2"), ["HN2"]);
    assert!(hydrogens(&output, "H", "C3").is_empty());
    let mut seen = std::collections::HashSet::new();
    assert!(output.atom_paths().into_iter().all(|p| seen.insert(output.atom(p).name.trim())));
    let h = output.atom(atom(&output, "H", "HN2")).xyz;
    let n = output.atom(atom(&output, "H", "N2")).xyz;
    assert!((h.dist(n) - 0.86).abs() < 1e-12);
    let a = output.atom(atom(&output, "H", "C1")).xyz - n;
    let b = output.atom(atom(&output, "H", "CA2")).xyz - n;
    assert!((h - n).dot(a.cross(b)).abs() < 1e-10);
}

#[test]
fn shared_dictionary_does_not_charge_a_neutral_copy() {
    let Some(library) = library() else { return };
    let mut input = fixture();
    let mut neutral = input.models[0].chains[0].clone();
    neutral.id = "J".into();
    neutral.residue_groups[0].resseq = "701".into();
    for group in &mut neutral.residue_groups[0].atom_groups {
        for atom in &mut group.atoms { atom.xyz += v3(50., 0., 0.); atom.src = Atom::NEW; }
    }
    input.models[0].chains.push(neutral);
    let output = pipeline::run(input, &library, &Params {
        stop_on_any_missing_hydrogen: true, ..Params::default()
    }).unwrap().structure;
    assert_eq!(hydrogens(&output, "H", "N").len(), 3);
    assert_eq!(hydrogens(&output, "J", "N").len(), 1);
    assert_eq!(hydrogens(&output, "J", "N2"), ["HN2"]);
}

#[test]
fn modified_terminal_policy_and_compat_mode_keep_their_meaning() {
    let Some(library) = library() else { return };
    for (seq, policy, count) in [(1, NTermCharge::NoCharge, 1), (701, NTermCharge::ResidueOne, 1),
        (701, NTermCharge::FirstInChain, 3)] {
        let mut input = fixture();
        input.models[0].chains[0].residue_groups[0].resseq = seq.to_string();
        let output = pipeline::run(input, &library, &Params {
            n_terminal_charge: policy, stop_on_any_missing_hydrogen: true, ..Params::default()
        }).unwrap().structure;
        assert_eq!(hydrogens(&output, "H", "N").len(), count);
        assert_eq!(hydrogens(&output, "H", "N2"), ["HN2"]);
    }
    assert!(pipeline::run(fixture(), &library, &Params {
        compat: true, stop_on_any_missing_hydrogen: true, ..Params::default()
    }).is_err());
}

#[test]
fn canonical_terminal_names_remain_unchanged() {
    let Some(library) = library() else { return };
    let mut input = fixture();
    let group = &mut input.models[0].chains[0].residue_groups[0].atom_groups[0];
    group.resname = "GLY".into();
    group.atoms.retain(|atom| ["N", "CA", "C", "O"].contains(&atom.name.trim()));
    let original = sites(&input);
    let output = pipeline::run(input, &library, &Params {
        stop_on_any_missing_hydrogen: true, ..Params::default()
    }).unwrap().structure;
    assert_eq!(sites(&output), original);
    assert_eq!(hydrogens(&output, "H", "N"), ["H1", "H2", "H3"]);
}
