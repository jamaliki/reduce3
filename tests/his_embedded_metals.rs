use reduce3::{mmcif, pipeline, MonLib, Params};

#[test]
fn long_embedded_magnesium_contacts_keep_the_other_histidine_ring_h() {
    let library = MonLib::load(&MonLib::locate(None).expect("set REDUCE3_CHEM_DATA for the exact histidine test")).unwrap();
    let input = mmcif::read_mmcif(include_str!("data/his_embedded_7vzg.cif")).unwrap();
    let required = input.atom_paths().into_iter().filter(|&p| input.atom_group(p).resname == "HIS")
        .map(|p| input.atom(p).src).collect();
    let original: Vec<_> = input.atom_paths().into_iter().map(|p| (
        input.chain(p).id.clone(), input.residue_group(p).resseq.clone(), input.atom(p).name.clone(), input.atom(p).xyz,
    )).collect();
    let output = pipeline::run(input, &library, &Params {
        stop_on_any_missing_hydrogen: true, required_hydrogen_parent_sources: Some(required), ..Params::default()
    }).unwrap();
    let atoms = output.structure.atom_paths();
    for (chain, sequence) in [("A", "362"), ("a", "775")] {
        let has = |name| atoms.iter().any(|&p| output.structure.chain(p).id == chain
            && output.structure.residue_group(p).resseq.trim() == sequence && output.structure.atom(p).name.trim() == name);
        assert!(has("HD1"));
        assert!(!has("HE2"));
    }
    for (chain, sequence, name, xyz) in original {
        let p = atoms.iter().find(|&&p| output.structure.chain(p).id == chain
            && output.structure.residue_group(p).resseq == sequence && output.structure.atom(p).name == name).unwrap();
        assert_eq!(output.structure.atom(*p).xyz, xyz);
    }
}
