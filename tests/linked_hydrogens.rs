use reduce3::geom::Vec3;
use reduce3::{hplace, mmcif, pipeline, MonLib, Params};
use reduce3::model::{Atom, Structure};

const PEPTIDE: &str = include_str!("data/linked_8zje.cif");
const ADDUCT: &str = include_str!("data/linked_9hj3.cif");

fn library() -> Option<MonLib> {
    let Some(root) = MonLib::locate(None) else {
        eprintln!("skipped: chem_data not found (set REDUCE3_CHEM_DATA)");
        return None;
    };
    Some(MonLib::load(&root).unwrap())
}

fn atoms<'a>(structure: &'a Structure, component: &str) -> Vec<&'a Atom> {
    structure.models.iter().flat_map(|m| &m.chains).flat_map(|c| &c.residue_groups)
        .flat_map(|r| &r.atom_groups).filter(|g| g.resname == component).flat_map(|g| &g.atoms).collect()
}

fn volume(center: Vec3, a: Vec3, b: Vec3, h: Vec3) -> f64 {
    (a - center).dot((b - center).cross(h - center))
}

#[test]
fn linked_modified_peptide_keeps_its_valid_amide_hydrogens() {
    let Some(library) = library() else { return };
    let out = pipeline::run(mmcif::read_mmcif(PEPTIDE).unwrap(), &library, &Params {
        stop_on_any_missing_hydrogen: true, ..Params::default()
    }).unwrap();
    assert!(!atoms(&out.structure, "HYP").iter().any(|a| a.name.trim() == "H"));
    assert_eq!(atoms(&out.structure, "XZA").iter().filter(|a| ["HN1", "HN2"].contains(&a.name.trim())).count(), 1);
    assert_eq!(atoms(&out.structure, "LEU").iter().filter(|a| a.name.trim() == "H").count(), 1);
    let mut input = mmcif::read_mmcif(PEPTIDE).unwrap();
    let placed = hplace::place_hydrogens(&mut input, &library, &hplace::HPlaceParams {
        neutron: false, n_terminal_charge: hplace::NTermCharge::ResidueOne,
        exclude_water: true, keep_existing_h: false, adp_scale: 1.2, compat: false, cell: None,
    });
    assert!(placed.site_labels_no_para.is_empty());
    assert!(placed.removed_on_links.iter().any(|(atom, reason)| atom.contains("HYP") && reason == "a tertiary amide"));
}

#[test]
fn explicit_cysteine_adduct_keeps_named_methylene_hydrogens() {
    let Some(library) = library() else { return };
    let input = mmcif::read_mmcif(ADDUCT).unwrap();
    assert_eq!(input.links.len(), 1);
    let out = pipeline::run(input, &library, &Params {
        stop_on_any_missing_hydrogen: true, ..Params::default()
    }).unwrap();
    let ligand = atoms(&out.structure, "Z41");
    let site = |name| ligand.iter().find(|a| a.name.trim() == name).unwrap().xyz;
    let sulfur = atoms(&out.structure, "CYS").iter().find(|a| a.name.trim() == "SG").unwrap().xyz;
    for name in ["H351", "H352"] {
        assert!((site(name).dist(site("C35")) - 0.970).abs() < 1e-10);
        assert!(site(name).as_array().iter().all(|x| x.is_finite()));
    }
    let ccd = library.ccd("Z41").unwrap();
    let reference = |name| Vec3::from_array(ccd.atoms.iter().find(|a| a.0 == name).unwrap().3.unwrap());
    let ideal = volume(reference("C35"), reference("C18"), reference("O5"), reference("H351"));
    let observed = volume(site("C35"), site("C18"), sulfur, site("H351"));
    assert_eq!(ideal.is_sign_positive(), observed.is_sign_positive());
    assert!(!atoms(&out.structure, "CYS").iter().any(|a| a.name.trim() == "HG"));
}

#[test]
fn absent_or_ambiguous_adduct_annotation_still_fails() {
    let Some(library) = library() else { return };
    for invalid in 0..4 {
        let mut input = mmcif::read_mmcif(ADDUCT).unwrap();
        match invalid {
            0 => input.links.clear(),
            1 => input.links.push(input.links[0].clone()),
            2 => input.links[0].atom2.resseq = "603".into(),
            3 => input.links[0].atom1.chain = "X".into(),
            _ => unreachable!(),
        }
        let error = pipeline::run(input, &library, &Params {
            stop_on_any_missing_hydrogen: true, ..Params::default()
        }).err().unwrap();
        assert!(error.contains("H351") && error.contains("H352"), "{error}");
    }
}

#[test]
fn connections_keep_auth_ids_and_reject_symmetry_or_metal_links() {
    let original = mmcif::read_mmcif(ADDUCT).unwrap();
    assert_eq!(original.links[0].atom1.chain, "G");
    assert_eq!(original.links[0].atom2.resseq, "602");
    assert_eq!(original.links[0].atom2.name, "C35");
    for (field, value) in [("ptnr1_symmetry", "2_555"), ("ptnr2_symmetry", "2_555"), ("conn_type_id", "metalc")] {
        let mut document = reduce3::cif::parse(ADDUCT);
        let table = document.blocks[0].categories.iter_mut().find(|c| c.name == "_struct_conn").unwrap();
        let column = table.col(field).unwrap();
        table.values[column] = value;
        assert!(mmcif::structure_from_cif(&document.blocks[0]).unwrap().links.is_empty());
    }
}

#[test]
fn compat_keeps_the_prior_missing_restraint_error() {
    let Some(library) = library() else { return };
    let error = pipeline::run(mmcif::read_mmcif(ADDUCT).unwrap(), &library, &Params {
        compat: true, stop_on_any_missing_hydrogen: true, ..Params::default()
    }).err().unwrap();
    assert!(error.contains("H351") && error.contains("H352"), "{error}");
}
