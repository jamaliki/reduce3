#[path = "common/chem_data.rs"]
mod chem_data;

use reduce3::{hplace, mmcif, MonLib};
use reduce3::model::Structure;

fn library() -> Option<MonLib> {
    let root = chem_data::locate()?;
    Some(MonLib::load(&root).unwrap())
}

fn place(st: &mut Structure, library: &MonLib, compat: bool) -> hplace::Placed {
    hplace::place_hydrogens(st, library, &hplace::HPlaceParams {
        neutron: false, n_terminal_charge: hplace::NTermCharge::ResidueOne,
        exclude_water: true, keep_existing_h: false, adp_scale: 1.0, compat, cell: None,
    })
}

fn heavy_sites(st: &Structure) -> Vec<(String, String, String, String, reduce3::geom::Vec3, String)> {
    let mut sites: Vec<_> = st.atom_paths().into_iter().filter(|&p| !st.atom(p).is_hydrogen()).map(|p| (
        st.chain(p).id.clone(), st.residue_group(p).resseq.trim().into(),
        st.atom_group(p).resname.trim().into(), st.atom(p).name.trim().into(), st.atom(p).xyz, st.atom(p).charge.to_string(),
    )).collect();
    sites.sort_by(|a, b| (&a.0, &a.1, &a.2, &a.3).cmp(&(&b.0, &b.1, &b.2, &b.3)));
    sites
}

fn has(st: &Structure, chain: &str, seq: &str, name: &str) -> bool {
    st.atom_paths().into_iter().any(|p| st.chain(p).id == chain
        && st.residue_group(p).resseq.trim() == seq && st.atom(p).name.trim() == name)
}

#[test]
fn declared_two_sulfur_product_removes_both_cysteine_thiol_h() {
    let Some(library) = library() else { return };
    let mut st = mmcif::read_mmcif(include_str!("data/covalent_9b5a.cif")).unwrap();
    let original = heavy_sites(&st);
    let placed = place(&mut st, &library, false);
    assert!(placed.site_labels_no_para.is_empty(), "{:?}", placed.site_labels_no_para);
    assert_eq!(heavy_sites(&st), original);
    assert!(!has(&st, "A", "85", "HG"));
    assert!(!has(&st, "C", "639", "HG"));
    assert!(!has(&st, "B", "101", "H031") && !has(&st, "B", "101", "H032"));
    let paths = st.atom_paths();
    let c = paths.iter().position(|&p| st.atom_group(p).resname == "A1AIV" && st.atom(p).name.trim() == "C03").unwrap() as u32;
    let sulfur_bonds = placed.bonds.iter().filter(|&&(a, b, _)| {
        let other = if a == c { b } else if b == c { a } else { return false };
        st.atom(paths[other as usize]).name.trim() == "SG"
    }).count();
    assert_eq!(sulfur_bonds, 2);
}

#[test]
fn declared_partial_glycine_cap_links_ca_to_cysteine_sg() {
    let Some(library) = library() else { return };
    let mut st = mmcif::read_mmcif(include_str!("data/covalent_6wgn.cif")).unwrap();
    let original = heavy_sites(&st);
    let placed = place(&mut st, &library, false);
    assert_eq!(heavy_sites(&st), original);
    assert!(!has(&st, "G", "13", "HG"));
    assert!(!placed.site_labels_no_para.iter().any(|label| label.contains("CYS")), "{:?}", placed.site_labels_no_para);
    let paths = st.atom_paths();
    let cap = paths.iter().position(|&p| st.residue_group(p).resseq.trim() == "0" && st.atom(p).name.trim() == "CA").unwrap() as u32;
    let sulfur = paths.iter().position(|&p| st.residue_group(p).resseq.trim() == "13" && st.atom(p).name.trim() == "SG").unwrap() as u32;
    assert!(placed.bonds.iter().any(|&(a, b, _)| (a, b) == (cap.min(sulfur), cap.max(sulfur))));
}

#[test]
fn matching_wpn_product_places_both_methylene_hydrogens_from_source_angles() {
    let Some(library) = library() else { return };
    let mut st = mmcif::read_mmcif(include_str!("data/covalent_8qsc.cif")).unwrap();
    let original = heavy_sites(&st);
    let placed = place(&mut st, &library, false);
    assert!(placed.site_labels_no_para.is_empty(), "{:?}", placed.site_labels_no_para);
    assert_eq!(heavy_sites(&st), original);
    for chain in ["A", "J"] {
        assert!(!has(&st, chain, "38", "HG"));
        let paths = st.atom_paths();
        let site = |name| paths.iter().find(|&&p| st.chain(p).id == chain
            && st.atom_group(p).resname == "WPN" && st.atom(p).name.trim() == name).map(|&p| st.atom(p).xyz).unwrap();
        for name in ["H1", "H2"] {
            assert!((site(name).dist(site("C2")) - 0.97).abs() < 1e-12);
            assert!(site(name).as_array().iter().all(|v| v.is_finite()));
        }
    }
}

#[test]
fn wpn_without_a_valid_source_product_keeps_its_restraint_error() {
    let Some(library) = library() else { return };
    for invalid in 0..4 {
        let mut st = mmcif::read_mmcif(include_str!("data/covalent_8qsc.cif")).unwrap();
        match invalid {
            0 => st.links.clear(),
            1 => st.links[0].atom2.resseq = "999".into(),
            2 => st.links[0].order = Some("doub".into()),
            3 => { let original = st.links[0].clone(); st.links.push(original); }
            _ => unreachable!(),
        }
        let placed = place(&mut st, &library, false);
        assert!(placed.site_labels_no_para.iter().any(|label| label.contains("WPN") && label.contains(" J")), "{:?}", placed.site_labels_no_para);
    }
}

#[test]
fn cif_missing_order_is_absent_and_explicit_order_stays_available() {
    let source = include_str!("data/covalent_6wgn.cif");
    for value in ["?", ".", "sing", "doub", "unsupported"] {
        let text = source.replace("sing", value);
        let st = mmcif::read_mmcif(&text).unwrap();
        let link = st.links.iter().find(|l| l.atom1.name.trim() == "CA").unwrap();
        assert_eq!(link.order.as_deref(), if value == "?" || value == "." { None } else { Some(value) });
    }
}
