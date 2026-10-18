use reduce3::{cif, hplace, mmcif, pipeline, MonLib, Params};
use reduce3::model::{AtomPath, Structure};
use reduce3::geom::v3;

fn library() -> MonLib {
    let root = MonLib::locate(None).expect("These source tests require REDUCE3_CHEM_DATA");
    MonLib::load(&root).unwrap()
}

fn fixtures() -> Vec<Structure> {
    let document = cif::parse(include_str!("data/symmetric_covalent_anchors.cif"));
    assert_eq!(document.blocks.len(), 4);
    document.blocks.iter().map(|block| mmcif::structure_from_cif(block).unwrap()).collect()
}

fn atom(st: &Structure, component: &str, name: &str) -> AtomPath {
    st.atom_paths().into_iter().find(|&p| st.atom_group(p).resname == component
        && st.atom(p).name.trim() == name).unwrap()
}

fn heavy(st: &Structure) -> Vec<(u32, String, [u32; 3], u64, String)> {
    let mut rows: Vec<_> = st.atom_paths().into_iter().filter(|&p| !st.atom(p).is_hydrogen()).map(|p| {
        let a = st.atom(p);
        (a.src, a.name.to_string(), a.xyz.as_array().map(|v| (v as f32).to_bits()), a.b.to_bits(), a.charge.to_string())
    }).collect();
    rows.sort_by_key(|row| row.0);
    rows
}

fn place(st: &mut Structure, library: &MonLib, compat: bool) -> hplace::Placed {
    hplace::place_hydrogens(st, library, &hplace::HPlaceParams {
        neutron: false, n_terminal_charge: hplace::NTermCharge::ResidueOne,
        exclude_water: true, keep_existing_h: false, adp_scale: 1.2, compat, cell: None,
    })
}

fn names<'a>(st: &'a Structure, component: &str) -> Vec<&'a str> {
    st.atom_paths().into_iter().filter(|&p| st.atom_group(p).resname == component)
        .map(|p| st.atom(p).name.trim()).collect()
}

#[test]
fn exact_source_links_keep_names_heavy_bits_charges_and_hydrogen_states() {
    let library = library();
    for input in fixtures() {
        assert_eq!(input.models.len(), 1);
        let model_id = input.models[0].id.clone();
        assert_eq!(input.links.len(), 1);
        let original = heavy(&input);
        let ring: Vec<_> = ["CD1", "CD2", "CE1", "CE2"].into_iter()
            .map(|name| (name, input.atom(atom(&input, "TYR", name)).xyz)).collect();
        let link = input.links[0].clone();
        let output = pipeline::run(input, &library, &Params {
            model_id: Some(1), stop_on_any_missing_hydrogen: true, ..Params::default()
        }).unwrap().structure;
        assert_eq!(heavy(&output), original);
        for (name, xyz) in ring { assert_eq!(output.atom(atom(&output, "TYR", name)).xyz, xyz); }
        assert_eq!(output.models[0].id, model_id);
        assert_eq!(output.links.len(), 1);
        assert_eq!(output.links[0].atom1, link.atom1);
        assert_eq!(output.links[0].atom2, link.atom2);
        assert_eq!(output.links[0].distance, link.distance);
        assert_eq!(output.links[0].order, link.order);
        assert_eq!(output.links[0].kind, link.kind);
        let cys = names(&output, "CYS");
        let tyr = names(&output, "TYR");
        assert!(!cys.contains(&"HG"));
        assert!(!tyr.contains(&"HE2"));
        for name in ["HD1", "HD2", "HE1", "HH"] { assert!(tyr.contains(&name), "{name}"); }
        assert!(output.atom_paths().into_iter().all(|p| output.atom(p).xyz.as_array().iter().all(|v| v.is_finite())));
    }
}

#[test]
fn ordinary_unlinked_tyr_still_normalizes_its_symmetric_names() {
    let library = library();
    let mut input = fixtures().remove(0);
    let (ce1, ce2) = (input.atom(atom(&input, "TYR", "CE1")).xyz, input.atom(atom(&input, "TYR", "CE2")).xyz);
    input.links.clear();
    let placed = place(&mut input, &library, false);
    assert!(placed.site_labels_no_para.is_empty());
    assert_eq!(input.atom(atom(&input, "TYR", "CE1")).xyz, ce2);
    assert_eq!(input.atom(atom(&input, "TYR", "CE2")).xyz, ce1);
}

#[test]
fn compat_keeps_its_prior_symmetric_name_behavior() {
    let library = library();
    let mut input = fixtures().remove(0);
    let original_ce1 = input.atom(atom(&input, "TYR", "CE1")).xyz;
    let original_ce2 = input.atom(atom(&input, "TYR", "CE2")).xyz;
    place(&mut input, &library, true);
    assert_eq!(input.atom(atom(&input, "TYR", "CE1")).xyz, original_ce2);
    assert_eq!(input.atom(atom(&input, "TYR", "CE2")).xyz, original_ce1);
}

#[test]
fn a_backbone_link_does_not_protect_symmetric_sidechain_names() {
    let library = library();
    let mut input = fixtures().remove(0);
    let original_ce1 = input.atom(atom(&input, "TYR", "CE1")).xyz;
    input.links[0].atom2.name = "N".into();
    place(&mut input, &library, false);
    assert_eq!(input.atom(atom(&input, "TYR", "CE2")).xyz, original_ce1);
}

#[test]
fn stale_geometry_keeps_anchor_names_without_creating_a_bond() {
    let library = library();
    let mut input = fixtures().remove(0);
    for p in input.atom_paths() {
        if input.atom_group(p).resname == "CYS" { input.atom_mut(p).xyz += v3(30., 0., 0.); }
    }
    let original = heavy(&input);
    let placed = place(&mut input, &library, false);
    assert!(placed.site_labels_no_para.is_empty());
    assert_eq!(heavy(&input), original);
    assert!(names(&input, "CYS").contains(&"HG"));
    assert!(names(&input, "TYR").contains(&"HE2"));
}

#[test]
fn absent_ambiguous_and_incompatible_endpoints_do_not_choose_an_anchor() {
    let library = library();
    for invalid in 0..5 {
        let mut input = fixtures().remove(0);
        let original_ce1 = input.atom(atom(&input, "TYR", "CE1")).xyz;
        match invalid {
            0 => input.links[0].atom2.resseq = "158".into(),
            1 => input.links[0].atom1.chain = "X".into(),
            2 => {
                let p = atom(&input, "TYR", "CE2");
                let copy = input.atom_group(p).clone();
                input.models[0].chains[p.chain as usize].residue_groups[p.rg as usize].atom_groups.push(copy);
            }
            3 => {
                for group in &mut input.models[0].chains[0].residue_groups {
                    for ag in &mut group.atom_groups { ag.altloc = if ag.resname == "CYS" { "A" } else { "B" }.into(); }
                }
            }
            4 => input.links[0].kind = "metalc".into(),
            _ => unreachable!(),
        }
        place(&mut input, &library, false);
        assert_eq!(input.atom(atom(&input, "TYR", "CE2")).xyz, original_ce1, "case {invalid}");
    }
}

#[test]
fn repeated_annotations_do_not_change_the_linked_hydrogen_states() {
    let library = library();
    let mut input = fixtures().remove(0);
    input.links.push(input.links[0].clone());
    let original = heavy(&input);
    let placed = place(&mut input, &library, false);
    assert!(placed.site_labels_no_para.is_empty());
    assert_eq!(heavy(&input), original);
    assert!(!names(&input, "CYS").contains(&"HG"));
    assert!(!names(&input, "TYR").contains(&"HE2"));
}
