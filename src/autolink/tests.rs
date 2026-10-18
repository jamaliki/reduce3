use super::*;
use crate::monlib::{Comp, CompAtom, CompBond, LinkBond};
use std::sync::atomic::{AtomicUsize, Ordering};

fn library() -> MonLib {
    static NEXT: AtomicUsize = AtomicUsize::new(0);
    let root = std::env::temp_dir().join(format!("reduce3-autolink-tests-{}-{}", std::process::id(), NEXT.fetch_add(1, Ordering::Relaxed)));
    std::fs::create_dir_all(root.join("mon_lib/list")).unwrap();
    std::fs::create_dir_all(root.join("geostd")).unwrap();
    for file in ["mon_lib/list/mon_lib_list.cif", "geostd/ener_lib.cif"] {
        std::fs::write(root.join(file), "data_empty\n").unwrap();
    }
    let library = MonLib::load(&root).unwrap();
    std::fs::remove_dir_all(&root).unwrap();
    library
}

fn carbon(orders: [&str; 2]) -> (Structure, FlatAtoms, Interp) {
    let st = crate::mmcif::read_mmcif(concat!(
        "data_test\nloop_\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_alt_id\n",
        "_atom_site.auth_comp_id\n_atom_site.auth_asym_id\n_atom_site.auth_seq_id\n",
        "_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n_atom_site.occupancy\n_atom_site.B_iso_or_equiv\n",
        "C CHB . BCL A 101 0 0 0 1 10\nC C1 . BCL A 101 1.4 0 0 1 10\n",
        "C C2 . BCL A 101 0 1.4 0 1 10\nH HHB . BCL A 101 -1 -1 0 1 10\n",
        "N N . HOA B 201 0 0 1.5 1 10\n",
    )).unwrap();
    let flat = FlatAtoms::from_structure(&st);
    let mut comp = Comp { id: "BCL".into(), ..Comp::default() };
    for (id, element) in [("CHB", "C"), ("C1", "C"), ("C2", "C"), ("HHB", "H")] {
        comp.atoms.push(CompAtom { id: id.into(), type_symbol: element.into(), type_energy: None });
    }
    let mut it = Interp::empty(flat.pos.len());
    for (i, order) in [orders[0], orders[1], "single"].into_iter().enumerate() {
        comp.bonds.push(CompBond { a1: "CHB".into(), a2: flat.name[i + 1].trim().into(), type_: order.into(), ..Default::default() });
        it.add_bond(0, (i + 1) as u32, 1.4, 0);
    }
    let comp = Arc::new(comp);
    for i in 0..4 {
        it.dict[i] = crate::interp::AtomDictInfo { comp: Some(comp.clone()), dict_id: Some(flat.name[i].trim().into()), instance: 0 };
    }
    (st, flat, it)
}

fn initial(it: &Interp) -> Initial {
    Initial::from_lists(&lists(it))
}

fn lists(it: &Interp) -> Vec<Vec<u32>> {
    let mut result = vec![Vec::new(); it.dict.len()];
    for bond in &it.bonds {
        result[bond.i as usize].push(bond.j);
        result[bond.j as usize].push(bond.i);
    }
    result
}

fn apply(st: &Structure, flat: &FlatAtoms, it: &mut Interp, compat: bool) {
    auto_link(it, st, flat, &library(), &AutoLinkParams { compat }, &mut String::new());
}

fn edit_source(it: &mut Interp, edit: impl FnOnce(&mut Comp)) {
    let mut comp = (**it.dict[0].comp.as_ref().unwrap()).clone();
    edit(&mut comp);
    let comp = Arc::new(comp);
    for info in &mut it.dict[..4] { info.comp = Some(comp.clone()); }
}

fn declare(st: &mut Structure, flat: &FlatAtoms, first: u32, second: u32) {
    let label = |i: u32| {
        let p = flat.path[i as usize];
        AtomLabel { name: flat.name[i as usize].trim().into(), resname: st.atom_group(p).resname.clone(),
            chain: st.chain(p).id.clone(), resseq: st.residue_group(p).resseq.trim().into(),
            altloc: flat.altloc[i as usize].to_string(), icode: st.residue_group(p).icode.clone() }
    };
    st.links.push(LinkRecord { atom1: label(first), atom2: label(second), distance: None, order: None, kind: "covale".into() });
}

#[test]
fn complete_source_carbon_rejects_custom_contacts_only_in_fixed_mode() {
    for orders in [["single", "double"], ["deloc", "deloc"], ["aromatic", "aromatic"]] {
        for compat in [false, true] {
            let (st, flat, mut it) = carbon(orders);
            let original = flat.pos.clone();
            apply(&st, &flat, &mut it, compat);
            assert_eq!(it.has_bond(0, 4), compat);
            assert_eq!(flat.pos, original);
        }
    }
}

#[test]
fn missing_unknown_or_ambiguous_source_bonds_do_not_prove_saturation() {
    for invalid in 0..6 {
        let (st, flat, mut it) = carbon(["single", "double"]);
        match invalid {
            0 => { it.bonds.retain(|bond| bond.j != 2); }
            1 => edit_source(&mut it, |comp| comp.bonds[0].type_ = "unknown".into()),
            2 => edit_source(&mut it, |comp| comp.bonds.push(comp.bonds[0].clone())),
            3 => it.dict[2].dict_id = None,
            4 => it.dict[2].comp = None,
            5 => it.dict[0].comp = None,
            _ => unreachable!(),
        }
        let ml = library();
        let params = AutoLinkParams { compat: false };
        let cx = Ctx::new(&st, &flat, &ml, &params);
        assert!(cx.carbon_valence(&it, &initial(&it), 0, 4).is_none_or(|value| value < 8));
        apply(&st, &flat, &mut it, false);
        assert!(it.has_bond(0, 4), "case {invalid}");
    }
}

#[test]
fn compatible_alternates_count_once_and_same_alt_duplicates_are_ambiguous() {
    let (mut st, _, mut it) = carbon(["single", "double"]);
    let mut copy = st.models[0].chains[0].residue_groups[0].atom_groups[0].clone();
    copy.atoms = vec![copy.atoms[1].clone()];
    copy.altloc = "B".into();
    st.models[0].chains[0].residue_groups[0].atom_groups.push(copy);
    let flat = FlatAtoms::from_structure(&st);
    it.dict.insert(4, it.dict[1].clone());
    it.add_bond(0, 4, 1.4, 0);
    let params = AutoLinkParams { compat: false };
    let ml = library();
    let cx = Ctx::new(&st, &flat, &ml, &params);
    assert_eq!(cx.carbon_valence(&it, &initial(&it), 0, 5), Some(8));
    let mut neighbors = lists(&it);
    neighbors[0].push(1);
    assert_eq!(cx.carbon_valence(&it, &Initial::from_lists(&neighbors), 0, 5), None);
}

#[test]
fn exact_single_covalent_record_bypasses_saturation() {
    for order in [None, Some("sing"), Some("single")] {
        let (mut st, flat, mut it) = carbon(["deloc", "deloc"]);
        declare(&mut st, &flat, 0, 4);
        st.links[0].order = order.map(String::from);
        apply(&st, &flat, &mut it, false);
        assert!(it.has_bond(0, 4));
    }
    for invalid in 0..5 {
        let (mut st, flat, mut it) = carbon(["deloc", "deloc"]);
        declare(&mut st, &flat, 0, 4);
        match invalid {
            0 => st.links.push(st.links[0].clone()),
            1 => st.links[0].order = Some("unknown".into()),
            2 => st.links[0].order = Some("double".into()),
            3 => st.links[0].atom2.chain = "X".into(),
            4 => st.links[0].kind = "metalc".into(),
            _ => unreachable!(),
        }
        apply(&st, &flat, &mut it, false);
        assert!(!it.has_bond(0, 4), "case {invalid}");
    }
}

fn metal() -> (Structure, FlatAtoms, Interp) {
    let st = crate::mmcif::read_mmcif(concat!(
        "data_test\nloop_\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_alt_id\n",
        "_atom_site.auth_comp_id\n_atom_site.auth_asym_id\n_atom_site.auth_seq_id\n",
        "_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n_atom_site.occupancy\n_atom_site.B_iso_or_equiv\n",
        "N N . HOA A 610 0 0 0 1 10\nO O . HOA A 610 -1.4 0 0 1 10\n",
        "FE FE . HEC B 601 1.936 0 0 1 10\nC C . HEC B 601 3.4 0 0 1 10\n",
    )).unwrap();
    let flat = FlatAtoms::from_structure(&st);
    let mut it = Interp::empty(4);
    it.add_bond(0, 1, 1.4, 0);
    it.add_bond(2, 3, 1.5, 0);
    (st, flat, it)
}

#[test]
fn embedded_metal_contacts_are_not_custom_covalent_bonds() {
    for compat in [false, true] {
        let (st, flat, mut it) = metal();
        apply(&st, &flat, &mut it, compat);
        assert_eq!(it.has_bond(0, 2), compat);
    }
    let (mut st, flat, mut it) = metal();
    declare(&mut st, &flat, 0, 2);
    apply(&st, &flat, &mut it, false);
    assert!(it.has_bond(0, 2));
}

#[test]
fn single_atom_metals_keep_the_existing_exclusion() {
    for compat in [false, true] {
        let (mut st, _, mut it) = metal();
        st.models[0].chains[1].residue_groups[0].atom_groups[0].atoms.pop();
        let flat = FlatAtoms::from_structure(&st);
        it.dict.pop();
        it.bonds.retain(|bond| bond.j != 3);
        apply(&st, &flat, &mut it, compat);
        assert!(!it.has_bond(0, 2));
    }
}

#[test]
fn only_histidine_ring_nitrogens_get_embedded_metal_coordination() {
    for name in ["ND1", "NE2", "N", "CA"] {
        for compat in [false, true] {
            let (mut st, _, mut it) = metal();
            let group = &mut st.models[0].chains[0].residue_groups[0].atom_groups[0];
            group.resname = "HIS".into();
            group.atoms[0].name = name.into();
            let flat = FlatAtoms::from_structure(&st);
            apply(&st, &flat, &mut it, compat);
            let ring = matches!(name, "ND1" | "NE2");
            assert_eq!(it.has_bond(0, 2), ring || compat);
            if it.has_bond(0, 2) {
                assert_eq!(it.bonds.last().unwrap().origin, if ring && !compat { ORIGIN_METAL } else { ORIGIN_MISC });
            }
        }
    }
}

#[test]
fn incompatible_alternates_and_other_models_never_link() {
    for alternate in [false, true] {
        let (mut st, _, mut it) = carbon(["single", "single"]);
        if alternate {
            st.models[0].chains[0].residue_groups[0].atom_groups[0].altloc = "A".into();
            st.models[0].chains[1].residue_groups[0].atom_groups[0].altloc = "B".into();
        } else {
            let mut model = st.models[0].clone();
            model.id = "2".into();
            model.chains = st.models[0].chains.split_off(1);
            st.models.push(model);
        }
        let flat = FlatAtoms::from_structure(&st);
        apply(&st, &flat, &mut it, false);
        assert!(!it.has_bond(0, 4));
    }
}

#[test]
fn named_dictionary_links_bypass_the_custom_guards() {
    let (mut st, _, mut it) = carbon(["deloc", "deloc"]);
    st.models[0].chains[0].residue_groups[0].atom_groups[0].resname = "MPR".into();
    let other = &mut st.models[0].chains[1].residue_groups[0].atom_groups[0];
    other.resname = "CYS".into();
    other.atoms[0].name = " SG ".into();
    other.atoms[0].element = "S".into();
    let flat = FlatAtoms::from_structure(&st);
    let mut library = library();
    let id = "CYS-MPR";
    library.link_index.insert(id.into(), library.links.len());
    library.links.push(ChemLink { id: id.into(), bonds: vec![LinkBond {
        c1: 1, a1: "SG".into(), c2: 2, a2: "CHB".into(), value_dist: Some(1.8), esd: Some(0.02),
    }], ..Default::default() });
    auto_link(&mut it, &st, &flat, &library, &AutoLinkParams { compat: false }, &mut String::new());
    assert!(it.has_bond(0, 4));
    assert_eq!(it.bonds.last().unwrap().origin, link_origin(id).unwrap());
}

#[test]
fn named_iron_cysteine_link_keeps_its_library_origin() {
    let (mut st, _, mut it) = metal();
    let sulfur = &mut st.models[0].chains[0].residue_groups[0].atom_groups[0];
    sulfur.resname = "CYS".into();
    sulfur.atoms[0].name = " SG ".into();
    sulfur.atoms[0].element = "S".into();
    st.models[0].chains[1].residue_groups[0].atom_groups[0].resname = "FE".into();
    let flat = FlatAtoms::from_structure(&st);
    let mut library = library();
    let id = "FE-CYS";
    library.link_index.insert(id.into(), library.links.len());
    library.links.push(ChemLink { id: id.into(), bonds: vec![LinkBond {
        c1: 1, a1: "FE".into(), c2: 2, a2: "SG".into(), value_dist: Some(2.3), esd: Some(0.02),
    }], ..Default::default() });
    auto_link(&mut it, &st, &flat, &library, &AutoLinkParams { compat: false }, &mut String::new());
    assert!(it.has_bond(0, 2));
    assert_eq!(it.bonds.last().unwrap().origin, link_origin(id).unwrap());
}

#[test]
fn rejected_carbon_contact_does_not_consume_the_other_endpoint() {
    let (mut st, _, mut it) = carbon(["deloc", "deloc"]);
    let mut chain = st.models[0].chains[1].clone();
    chain.id = "C".into();
    let atom = &mut chain.residue_groups[0].atom_groups[0].atoms[0];
    atom.name = " O  ".into();
    atom.element = "O".into();
    atom.xyz.z = 3.2;
    st.models[0].chains.push(chain);
    let flat = FlatAtoms::from_structure(&st);
    it.dict.push(Default::default());
    apply(&st, &flat, &mut it, false);
    assert!(!it.has_bond(0, 4));
    assert!(it.has_bond(4, 5));
}
