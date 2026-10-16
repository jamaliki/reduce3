use super::*;
use crate::monlib::{CompAtom, CompBond};

fn add_thioethers(it: &mut Interp, st: &Structure, flat: &FlatAtoms) {
    super::add_thioethers(it, st, flat, &|_, bond| order(&bond.type_));
}

fn heavy_valence(it: &Interp, flat: &FlatAtoms, parent: u32) -> Option<u8> {
    super::heavy_valence(it, flat, parent, &|_, bond| order(&bond.type_))
}

fn fixture() -> (Structure, FlatAtoms, Interp) {
    let mut st = crate::mmcif::read_mmcif(concat!(
        "data_test\nloop_\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_alt_id\n",
        "_atom_site.auth_comp_id\n_atom_site.auth_asym_id\n_atom_site.auth_seq_id\n",
        "_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n_atom_site.occupancy\n_atom_site.B_iso_or_equiv\n",
        "S SG . CYS A 85 1.8 0 0 1 10\nC CB . CYS A 85 3.6 0 0 1 10\n",
        "C C03 . A1AIV B 101 0 0 0 1 10\nC C02 . A1AIV B 101 0 1.5 0 1 10\nC C04 . A1AIV B 101 0 -1.5 0 1 10\n",
        "S SG . CYS C 639 -1.8 0 0 1 10\nC CB . CYS C 639 -3.6 0 0 1 10\n",
    )).unwrap();
    let label = |name: &str, comp: &str, chain: &str, seq: &str| AtomLabel {
        name: name.into(), resname: comp.into(), chain: chain.into(), resseq: seq.into(),
        altloc: String::new(), icode: String::new(),
    };
    for (chain, seq) in [("A", "85"), ("C", "639")] {
        st.links.push(LinkRecord { atom1: label("SG", "CYS", chain, seq),
            atom2: label("C03", "A1AIV", "B", "101"), order: Some("sing".into()),
            distance: Some(1.8), kind: "covale".into() });
    }
    let comp = |id: &str, atoms: &[(&str, &str)], pairs: &[(&str, &str)]| {
        Arc::new(Comp { id: id.into(),
            atoms: atoms.iter().map(|&(id, element)| CompAtom { id: id.into(), type_symbol: element.into(), type_energy: None }).collect(),
            bonds: pairs.iter().map(|&(a, b)| CompBond { a1: a.into(), a2: b.into(), type_: "single".into(), ..Default::default() }).collect(),
            ..Default::default() })
    };
    let sulfur = comp("CYS", &[("SG", "S"), ("CB", "C")], &[("SG", "CB")]);
    let carbon = comp("A1AIV", &[("C03", "C"), ("C02", "C"), ("C04", "C")], &[("C03", "C02"), ("C03", "C04")]);
    let flat = FlatAtoms::from_structure(&st);
    let mut it = Interp::empty(flat.pos.len());
    for (i, source) in [(0, &sulfur), (1, &sulfur), (2, &carbon), (3, &carbon), (4, &carbon), (5, &sulfur), (6, &sulfur)] {
        it.dict[i] = AtomDictInfo { comp: Some(source.clone()), dict_id: Some(flat.name[i].trim().into()) };
    }
    for (a, b) in [(0, 1), (2, 3), (2, 4), (5, 6)] { it.add_bond(a, b, 1.5, ORIGIN_COVALENT); }
    (st, flat, it)
}

#[test]
fn two_declared_thioethers_share_a_tetravalent_carbon() {
    let (st, flat, mut it) = fixture();
    let positions = flat.pos.clone();
    add_thioethers(&mut it, &st, &flat);
    assert!(it.has_bond(0, 2) && it.has_bond(5, 2));
    assert_eq!(heavy_valence(&it, &flat, 2), Some(8));
    assert_eq!(heavy_valence(&it, &flat, 0), Some(4));
    assert_eq!(flat.pos, positions);
    for (i, &p) in flat.path.iter().enumerate() { assert_eq!(st.atom(p).xyz, positions[i]); }
}

#[test]
fn declared_ca_sg_link_does_not_require_a_peptide_n_c_pair() {
    let (mut st, mut flat, mut it) = fixture();
    st.links.truncate(1);
    let p = flat.path[2];
    st.models[p.model as usize].chains[p.chain as usize].residue_groups[p.rg as usize].atom_groups[p.ag as usize].resname = "GLY".into();
    st.atom_mut(p).name = "CA".into();
    st.links[0].atom2.resname = "GLY".into();
    st.links[0].atom2.name = "CA".into();
    flat.name[2] = "CA".into();
    let dictionary = Arc::make_mut(it.dict[2].comp.as_mut().unwrap());
    dictionary.id = "GLY".into();
    dictionary.atoms[0].id = "CA".into();
    for b in &mut dictionary.bonds { b.a1 = "CA".into(); }
    it.dict[2].dict_id = Some("CA".into());
    add_thioethers(&mut it, &st, &flat);
    assert!(it.has_bond(0, 2));
}

#[test]
fn invalid_declared_bonds_do_not_enter_the_graph() {
    for invalid in 0..12 {
        let (mut st, mut flat, mut it) = fixture();
        st.links.truncate(1);
        match invalid {
            0 => st.links[0].atom2.chain = "X".into(),
            1 => st.links[0].atom2.icode = "A".into(),
            2 => { let copy = st.models[0].chains[0].clone(); st.models[0].chains.push(copy); flat = FlatAtoms::from_structure(&st); }
            3 => flat.pos[0] = Vec3::from_array([2.31, 0.0, 0.0]),
            4 => flat.pos[0] = Vec3::from_array([1.49, 0.0, 0.0]),
            5 => flat.element[2] = "N".into(),
            6 => st.links[0].order = Some("doub".into()),
            7 => st.links[0].order = Some("unknown".into()),
            8 => { flat.altloc[0] = "A".into(); flat.altloc[2] = "B".into(); }
            9 => { st.models.push(st.models[0].clone()); flat.path[2].model += 1; }
            10 => it.dict[2].comp = None,
            11 => Arc::make_mut(it.dict[2].comp.as_mut().unwrap()).bonds[0].type_ = "unknown".into(),
            _ => unreachable!(),
        }
        add_thioethers(&mut it, &st, &flat);
        assert!(!it.has_bond(0, 2), "case {invalid}");
    }
}

#[test]
fn repeated_annotations_do_not_use_extra_valence() {
    let (mut st, flat, mut it) = fixture();
    st.links.push(st.links[0].clone());
    st.links[0].order = None;
    add_thioethers(&mut it, &st, &flat);
    assert!(it.has_bond(0, 2) && it.has_bond(5, 2));
    assert_eq!(it.bonds.len(), 6);
}

#[test]
fn ambiguous_sulfur_or_overvalent_carbon_do_not_add_links() {
    for sulfur in [false, true] {
        let (mut st, mut flat, mut it) = fixture();
        if sulfur {
            let mut extra = st.links[0].clone();
            extra.atom2.name = "C02".into();
            flat.pos[3] = Vec3::from_array([0.0, 1.0, 0.0]);
            st.links.push(extra);
        } else {
            Arc::make_mut(it.dict[2].comp.as_mut().unwrap()).bonds[0].type_ = "double".into();
        }
        add_thioethers(&mut it, &st, &flat);
        assert!(!it.has_bond(0, 2));
        assert_eq!(it.has_bond(5, 2), sulfur);
    }
}
