use super::*;

fn add_replacement_angles(it: &mut Interp, st: &Structure, flat: &FlatAtoms) -> Vec<(u32, u32, String)> {
    super::linked::add_replacement_angles(it, st, flat, &|_, bond, _|
        (bond.type_.eq_ignore_ascii_case("single") || bond.type_.eq_ignore_ascii_case("sing")).then_some(2))
}

fn fixture() -> (Structure, FlatAtoms, Interp) {
    let mut st = crate::mmcif::read_mmcif(concat!(
        "data_test\nloop_\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_alt_id\n",
        "_atom_site.auth_comp_id\n_atom_site.auth_asym_id\n_atom_site.auth_seq_id\n",
        "_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n_atom_site.occupancy\n_atom_site.B_iso_or_equiv\n",
        "C C35 . Z41 G 602 0 0 0 1 10\nC C18 . Z41 G 602 1.5 0 0 1 10\n",
        "H H351 . Z41 G 602 0 0 1 1 10\nH H352 . Z41 G 602 0 0 -1 1 10\n",
        "C EXTRA . Z41 G 602 1 1 1 1 10\nS SG . CYS G 18 0 1.8 0 1 10\nC CB . CYS G 18 0 3.6 0 1 10\n",
    )).unwrap();
    let atom = |name: &str, resname: &str, resseq: &str| AtomLabel {
        name: name.into(), resname: resname.into(), resseq: resseq.into(), chain: "G".into(), altloc: "".into(), icode: "".into(),
    };
    st.links.push(LinkRecord { atom1: atom("C35", "Z41", "602"), atom2: atom("SG", "CYS", "18"), distance: Some(1.8), order: None, kind: "covale".into() });
    let flat = FlatAtoms::from_structure(&st);
    let mut it = Interp::empty(flat.pos.len());
    let mut comp = Comp { id: "Z41".into(), ..Comp::default() };
    for (id, element) in [("C35", "C"), ("C18", "C"), ("O5", "O"), ("H351", "H"), ("H352", "H")] {
        comp.atoms.push(crate::monlib::CompAtom { id: id.into(), type_symbol: element.into(), type_energy: None });
    }
    for name in ["C18", "O5", "H351", "H352"] {
        comp.bonds.push(crate::monlib::CompBond { a1: "C35".into(), a2: name.into(), type_: "single".into(), ..Default::default() });
    }
    for (name, value) in [("C18", 108.21), ("H351", 112.03), ("H352", 112.07)] {
        comp.angles.push(crate::monlib::CompAngle { a1: name.into(), a2: "C35".into(), a3: "O5".into(), value: Some(value), esd: Some(3.0) });
    }
    let comp = Arc::new(comp);
    for (i, id) in ["C35", "C18", "H351", "H352"].iter().enumerate() {
        it.dict[i] = interp::AtomDictInfo { comp: Some(comp.clone()), dict_id: Some((*id).into()), instance: 0 };
    }
    let sulfur = Arc::new(Comp { id: "CYS".into(), bonds: vec![crate::monlib::CompBond {
        a1: "SG".into(), a2: "CB".into(), type_: "single".into(), ..Default::default()
    }], ..Default::default() });
    for (i, id) in [(5, "SG"), (6, "CB")] {
        it.dict[i] = interp::AtomDictInfo { comp: Some(sulfur.clone()), dict_id: Some(id.into()), instance: 0 };
    }
    it.add_bond(5, 6, 1.8, 0);
    for (i, ideal, origin) in [(1, 1.5, 0), (2, 0.97, 0), (3, 0.97, 0), (5, 1.766, 10)] {
        it.add_bond(0, i, ideal, origin);
    }
    (st, flat, it)
}

fn unparameterized(it: &Interp, flat: &FlatAtoms) -> Vec<u32> {
    let occupancy = vec![1.0; flat.pos.len()];
    let residues = vec![0; flat.pos.len()];
    let input = riding::RidingAtoms {
        is_h: &flat.is_h, altloc: &flat.altloc, name: &flat.name,
        occ: &occupancy, rg_of: &residues, expected_heavy: &|_| None,
        dictionary_nh2_torsion: true, cis_partner: &|_| None, reroute_axial_reference: true,
        use_plane_restraints: true,
    };
    riding::riding(it, &mut flat.pos.clone(), &input, true).unparameterized
}

fn paired_fixture() -> (Structure, FlatAtoms, Interp) {
    let (mut st, flat, old) = fixture();
    let mut copy = st.models[0].chains[0].clone();
    copy.id = "J".into();
    for residue in &mut copy.residue_groups {
        for group in &mut residue.atom_groups {
            for atom in &mut group.atoms { atom.xyz.x += 20.0; }
        }
    }
    st.models[0].chains.push(copy);
    st.reset_i_seq();
    let mut it = Interp::empty(flat.pos.len() * 2);
    it.dict = [old.dict.clone(), old.dict].concat();
    for bond in old.bonds {
        for offset in [0, flat.pos.len() as u32] {
            it.add_bond(bond.i + offset, bond.j + offset, bond.ideal, bond.origin);
        }
    }
    let flat = FlatAtoms::from_structure(&st);
    (st, flat, it)
}

#[test]
fn matching_unannotated_product_uses_declared_source_angles() {
    let (st, flat, mut it) = paired_fixture();
    let positions = flat.pos.clone();
    assert_eq!(add_replacement_angles(&mut it, &st, &flat), [(0, 5, "O5".into()), (7, 12, "O5".into())]);
    assert_eq!(it.angles.iter().map(|a| a.ideal).collect::<Vec<_>>(), [108.21, 112.03, 112.07, 108.21, 112.03, 112.07]);
    assert_eq!(flat.pos, positions);
}

#[test]
fn product_precedent_requires_matching_source_site_and_unique_valid_bonds() {
    for invalid in 0..10 {
        let (mut st, mut flat, mut it) = paired_fixture();
        match invalid {
            0 => Arc::make_mut(it.dict[7].comp.as_mut().unwrap()).id = "OTHER".into(),
            1 => it.dict[7].dict_id = Some("OTHER".into()),
            2 => Arc::make_mut(it.dict[12].comp.as_mut().unwrap()).id = "OTHER".into(),
            3 => it.dict[12].dict_id = Some("OTHER".into()),
            4 => { it.add_bond(7, 6, 1.8, interp::ORIGIN_MISC); }
            5 => flat.pos[12].y = 8.0,
            6 => { Arc::make_mut(it.dict[7].comp.as_mut().unwrap()).bonds.push(crate::monlib::CompBond {
                a1: "O5".into(), a2: "C18".into(), type_: "single".into(), ..Default::default()
            }); },
            7 => { let mut bad = st.links[0].clone(); bad.atom1.chain = "J".into(); bad.atom2.chain = "J".into();
                bad.atom2.resseq = "999".into(); st.links.push(bad); },
            8 => { let mut bad = st.links[0].clone(); bad.atom1.chain = "J".into(); bad.atom2.chain = "J".into();
                bad.order = Some("unknown".into()); st.links.push(bad); },
            9 => { it.add_bond(12, 4, 1.8, interp::ORIGIN_MISC); }
            _ => unreachable!(),
        };
        assert_eq!(add_replacement_angles(&mut it, &st, &flat), [(0, 5, "O5".into())], "case {invalid}");
        assert_eq!(it.angles.len(), 3);
    }
}

#[test]
fn absent_ambiguous_or_distant_source_product_gives_no_precedent() {
    for invalid in 0..5 {
        let (mut st, mut flat, mut it) = paired_fixture();
        match invalid {
            0 => st.links.clear(),
            1 => st.links.push(st.links[0].clone()),
            2 => flat.pos[5].y = 8.0,
            3 => st.links[0].order = Some("doub".into()),
            4 => st.links[0].atom1.chain = "OTHER".into(),
            _ => unreachable!(),
        }
        assert!(add_replacement_angles(&mut it, &st, &flat).is_empty(), "case {invalid}");
        assert!(it.angles.is_empty());
    }
}

#[test]
fn replacement_adds_only_dictionary_angles_and_keeps_bonds() {
    let (st, flat, mut it) = fixture();
    let bonds = it.bonds.iter().map(|b| (b.i, b.j, b.ideal)).collect::<Vec<_>>();
    assert_eq!(add_replacement_angles(&mut it, &st, &flat), [(0, 5, "O5".into())]);
    assert_eq!(it.angles.iter().map(|a| a.ideal).collect::<Vec<_>>(), [108.21, 112.03, 112.07]);
    assert_eq!(bonds, it.bonds.iter().map(|b| (b.i, b.j, b.ideal)).collect::<Vec<_>>());
}

#[test]
fn overvalent_or_unlinked_parents_get_no_replacement_restraints() {
    for invalid in 0..4 {
        let (mut st, flat, mut it) = fixture();
        match invalid {
            0 => { it.add_bond(0, 4, 1.5, 10); }
            1 => { st.links.clear(); }
            2 => { st.links.push(st.links[0].clone()); }
            3 => { Arc::make_mut(it.dict[0].comp.as_mut().unwrap()).bonds[0].type_ = "double".into(); }
            _ => unreachable!(),
        }
        assert!(add_replacement_angles(&mut it, &st, &flat).is_empty());
        assert!(it.angles.is_empty());
        assert_eq!(unparameterized(&it, &flat), [2, 3]);
    }
}

#[test]
fn nonterminal_or_ambiguous_missing_atoms_get_no_replacement_restraints() {
    for ambiguous in [false, true] {
        let (st, flat, mut it) = fixture();
        let comp = Arc::make_mut(it.dict[0].comp.as_mut().unwrap());
        if ambiguous {
            comp.bonds.push(crate::monlib::CompBond { a1: "C35".into(), a2: "OX".into(), type_: "single".into(), ..Default::default() });
        } else {
            comp.bonds.push(crate::monlib::CompBond { a1: "O5".into(), a2: "C18".into(), type_: "single".into(), ..Default::default() });
        }
        assert!(add_replacement_angles(&mut it, &st, &flat).is_empty());
        assert!(it.angles.is_empty());
        assert_eq!(unparameterized(&it, &flat), [2, 3]);
    }
}

#[test]
fn actual_dt_replacement_keeps_source_angles_bonds_and_coordinates() {
    let ml = crate::monlib::MonLib::load(&crate::monlib::MonLib::locate(None).unwrap()).unwrap();
    let mut st = crate::mmcif::read_mmcif(include_str!("../../tests/data/ptr_dt_3px7.cif")).unwrap();
    let params = HPlaceParams { neutron: false, n_terminal_charge: NTermCharge::ResidueOne,
        exclude_water: true, keep_existing_h: false, adp_scale: 1.0, compat: false, cell: None };
    let mut auto_comps = FxHashMap::default();
    add_missing_h(&mut st, &ml, &params, &mut Vec::new(), &mut auto_comps);
    st.sort_atoms_in_place();
    let flat = FlatAtoms::from_structure(&st);
    let mut it = interp::interpret(&st, &flat, &ml, &InterpParams {
        neutron: false, link_distance_cutoff: 3.0, compat: false, auto_comps,
    });
    let bonds = it.bonds.iter().map(|b| (b.i, b.j, b.ideal, b.origin)).collect::<Vec<_>>();
    let positions = flat.pos.clone();
    let initial_angles = it.angles.len();
    let replacements = super::linked::add_replacement_angles(&mut it, &st, &flat,
        &|comp, bond, actual| interp::source_bond_order(&ml, comp, bond, actual));
    let [(parent, partner, leaving)] = &replacements[..] else { panic!("expected the declared DT replacement") };
    assert_eq!(flat.name[*parent as usize].trim(), "C5'");
    assert_eq!(flat.name[*partner as usize].trim(), "O3P");
    assert_eq!(leaving, "O5'");
    let comp = it.dict[*parent as usize].comp.as_ref().unwrap();
    assert_eq!(comp.id, "Td");
    assert_eq!(comp.applied_mods, ["5*END"]);
    assert_eq!(it.angles[initial_angles..].iter().map(|a| (a.ideal, a.esd)).collect::<Vec<_>>(),
        [(109.4, 0.8), (109.0, 1.5), (109.0, 1.5)]);
    assert!(it.angles[initial_angles..].iter().all(|a| a.i[1] == *parent && a.i.contains(partner)));
    assert_eq!(bonds, it.bonds.iter().map(|b| (b.i, b.j, b.ideal, b.origin)).collect::<Vec<_>>());
    assert_eq!(flat.pos, positions);
}
