use super::linked::*;
use super::*;

fn fixture() -> (Structure, FlatAtoms, Interp) {
    let mut st = crate::mmcif::read_mmcif(concat!(
        "data_test\nloop_\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_alt_id\n",
        "_atom_site.auth_comp_id\n_atom_site.auth_asym_id\n_atom_site.auth_seq_id\n",
        "_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n_atom_site.occupancy\n_atom_site.B_iso_or_equiv\n",
        "C C35 . Z41 G 602 0 0 0 1 10\nC C18 . Z41 G 602 1.5 0 0 1 10\n",
        "H H351 . Z41 G 602 0 0 1 1 10\nH H352 . Z41 G 602 0 0 -1 1 10\n",
        "C EXTRA . Z41 G 602 1 1 1 1 10\nS SG . CYS G 18 0 1.8 0 1 10\n",
    )).unwrap();
    let atom = |name: &str, resname: &str, resseq: &str| AtomLabel {
        name: name.into(), resname: resname.into(), resseq: resseq.into(), chain: "G".into(), altloc: "".into(), icode: "".into(),
    };
    st.links.push(LinkRecord { atom1: atom("C35", "Z41", "602"), atom2: atom("SG", "CYS", "18"), distance: Some(1.8), kind: "covale".into() });
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
        it.dict[i] = interp::AtomDictInfo { comp: Some(comp.clone()), dict_id: Some((*id).into()) };
    }
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
        dictionary_nh2_torsion: true, reroute_axial_reference: true,
        use_plane_restraints: true,
    };
    riding::riding(it, &mut flat.pos.clone(), &input, true).unparameterized
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
