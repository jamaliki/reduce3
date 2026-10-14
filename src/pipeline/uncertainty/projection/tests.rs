use super::*;
use crate::geom::{v3, Vec3};

fn fixture(observed: bool) -> (Structure, RidingCoef, Vec<u32>, Vec<bool>, Vec<(u32, u32, u16)>) {
    let mut st = crate::mmcif::read_mmcif(concat!(
        "data_frame\nloop_\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_alt_id\n",
        "_atom_site.auth_comp_id\n_atom_site.auth_asym_id\n_atom_site.auth_seq_id\n",
        "_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n",
        "_atom_site.occupancy\n_atom_site.B_iso_or_equiv\n",
        "O O . LIG A 1 0 0 0 1 10\nC C1 . LIG A 1 0 0 1 1 10\n",
        "C C2 . LIG A 1 1 0 1 1 10\nH HO . LIG A 1 0 1 0 1 10\nH HC . LIG A 1 -1 0 1 1 10\n",
    )).unwrap();
    let paths = st.atom_paths();
    let coef = RidingCoef { htype: HType::Alg1b, ih: 3, a0: 0, a1: 1, a2: 4, a3: -1,
        a: 1.864, b: -0.7, h: 0., n: 0, disth: 0.85 };
    let sites = paths.iter().map(|&p| st.atom(p).xyz).collect::<Vec<_>>();
    st.atom_mut(paths[3]).xyz = compute_h_position(&coef, &sites).unwrap();
    if !observed { st.atom_mut(paths[3]).src = Atom::NEW; }
    st.atom_mut(paths[3]).occ = 0.43;
    st.atom_mut(paths[3]).b = 17.8;
    st.retain_atoms(|atom| atom.name.trim() != "HC");
    (st, coef, vec![0, 1, 2, 3, u32::MAX], vec![false, false, false, true, true],
        vec![(0, 1, interp::ORIGIN_COVALENT), (1, 2, interp::ORIGIN_COVALENT), (0, 3, interp::ORIGIN_COVALENT)])
}

fn project(st: &Structure, coef: RidingCoef, map: &[u32], old_h: &[bool], bonds: &[(u32, u32, u16)]) -> Result<RidingCoef, String> {
    let paths = st.atom_paths();
    let sites = paths.iter().map(|&p| st.atom(p).xyz).collect::<Vec<Vec3>>();
    projected_riding(coef, st, &paths, &sites, map, old_h, bonds)
}

#[test]
fn removed_and_sentinel_h_refs_rebase_without_changing_observed_or_generated_sites() {
    for element in ["N", "O", "S"] {
        for (observed, sentinel) in [(false, false), (false, true), (true, false), (true, true)] {
            let (mut st, mut coef, map, old_h, bonds) = fixture(observed);
            st.atom_mut(st.atom_paths()[0]).element = element.into();
            if sentinel { coef.a2 = u32::MAX as i64; }
            let before = format!("{st:?}");
            let actual = project(&st, coef, &map, &old_h, &bonds).unwrap();
            assert_eq!(actual.a2, 2);
            let sites = st.atom_paths().iter().map(|&p| st.atom(p).xyz).collect::<Vec<_>>();
            assert!(compute_h_position(&actual, &sites).unwrap().dist(sites[3]) < 1e-12);
            assert_eq!(format!("{st:?}"), before);
        }
    }
}

#[test]
fn retained_references_keep_their_coefficients() {
    let (st, mut coef, map, old_h, bonds) = fixture(false);
    coef.a2 = 2;
    assert_eq!(format!("{:?}", project(&st, coef, &map, &old_h, &bonds).unwrap()), format!("{coef:?}"));
}

#[test]
fn absent_bonds_or_changed_geometry_cannot_supply_a_frame() {
    for fault in 0..5 {
        let (mut st, coef, map, old_h, mut bonds) = fixture(true);
        match fault {
            0 => bonds.retain(|&(a, b, _)| (a, b) != (1, 2)),
            1 => bonds.retain(|&(a, b, _)| (a, b) != (0, 1)),
            2 => st.atom_mut(st.atom_paths()[2]).xyz = v3(0., 0., 2.),
            3 => st.atom_mut(st.atom_paths()[3]).xyz += v3(0.2, 0.1, 0.),
            _ => st.atom_mut(st.atom_paths()[0]).element = "C".into(),
        }
        assert!(project(&st, coef, &map, &old_h, &bonds).is_err(), "fault {fault}");
    }
}

#[test]
fn other_types_removed_heavy_refs_and_invalid_indices_keep_the_guard() {
    for fault in 0..4 {
        let (st, mut coef, map, mut old_h, bonds) = fixture(false);
        match fault {
            0 => coef.htype = HType::Prop,
            1 => old_h[4] = false,
            2 => coef.a2 = 99,
            _ => coef.a0 = 4,
        }
        assert!(project(&st, coef, &map, &old_h, &bonds).is_err(), "fault {fault}");
    }
}

#[test]
fn source_rows_models_and_altlocs_must_share_an_accepted_frame() {
    for fault in 0..6 {
        let (mut st, coef, mut map, old_h, bonds) = fixture(false);
        if fault < 4 {
            let index = if fault < 2 { 2 } else { 1 };
            let group = &mut st.models[0].chains[0].residue_groups[0].atom_groups[0];
            group.altloc = "A".into();
            let atom = group.atoms.remove(index);
            let other = AtomGroup { atoms: vec![atom], resname: "LIG".into(), altloc: "B".into() };
            if fault % 2 == 0 {
                st.models[0].chains[0].residue_groups[0].atom_groups.push(other);
            } else {
                let mut model = st.models[0].clone();
                model.id = "2".into();
                model.chains[0].residue_groups[0].atom_groups = vec![other];
                st.models.push(model);
            }
            let paths = st.atom_paths();
            for (old, name) in ["O", "C1", "C2", "HO"].into_iter().enumerate() {
                map[old] = paths.iter().position(|&p| st.atom(p).name.trim() == name).unwrap() as u32;
            }
        } else {
            st.atom_mut(st.atom_paths()[if fault == 4 { 2 } else { 1 }]).src = Atom::NEW;
        }
        assert!(project(&st, coef, &map, &old_h, &bonds).is_err(), "fault {fault}");
    }
}

#[test]
fn exact_bjm_imine_and_hse_hydroxyl_keep_placed_sites_and_dictionary_restraints() {
    let ml = MonLib::load(&MonLib::locate(None).expect("Set REDUCE3_CHEM_DATA")).unwrap();
    for (cif, component, hydrogen, parent, reference) in [
        (include_str!("../../../../tests/data/frame_3bjm.cif"), "BJM", "H19", "N23", "C18"),
        (include_str!("../../../../tests/data/frame_1waw.cif"), "HSE", "HO3", "O3", "C3"),
    ] {
        let mut st = crate::mmcif::read_mmcif(cif).unwrap();
        let policy = Policy::read(&st, &ml, &crate::Params { stop_on_any_missing_hydrogen: true,
            allow_unresolved_bond_targets: true, ..crate::Params::default() }).unwrap();
        let mut placed = crate::hplace::place_hydrogens(&mut st, &ml, &crate::hplace::HPlaceParams {
            neutron: false, n_terminal_charge: crate::hplace::NTermCharge::ResidueOne,
            exclude_water: true, keep_existing_h: false, adp_scale: 1.0, compat: false, cell: None,
        });
        let paths = st.atom_paths();
        let hi = paths.iter().position(|&p| st.atom_group(p).resname == component && st.atom(p).name.trim() == hydrogen).unwrap();
        let atom = st.atom(paths[hi]).clone();
        let coef = placed.riding[hi].unwrap();
        assert_eq!(coef.htype, HType::Alg1b);
        assert_eq!(st.atom(paths[coef.a0 as usize]).name.trim(), parent);
        let comp = placed.dict[hi].comp.as_ref().unwrap();
        let dictionary = format!("{comp:?}");
        if component == "BJM" {
            assert_eq!(coef.a2, u32::MAX as i64);
            assert!(comp.bonds.iter().any(|bond| [bond.a1.as_str(), bond.a2.as_str()] == ["N23", "C22"] && bond.type_ == "double"));
            assert!(comp.tors.iter().any(|tor| tor.a == ["H19", "N23", "C22", "C18"]));
        }
        let reports = interp::target_uncertainty::inspect(&st, &placed.bonds);
        assert!(!reports.is_empty());
        assert!(prepare(&policy, &mut st, &mut placed, &reports, &ml).unwrap().is_empty());
        let paths = st.atom_paths();
        let hi = paths.iter().position(|&p| st.atom_group(p).resname == component && st.atom(p).name.trim() == hydrogen).unwrap();
        let rebased = placed.riding[hi].unwrap();
        assert_eq!(st.atom(paths[rebased.a2 as usize]).name.trim(), reference);
        assert_eq!(st.atom(paths[hi]).xyz, atom.xyz);
        assert_eq!(st.atom(paths[hi]).src, atom.src);
        assert_eq!((rebased.htype, rebased.a, rebased.h, rebased.n, rebased.disth),
            (coef.htype, coef.a, coef.h, coef.n, coef.disth));
        assert_eq!(format!("{:?}", placed.dict[hi].comp.as_ref().unwrap()), dictionary);
        let sites = paths.iter().map(|&p| st.atom(p).xyz).collect::<Vec<_>>();
        assert!(compute_h_position(&rebased, &sites).unwrap().dist(atom.xyz) < 1e-9);
    }
}
