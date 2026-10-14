use super::*;
use crate::hplace::{self, HPlaceParams, NTermCharge, Placed};

fn library() -> MonLib {
    MonLib::load(&MonLib::locate(None).expect("Set REDUCE3_CHEM_DATA")).unwrap()
}

fn input() -> Structure {
    crate::mmcif::read_mmcif(include_str!("../../../../../tests/data/frame_3e6p.cif")).unwrap()
}

fn place(st: &mut Structure, ml: &MonLib, compat: bool) -> Placed {
    hplace::place_hydrogens(st, ml, &HPlaceParams { neutron: false, n_terminal_charge: NTermCharge::ResidueOne,
        exclude_water: true, keep_existing_h: false, adp_scale: 1.0, compat, cell: None })
}

fn site(st: &Structure, name: &str) -> Option<usize> {
    st.atom_paths().into_iter().position(|p| st.atom_group(p).resname == "DFK" && st.atom(p).name.trim() == name)
}

#[test]
fn exact_3e6p_product_keeps_both_fixed_methylene_h_sites_and_accepted_graph() {
    let ml = library();
    let mut st = input();
    assert!(st.atom_paths().into_iter().all(|p| !st.atom(p).is_hydrogen()));
    let original_heavy = st.atom_paths().into_iter().map(|p| (st.atom(p).src, st.atom(p).xyz)).collect::<Vec<_>>();
    let placed = place(&mut st, &ml, false);
    let paths = st.atom_paths();
    for (source, xyz) in original_heavy {
        assert_eq!(st.atom(*paths.iter().find(|&&p| st.atom(p).src == source).unwrap()).xyz, xyz);
    }
    for (name, xyz) in [
        ("HH21", v3(-49.62646542262007, -12.834117340010682, -20.747742889998275)),
        ("HH22", v3(-48.18460125326361, -12.337171219872829, -20.319839699919214)),
    ] {
        let index = site(&st, name).unwrap();
        assert!(placed.riding[index].is_none());
        assert_eq!(st.atom(paths[index]).src, Atom::NEW);
        assert!(st.atom(paths[index]).xyz.dist(xyz) < 1e-12);
        assert!(placed.bonds.iter().any(|&(a, b, _)| (a as usize == index && b as usize == site(&st, "CH2").unwrap())
            || (b as usize == index && a as usize == site(&st, "CH2").unwrap())));
    }
    assert!(site(&st, "HH23").is_none() && site(&st, "H3").is_none());
    let parent = site(&st, "CH2").unwrap() as u32;
    let heavy = placed.bonds.iter().filter_map(|&(a, b, origin)| {
        let other = if a == parent { b } else if b == parent { a } else { return None };
        (!st.atom(paths[other as usize]).is_hydrogen()).then_some((st.atom(paths[other as usize]).name.trim(), origin))
    }).collect::<Vec<_>>();
    assert_eq!(heavy, [("C3", interp::ORIGIN_COVALENT), ("NE2", interp::ORIGIN_MISC)]);
    let policy = Policy::read(&input(), &ml, &crate::Params { stop_on_any_missing_hydrogen: true,
        allow_unresolved_bond_targets: true, ..crate::Params::default() }).unwrap();
    let mut placed = placed;
    let reports = interp::target_uncertainty::inspect(&st, &placed.bonds);
    let before = format!("{st:?}");
    let bonds = placed.bonds.clone();
    assert!(prepare(&policy, &mut st, &mut placed, &reports, &ml).unwrap().is_empty());
    assert_eq!(format!("{st:?}"), before);
    assert_eq!(placed.bonds, bonds);
}

#[test]
fn absent_or_ambiguous_product_links_keep_three_h_prop_recipes() {
    let ml = library();
    for invalid in 0..3 {
        let mut st = input();
        let index = st.links.iter().position(|link| link.atom1.name == "CH2" || link.atom2.name == "CH2").unwrap();
        match invalid {
            0 => { st.links.remove(index); },
            1 => { let link = st.links[index].clone(); st.links.push(link); },
            _ => st.links[index].atom2.resseq = "999".into(),
        }
        let placed = place(&mut st, &ml, false);
        for name in ["HH21", "HH22", "HH23"] {
            let coef = placed.riding[site(&st, name).unwrap()].unwrap();
            assert_eq!(coef.htype, HType::Prop);
        }
    }
}

#[test]
fn compatibility_mode_keeps_its_existing_product_recipes() {
    let ml = library();
    let mut st = input();
    let placed = place(&mut st, &ml, true);
    for name in ["HH21", "HH22"] {
        assert_eq!(placed.riding[site(&st, name).unwrap()].unwrap().htype, HType::Prop);
    }
}
