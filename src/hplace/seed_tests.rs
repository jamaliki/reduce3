use super::*;
use crate::monlib::{CompAtom, CompBond};

fn serine() -> (Comp, Vec<Atom>) {
    let sites = [
        ("N", "N", [-9.168000221252441, 3.950000047683716, 23.77199935913086]),
        ("CA", "C", [-8.21399974822998, 2.8529999256134033, 24.030000686645508]),
        ("C", "C", [-7.635000228881836, 2.2309999465942383, 22.75200080871582]),
        ("O", "O", [-8.364999771118164, 2.109999895095825, 21.764999389648438]),
        ("CB", "C", [-8.914999961853027, 1.7419999837875366, 24.819000244140625]),
        ("OG", "O", [-9.987000465393066, 1.2289999723434448, 24.047000885009766]),
    ];
    let atoms = sites.iter().map(|(name, element, xyz)| Atom::new(name, element, v3(xyz[0], xyz[1], xyz[2]))).collect();
    let comp = Comp {
        atoms: sites.iter().map(|(id, symbol, _)| CompAtom {
            id: id.to_string(), type_symbol: symbol.to_string(), ..CompAtom::default()
        }).chain([CompAtom { id: "HA".into(), type_symbol: "H".into(), ..CompAtom::default() }]).collect(),
        bonds: [("CA", "N"), ("CA", "C"), ("CA", "CB"), ("CA", "HA")].into_iter()
            .map(|(a1, a2)| CompBond { a1: a1.into(), a2: a2.into(), ..CompBond::default() }).collect(),
        ..Comp::default()
    };
    (comp, atoms)
}

#[test]
fn afdb_ser200_placeholder_near_parent_gets_a_valid_methine_seed() {
    let (comp, atoms) = serine();
    let bogus = atoms.iter().fold(Vec3::ZERO, |sum, atom| sum + atom.xyz) / atoms.len() as f64 + v3(0.5, 0.5, 0.5);
    let ca = atoms.iter().find(|atom| atom.name.trim() == "CA").unwrap().xyz;
    assert!(bogus.dist(ca) < 0.001);
    let fixed = methine_seed(&comp, &atoms, "HA", bogus);
    assert!((fixed.dist(ca) - 1.).abs() < 1e-12);
    assert_eq!(atoms[1].xyz, ca);
}

#[test]
fn valid_seeds_and_other_hydrogens_do_not_change() {
    let (comp, atoms) = serine();
    let seed = v3(1., 2., 3.);
    assert_eq!(methine_seed(&comp, &atoms, "HA", seed), seed);
    assert_eq!(methine_seed(&comp, &atoms, "HG", atoms[1].xyz), atoms[1].xyz);
    let mut incomplete = comp.clone();
    incomplete.bonds.retain(|bond| bond.a2 != "CB");
    assert_eq!(methine_seed(&incomplete, &atoms, "HA", atoms[1].xyz), atoms[1].xyz);
}
