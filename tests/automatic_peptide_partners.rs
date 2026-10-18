#[path = "common/validation_fixture.rs"]
mod fixture;

use reduce3::autolink::{auto_link, AutoLinkParams};
use reduce3::geom::v3;
use reduce3::interp::{FlatAtoms, Interp, ORIGIN_COVALENT, ORIGIN_MISC};
use reduce3::model::{Atom, AtomGroup, AtomLabel, Chain, LinkRecord, Model, ResidueGroup, Structure};
use reduce3::MonLib;

fn chain(id: &str, residues: &[(&str, &str, Vec<Atom>)]) -> Chain {
    Chain { id: id.into(), residue_groups: residues.iter().map(|(seq, name, atoms)| ResidueGroup {
        resseq: (*seq).into(), icode: String::new(), link_to_previous: false,
        atom_groups: vec![AtomGroup { resname: (*name).into(), altloc: String::new(), atoms: atoms.clone() }],
    }).collect() }
}

fn atom(name: &str, element: &str, x: f64) -> Atom { Atom::new(&format!(" {name:<3}"), element, v3(x, 0., 0.)) }

fn structure(chains: Vec<Chain>) -> Structure {
    Structure { models: vec![Model { id: "1".into(), chains }], ..Default::default() }
}

fn run(st: &Structure, it: &mut Interp, compat: bool) {
    let files = fixture::library();
    let mut ml = MonLib::load(&files.0).unwrap();
    // Automatic TRANS links use N/C group order; this fixture uses C/N for polymer links.
    let link = &mut ml.links[ml.link_index["TRANS"]];
    link.bonds[0].c1 = 2;
    link.bonds[0].c2 = 1;
    let flat = FlatAtoms::from_structure(st);
    let before = format!("{st:?}");
    auto_link(it, st, &flat, &ml, &AutoLinkParams { compat }, &mut String::new());
    assert_eq!(format!("{st:?}"), before);
}

fn declare(st: &mut Structure, a: u32, b: u32) {
    let flat = FlatAtoms::from_structure(st);
    let label = |i: u32| {
        let p = flat.path[i as usize];
        AtomLabel { name: st.atom(p).name.to_string(), resname: st.atom_group(p).resname.clone(),
            chain: st.chain(p).id.clone(), resseq: st.residue_group(p).resseq.clone(),
            altloc: st.atom_group(p).altloc.clone(), icode: st.residue_group(p).icode.clone() }
    };
    st.links.push(LinkRecord { atom1: label(a), atom2: label(b), kind: "covale".into(), order: None, distance: None });
}

#[test]
fn initial_peptide_at_either_endpoint_blocks_only_fixed_automatic_duplicates() {
    for occupied_n in [false, true] {
        for origin in [ORIGIN_COVALENT, 35, 63, 65, 66, 70, 79] {
            for compat in [false, true] {
                let other = if occupied_n { atom("C", "C", -1.2) } else { atom("N", "N", 2.5) };
                let st = structure(vec![
                    chain("A", &[("1", "ALA", vec![atom("N", "N", 0.)])]),
                    chain("B", &[("2", "ALA", vec![atom("C", "C", 1.3)])]),
                    chain("C", &[("3", "ALA", vec![other])]),
                ]);
                let mut it = Interp::empty(3);
                let occupied = if occupied_n { 0 } else { 1 };
                it.add_bond(occupied, 2, 1.33, origin);
                run(&st, &mut it, compat);
                assert_eq!(it.has_bond(0, 1), compat, "origin={origin}, occupied_n={occupied_n}");
                assert!(it.has_bond(occupied, 2));
            }
        }
    }
}

#[test]
fn accepted_automatic_peptide_also_blocks_a_later_candidate() {
    for occupied_n in [false, true] {
        let mut st = structure(vec![
            chain("A", &[("1", "ALA", vec![atom("N", "N", 0.)])]),
            chain("B", &[("2", "ALA", vec![atom("C", "C", 1.2)])]),
            chain("C", &[("3", "ALA", vec![atom("C", "C", -1.4)])]),
        ]);
        if !occupied_n {
            for p in st.atom_paths() {
                let atom = st.atom_mut(p);
                let name = if atom.name.trim() == "N" { "C" } else { "N" };
                atom.name = format!(" {name:<3}").into();
                atom.element = name.into();
            }
        }
        let mut it = Interp::empty(3);
        run(&st, &mut it, false);
        assert!(it.has_bond(0, 1));
        assert!(!it.has_bond(0, 2));
        assert_eq!(it.bonds.len(), 1);
    }
}

#[test]
fn explicit_peptide_records_keep_the_existing_link_behavior() {
    for repeats in [1, 2] {
        let mut st = structure(vec![
            chain("A", &[("1", "ALA", vec![atom("N", "N", 0.)])]),
            chain("B", &[("2", "ALA", vec![atom("C", "C", 1.3)])]),
            chain("C", &[("3", "ALA", vec![atom("C", "C", -1.2)])]),
        ]);
        for _ in 0..repeats { declare(&mut st, 0, 1); }
        let mut it = Interp::empty(3);
        it.add_bond(0, 2, 1.33, ORIGIN_COVALENT);
        run(&st, &mut it, false);
        assert!(it.has_bond(0, 1));
        assert!(it.has_bond(0, 2));
    }
}

#[test]
fn stale_or_ambiguous_records_do_not_authorize_a_duplicate_peptide() {
    for ambiguous in [false, true] {
        let mut st = structure(vec![
            chain("A", &[("1", "ALA", vec![atom("N", "N", 0.)])]),
            chain("B", &[("2", "ALA", vec![atom("C", "C", 1.3)])]),
            chain("C", &[("3", "ALA", vec![atom("C", "C", -1.2)])]),
        ]);
        declare(&mut st, 0, 1);
        if ambiguous {
            st.models[0].chains[1].residue_groups[0].atom_groups[0].atoms.push(atom("C", "C", 10.));
        } else { st.links[0].atom2.chain = "missing".into(); }
        let flat = FlatAtoms::from_structure(&st);
        let previous = flat.path.iter().position(|&p| st.chain(p).id == "C").unwrap() as u32;
        let mut it = Interp::empty(flat.pos.len());
        it.add_bond(0, previous, 1.33, ORIGIN_COVALENT);
        run(&st, &mut it, false);
        assert!(!it.has_bond(0, 1));
    }
}

#[test]
fn a_rejected_duplicate_leaves_its_free_carbon_available() {
    let st = structure(vec![
        chain("A", &[("1", "ALA", vec![atom("N", "N", 0.)])]),
        chain("B", &[("2", "ALA", vec![atom("C", "C", 1.2)])]),
        chain("C", &[("3", "ALA", vec![atom("C", "C", -1.3)])]),
        chain("D", &[("4", "ALA", vec![atom("N", "N", 2.5)])]),
    ]);
    let mut it = Interp::empty(4);
    it.add_bond(0, 2, 1.33, ORIGIN_COVALENT);
    run(&st, &mut it, false);
    assert!(!it.has_bond(0, 1));
    assert!(it.has_bond(1, 3));
}

#[test]
fn free_chain_ends_and_cyclic_peptides_still_join() {
    for same_chain in [false, true] {
        let residues = [("1", "ALA", vec![atom("N", "N", 0.)]),
            ("2", "ALA", vec![atom("C", "C", 1.3)])];
        let chains = if same_chain { vec![chain("A", &residues)] }
            else { vec![chain("A", &residues[..1]), chain("B", &residues[1..])] };
        let st = structure(chains);
        let mut it = Interp::empty(2);
        run(&st, &mut it, false);
        assert!(it.has_bond(0, 1));
    }
}

#[test]
fn incompatible_alternate_partner_does_not_block_a_free_conformer() {
    let mut st = structure(vec![
        chain("A", &[("1", "ALA", vec![atom("N", "N", 0.)])]),
        chain("B", &[("2", "ALA", vec![atom("C", "C", 1.3)])]),
        chain("C", &[("3", "ALA", vec![atom("C", "C", -1.2)])]),
    ]);
    st.models[0].chains[1].residue_groups[0].atom_groups[0].altloc = "B".into();
    st.models[0].chains[2].residue_groups[0].atom_groups[0].altloc = "A".into();
    let mut it = Interp::empty(3);
    it.add_bond(0, 2, 1.33, ORIGIN_COVALENT);
    run(&st, &mut it, false);
    assert!(it.has_bond(0, 1));
}

#[test]
fn sidechain_custom_products_and_nonpeptide_origins_do_not_prove_a_peptide_partner() {
    for backbone in [false, true] {
        let st = structure(vec![
            chain("A", &[("1", "ALA", vec![atom("N", "N", 0.)])]),
            chain("B", &[("2", "ALA", vec![atom("C", "C", 1.3)])]),
            chain("C", &[("3", "ASP", vec![atom(if backbone { "C" } else { "CG" }, "C", -1.2)])]),
        ]);
        let mut it = Interp::empty(3);
        it.add_bond(0, 2, 1.33, if backbone { ORIGIN_MISC } else { ORIGIN_COVALENT });
        run(&st, &mut it, false);
        assert!(it.has_bond(0, 1));
        assert!(it.has_bond(0, 2));
    }
    let st = structure(vec![
        chain("A", &[("1", "ALA", vec![atom("N", "N", 0.)])]),
        chain("B", &[("2", "ASP", vec![atom("CG", "C", 1.3)])]),
        chain("C", &[("3", "ALA", vec![atom("C", "C", -1.2)])]),
    ]);
    let mut it = Interp::empty(3);
    it.add_bond(0, 2, 1.33, ORIGIN_COVALENT);
    run(&st, &mut it, false);
    assert!(it.has_bond(0, 1));
}

#[test]
fn exact_4u90_assembly_contacts_keep_the_four_deposited_peptides() {
    // owner-19/train-001139-0031: original atom-site rows and deposited coordinates.
    let source = [
        ("D", 3139, [15.967, 113.711, 109.913], 3148, [15.172, 113.886, 108.881], 3149, [15.705, 114.099, 107.551]),
        ("E", 3217, [15.773, 84.529, 118.942], 3226, [15.021, 84.637, 120.023], 3227, [15.591, 84.426, 121.338]),
        ("D-2", 6715, [15.967, 84.769, 118.899], 6724, [15.172, 84.594, 119.931], 6725, [15.705, 84.381, 121.261]),
        ("E-2", 6793, [15.773, 113.951, 109.870], 6802, [15.021, 113.843, 108.789], 6803, [15.591, 114.054, 107.474]),
    ];
    let chains = source.iter().map(|(id, c_row, c_xyz, n_row, n_xyz, ca_row, ca_xyz)| {
        let make = |name, element, row, xyz: &[f64; 3]| {
            let mut atom = Atom::new(&format!(" {name:<3}"), element, v3(xyz[0], xyz[1], xyz[2]));
            atom.src = row;
            atom
        };
        let mut chain = chain(id, &[("368", "PHE", vec![make("C", "C", *c_row, c_xyz)]),
            ("369", "ASN", vec![make("N", "N", *n_row, n_xyz), make("CA", "C", *ca_row, ca_xyz)])]);
        let residue = &mut chain.residue_groups[1];
        let mut alternate = residue.atom_groups[0].clone();
        alternate.altloc = "A".into();
        alternate.atoms.remove(0);
        residue.atom_groups[0].atoms.pop();
        residue.atom_groups.push(alternate);
        chain
    }).collect();
    let st = structure(chains);
    let mut it = Interp::empty(12);
    for c in [0, 3, 6, 9] {
        it.add_bond(c, c + 1, 1.33, ORIGIN_COVALENT);
        it.add_bond(c + 1, c + 2, 1.458, ORIGIN_COVALENT);
    }
    let before = it.bonds.iter().map(|b| (b.i, b.j, b.origin, b.ideal.to_bits())).collect::<Vec<_>>();
    run(&st, &mut it, false);
    let after = it.bonds.iter().map(|b| (b.i, b.j, b.origin, b.ideal.to_bits())).collect::<Vec<_>>();
    assert_eq!(&after[..before.len()], &before);
    let flat = FlatAtoms::from_structure(&st);
    let peptide = |i: u32, j: u32| matches!((flat.name[i as usize].trim(), flat.name[j as usize].trim()),
        ("C", "N") | ("N", "C"));
    assert_eq!(after.iter().filter(|&&(i, j, _, _)| peptide(i, j)).collect::<Vec<_>>(),
        before.iter().filter(|&&(i, j, _, _)| peptide(i, j)).collect::<Vec<_>>());
    assert_eq!(st.models[0].id, "1");
}
