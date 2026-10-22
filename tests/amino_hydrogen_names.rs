//! Fixed mode gives the hydrogens of amino groups their IUPAC names whether
//! riding or a flip Mover places them (fix 5 in the README). HD22, HE22, HH11
//! and HH21 of Asn, Gln and Arg are cis to OD1, OE1 and NE with and without
//! flips. H21 of G and DG is cis to N1, as H61 of A and H41 of C are.
//! `--compat` keeps Reduce2's names, which tools/compare_outputs.sh checks.

use reduce3::geom::{dihedral_deg, Vec3};
use reduce3::model::Structure;
use reduce3::{mmcif, pipeline, MonLib, Params};
use std::collections::BTreeMap;

/// Every residue of 1A28 chain A within 4 A of Asn 689, Gln 725 and His 881,
/// with Arg 724 and 766 and Asn 828 among them.
const PROTEIN: &str = include_str!("data/flips_1a28.cif");
/// G, DG, A, DA, C and DC from the ideal coordinates of their CCD entries.
const BASES: &str = include_str!("data/amino_bases_ccd.cif");

/// Each amino group: the reference atom, the two atoms of the bond, the H cis
/// to the reference atom and the H trans to it.
type Group = (&'static str, &'static str, &'static str, &'static str, &'static str);

fn groups(resname: &str) -> &'static [Group] {
    match resname {
        "ASN" => &[("OD1", "CG", "ND2", "HD22", "HD21")],
        "GLN" => &[("OE1", "CD", "NE2", "HE22", "HE21")],
        "ARG" => &[("NE", "CZ", "NH1", "HH11", "HH12"), ("NE", "CZ", "NH2", "HH21", "HH22")],
        "G" | "DG" => &[("N1", "C2", "N2", "H21", "H22")],
        "A" | "DA" => &[("N1", "C6", "N6", "H61", "H62")],
        "C" | "DC" => &[("N3", "C4", "N4", "H41", "H42")],
        _ => &[],
    }
}

fn library() -> MonLib {
    MonLib::load(&MonLib::locate(None).expect("set REDUCE3_CHEM_DATA for the amino hydrogen tests")).unwrap()
}

fn run(library: &MonLib, model: &str, flips: bool) -> pipeline::Output {
    let mut params = Params::default();
    params.opt.add_flip_movers = flips;
    pipeline::run(mmcif::read_mmcif(model).unwrap(), library, &params).unwrap()
}

/// The atoms of each residue by chain, number and name.
fn residues(st: &Structure) -> BTreeMap<(String, String, String), BTreeMap<String, Vec3>> {
    let mut out: BTreeMap<_, BTreeMap<_, _>> = BTreeMap::new();
    for p in st.atom_paths() {
        let key = (st.chain(p).id.clone(), st.residue_group(p).resseq.trim().to_string(), st.atom_group(p).resname.trim().to_string());
        out.entry(key).or_default().insert(st.atom(p).name.trim().to_string(), st.atom(p).xyz);
    }
    out
}

/// Check every amino group of `st`, returning the residues checked.
fn check_names(st: &Structure) -> Vec<String> {
    let mut checked = Vec::new();
    for ((chain, resseq, resname), atoms) in residues(st) {
        for &(reference, a1, a0, cis, trans) in groups(&resname) {
            let torsion = |h: &str| {
                let at = |n: &str| atoms.get(n).copied().unwrap_or_else(|| panic!("{resname} {resseq} has no {n}"));
                dihedral_deg(at(reference), at(a1), at(a0), at(h)).unwrap()
            };
            let label = format!("{chain} {resname} {resseq}");
            assert!(torsion(cis).abs() < 30.0, "{label}: {reference}-{a1}-{a0}-{cis} is {:.1}", torsion(cis));
            assert!(torsion(trans).abs() > 150.0, "{label}: {reference}-{a1}-{a0}-{trans} is {:.1}", torsion(trans));
        }
        if !groups(&resname).is_empty() {
            checked.push(format!("{resname} {resseq}"));
        }
    }
    checked
}

fn flipped(output: &pipeline::Output, mover: &str) -> bool {
    output.description.lines().any(|l| l.contains(mover) && l.contains("pose Flipped"))
}

#[test]
fn amide_and_guanidinium_names_do_not_depend_on_flips() {
    let library = library();
    for flips in [true, false] {
        let output = run(&library, PROTEIN, flips);
        assert_eq!(check_names(&output.structure), ["ASN 689", "ARG 724", "GLN 725", "ARG 766", "ASN 828"]);
        // the Mover names both of its poses
        assert_eq!(flipped(&output, "AmideFlip at chain A ASN 689 ND2"), flips);
        assert_eq!(flipped(&output, "AmideFlip at chain A GLN 725 NE2"), flips);
        assert!(!flipped(&output, "ASN 828"));
    }
}

#[test]
fn base_amino_names_match_in_rna_and_dna() {
    let output = run(&library(), BASES, true);
    assert_eq!(check_names(&output.structure), ["G 1", "DG 1", "A 1", "DA 1", "C 1", "DC 1"]);
}
