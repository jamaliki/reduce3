//! Atom classification tables from mmtbx/probe/AtomTypes.py.

const AROMATIC_ACCEPTORS: &[(&[&str], &[&str])] = &[
    (&["HIS"], &["ND1", "NE2"]),
    (&["TRP"], &["CH2", "CZ3", "CZ2", "CE3", "CE2", "CD2"]),
    (&["ADE", "A"], &["N1", "N3", "N7", "C2", "C4", "C5", "C6", "C8", "N9"]),
    (&["CYT", "C"], &["N3", "N1", "C2", "C4", "C5", "C6"]),
    (&["GUA", "G"], &["N3", "N7", "N1", "C2", "C4", "C5", "C6", "C8", "N9"]),
    (&["THY", "T"], &["N1", "C2", "N3", "C4", "C5", "C6"]),
    (&["URA", "U"], &["N1", "C2", "N3", "C4", "C5", "C6"]),
    (&["DA"], &["N1", "N3", "N7", "C2", "C4", "C5", "C6", "C8", "N9"]),
    (&["DC"], &["N3", "N1", "C2", "C4", "C5", "C6"]),
    (&["DG"], &["N3", "N7", "N1", "C2", "C4", "C5", "C6", "C8", "N9"]),
    (&["DT"], &["N1", "C2", "N3", "C4", "C5", "C6"]),
    (&["HEM"], &["N A", "N B", "N C", "N D"]),
    (
        &["HEM"],
        &[
            "C1A", "C2A", "C3A", "C4A", "C1B", "C2B", "C3B", "C4B", "C1C", "C2C", "C3C", "C4C", "C1D", "C2D",
            "C3D", "C4D",
        ],
    ),
    (&["PHE"], &["CZ", "CE2", "CE1", "CD2", "CD1", "CG"]),
    (&["TYR"], &["CZ", "CE2", "CE1", "CD2", "CD1", "CG"]),
];

/// `IsAromaticAcceptor(resName, atomName)` (both compared stripped).
pub fn is_aromatic_acceptor(resname: &str, atom_name: &str) -> bool {
    let r = resname.trim();
    let a = atom_name.trim();
    AROMATIC_ACCEPTORS.iter().any(|(rs, ats)| rs.contains(&r) && ats.contains(&a))
}

/// `IsSpecialAminoAcidCarbonyl`: CG of ASP/ASN/ASX and CD of GLU/GLN/GLX.
///
/// The original compares the padded PDB name (`" CG"`), so it never matches
/// names read from mmCIF; `padded_name_only` reproduces that.
pub fn is_special_amino_acid_carbonyl(resname: &str, raw_atom_name: &str, padded_name_only: bool) -> bool {
    let n = if padded_name_only { raw_atom_name.trim_end().to_string() } else { format!(" {}", raw_atom_name.trim()) };
    if n == " CG" {
        return matches!(resname, "ASP" | "ASN" | "ASX");
    }
    if n == " CD" {
        return matches!(resname, "GLU" | "GLN" | "GLX");
    }
    false
}
