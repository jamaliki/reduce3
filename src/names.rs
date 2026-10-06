//! Atom-name interpretation: maps model atom names to monomer-library atom ids
//! (port of iotbx.pdb.atom_name_interpretation and the mapping logic of
//! mmtbx.monomer_library.pdb_interpretation.monomer_mapping).

use crate::monlib::Comp;
use rustc_hash::FxHashMap;
use std::sync::OnceLock;

struct Interpreter {
    /// name (H/D expanded) -> expected pattern
    expected: FxHashMap<String, String>,
    /// synonym name -> expected name (H/D expanded)
    synonyms: FxHashMap<String, String>,
    meps: Vec<[&'static str; 3]>,
}

fn alternative_hydrogen_pattern(p: &str) -> Option<String> {
    let b = p.as_bytes();
    if b.len() > 1 && b[1] == b'h' && (b'1'..=b'9').contains(&b[0]) {
        Some(format!("{}{}", &p[1..], &p[..1]))
    } else {
        None
    }
}

fn build(expected_patterns: &[&str], synonym_patterns: &[(&str, &str)], meps: &[[&'static str; 3]]) -> Interpreter {
    let mut expected = FxHashMap::default();
    for &ep in expected_patterns {
        for h in ["H", "D"] {
            let name = ep.replace('h', h);
            let same = name == ep;
            expected.insert(name, ep.to_string());
            if same {
                break;
            }
        }
    }
    let mut syn: Vec<(String, String)> = synonym_patterns.iter().map(|(a, b)| (a.to_string(), b.to_string())).collect();
    for &ep in expected_patterns {
        if let Some(alt) = alternative_hydrogen_pattern(ep) {
            syn.push((alt, ep.to_string()));
        }
    }
    let snapshot = syn.clone();
    for (sp, ep) in snapshot {
        if let Some(alt) = alternative_hydrogen_pattern(&sp) {
            syn.push((alt, ep));
        }
    }
    let mut synonyms = FxHashMap::default();
    for (sp, ep) in &syn {
        for h in ["H", "D"] {
            let name = sp.replace('h', h);
            let same = name == *sp;
            synonyms.insert(name, ep.replace('h', h));
            if same {
                break;
            }
        }
    }
    Interpreter { expected, synonyms, meps: meps.to_vec() }
}

const PEP: &[&str] = &["N", "h", "1h", "2h", "3h", "CA", "C", "O", "OXT", "hXT"];
const PEP_SYN: &[(&str, &str)] = &[
    ("OT1", "O"),
    ("OT2", "OXT"),
    ("OC", "OXT"),
    ("hC", "hXT"),
    ("hN", "h"),
    ("1hN", "1h"),
    ("2hN", "2h"),
    ("3hN", "3h"),
    ("1hT", "1h"),
    ("2hT", "2h"),
    ("3hT", "3h"),
    ("h0A", "1h"),
    ("h0B", "2h"),
    ("h0C", "3h"),
];
const MB: [&str; 3] = ["1hB", "2hB", "3hB"];
const MG: [&str; 3] = ["1hG", "2hG", "3hG"];
const MD: [&str; 3] = ["1hD", "2hD", "3hD"];
const ME: [&str; 3] = ["1hE", "2hE", "3hE"];

fn interpreters() -> &'static FxHashMap<&'static str, Interpreter> {
    static I: OnceLock<FxHashMap<&'static str, Interpreter>> = OnceLock::new();
    I.get_or_init(|| {
        let mut m = FxHashMap::default();
        let mk = |extra: &[&str], syn: &[(&str, &str)], meps: &[[&'static str; 3]]| {
            let mut e: Vec<&str> = PEP.to_vec();
            e.extend_from_slice(extra);
            let mut s: Vec<(&str, &str)> = PEP_SYN.to_vec();
            s.extend_from_slice(syn);
            build(&e, &s, meps)
        };
        m.insert("GLY", mk(&["1hA", "2hA", "3hA"], &[], &[["1hA", "2hA", "3hA"]]));
        m.insert("ALA", mk(&["hA", "CB", "1hB", "2hB", "3hB"], &[], &[]));
        m.insert(
            "VAL",
            mk(&["hA", "CB", "hB", "CG1", "1hG1", "2hG1", "3hG1", "CG2", "1hG2", "2hG2", "3hG2"], &[], &[]),
        );
        m.insert(
            "LEU",
            mk(
                &["hA", "CB", "1hB", "2hB", "3hB", "CG", "hG", "CD1", "1hD1", "2hD1", "3hD1", "CD2", "1hD2", "2hD2", "3hD2"],
                &[("1hG", "hG")],
                &[MB],
            ),
        );
        m.insert(
            "ILE",
            mk(
                &[
                    "hA", "CB", "hB", "CG1", "1hG1", "2hG1", "3hG1", "CG2", "1hG2", "2hG2", "3hG2", "CD1", "1hD1", "2hD1", "3hD1",
                ],
                &[("CD", "CD1"), ("1hD", "1hD1"), ("2hD", "2hD1"), ("3hD", "3hD1")],
                &[["1hG1", "2hG1", "3hG1"]],
            ),
        );
        m.insert(
            "MET",
            mk(&["hA", "CB", "1hB", "2hB", "3hB", "CG", "1hG", "2hG", "3hG", "SD", "CE", "1hE", "2hE", "3hE"], &[], &[MB, MG]),
        );
        m.insert(
            "MSE",
            mk(
                &["hA", "CB", "1hB", "2hB", "3hB", "CG", "1hG", "2hG", "3hG", "SE", "CE", "1hE", "2hE", "3hE"],
                &[("SED", "SE")],
                &[MB, MG],
            ),
        );
        m.insert(
            "PRO",
            mk(&["hA", "CB", "1hB", "2hB", "3hB", "CG", "1hG", "2hG", "3hG", "CD", "1hD", "2hD", "3hD"], &[], &[MB, MG, MD]),
        );
        m.insert(
            "PHE",
            mk(
                &["hA", "CB", "1hB", "2hB", "3hB", "CG", "CD1", "hD1", "CD2", "hD2", "CE1", "hE1", "CE2", "hE2", "CZ", "hZ"],
                &[("1hZ", "hZ")],
                &[MB],
            ),
        );
        m.insert(
            "TRP",
            mk(
                &[
                    "hA", "CB", "1hB", "2hB", "3hB", "CG", "CD1", "hD1", "CD2", "NE1", "hE1", "CE2", "CE3", "hE3", "CZ2", "hZ2", "CZ3",
                    "hZ3", "CH2", "hH2",
                ],
                &[],
                &[MB],
            ),
        );
        m.insert("SER", mk(&["hA", "CB", "1hB", "2hB", "3hB", "OG", "hG"], &[], &[MB]));
        m.insert("THR", mk(&["hA", "CB", "hB", "OG1", "hG1", "CG2", "1hG2", "2hG2", "3hG2"], &[], &[]));
        m.insert("ASN", mk(&["hA", "CB", "1hB", "2hB", "3hB", "CG", "OD1", "ND2", "1hD2", "2hD2"], &[], &[MB]));
        m.insert(
            "GLN",
            mk(&["hA", "CB", "1hB", "2hB", "3hB", "CG", "1hG", "2hG", "3hG", "CD", "OE1", "NE2", "1hE2", "2hE2"], &[], &[MB, MG]),
        );
        m.insert(
            "TYR",
            mk(
                &["hA", "CB", "1hB", "2hB", "3hB", "CG", "CD1", "hD1", "CD2", "hD2", "CE1", "hE1", "CE2", "hE2", "CZ", "OH", "hH"],
                &[],
                &[MB],
            ),
        );
        m.insert("CYS", mk(&["hA", "CB", "1hB", "2hB", "3hB", "SG", "hG"], &[("1hG", "hG")], &[MB]));
        m.insert(
            "LYS",
            mk(
                &[
                    "hA", "CB", "1hB", "2hB", "3hB", "CG", "1hG", "2hG", "3hG", "CD", "1hD", "2hD", "3hD", "CE", "1hE", "2hE", "3hE", "NZ",
                    "1hZ", "2hZ", "3hZ",
                ],
                &[],
                &[MB, MG, MD, ME],
            ),
        );
        m.insert(
            "ARG",
            mk(
                &[
                    "hA", "CB", "1hB", "2hB", "3hB", "CG", "1hG", "2hG", "3hG", "CD", "1hD", "2hD", "3hD", "NE", "hE", "CZ", "NH1", "1hH1",
                    "2hH1", "NH2", "1hH2", "2hH2",
                ],
                &[],
                &[MB, MG, MD],
            ),
        );
        m.insert(
            "HIS",
            mk(
                &["hA", "CB", "1hB", "2hB", "3hB", "CG", "ND1", "hD1", "CD2", "hD2", "CE1", "hE1", "NE2", "hE2"],
                &[],
                &[MB],
            ),
        );
        m.insert(
            "ASP",
            mk(&["hA", "CB", "1hB", "2hB", "3hB", "CG", "OD1", "hD1", "OD2", "hD2"], &[], &[MB, ["hD1", "hD2", "hD2"]]),
        );
        m.insert(
            "GLU",
            mk(&["hA", "CB", "1hB", "2hB", "3hB", "CG", "1hG", "2hG", "3hG", "CD", "OE1", "hE1", "OE2", "hE2"], &[], &[MB, MG]),
        );
        m
    })
}

pub fn has_protein_interpreter(resname: &str) -> bool {
    interpreters().contains_key(resname)
}

/// `interpreter.match_atom_names(names).mon_lib_names()` for a protein residue.
pub fn protein_mon_lib_names(resname: &str, names: &[String]) -> Option<Vec<Option<String>>> {
    let it = interpreters().get(resname)?;
    let mut matched: FxHashMap<String, Vec<usize>> = FxHashMap::default();
    for (i, n) in names.iter().enumerate() {
        let name = n.trim().to_ascii_uppercase();
        let key = it.synonyms.get(&name).cloned().unwrap_or(name);
        if let Some(pat) = it.expected.get(&key) {
            matched.entry(pat.clone()).or_default().push(i);
        }
    }
    let mut transl: FxHashMap<&str, &str> = FxHashMap::default();
    for mep in &it.meps {
        if matched.contains_key(mep[2]) {
            transl.insert(mep[1], mep[0]);
            transl.insert(mep[2], mep[1]);
        } else {
            transl.insert(mep[0], mep[0]);
            transl.insert(mep[1], mep[1]);
        }
    }
    let mut out = vec![None; names.len()];
    for (pat, idxs) in &matched {
        let p = transl.get(pat.as_str()).copied().unwrap_or(pat.as_str());
        let mut ml = p.to_ascii_uppercase();
        if ml.as_bytes()[0].is_ascii_digit() && ml.as_bytes()[0] != b'0' {
            ml = format!("{}{}", &ml[1..], &ml[..1]);
        }
        for &i in idxs {
            out[i] = Some(ml.clone());
        }
    }
    Some(out)
}

/// The D-amino-acid names that are interpreted with the L dictionary.
pub fn l_given_d(resname: &str) -> Option<&'static str> {
    Some(match resname {
        "DAL" => "ALA",
        "DAR" => "ARG",
        "DAS" => "ASP",
        "DCY" => "CYS",
        "DGL" => "GLU",
        "DGN" => "GLN",
        "DHI" => "HIS",
        "DIL" => "ILE",
        "DLE" => "LEU",
        "DLY" => "LYS",
        "DPN" => "PHE",
        "DPR" => "PRO",
        "DSG" => "ASN",
        "DSN" => "SER",
        "DTH" => "THR",
        "DTR" => "TRP",
        "DTY" => "TYR",
        "DVA" => "VAL",
        "MED" => "MET",
        _ => return None,
    })
}

/// RNA/DNA reference residue name for interpretation, if this is a standard
/// nucleotide name. `has_o2` decides RNA vs DNA for A/C/G.
pub fn rna_dna_mon_lib_name(resname: &str, has_o2: bool) -> Option<&'static str> {
    let r = resname.trim().trim_start_matches('+');
    Some(match r {
        "A" | "ADE" => {
            if has_o2 { "A" } else { "AD" }
        }
        "C" | "CYT" => {
            if has_o2 { "C" } else { "CD" }
        }
        "G" | "GUA" => {
            if has_o2 { "G" } else { "GD" }
        }
        "U" | "URI" => "U",
        "T" | "THY" | "DT" | "TD" => "TD",
        "DA" | "AD" => "AD",
        "DC" | "CD" => "CD",
        "DG" | "GD" => "GD",
        _ => return None,
    })
}

/// Normalize a nucleotide atom name to the v3 reference name.
pub fn na_reference_name(name: &str) -> String {
    let n = name.trim().to_ascii_uppercase();
    let fixed = match n.as_str() {
        "O1P" => "OP1".to_string(),
        "O2P" => "OP2".to_string(),
        "O3P" => "OP3".to_string(),
        "O3T" => "OP3".to_string(),
        "1H5*" | "H5*1" | "H5'1" | "1H5'" => "H5'".to_string(),
        "2H5*" | "H5*2" | "H5'2" | "2H5'" => "H5''".to_string(),
        "1H2*" | "H2*1" | "H2'1" | "1H2'" => "H2'".to_string(),
        "2H2*" | "H2*2" | "H2'2" | "2H2'" => "H2''".to_string(),
        "2HO*" | "HO2*" | "2HO'" => "HO2'".to_string(),
        "3HO*" | "HO3*" | "H3T" => "HO3'".to_string(),
        "5HO*" | "HO5*" | "H5T" => "HO5'".to_string(),
        _ => n.replace('*', "'"),
    };
    fixed
}

/// `iotbx.pdb.rna_dna_atom_names_backbone_aliases` (tools/gen_na_aliases.py):
/// backbone atom name in any spelling -> reference name (stripped).
const NA_BACKBONE_ALIASES: &[(&str, &str)] = &[
    ("1D2'", "H2'"),
    ("1D2*", "H2'"),
    ("1D5'", "H5'"),
    ("1D5*", "H5'"),
    ("1H2'", "H2'"),
    ("1H2*", "H2'"),
    ("1H5'", "H5'"),
    ("1H5*", "H5'"),
    ("2D2'", "H2''"),
    ("2D2*", "H2''"),
    ("2D5'", "H5''"),
    ("2D5*", "H5''"),
    ("2DO'", "HO2'"),
    ("2DO*", "HO2'"),
    ("2H2'", "H2''"),
    ("2H2*", "H2''"),
    ("2H5'", "H5''"),
    ("2H5*", "H5''"),
    ("2HO'", "HO2'"),
    ("2HO*", "HO2'"),
    ("3DOP", "HOP3"),
    ("3HOP", "HOP3"),
    ("C1'", "C1'"),
    ("C1*", "C1'"),
    ("C2'", "C2'"),
    ("C2*", "C2'"),
    ("C3'", "C3'"),
    ("C3*", "C3'"),
    ("C4'", "C4'"),
    ("C4*", "C4'"),
    ("C5'", "C5'"),
    ("C5*", "C5'"),
    ("D1'", "H1'"),
    ("D1*", "H1'"),
    ("D2'", "H2'"),
    ("D2''", "H2''"),
    ("D2'1", "H2'"),
    ("D2'2", "H2''"),
    ("D2*", "H2'"),
    ("D2*1", "H2'"),
    ("D2*2", "H2''"),
    ("D3'", "H3'"),
    ("D3*", "H3'"),
    ("D3T", "HO3'"),
    ("D4'", "H4'"),
    ("D4*", "H4'"),
    ("D5'", "H5'"),
    ("D5''", "H5''"),
    ("D5'1", "H5'"),
    ("D5'2", "H5''"),
    ("D5*", "HO5'"),
    ("D5*1", "H5'"),
    ("D5*2", "H5''"),
    ("D5T", "HO5'"),
    ("DO2'", "HO2'"),
    ("DO2*", "HO2'"),
    ("H1'", "H1'"),
    ("H1*", "H1'"),
    ("H2'", "H2'"),
    ("H2''", "H2''"),
    ("H2'1", "H2'"),
    ("H2'2", "H2''"),
    ("H2*", "H2'"),
    ("H2*1", "H2'"),
    ("H2*2", "H2''"),
    ("H3'", "H3'"),
    ("H3*", "H3'"),
    ("H3T", "HO3'"),
    ("H4'", "H4'"),
    ("H4*", "H4'"),
    ("H5'", "H5'"),
    ("H5''", "H5''"),
    ("H5'1", "H5'"),
    ("H5'2", "H5''"),
    ("H5*", "HO5'"),
    ("H5*1", "H5'"),
    ("H5*2", "H5''"),
    ("H5T", "HO5'"),
    ("HO2'", "HO2'"),
    ("HO2*", "HO2'"),
    ("HO3'", "HO3'"),
    ("HO3*", "HO3'"),
    ("HO5'", "HO5'"),
    ("HO5*", "HO5'"),
    ("HOP3", "HOP3"),
    ("O1P", "OP1"),
    ("O2'", "O2'"),
    ("O2*", "O2'"),
    ("O2P", "OP2"),
    ("O3'", "O3'"),
    ("O3*", "O3'"),
    ("O3P", "OP3"),
    ("O3T", "OP3"),
    ("O4'", "O4'"),
    ("O4*", "O4'"),
    ("O5'", "O5'"),
    ("O5*", "O5'"),
    ("O5T", "OP3"),
    ("OP1", "OP1"),
    ("OP2", "OP2"),
    ("OP3", "OP3"),
    ("P", "P"),
];

pub fn na_backbone_alias(name: &str) -> Option<&'static str> {
    NA_BACKBONE_ALIASES.iter().find(|(a, _)| *a == name).map(|&(_, r)| r)
}

/// Map the atoms of one residue to dictionary ids.
///
/// Returns, per atom, `Some(id)` when the atom is an expected dictionary atom
/// (first occurrence), or `None` for unexpected/duplicate atoms.
pub struct Mapping {
    pub ids: Vec<Option<String>>,
    /// mapped names that were not found in the dictionary (unexpected)
    pub unexpected: Vec<String>,
}

pub fn map_atoms(comp: &Comp, model_names: &[String], ani: Option<&[Option<String>]>, is_na: bool, atom_synonyms: Option<&FxHashMap<String, String>>) -> Mapping {
    let mut ids = vec![None; model_names.len()];
    let mut unexpected = Vec::new();
    let mut seen: rustc_hash::FxHashSet<String> = rustc_hash::FxHashSet::default();
    let all_upper = comp.atoms.iter().all(|a| a.id == a.id.to_ascii_uppercase());
    let replace_primes = ani.is_none()
        && (is_na || {
            let np = model_names.iter().filter(|n| n.contains('\'')).count();
            let ns = model_names.iter().filter(|n| n.contains('*')).count();
            np > 0 && ns == 0
        });
    for (i, given0) in model_names.iter().enumerate() {
        let mut given: String = given0.chars().filter(|c| *c != ' ').collect();
        if all_upper {
            given = given.to_ascii_uppercase();
        }
        let mut name = ani.and_then(|a| a[i].clone()).unwrap_or_else(|| given.clone());
        if !comp.has_atom(&name) {
            let mut cands: Vec<String> = Vec::new();
            let rot = |s: &str| -> Option<String> {
                let b = s.as_bytes();
                if b.len() > 1 && b[0].is_ascii_digit() {
                    Some(format!("{}{}", &s[1..], &s[..1]))
                } else if b.len() > 1 && b[b.len() - 1].is_ascii_digit() {
                    Some(format!("{}{}", &s[s.len() - 1..], &s[..s.len() - 1]))
                } else {
                    None
                }
            };
            if let Some(r) = rot(&name) {
                cands.push(r);
            }
            if replace_primes {
                let p = name.replace('\'', "*");
                if let Some(r) = rot(&p) {
                    cands.push(r);
                }
                cands.push(p);
            }
            let mut found = cands.into_iter().find(|c| comp.has_atom(c));
            if found.is_none() {
                if comp.has_atom(&given) {
                    found = Some(given.clone());
                } else if let Some(h) = comp.hd_alias(&given) {
                    found = Some(h);
                } else if let Some(s) = atom_synonyms.and_then(|m| m.get(&given)) {
                    found = Some(s.clone());
                }
            }
            if let Some(f) = found {
                name = f;
            }
            if !comp.has_atom(&name) && is_na {
                // cctbx: reference name of the given spelling -> the
                // dictionary's spelling of that reference (last one wins)
                if let Some(r) = na_backbone_alias(&name) {
                    if let Some(c) = comp.atoms.iter().rev().find(|a| na_backbone_alias(&a.id) == Some(r)) {
                        name = c.id.clone();
                    }
                }
            }
            if !comp.has_atom(&name) && is_na {
                // backbone v2/v3 aliases (and the thymine methyl)
                let r = na_reference_name(&given);
                let alts: Vec<String> = match r.as_str() {
                    "OP3" => vec!["OP3".into(), "O3T".into()],
                    "C7" => vec!["C7".into(), "C5M".into()],
                    "H71" => vec!["H71".into(), "H5M1".into()],
                    "H72" => vec!["H72".into(), "H5M2".into()],
                    "H73" => vec!["H73".into(), "H5M3".into()],
                    _ => vec![r.clone()],
                };
                if let Some(a) = alts.into_iter().find(|a| comp.has_atom(a)) {
                    name = a;
                } else if r == "OP3" {
                    name = "O3T".into();
                }
            }
        }
        if comp.has_atom(&name) && !seen.contains(&name) {
            seen.insert(name.clone());
            ids[i] = Some(name);
        } else if !seen.contains(&name) {
            seen.insert(name.clone());
            unexpected.push(name);
        }
    }
    Mapping { ids, unexpected }
}
