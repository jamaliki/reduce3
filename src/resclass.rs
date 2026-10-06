//! Residue classification identical to `iotbx.pdb.common_residue_names_get_class`.

use crate::resclass_data as d;
use rustc_hash::FxHashSet;
use std::sync::OnceLock;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum ResClass {
    CommonAminoAcid,
    DAminoAcid,
    ModifiedAminoAcid,
    CommonRnaDna,
    ModifiedRnaDna,
    Ccp4MonLibRnaDna,
    CommonWater,
    CommonSmallMolecule,
    CommonSaccharide,
    CommonElement,
    Other,
}

impl ResClass {
    pub fn name(self) -> &'static str {
        match self {
            ResClass::CommonAminoAcid => "common_amino_acid",
            ResClass::DAminoAcid => "d_amino_acid",
            ResClass::ModifiedAminoAcid => "modified_amino_acid",
            ResClass::CommonRnaDna => "common_rna_dna",
            ResClass::ModifiedRnaDna => "modified_rna_dna",
            ResClass::Ccp4MonLibRnaDna => "ccp4_mon_lib_rna_dna",
            ResClass::CommonWater => "common_water",
            ResClass::CommonSmallMolecule => "common_small_molecule",
            ResClass::CommonSaccharide => "common_saccharide",
            ResClass::CommonElement => "common_element",
            ResClass::Other => "other",
        }
    }
    pub fn is_amino_acid(self) -> bool {
        matches!(self, ResClass::CommonAminoAcid | ResClass::ModifiedAminoAcid | ResClass::DAminoAcid)
    }
    pub fn is_rna_dna(self) -> bool {
        matches!(self, ResClass::CommonRnaDna | ResClass::ModifiedRnaDna)
    }
}

struct Sets {
    aa: FxHashSet<&'static str>,
    daa: FxHashSet<&'static str>,
    maa: FxHashSet<&'static str>,
    na: FxHashSet<&'static str>,
    mna: FxHashSet<&'static str>,
    mlna: FxHashSet<&'static str>,
    water: FxHashSet<&'static str>,
    small: FxHashSet<&'static str>,
    sacch: FxHashSet<&'static str>,
    elem: FxHashSet<&'static str>,
}

fn sets() -> &'static Sets {
    static S: OnceLock<Sets> = OnceLock::new();
    S.get_or_init(|| {
        let mk = |a: &'static [&'static str]| a.iter().copied().collect::<FxHashSet<_>>();
        Sets {
            aa: mk(d::AMINO_ACID),
            daa: mk(d::D_AMINO_ACID),
            maa: mk(d::MODIFIED_AMINO_ACID),
            na: mk(d::RNA_DNA),
            mna: mk(d::MODIFIED_RNA_DNA),
            mlna: mk(d::CCP4_MON_LIB_RNA_DNA),
            water: mk(d::WATER),
            small: mk(d::SMALL_MOLECULE),
            sacch: mk(d::COMMON_SACCHARIDE),
            elem: mk(d::ELEMENT),
        }
    })
}

/// Classify a residue name. Names shorter than three characters are padded on
/// the left, exactly like the C++ implementation.
pub fn get_class(name: &str) -> ResClass {
    get_class_ext(name, false)
}

pub fn get_class_ext(name: &str, consider_ccp4_mon_lib_rna_dna: bool) -> ResClass {
    let mut buf = String::with_capacity(4);
    if name.len() < 3 {
        for _ in 0..(3 - name.len()) {
            buf.push(' ');
        }
    }
    buf.push_str(name);
    let p = buf.as_str();
    let s = sets();
    if s.aa.contains(p) {
        ResClass::CommonAminoAcid
    } else if s.daa.contains(p) {
        ResClass::DAminoAcid
    } else if s.maa.contains(p) {
        ResClass::ModifiedAminoAcid
    } else if s.na.contains(p) {
        ResClass::CommonRnaDna
    } else if s.mna.contains(p) {
        ResClass::ModifiedRnaDna
    } else if consider_ccp4_mon_lib_rna_dna && s.mlna.contains(p) {
        ResClass::Ccp4MonLibRnaDna
    } else if s.water.contains(p) {
        ResClass::CommonWater
    } else if s.small.contains(p) {
        ResClass::CommonSmallMolecule
    } else if s.sacch.contains(p) {
        ResClass::CommonSaccharide
    } else if s.elem.contains(p) {
        ResClass::CommonElement
    } else {
        ResClass::Other
    }
}

/// The 20 standard amino acids plus MSE etc. as used by
/// `iotbx.pdb.amino_acid_codes.one_letter_given_three_letter`.
pub const ONE_LETTER_GIVEN_THREE_LETTER: &[(&str, char)] = &[
    ("ALA", 'A'), ("ARG", 'R'), ("ASN", 'N'), ("ASP", 'D'), ("CYS", 'C'), ("GLN", 'Q'),
    ("GLU", 'E'), ("GLY", 'G'), ("HIS", 'H'), ("ILE", 'I'), ("LEU", 'L'), ("LYS", 'K'),
    ("MET", 'M'), ("MSE", 'M'), ("PHE", 'F'), ("PRO", 'P'), ("PYL", 'O'), ("SEC", 'U'),
    ("SER", 'S'), ("THR", 'T'), ("TRP", 'W'), ("TYR", 'Y'), ("VAL", 'V'), ("UNK", 'X'),
];

/// RNA/DNA names from `iotbx.pdb.nucleic_acid_codes` (rna_one_letter_code_dict
/// union dna_one_letter_code_dict), used by Reduce2's kinemage helpers.
pub const NUCLEIC_ACID_RESNAMES: &[&str] =
    &["A", "ADE", "C", "CYT", "G", "GUA", "U", "URI", "DA", "DC", "DG", "DT", "T", "THY"];

pub fn is_standard_amino_acid_code(resname: &str) -> bool {
    ONE_LETTER_GIVEN_THREE_LETTER.iter().any(|(n, _)| *n == resname)
}

/// Elements considered positive ions by iotbx (`atom::element_is_positive_ion`).
const POSITIVE_IONS: &[&str] = &[
    "LI", "NA", "AL", "K", "MG", "CA", "MN", "FE", "CO", "NI", "CU", "ZN", "RB", "SR", "MO", "AG",
    "CD", "IN", "CS", "BA", "AU", "HG", "TL", "PB", "V", "CR", "TE", "SM", "GD", "YB", "W", "PT",
    "U", "BE", "SI", "SC", "TI", "FA", "GE", "Y", "ZR", "SN", "SB", "LA", "CE", "FR", "RA", "TH",
    "NB", "TC", "RU", "RH", "PD", "PR", "ND", "PM", "EU", "TB", "DY", "HO", "ER", "TM", "LU", "HF",
    "TA", "RE", "OS", "IR", "BI", "PO", "AT", "AC", "PA", "NP", "PU", "AM", "CM", "BK", "CF", "ES",
    "FM", "MD", "NO", "B",
];
const NEGATIVE_IONS: &[&str] = &["F", "CL", "BR", "I"];

/// `element` is the stripped, upper-case element symbol.
pub fn element_is_positive_ion(element: &str) -> bool {
    POSITIVE_IONS.contains(&element)
}
pub fn element_is_negative_ion(element: &str) -> bool {
    NEGATIVE_IONS.contains(&element)
}
pub fn element_is_ion(element: &str) -> bool {
    element_is_positive_ion(element) || element_is_negative_ion(element)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn classes() {
        assert_eq!(get_class("ALA"), ResClass::CommonAminoAcid);
        assert_eq!(get_class("HOH"), ResClass::CommonWater);
        assert_eq!(get_class("A"), ResClass::CommonRnaDna);
        assert_eq!(get_class("DA"), ResClass::CommonRnaDna);
        assert_eq!(get_class("ZN"), ResClass::CommonElement);
        assert_eq!(get_class("STR"), ResClass::Other);
    }
}
