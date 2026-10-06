//! Macromolecular hierarchy (model → chain → residue group → atom group → atom)
//! with the semantics of `iotbx.pdb.hierarchy` that Reduce2 relies on:
//! atom-group ordering, blank-altloc handling, conformers and atom sorting.

use crate::geom::Vec3;
use crate::resclass::{self, ResClass};

#[derive(Clone, Debug)]
pub struct Atom {
    /// Atom name as stored by iotbx: the 4-column PDB field for PDB input
    /// (e.g. `" CA "`), or the bare name for mmCIF input.
    pub name: String,
    /// Upper-case element symbol, right-justified in two columns (e.g. `" C"`).
    pub element: String,
    /// Two-column charge field such as `"1+"` or `"  "`.
    pub charge: String,
    pub serial: String,
    pub xyz: Vec3,
    pub occ: f64,
    pub b: f64,
    pub segid: String,
    pub hetero: bool,
    pub uij: Option<[f64; 6]>,
    /// Index in the flattened atom array (set by `Structure::reset_i_seq`).
    pub i_seq: usize,
    /// Where the atom came from: its row in the source `_atom_site` loop, or
    /// its ATOM/HETATM record number in a PDB file; [`Atom::NEW`] if added.
    pub src: u32,
}

impl Atom {
    pub fn new(name: &str, element: &str, xyz: Vec3) -> Atom {
        Atom {
            name: name.to_string(),
            element: format!("{:>2}", element.trim().to_ascii_uppercase()),
            charge: "  ".to_string(),
            serial: String::new(),
            xyz,
            occ: 1.0,
            b: 0.0,
            segid: String::new(),
            hetero: false,
            uij: None,
            i_seq: 0,
            src: Atom::NEW,
        }
    }

    /// `src` of an atom that is not from the input.
    pub const NEW: u32 = u32::MAX;
    #[inline]
    pub fn name_trim(&self) -> &str {
        self.name.trim()
    }
    /// Stripped upper-case element symbol.
    #[inline]
    pub fn elem(&self) -> &str {
        self.element.trim()
    }
    /// `atom.element_is_hydrogen()` from iotbx: H or D.
    #[inline]
    pub fn is_hydrogen(&self) -> bool {
        let e = self.element.as_bytes();
        match e.len() {
            0 => false,
            1 => e[0] == b'H' || e[0] == b'D',
            _ => {
                (e[0] == b' ' && (e[1] == b'H' || e[1] == b'D'))
                    || ((e[0] == b'H' || e[0] == b'D') && e[1] == b' ')
            }
        }
    }
    pub fn is_positive_ion(&self) -> bool {
        resclass::element_is_positive_ion(self.elem())
    }
    pub fn is_ion(&self) -> bool {
        resclass::element_is_ion(self.elem())
    }
    /// Integer formal charge from the charge field (`atom_charge` in probe).
    pub fn charge_value(&self) -> i32 {
        match charge_tidy(&self.charge) {
            Some(s) if s.len() == 2 => {
                let b = s.as_bytes();
                let mag = (b[0] as i32) - ('0' as i32);
                if b[1] == b'-' { -mag } else { mag }
            }
            _ => 0,
        }
    }
}

/// iotbx `atom::charge_tidy(strip=true)`.
pub fn charge_tidy(charge: &str) -> Option<String> {
    let b = charge.as_bytes();
    let mut c = [b' ', b' '];
    if b.is_empty() {
        return Some(String::new());
    }
    c[0] = b[0];
    if b.len() > 1 {
        c[1] = b[1];
    }
    if (c[0] == b' ' || c[0] == b'0') && (c[1] == b' ' || c[1] == b'0') {
        return Some(String::new());
    }
    if c[0] == b'+' || c[0] == b'-' {
        if c[1] == c[0] {
            return Some(format!("2{}", c[0] as char));
        }
        c.swap(0, 1);
    }
    if c[1] == b'+' || c[1] == b'-' {
        if c[0] == b' ' {
            return Some(format!("1{}", c[1] as char));
        }
        if c[0].is_ascii_digit() {
            return Some(format!("{}{}", c[0] as char, c[1] as char));
        }
    }
    None
}

#[derive(Clone, Debug)]
pub struct AtomGroup {
    /// "" for the main conformation; otherwise a single character.
    pub altloc: String,
    pub resname: String,
    pub atoms: Vec<Atom>,
}

impl AtomGroup {
    pub fn get_atom(&self, name: &str) -> Option<&Atom> {
        self.atoms.iter().find(|a| a.name_trim() == name)
    }
    pub fn resname_trim(&self) -> &str {
        self.resname.trim()
    }
}

#[derive(Clone, Debug)]
pub struct ResidueGroup {
    /// Right-justified 4-column residue number (may be hybrid-36).
    pub resseq: String,
    pub icode: String,
    pub link_to_previous: bool,
    pub atom_groups: Vec<AtomGroup>,
}

impl ResidueGroup {
    pub fn resseq_as_int(&self) -> i32 {
        hy36_decode(&self.resseq).unwrap_or(0)
    }
    pub fn icode_trim(&self) -> &str {
        self.icode.trim()
    }
    /// Residue id string as iotbx `rg.resid()`: resseq + icode (5 columns).
    pub fn resid(&self) -> String {
        format!("{:>4}{}", self.resseq, if self.icode.is_empty() { " " } else { &self.icode })
    }
    pub fn atoms(&self) -> impl Iterator<Item = &Atom> {
        self.atom_groups.iter().flat_map(|ag| ag.atoms.iter())
    }
}

#[derive(Clone, Debug)]
pub struct Chain {
    pub id: String,
    pub residue_groups: Vec<ResidueGroup>,
}

impl Chain {
    /// `chain.is_polymer_chain()`: decides whether a TER card follows the chain.
    pub fn is_polymer_chain(&self) -> bool {
        let (mut n_poly, mut n_unk, mut n_non) = (0, 0, 0);
        for rg in &self.residue_groups {
            let Some(ag) = rg.atom_groups.first() else { continue };
            let cls = resclass::get_class_ext(&ag.resname, true);
            match cls {
                ResClass::CommonAminoAcid
                | ResClass::DAminoAcid
                | ResClass::ModifiedAminoAcid
                | ResClass::CommonRnaDna
                | ResClass::ModifiedRnaDna
                | ResClass::Ccp4MonLibRnaDna => n_poly += 1,
                ResClass::Other | ResClass::CommonElement => n_non += 1,
                _ => {}
            }
            if ag.resname == "UNK" {
                n_unk += 1;
            }
        }
        n_poly > n_non || n_unk > n_non
    }

    /// Altlocs of the chain's conformers in order of first appearance; `[""]`
    /// when the chain has no alternate conformations (iotbx `chain.conformers()`).
    pub fn conformer_altlocs(&self) -> Vec<String> {
        let mut alts: Vec<String> = Vec::new();
        let mut any = false;
        for rg in &self.residue_groups {
            for ag in &rg.atom_groups {
                any = true;
                if ag.altloc.is_empty() {
                    continue;
                }
                if !alts.contains(&ag.altloc) {
                    alts.push(ag.altloc.clone());
                }
            }
        }
        if !any {
            return vec![];
        }
        if alts.is_empty() {
            alts.push(String::new());
        }
        alts
    }
}

#[derive(Clone, Debug)]
pub struct Model {
    pub id: String,
    pub chains: Vec<Chain>,
}

#[derive(Clone, Debug, Default)]
pub struct CrystalSymmetry {
    pub cell: [f64; 6],
    pub space_group: String,
}

/// LINK / SSBOND / struct_conn style covalent annotation.
#[derive(Clone, Debug)]
pub struct LinkRecord {
    pub atom1: AtomLabel,
    pub atom2: AtomLabel,
    pub distance: Option<f64>,
    pub kind: String,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AtomLabel {
    pub name: String,
    pub altloc: String,
    pub resname: String,
    pub chain: String,
    pub resseq: String,
    pub icode: String,
}

/// Records of an input PDB file that refer to atoms, which fixed mode writes
/// back (Reduce2 drops them).
#[derive(Clone, Debug, Default)]
pub struct PdbRecords {
    /// SSBOND and LINK lines, verbatim, in input order.
    pub links: Vec<String>,
    /// CONECT lines, verbatim.
    pub conect: Vec<String>,
    /// The input serial of each ATOM/HETATM record (indexed by `Atom::src`).
    pub serials: Vec<String>,
}

#[derive(Clone, Debug, Default)]
pub struct Structure {
    pub models: Vec<Model>,
    pub crystal: Option<CrystalSymmetry>,
    /// True when the input carried a crystal symmetry record (even a dummy one).
    pub had_cell_record: bool,
    pub links: Vec<LinkRecord>,
    pub pdb_records: PdbRecords,
}

/// Index path of an atom inside a structure.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct AtomPath {
    pub model: u32,
    pub chain: u32,
    pub rg: u32,
    pub ag: u32,
    pub atom: u32,
}

impl Structure {
    pub fn atom(&self, p: AtomPath) -> &Atom {
        &self.models[p.model as usize].chains[p.chain as usize].residue_groups[p.rg as usize]
            .atom_groups[p.ag as usize]
            .atoms[p.atom as usize]
    }
    pub fn atom_mut(&mut self, p: AtomPath) -> &mut Atom {
        &mut self.models[p.model as usize].chains[p.chain as usize].residue_groups[p.rg as usize]
            .atom_groups[p.ag as usize]
            .atoms[p.atom as usize]
    }
    pub fn atom_group(&self, p: AtomPath) -> &AtomGroup {
        &self.models[p.model as usize].chains[p.chain as usize].residue_groups[p.rg as usize]
            .atom_groups[p.ag as usize]
    }
    pub fn residue_group(&self, p: AtomPath) -> &ResidueGroup {
        &self.models[p.model as usize].chains[p.chain as usize].residue_groups[p.rg as usize]
    }
    pub fn chain(&self, p: AtomPath) -> &Chain {
        &self.models[p.model as usize].chains[p.chain as usize]
    }

    /// All atom paths in hierarchy (i_seq) order.
    pub fn atom_paths(&self) -> Vec<AtomPath> {
        let mut out = Vec::new();
        for (mi, m) in self.models.iter().enumerate() {
            for (ci, c) in m.chains.iter().enumerate() {
                for (ri, rg) in c.residue_groups.iter().enumerate() {
                    for (gi, ag) in rg.atom_groups.iter().enumerate() {
                        for ai in 0..ag.atoms.len() {
                            out.push(AtomPath {
                                model: mi as u32,
                                chain: ci as u32,
                                rg: ri as u32,
                                ag: gi as u32,
                                atom: ai as u32,
                            });
                        }
                    }
                }
            }
        }
        out
    }

    pub fn atoms_size(&self) -> usize {
        self.models
            .iter()
            .flat_map(|m| m.chains.iter())
            .flat_map(|c| c.residue_groups.iter())
            .flat_map(|r| r.atom_groups.iter())
            .map(|g| g.atoms.len())
            .sum()
    }

    /// Assign `i_seq` in hierarchy order.
    pub fn reset_i_seq(&mut self) {
        let mut i = 0usize;
        for m in &mut self.models {
            for c in &mut m.chains {
                for rg in &mut c.residue_groups {
                    for ag in &mut rg.atom_groups {
                        for a in &mut ag.atoms {
                            a.i_seq = i;
                            i += 1;
                        }
                    }
                }
            }
        }
    }

    /// `atoms().reset_serial()`: serial numbers 1..N in hybrid-36.
    pub fn reset_serial(&mut self) {
        let mut i = 1i64;
        for m in &mut self.models {
            for c in &mut m.chains {
                for rg in &mut c.residue_groups {
                    for ag in &mut rg.atom_groups {
                        for a in &mut ag.atoms {
                            a.serial = hy36_encode(5, i).unwrap_or_else(|| "*****".into());
                            i += 1;
                        }
                    }
                }
            }
        }
    }

    /// iotbx `hierarchy.sort_atoms_in_place()`.
    pub fn sort_atoms_in_place(&mut self) {
        for m in &mut self.models {
            for c in &mut m.chains {
                for rg in &mut c.residue_groups {
                    for ag in &mut rg.atom_groups {
                        sort_atom_group(ag);
                    }
                }
            }
        }
    }

    /// Remove atoms for which `keep` returns false, dropping empty groups.
    pub fn retain_atoms<F: FnMut(&Atom) -> bool>(&mut self, mut keep: F) {
        for m in &mut self.models {
            for c in &mut m.chains {
                for rg in &mut c.residue_groups {
                    for ag in &mut rg.atom_groups {
                        ag.atoms.retain(|a| keep(a));
                    }
                    rg.atom_groups.retain(|ag| !ag.atoms.is_empty());
                }
                c.residue_groups.retain(|rg| !rg.atom_groups.is_empty());
            }
            m.chains.retain(|c| !c.residue_groups.is_empty());
        }
    }

    pub fn has_hydrogens(&self) -> bool {
        self.models
            .iter()
            .flat_map(|m| m.chains.iter())
            .flat_map(|c| c.residue_groups.iter())
            .flat_map(|r| r.atom_groups.iter())
            .flat_map(|g| g.atoms.iter())
            .any(|a| a.is_hydrogen())
    }
}

/// Atom-name orderings used by iotbx's `sort_atoms_in_place`.
const AA_ORDER: &[&str] = &[
    "N", "CA", "C", "O", "CB", "NB", "OB", "SB", "CB1", "NB1", "OB1", "SB1", "CB2", "NB2", "OB2",
    "SB2", "CB3", "NB3", "OB3", "SB3", "CG", "NG", "OG", "SG", "CG1", "NG1", "OG1", "SG1", "CG2",
    "NG2", "OG2", "SG2", "CG3", "NG3", "OG3", "SG3", "CD", "ND", "OD", "SD", "CD1", "ND1", "OD1",
    "SD1", "CD2", "ND2", "OD2", "SD2", "CD3", "ND3", "OD3", "SD3", "SE", "CE", "NE", "OE", "SE",
    "CE1", "NE1", "OE1", "SE1", "CE2", "NE2", "OE2", "SE2", "CE3", "NE3", "OE3", "SE3", "CZ", "NZ",
    "OZ", "SZ", "CZ1", "NZ1", "OZ1", "SZ1", "CZ2", "NZ2", "OZ2", "SZ2", "CZ3", "NZ3", "OZ3", "SZ3",
    "CH", "NH", "OH", "SH", "CH1", "NH1", "OH1", "SH1", "CH2", "NH2", "OH2", "SH2", "CH3", "NH3",
    "OH3", "SH3", "OXT", "H", "H1", "1H", "H2", "2H", "H3", "3H", "HA", "HA1", "1HA", "HA11",
    "1HA1", "HA12", "2HA1", "HA13", "3HA1", "HA2", "2HA", "HA21", "1HA2", "HA22", "2HA2", "HA23",
    "3HA2", "HA3", "3HA", "HA31", "1HA3", "HA32", "2HA3", "HA33", "3HA3", "HB", "HB1", "1HB",
    "HB11", "1HB1", "HB12", "2HB1", "HB13", "3HB1", "HB2", "2HB", "HB21", "1HB2", "HB22", "2HB2",
    "HB23", "3HB2", "HB3", "3HB", "HB31", "1HB3", "HB32", "2HB3", "HB33", "3HB3", "HG", "HG1",
    "1HG", "HG11", "1HG1", "HG12", "2HG1", "HG13", "3HG1", "HG2", "2HG", "HG21", "1HG2", "HG22",
    "2HG2", "HG23", "3HG2", "HG3", "3HG", "HG31", "1HG3", "HG32", "2HG3", "HG33", "3HG3", "HD",
    "HD1", "1HD", "HD11", "1HD1", "HD12", "2HD1", "HD13", "3HD1", "HD2", "2HD", "HD21", "1HD2",
    "HD22", "2HD2", "HD23", "3HD2", "HD3", "3HD", "HD31", "1HD3", "HD32", "2HD3", "HD33", "3HD3",
    "HE", "HE1", "1HE", "HE11", "1HE1", "HE12", "2HE1", "HE13", "3HE1", "HE2", "2HE", "HE21",
    "1HE2", "HE22", "2HE2", "HE23", "3HE2", "HE3", "3HE", "HE31", "1HE3", "HE32", "2HE3", "HE33",
    "3HE3", "HZ", "HZ1", "1HZ", "HZ11", "1HZ1", "HZ12", "2HZ1", "HZ13", "3HZ1", "HZ2", "2HZ",
    "HZ21", "1HZ2", "HZ22", "2HZ2", "HZ23", "3HZ2", "HZ3", "3HZ", "HZ31", "1HZ3", "HZ32", "2HZ3",
    "HZ33", "3HZ3", "HH", "HH1", "1HH", "HH11", "1HH1", "HH12", "2HH1", "HH13", "3HH1", "HH2",
    "2HH", "HH21", "1HH2", "HH22", "2HH2", "HH23", "3HH2", "HH3", "3HH", "HH31", "1HH3", "HH32",
    "2HH3", "HH33", "3HH3",
];

const SMALL_NA_ORDER: &[&str] = &[
    "P", "OP1", "O1P", "OP2", "O2P", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'",
    "N1", "C2", "O2", "N2", "N3", "C4", "O4", "N4", "C5", "C7", "C6", "N6", "O6", "N7", "C8", "N9",
    "H5'", "H5''", "H4'", "H3'", "HO3'", "H2'", "H2''", "HO2'", "H1'", "H8", "H41", "H42", "H5",
    "H5'1", "H5'2", "H3", "H5M1", "H5M2", "H5M3", "H71", "H72", "H73", "H6", "H61", "H62", "H1",
    "H2", "H21", "H22", "HO5'", "H2'1", "H2'2",
];

const BIG_NA_ORDER: &[&str] = &[
    "P", "OP1", "O1P", "OP2", "O2P", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'",
    "N9", "C8", "N7", "C7", "C5", "C6", "N6", "O6", "N1", "C2", "O2", "N2", "N3", "C4", "O4", "N4",
    "H5'", "H5''", "H4'", "H3'", "HO3'", "H2'", "H2''", "HO2'", "H1'", "H8", "H41", "H42", "H5",
    "H5'1", "H5'2", "H3", "H5M1", "H5M2", "H5M3", "H71", "H72", "H73", "H6", "H61", "H62", "H1",
    "H2", "H21", "H22", "HO5'", "H2'1", "H2'2",
];

fn sort_atom_group(ag: &mut AtomGroup) {
    let cls = resclass::get_class(&ag.resname);
    let order: &[&str] = if cls == ResClass::CommonRnaDna || cls == ResClass::ModifiedRnaDna {
        if ag.get_atom("N9").is_none() { SMALL_NA_ORDER } else { BIG_NA_ORDER }
    } else {
        AA_ORDER
    };
    // Precompute keys: (position in list or len, is_hydrogen-by-name, padded name).
    let n = order.len();
    let key = |a: &Atom| {
        let stripped = a.name.trim().replace('*', "'");
        let pos = order.iter().position(|x| *x == stripped).unwrap_or(n);
        let is_h = stripped.as_bytes().first() == Some(&b'H');
        (pos, is_h)
    };
    // The iotbx comparator: atoms not in the list compare by (H last, then
    // padded name); otherwise non-H before H, then by list position.
    ag.atoms.sort_by(|a1, a2| {
        let (p1, h1) = key(a1);
        let (p2, h2) = key(a2);
        use std::cmp::Ordering::*;
        if p1 == p2 {
            if h1 && !h2 {
                return Greater;
            }
            if h2 && !h1 {
                return Less;
            }
            return a1.name.cmp(&a2.name);
        }
        if h1 && !h2 {
            return Greater;
        }
        if !h1 && h2 {
            return Less;
        }
        p1.cmp(&p2)
    });
}

// ----------------------------------------------------------------------------
// Hybrid-36 encoding (PDB serial numbers and residue numbers beyond 4/5 digits).

const DIGITS_UPPER: &[u8] = b"0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
const DIGITS_LOWER: &[u8] = b"0123456789abcdefghijklmnopqrstuvwxyz";

fn encode_pure(digits: &[u8], width: usize, mut value: i64) -> String {
    let neg = value < 0;
    if neg {
        value = -value;
    }
    let base = digits.len() as i64;
    let mut buf = Vec::new();
    loop {
        buf.push(digits[(value % base) as usize]);
        value /= base;
        if value == 0 {
            break;
        }
    }
    if neg {
        buf.push(b'-');
    }
    buf.reverse();
    let s = String::from_utf8(buf).unwrap();
    format!("{:>width$}", s, width = width)
}

pub fn hy36_encode(width: usize, value: i64) -> Option<String> {
    let w = width as u32;
    let min = -(10i64.pow(w - 1)) + 1;
    let max = 10i64.pow(w) - 1;
    if value >= min && value <= max {
        return Some(format!("{:>width$}", value, width = width));
    }
    let lim = 26 * 36i64.pow(w - 1);
    let mut v = value - (max + 1);
    if v >= 0 && v < lim {
        v += 10 * 36i64.pow(w - 1);
        return Some(encode_pure(DIGITS_UPPER, width, v));
    }
    v -= lim;
    if v >= 0 && v < lim {
        v += 10 * 36i64.pow(w - 1);
        return Some(encode_pure(DIGITS_LOWER, width, v));
    }
    None
}

pub fn hy36_decode(s: &str) -> Option<i32> {
    let t = s.trim();
    if t.is_empty() {
        return Some(0);
    }
    if let Ok(v) = t.parse::<i32>() {
        return Some(v);
    }
    let width = s.len() as u32;
    let b = s.as_bytes();
    let first = b[0];
    let decode = |digits: &[u8]| -> Option<i64> {
        let mut v: i64 = 0;
        for &c in b {
            let d = digits.iter().position(|&x| x == c)? as i64;
            v = v * 36 + d;
        }
        Some(v)
    };
    if first.is_ascii_uppercase() {
        let v = decode(DIGITS_UPPER)?;
        Some((v - 10 * 36i64.pow(width - 1) + 10i64.pow(width)) as i32)
    } else if first.is_ascii_lowercase() {
        let v = decode(DIGITS_LOWER)?;
        Some((v + 16 * 36i64.pow(width - 1) + 10i64.pow(width)) as i32)
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use crate::geom::v3;

    fn atom(name: &str, element: &str) -> Atom {
        Atom::new(name, element, v3(0.0, 0.0, 0.0))
    }

    fn residue(resseq: &str, groups: Vec<(&str, &str, Vec<Atom>)>) -> ResidueGroup {
        ResidueGroup {
            resseq: resseq.to_string(),
            icode: " ".into(),
            link_to_previous: true,
            atom_groups: groups
                .into_iter()
                .map(|(altloc, resname, atoms)| AtomGroup { altloc: altloc.into(), resname: resname.into(), atoms })
                .collect(),
        }
    }

    fn layout(st: &Structure) -> Vec<String> {
        let mut out = Vec::new();
        for m in &st.models {
            for c in &m.chains {
                for rg in &c.residue_groups {
                    for ag in &rg.atom_groups {
                        for a in &ag.atoms {
                            out.push(format!("{}/{}/{}{}/{}", c.id, rg.resseq.trim(), ag.altloc, ag.resname.trim(), a.name.trim()));
                        }
                    }
                }
            }
        }
        out
    }

    #[test]
    fn held_atoms_go_back_where_they_were() {
        let chain_a = Chain {
            id: "A".into(),
            residue_groups: vec![
                residue("   1", vec![("", "ALA", vec![atom(" N  ", "N"), atom(" CA ", "C")])]),
                residue("   2", vec![("", "UNX", vec![atom(" UNK", "X")])]),
                residue("   3", vec![("", "GLY", vec![atom(" N  ", "N"), atom(" X1 ", "X")]), ("A", "GLY", vec![atom(" CA ", "C")])]),
            ],
        };
        let chain_x = Chain { id: "X".into(), residue_groups: vec![residue("   9", vec![("", "UNX", vec![atom(" UNK", "X")])])] };
        let mut st = Structure { models: vec![Model { id: String::new(), chains: vec![chain_a, chain_x] }], ..Default::default() };
        let before = layout(&st);
        let held = st.take_atoms(|_, a| a.elem() == "X");
        assert_eq!(held.len(), 3);
        assert_eq!(layout(&st), ["A/1/ALA/N", "A/1/ALA/CA", "A/3/GLY/N", "A/3/AGLY/CA"]);
        // processing adds a hydrogen meanwhile
        st.models[0].chains[0].residue_groups[0].atom_groups[0].atoms.push(atom(" H  ", "H"));
        st.restore_atoms(held);
        let mut expected = before;
        expected.insert(2, "A/1/ALA/H".into());
        assert_eq!(layout(&st), expected);
    }

    use super::*;
    #[test]
    fn hy36_roundtrip() {
        assert_eq!(hy36_encode(5, 99999).unwrap(), "99999");
        assert_eq!(hy36_encode(5, 100000).unwrap(), "A0000");
        assert_eq!(hy36_decode("A0000"), Some(100000));
        assert_eq!(hy36_encode(4, 10000).unwrap(), "A000");
        assert_eq!(hy36_decode("A000"), Some(10000));
        assert_eq!(hy36_decode("   1"), Some(1));
    }
    #[test]
    fn charge() {
        assert_eq!(charge_tidy("1+").as_deref(), Some("1+"));
        assert_eq!(charge_tidy("+1").as_deref(), Some("1+"));
        assert_eq!(charge_tidy("  ").as_deref(), Some(""));
        assert_eq!(charge_tidy("--").as_deref(), Some("2-"));
    }
}

/// Developer aid: with REDUCE3_MEM set, print elapsed time and resident memory.
pub fn mem_checkpoint(tag: &str) {
    if std::env::var_os("REDUCE3_MEM").is_none() {
        return;
    }
    static START: std::sync::OnceLock<std::time::Instant> = std::sync::OnceLock::new();
    let t0 = START.get_or_init(std::time::Instant::now);
    let pid = std::process::id().to_string();
    let rss = std::process::Command::new("ps")
        .args(["-o", "rss=", "-p", &pid])
        .output()
        .ok()
        .and_then(|o| String::from_utf8(o.stdout).ok())
        .and_then(|s| s.trim().parse::<u64>().ok())
        .unwrap_or(0);
    eprintln!("[mem] {:>8.2}s {:>8} MB  {}", t0.elapsed().as_secs_f64(), rss / 1024, tag);
}

/// Atoms taken out of a structure and put back later, in their own chains,
/// residue groups and atom groups (fixed mode keeps atoms of unknown element,
/// which Reduce2 deletes).
#[derive(Clone, Debug, Default)]
pub struct HeldAtoms {
    groups: Vec<HeldGroup>,
}

/// A place in the hierarchy: an id and which occurrence of that id it is, so
/// that repeated chain ids and residue numbers stay apart.
type Place<K> = (K, usize);

#[derive(Clone, Debug)]
struct HeldGroup {
    model: (String, usize),
    chain: Place<String>,
    chain_before: Option<Place<String>>,
    residue: Place<(String, String)>,
    residue_before: Option<Place<(String, String)>>,
    link_to_previous: bool,
    altloc: String,
    resname: String,
    atoms: Vec<Atom>,
}

impl HeldAtoms {
    pub fn len(&self) -> usize {
        self.groups.iter().map(|g| g.atoms.len()).sum()
    }
    pub fn is_empty(&self) -> bool {
        self.groups.is_empty()
    }
    /// Add atoms held from the same structure.
    pub fn extend(&mut self, other: HeldAtoms) {
        self.groups.extend(other.groups);
    }
}

fn occurrence<T, K: PartialEq>(items: &[T], key: impl Fn(&T) -> K, want: &K, n: usize) -> Option<usize> {
    items.iter().enumerate().filter(|(_, x)| key(x) == *want).map(|(i, _)| i).nth(n)
}

impl Structure {
    /// Take out the atoms `pick` selects (given each atom's group), remembering
    /// where they were.
    pub fn take_atoms(&mut self, mut pick: impl FnMut(&AtomGroup, &Atom) -> bool) -> HeldAtoms {
        let mut held = HeldAtoms::default();
        for (mi, m) in self.models.iter().enumerate() {
            let mut chain_seen: Vec<String> = Vec::new();
            let mut chain_before: Option<Place<String>> = None;
            for c in &m.chains {
                let chain = (c.id.clone(), chain_seen.iter().filter(|x| **x == c.id).count());
                chain_seen.push(c.id.clone());
                let mut rg_seen: Vec<(String, String)> = Vec::new();
                let mut residue_before: Option<Place<(String, String)>> = None;
                for rg in &c.residue_groups {
                    let key = (rg.resseq.clone(), rg.icode.clone());
                    let residue = (key.clone(), rg_seen.iter().filter(|x| **x == key).count());
                    rg_seen.push(key);
                    for ag in &rg.atom_groups {
                        let atoms: Vec<Atom> = ag.atoms.iter().filter(|a| pick(ag, a)).cloned().collect();
                        if !atoms.is_empty() {
                            held.groups.push(HeldGroup {
                                model: (m.id.clone(), mi),
                                chain: chain.clone(),
                                chain_before: chain_before.clone(),
                                residue: residue.clone(),
                                residue_before: residue_before.clone(),
                                link_to_previous: rg.link_to_previous,
                                altloc: ag.altloc.clone(),
                                resname: ag.resname.clone(),
                                atoms,
                            });
                        }
                    }
                    residue_before = Some(residue);
                }
                chain_before = Some(chain);
            }
        }
        if !held.is_empty() {
            for m in &mut self.models {
                for c in &mut m.chains {
                    for rg in &mut c.residue_groups {
                        for ag in &mut rg.atom_groups {
                            let taken: Vec<bool> = ag.atoms.iter().map(|a| pick(ag, a)).collect();
                            let mut k = 0;
                            ag.atoms.retain(|_| {
                                k += 1;
                                !taken[k - 1]
                            });
                        }
                    }
                }
            }
            self.retain_atoms(|_| true);
        }
        held
    }

    /// Put held atoms back, recreating any group that no longer exists right
    /// after the group that preceded it.
    pub fn restore_atoms(&mut self, held: HeldAtoms) {
        for g in held.groups {
            let mi = match self.models.iter().position(|m| m.id == g.model.0) {
                Some(i) => i,
                None => {
                    let at = g.model.1.min(self.models.len());
                    self.models.insert(at, Model { id: g.model.0.clone(), chains: Vec::new() });
                    at
                }
            };
            let chains = &mut self.models[mi].chains;
            let ci = match occurrence(chains, |c| c.id.clone(), &g.chain.0, g.chain.1) {
                Some(i) => i,
                None => {
                    let at = g
                        .chain_before
                        .as_ref()
                        .and_then(|(id, n)| occurrence(chains, |c| c.id.clone(), id, *n))
                        .map_or(0, |i| i + 1);
                    chains.insert(at, Chain { id: g.chain.0.clone(), residue_groups: Vec::new() });
                    at
                }
            };
            let rgs = &mut chains[ci].residue_groups;
            let key = |r: &ResidueGroup| (r.resseq.clone(), r.icode.clone());
            let ri = match occurrence(rgs, key, &g.residue.0, g.residue.1) {
                Some(i) => i,
                None => {
                    let at = g.residue_before.as_ref().and_then(|(k, n)| occurrence(rgs, key, k, *n)).map_or(0, |i| i + 1);
                    rgs.insert(
                        at,
                        ResidueGroup {
                            resseq: g.residue.0 .0.clone(),
                            icode: g.residue.0 .1.clone(),
                            link_to_previous: g.link_to_previous,
                            atom_groups: Vec::new(),
                        },
                    );
                    at
                }
            };
            let ags = &mut rgs[ri].atom_groups;
            match ags.iter_mut().find(|ag| ag.altloc == g.altloc && ag.resname == g.resname) {
                Some(ag) => ag.atoms.extend(g.atoms),
                None => {
                    let ag = AtomGroup { altloc: g.altloc, resname: g.resname, atoms: g.atoms };
                    if ag.altloc.is_empty() {
                        ags.insert(0, ag);
                    } else {
                        ags.push(ag);
                    }
                }
            }
        }
    }
}
