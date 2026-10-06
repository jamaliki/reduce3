//! Automatic linking: covalent bonds between residues that are not polymer
//! neighbours (N- and O-glycosylation, glycosidic bonds, cyclic peptides,
//! ligand-protein bonds). Port of `process_nonbonded_for_links`
//! (mmtbx/monomer_library/linking_mixins.py) with the helpers from
//! linking_utils.py, linking_setup.py and glyco_utils.py, run with the
//! automatic-linking parameters Reduce2 uses (link_metals=Auto, which makes no
//! metal links here; link_residues, link_carbohydrates and link_ligands on).
//!
//! Compat mode keeps one quirk of the original: the per-residue-pair
//! bookkeeping is keyed by a slice of `atom.id_str()`, which for multi-model
//! files includes part of the model id and the atom name, so the
//! per-residue-pair limits act per atom pair there. Fixed mode keys it by
//! model and atom groups.

use crate::geom::*;
use crate::interp::{FlatAtoms, Interp, ORIGIN_GLYCO_CUSTOM, ORIGIN_METAL, ORIGIN_MISC};
use crate::linkdata;
use crate::model::*;
use crate::monlib::{ChemLink, MonLib};
use crate::resclass::{self, ONE_LETTER_GIVEN_THREE_LETTER};
use rustc_hash::FxHashMap;
use std::cell::RefCell;

const METALS: &[&str] = &[
    "ZN", "CA", "MG", "NA", "MN", "K", "FE", "CU", "CD", "HG", "NI", "CO", "SR", "CS", "PT", "BA", "TL", "PB", "SM",
    "AU", "RB", "YB", "LI", "MO", "LU", "CR", "OS", "GD", "TB", "LA", "AG", "HO", "GA", "CE", "W", "RU", "RE", "PR",
    "IR", "EU", "AL", "V", "PD", "U", "SB", "SE", "TE",
];
const NON_LINKING: &[&str] = &["H", "D", "F", "CL", "BR", "I", "AT", "HE", "NE", "AR", "KR", "XE"];
const FIRST_ROW: &[&str] = &["LI", "BE", "B", "C", "N", "O", "F"];
const SUGAR_TYPES: &[&str] = &[
    "SACCHARIDE",
    "D-SACCHARIDE",
    "L-SACCHARIDE",
    "D-SACCHARIDE, ALPHA LINKING",
    "D-SACCHARIDE, BETA LINKING",
    "L-SACCHARIDE, ALPHA LINKING",
    "L-SACCHARIDE, BETA LINKING",
];
const AMINO_TYPES: &[&str] = &["\"L-PEPTIDE LINKING\"", "\"D-PEPTIDE LINKING\"", "L-PEPTIDE LINKING", "D-PEPTIDE LINKING"];

/// `linking_class` origin ids of the named links (`link_<KEY>`).
const LINK_ORIGINS: &[(&str, u16)] = &[
    ("ACE_C-N", 16), ("AHT-ALA", 17), ("ALPHA1-2", 18), ("ALPHA1-3", 19), ("ALPHA1-4", 20), ("ALPHA1-6", 21),
    ("ALPHA2-3", 22), ("ALPHA2-6", 23), ("ASP_CG-ANY_N", 24), ("BETA1-2", 25), ("BETA1-3", 26), ("BETA1-4", 27),
    ("BETA1-6", 28), ("BETA2-3", 29), ("BOC_C-N", 30), ("BR-C5", 31), ("CH2-N2", 32), ("CH3-N1", 33),
    ("CH3-O2*", 34), ("CIS", 35), ("CYS-MPR", 36), ("DFO-NME", 37), ("DFO_C-N", 38), ("DFO_DFO", 39),
    ("DFO_N-C", 40), ("DFO_STA", 41), ("DM1-CH2", 42), ("FE-CYS", 43), ("FOR-LYZ", 44), ("FOR_C-C", 45),
    ("FOR_C-N", 46), ("ILG_CD-N", 47), ("ILG_CD-p", 48), ("IVA_C-N", 49), ("LINK_C-N", 50), ("LINK_CNp", 51),
    ("LINK_CpN", 52), ("MAN-ASN", 53), ("MAN-SER", 54), ("MAN-THR", 55), ("MG-O1P", 56), ("MG-O2P", 57),
    ("MPR-CYS", 58), ("NAG-ASN", 59), ("NAG-SER", 60), ("NAG-THR", 61), ("NH2_CTERM", 62), ("NMCIS", 63),
    ("NME_N-C", 64), ("NMTRANS", 65), ("PCIS", 66), ("PEPTIDE-PLANE", 67), ("POST-BETA-TRANS", 68),
    ("PRE-BETA-TRANS", 69), ("PTRANS", 70), ("SFN-TYR", 71), ("SS", 72), ("SSRAD", 73), ("STA-NME", 74),
    ("STA_C-N", 75), ("STA_DFO", 76), ("STA_N-C", 77), ("STA_STA", 78), ("TRANS", 79), ("XYS-ASN", 80),
    ("XYS-SER", 81), ("XYS-THR", 82), ("ZN-CYS", 83), ("gap", 84), ("p", 85), ("symmetry", 86),
];

fn link_origin(key: &str) -> Option<u16> {
    LINK_ORIGINS.iter().find(|(k, _)| *k == key).map(|&(_, o)| o)
}

// automatic_linking cutoffs as Reduce2's pdb_interpretation parameters set them
const AMINO_ACID_CUTOFF: f64 = 1.9;
const RNA_DNA_CUTOFF: f64 = 3.4;
const INTER_RESIDUE_CUTOFF: f64 = 2.2;
const CARBOHYDRATE_CUTOFF: f64 = 1.99;
const METAL_CUTOFF: f64 = 3.0;
const SULFUR_CUTOFF: f64 = 2.5;
const OTHER_CUTOFF: f64 = 2.0;
const SECOND_ROW_BUFFER: f64 = 0.5;
const MAX_BONDED_CUTOFF: f64 = 3.0;

/// `update_skip_if_longer`: squared distance limit for a sorted class pair.
fn skip_if_longer(a: &str, b: &str) -> Option<f64> {
    let v = match (a, b) {
        ("common_amino_acid", "common_amino_acid") => AMINO_ACID_CUTOFF * AMINO_ACID_CUTOFF,
        ("common_amino_acid", "other") => INTER_RESIDUE_CUTOFF * INTER_RESIDUE_CUTOFF,
        ("common_rna_dna", "common_rna_dna") => RNA_DNA_CUTOFF * RNA_DNA_CUTOFF,
        ("common_rna_dna", "metal") => METAL_CUTOFF * METAL_CUTOFF,
        ("common_rna_dna", "other") => OTHER_CUTOFF * OTHER_CUTOFF,
        ("ccp4_mon_lib_rna_dna", "other") => OTHER_CUTOFF * OTHER_CUTOFF,
        ("common_amino_acid", "common_saccharide") => CARBOHYDRATE_CUTOFF * CARBOHYDRATE_CUTOFF,
        ("common_saccharide", "common_saccharide") => CARBOHYDRATE_CUTOFF * CARBOHYDRATE_CUTOFF,
        ("common_element", "common_water") => METAL_CUTOFF * METAL_CUTOFF,
        ("other", "other") => OTHER_CUTOFF * OTHER_CUTOFF,
        ("sulfur", "sulfur") => SULFUR_CUTOFF * SULFUR_CUTOFF,
        ("common_amino_acid", "common_rna_dna") => AMINO_ACID_CUTOFF * AMINO_ACID_CUTOFF,
        ("common_rna_dna", "common_small_molecule") => OTHER_CUTOFF * OTHER_CUTOFF,
        ("common_amino_acid", "common_small_molecule") => AMINO_ACID_CUTOFF * OTHER_CUTOFF,
        ("common_small_molecule", "other") => OTHER_CUTOFF * OTHER_CUTOFF,
        ("metal", "metal") => METAL_CUTOFF * METAL_CUTOFF,
        ("metal", "other") | ("common_water", "metal") | ("common_amino_acid", "metal") => METAL_CUTOFF * METAL_CUTOFF,
        _ => return None,
    };
    Some(v)
}

fn maximum_inter_residue_links(a: &str, b: &str) -> usize {
    match (a, b) {
        ("common_element", "other") => 8,
        ("metal", "other") => 6,
        ("common_amino_acid", "metal") => 2,
        ("common_saccharide", "metal") => 3,
        ("common_rna_dna", "metal") => 2,
        ("common_rna_dna", "common_rna_dna") => 5,
        _ => 1,
    }
}

fn has_maximum_per_atom_links(c: &str) -> bool {
    matches!(c, "common_saccharide" | "common_rna_dna" | "common_amino_acid" | "other")
}

/// `linking_utils.get_classes` result.
#[derive(Clone, Debug)]
struct Classes {
    important: &'static str,
    flags: u16,
}

const CLASS_NAMES: [&str; 11] = [
    "common_saccharide",
    "common_water",
    "common_element",
    "common_small_molecule",
    "common_amino_acid",
    "common_rna_dna",
    "ccp4_mon_lib_rna_dna",
    "other",
    "uncommon_amino_acid",
    "unknown",
    "d_amino_acid",
];

fn class_bit(f: &str) -> u16 {
    CLASS_NAMES.iter().position(|x| *x == f).map(|k| 1u16 << k).unwrap_or(0)
}

impl Classes {
    fn has(&self, f: &str) -> bool {
        self.flags & class_bit(f) != 0
    }
}

#[derive(Clone, Debug)]
struct Ag {
    model: u32,
    chain: u32,
    rg: u32,
    atoms: Vec<u32>,
    /// All atoms of the residue group (every altloc).
    rg_atoms: Vec<u32>,
    resname: String,
    altloc: String,
}

pub struct AutoLinkParams {
    pub compat: bool,
}

struct Ctx<'a> {
    st: &'a Structure,
    flat: &'a FlatAtoms,
    ml: &'a MonLib,
    p: &'a AutoLinkParams,
    ag_of: Vec<u32>,
    ags: Vec<Ag>,
    /// 0 or 1 for atoms in the first or last residue group of a chain.
    first_last: Vec<i8>,
    type_cache: RefCell<FxHashMap<String, Option<String>>>,
    class_cache: RefCell<FxHashMap<u32, (&'static str, bool, bool)>>,
    /// interned residue-pair key part and raw atom name per atom (u32::MAX: not yet)
    key_id: RefCell<Vec<u32>>,
    name_id: RefCell<Vec<u32>>,
    interner: RefCell<FxHashMap<String, u32>>,
    interned: RefCell<Vec<String>>,
}

impl<'a> Ctx<'a> {
    fn new(st: &'a Structure, flat: &'a FlatAtoms, ml: &'a MonLib, p: &'a AutoLinkParams) -> Ctx<'a> {
        let n = flat.pos.len();
        let mut ag_of = vec![0u32; n];
        let mut ags: Vec<Ag> = Vec::new();
        let mut last: Option<(u32, u32, u32, u32)> = None;
        for (k, pth) in flat.path.iter().enumerate() {
            let key = (pth.model, pth.chain, pth.rg, pth.ag);
            if last != Some(key) {
                let g = st.atom_group(*pth);
                ags.push(Ag {
                    model: pth.model,
                    chain: pth.chain,
                    rg: pth.rg,
                    atoms: Vec::new(),
                    rg_atoms: Vec::new(),
                    resname: g.resname.clone(),
                    altloc: g.altloc.clone(),
                });
                last = Some(key);
            }
            let gi = ags.len() - 1;
            ags[gi].atoms.push(k as u32);
            ag_of[k] = gi as u32;
        }
        // residue-group atom lists
        let mut s = 0;
        while s < ags.len() {
            let mut e = s;
            while e < ags.len() && (ags[e].model, ags[e].chain, ags[e].rg) == (ags[s].model, ags[s].chain, ags[s].rg) {
                e += 1;
            }
            let all: Vec<u32> = ags[s..e].iter().flat_map(|g| g.atoms.iter().copied()).collect();
            for g in &mut ags[s..e] {
                g.rg_atoms = all.clone();
            }
            s = e;
        }
        let mut first_last = vec![-1i8; n];
        for (k, pth) in flat.path.iter().enumerate() {
            let c = st.chain(*pth);
            let nrg = c.residue_groups.len() as u32;
            if nrg == 0 {
                continue;
            }
            if pth.rg == nrg - 1 {
                first_last[k] = 1;
            } else if pth.rg == 0 {
                first_last[k] = 0;
            }
        }
        Ctx {
            st,
            flat,
            ml,
            p,
            ag_of,
            ags,
            first_last,
            type_cache: RefCell::new(FxHashMap::default()),
            class_cache: RefCell::new(FxHashMap::default()),
            key_id: RefCell::new(vec![u32::MAX; n]),
            name_id: RefCell::new(vec![u32::MAX; n]),
            interner: RefCell::new(FxHashMap::default()),
            interned: RefCell::new(Vec::new()),
        }
    }

    fn name(&self, a: u32) -> &str {
        &self.flat.name[a as usize]
    }
    fn elem(&self, a: u32) -> &str {
        &self.flat.element[a as usize]
    }
    fn ag(&self, a: u32) -> &Ag {
        &self.ags[self.ag_of[a as usize] as usize]
    }
    fn pos(&self, a: u32) -> Vec3 {
        self.flat.pos[a as usize]
    }
    fn d2(&self, a: u32, b: u32) -> f64 {
        let d = self.pos(a) - self.pos(b);
        d.x * d.x + d.y * d.y + d.z * d.z
    }

    /// `atom.id_str()`.
    fn id_str(&self, a: u32) -> String {
        let pth = self.flat.path[a as usize];
        let m = &self.st.models[pth.model as usize];
        let c = self.st.chain(pth);
        let rg = self.st.residue_group(pth);
        let g = self.st.atom_group(pth);
        let at = self.st.atom(pth);
        let mut s = String::with_capacity(48);
        if !m.id.is_empty() {
            s.push_str(&format!("model=\"{:>4}\" ", m.id));
        }
        s.push_str("pdb=\"");
        s.push_str(&format!("{:<4}", at.name));
        s.push_str(if g.altloc.is_empty() { " " } else { &g.altloc });
        s.push_str(&format!("{:>3}", g.resname));
        s.push_str(&format!("{:>2}", c.id));
        s.push_str(&format!("{:>4}", rg.resseq));
        s.push_str(if rg.icode.is_empty() { " " } else { &rg.icode });
        s.push('"');
        if !at.segid.trim().is_empty() {
            s.push_str(&format!(" segid=\"{:<4}\"", at.segid));
        }
        s
    }

    /// `atom_group.id_str()` (altloc, resname, chain, resseq, icode).
    fn ag_id_str(&self, a: u32) -> String {
        let pth = self.flat.path[a as usize];
        let c = self.st.chain(pth);
        let rg = self.st.residue_group(pth);
        let g = self.st.atom_group(pth);
        format!(
            "{}{:>3}{:>2}{:>4}{}",
            if g.altloc.is_empty() { " " } else { &g.altloc },
            g.resname,
            c.id,
            rg.resseq,
            if rg.icode.is_empty() { " " } else { &rg.icode }
        )
    }
    /// `residue_group.id_str()` (chain, resseq, icode).
    fn rg_id_str(&self, a: u32) -> String {
        let pth = self.flat.path[a as usize];
        let c = self.st.chain(pth);
        let rg = self.st.residue_group(pth);
        format!("{:>2}{:>4}{}", c.id, rg.resseq, if rg.icode.is_empty() { " " } else { &rg.icode })
    }

    /// Residue-pair key part of an atom.
    fn res_key(&self, a: u32) -> String {
        if self.p.compat {
            let s = self.id_str(a);
            s[9..s.len() - 1].to_string()
        } else {
            let pth = self.flat.path[a as usize];
            format!("{}:{}", pth.model, self.ag_of[a as usize])
        }
    }

    fn intern(&self, s: String) -> u32 {
        if let Some(&i) = self.interner.borrow().get(&s) {
            return i;
        }
        let mut v = self.interned.borrow_mut();
        let i = v.len() as u32;
        v.push(s.clone());
        self.interner.borrow_mut().insert(s, i);
        i
    }
    /// Interned `res_key` of an atom.
    fn res_key_id(&self, a: u32) -> u32 {
        let c = self.key_id.borrow()[a as usize];
        if c != u32::MAX {
            return c;
        }
        let i = self.intern(self.res_key(a));
        self.key_id.borrow_mut()[a as usize] = i;
        i
    }
    /// Interned raw atom name.
    fn name_id(&self, a: u32) -> u32 {
        let c = self.name_id.borrow()[a as usize];
        if c != u32::MAX {
            return c;
        }
        let i = self.intern(self.name(a).to_string());
        self.name_id.borrow_mut()[a as usize] = i;
        i
    }
    /// Order two interned strings as the strings compare.
    fn ordered(&self, a: u32, b: u32) -> (u32, u32) {
        let v = self.interned.borrow();
        if v[a as usize] <= v[b as usize] { (a, b) } else { (b, a) }
    }

    /// `mmtbx.chemical_components.get_type`, upper-cased.
    fn ccd_type(&self, resname: &str) -> Option<String> {
        let key = resname.trim().to_string();
        if key.is_empty() {
            return None;
        }
        if let Some(v) = self.type_cache.borrow().get(&key) {
            return v.clone();
        }
        let v = self.ml.ccd(&key).map(|e| e.type_.trim().to_ascii_uppercase());
        self.type_cache.borrow_mut().insert(key, v.clone());
        v
    }
    fn is_sugar_type(&self, resname: &str) -> bool {
        self.ccd_type(resname).map_or(false, |t| SUGAR_TYPES.contains(&t.as_str()))
    }

    /// `linking_utils.get_classes(atom)`.
    fn classes(&self, a: u32) -> Classes {
        let gi = self.ag_of[a as usize];
        let cached = self.class_cache.borrow().get(&gi).copied();
        let (gc, sacch, uncommon) = match cached {
            Some(v) => v,
            None => {
                let g = &self.ags[gi as usize];
                let consider = g.atoms.len() != 1;
                let mut gc = resclass::get_class_ext(&g.resname, consider).name();
                if g.resname == "UNK" {
                    gc = "common_amino_acid";
                }
                if gc == "modified_amino_acid" || gc == "modified_rna_dna" {
                    gc = "other";
                }
                let sacch = if ONE_LETTER_GIVEN_THREE_LETTER.iter().any(|(n, _)| *n == g.resname) || g.resname == "HOH" {
                    false
                } else if gc == "common_saccharide" {
                    true
                } else {
                    self.is_sugar_type(&g.resname)
                };
                let bb = g
                    .rg_atoms
                    .iter()
                    .filter(|&&x| matches!(self.name(x).trim(), "C" | "CA" | "N" | "O" | "OXT"))
                    .count();
                let v = (gc, sacch, bb >= 4);
                self.class_cache.borrow_mut().insert(gi, v);
                v
            }
        };
        let mut flags: u16 = 0;
        if sacch {
            flags |= class_bit("common_saccharide");
        }
        if gc != "common_saccharide" {
            flags |= class_bit(gc);
        }
        let important = if sacch {
            "common_saccharide"
        } else {
            let e = self.elem(a);
            if gc == "common_element" && METALS.contains(&e) {
                "metal"
            } else if gc == "other" && self.ag(a).rg_atoms.len() == 1 && METALS.contains(&e) {
                "metal"
            } else {
                gc
            }
        };
        if flags & class_bit("other") != 0 && uncommon {
            flags |= class_bit("uncommon_amino_acid");
        }
        Classes { important, flags }
    }

    /// `linking_utils.get_bonded(hierarchy, atom, bond_cutoff)`: atoms of the
    /// first atom group in the hierarchy with the atom's model, chain id,
    /// resseq and resname that lie within the cutoff (other names only).
    fn get_bonded(&self, a: u32, cutoff: f64) -> Vec<u32> {
        let pth = self.flat.path[a as usize];
        let model_id = &self.st.models[pth.model as usize].id;
        let chain_id = &self.st.chain(pth).id;
        let resseq = &self.st.residue_group(pth).resseq;
        let resname = &self.st.atom_group(pth).resname;
        let c2 = cutoff * cutoff;
        for g in &self.ags {
            let m = &self.st.models[g.model as usize];
            if &m.id != model_id {
                continue;
            }
            let c = &m.chains[g.chain as usize];
            if &c.id != chain_id {
                continue;
            }
            let rg = &c.residue_groups[g.rg as usize];
            if &rg.resseq != resseq || &g.resname != resname {
                continue;
            }
            let mut out = Vec::new();
            for &x in &g.atoms {
                if self.name(x) == self.name(a) {
                    continue;
                }
                if self.d2(a, x) <= c2 {
                    out.push(x);
                }
            }
            return out;
        }
        Vec::new()
    }

    fn check_valence(&self, a: u32) -> bool {
        if self.elem(a) != "O" {
            return true;
        }
        self.get_bonded(a, 1.8).len() != 2
    }

    /// `linking_utils.is_atom_pair_linked`.
    fn is_atom_pair_linked(&self, a1: u32, a2: u32, imp1: &'static str, imp2: &'static str, distance: f64, only_cutoff: bool) -> bool {
        let (e1, e2) = (self.elem(a1), self.elem(a2));
        if NON_LINKING.contains(&e1) || NON_LINKING.contains(&e2) {
            return false;
        }
        if e1 == "O" && e2 == "O" {
            return false;
        }
        let adjust = |e: &str, c: &'static str| -> &'static str {
            if c == "common_element" {
                if METALS.contains(&e) {
                    return "metal";
                } else if NON_LINKING.contains(&e) {
                    return "ion";
                }
            } else if c == "other" && METALS.contains(&e) {
                return "metal";
            }
            c
        };
        let mut c1 = adjust(e1, imp1);
        let mut c2 = adjust(e2, imp2);
        let sulfur = |e: &str, c: &str| (c == "common_amino_acid" || c == "other") && e == "S";
        if sulfur(e1, c1) && sulfur(e2, c2) {
            c1 = "sulfur";
            c2 = "sulfur";
        }
        let (l0, l1) = if c1 <= c2 { (c1, c2) } else { (c2, c1) };
        if matches!(
            (l0, l1),
            ("common_water", "common_water")
                | ("common_amino_acid", "common_water")
                | ("common_saccharide", "common_water")
                | ("common_water", "other")
        ) {
            return false;
        }
        let mut limit = skip_if_longer(l0, l1);
        if let Some(l) = limit.as_mut() {
            // elements are stripped by the time cctbx links (" C" -> "C")
            if !FIRST_ROW.contains(&e1) || !FIRST_ROW.contains(&e2) {
                *l += SECOND_ROW_BUFFER * SECOND_ROW_BUFFER;
            }
        }
        let d2 = distance * distance;
        if let Some(l) = limit {
            if l < d2 {
                return false;
            }
        }
        if only_cutoff {
            return true;
        }
        let has = |x: &str| l0 == x || l1 == x;
        if has("common_rna_dna") && has("common_amino_acid") {
            return true;
        }
        if has("common_rna_dna") && has("common_small_molecule") {
            return true;
        }
        if has("common_amino_acid") && has("common_small_molecule") {
            return true;
        }
        if has("other") && has("common_small_molecule") {
            return true;
        }
        if c1 == "sulfur" && c2 == "sulfur" {
            return true;
        }
        if has("common_saccharide") {
            let lim = if has("metal") { METAL_CUTOFF * METAL_CUTOFF } else { CARBOHYDRATE_CUTOFF * CARBOHYDRATE_CUTOFF };
            return d2 <= lim;
        }
        if has("common_element") {
            // a non-metal common element raises Sorry in cctbx
            return METALS.contains(&e1) || METALS.contains(&e2);
        }
        if has("metal") {
            for &x in &[a1, a2] {
                let rn = self.ag(x).resname.trim();
                let nm = self.name(x).trim();
                if rn == "HIS" && matches!(nm, "CE1" | "CD2" | "CB") {
                    return false;
                }
            }
            // link_metals = Auto (truthy)
            return true;
        }
        if c1 == "common_amino_acid" && c2 == "common_amino_acid" {
            if (e1 == "N" && e2 == "C") || (e1 == "C" && e2 == "N") {
                return true;
            }
        }
        if c1 == "d_amino_acid" || c2 == "d_amino_acid" {
            let d = if c1 == "d_amino_acid" { a1 } else { a2 };
            if self.name(d).trim() == "O" {
                return false;
            }
        }
        has("other")
    }

    fn possible_cyclic_peptide(&self, a1: u32, a2: u32) -> bool {
        let mut names = [self.name(a1), self.name(a2)];
        names.sort();
        if names != [" C  ", " N  "] {
            return false;
        }
        let fl = self.first_last[a1 as usize] as i32 + self.first_last[a2 as usize] as i32;
        if fl == 1 {
            let c1 = &self.st.chain(self.flat.path[a1 as usize]).id;
            let c2 = &self.st.chain(self.flat.path[a2 as usize]).id;
            return c1 == c2;
        }
        false
    }

    fn find_by_name(&self, g: u32, name: &str) -> Option<u32> {
        let n = name.trim();
        self.ags[g as usize].atoms.iter().copied().find(|&x| self.name(x).trim() == n)
    }

    /// `_apply_link_using_proxies`; returns the bonds made.
    fn apply_link(&self, it: &mut Interp, link: &ChemLink, g1: u32, g2: u32, origin: u16) -> Vec<(u32, u32)> {
        let (mut g1, mut g2) = (g1, g2);
        let side = |c: u8, g1: u32, g2: u32| if c == 1 { g1 } else { g2 };
        let mut bonds = Vec::new();
        for b in &link.bonds {
            let get = |g1: u32, g2: u32| -> Option<(u32, u32)> {
                Some((self.find_by_name(side(b.c1, g1, g2), &b.a1)?, self.find_by_name(side(b.c2, g1, g2), &b.a2)?))
            };
            let Some(mut ij) = get(g1, g2) else { continue };
            if self.d2(ij.0, ij.1) > 9.0 {
                std::mem::swap(&mut g1, &mut g2);
                match get(g1, g2) {
                    Some(x) => ij = x,
                    None => continue,
                }
            }
            let Some(v) = b.value_dist else { continue };
            it.upsert_bond(ij.0, ij.1, v, origin);
            bonds.push(ij);
        }
        for a in &link.angles {
            let ids: Option<Vec<u32>> = (0..3).map(|k| self.find_by_name(side(a.c[k], g1, g2), &a.a[k])).collect();
            let Some(ids) = ids else { continue };
            let (Some(v), Some(e)) = (a.value, a.esd) else { continue };
            it.add_angle(ids[0], ids[1], ids[2], v, e);
        }
        for t in &link.tors {
            let ids: Option<Vec<u32>> = (0..4).map(|k| self.find_by_name(side(t.c[k], g1, g2), &t.a[k])).collect();
            let Some(ids) = ids else { continue };
            let Some(v) = t.value else { continue };
            let ids = [ids[0], ids[1], ids[2], ids[3]];
            match t.esd {
                Some(e) if e > 0.0 => it.add_dihedral(ids, v, t.period, false),
                _ => it.add_dihedral(ids, v, 0, true),
            }
        }
        let mut planes: Vec<(String, Vec<u32>)> = Vec::new();
        for pl in &link.planes {
            let Some(x) = self.find_by_name(side(pl.c, g1, g2), &pl.atom) else { continue };
            match planes.iter_mut().find(|q| q.0 == pl.plane_id) {
                Some(q) => q.1.push(x),
                None => planes.push((pl.plane_id.clone(), vec![x])),
            }
        }
        for (_, atoms) in planes {
            it.add_plane(atoms);
        }
        bonds
    }

    /// `glyco_utils.get_glyco_link_atoms(atom_group1, atom_group2)`.
    fn glyco_link_atoms(&self, g1: u32, g2: u32, cutoff: f64) -> Option<GlycoAtoms> {
        // bonds by distance over the atoms of both residue groups
        let all: Vec<u32> =
            self.ags[g1 as usize].rg_atoms.iter().chain(self.ags[g2 as usize].rg_atoms.iter()).copied().collect();
        let c2 = cutoff * cutoff;
        let mut bonds: FxHashMap<u32, Vec<u32>> = FxHashMap::default();
        for (i, &x) in all.iter().enumerate() {
            for &y in &all[i + 1..] {
                let (ax, ay) = (&self.ag(x).altloc, &self.ag(y).altloc);
                if !(ax == ay || ax.is_empty() || ay.is_empty()) {
                    continue;
                }
                if self.d2(x, y) <= c2 {
                    bonds.entry(x).or_default().push(y);
                    bonds.entry(y).or_default().push(x);
                }
            }
        }
        let nb = |x: u32| bonds.get(&x).map(|v| v.as_slice()).unwrap_or(&[]);
        let el = |x: u32| self.elem(x);
        let order: Vec<u32> =
            self.ags[g1 as usize].atoms.iter().chain(self.ags[g2 as usize].atoms.iter()).copied().collect();
        // anomeric carbon: a C bonded to two O of different atom groups
        let mut found: Option<(u32, bool)> = None;
        for &x in &order {
            if el(x) != "C" {
                continue;
            }
            let ox: Vec<u32> = nb(x).iter().copied().filter(|&o| el(o) == "O").collect();
            if ox.len() == 2 && self.ag_id_str(ox[0]) != self.ag_id_str(ox[1]) {
                found = Some((x, true));
                break;
            }
        }
        if found.is_none() {
            // get_any_linking_carbon
            for &x in &order {
                if el(x) != "C" {
                    continue;
                }
                let mut nox = 0;
                let mut linking = false;
                for &o in nb(x) {
                    if el(o) == "O" {
                        if self.ag_id_str(x) != self.ag_id_str(o) {
                            linking = true;
                        }
                        nox += 1;
                    }
                }
                if nox == 2 {
                    return None; // Sorry in cctbx
                }
                if linking {
                    found = Some((x, true));
                    break;
                }
            }
        }
        if found.is_none() {
            // get_C1_carbon
            let c1s: Vec<u32> = order.iter().copied().filter(|&x| self.name(x).trim().starts_with("C1")).collect();
            let mut res = None;
            for &c1 in &c1s {
                let ox: Vec<u32> = order.iter().copied().filter(|&o| el(o) == "O" && self.d2(c1, o) < c2).collect();
                if ox.len() == 2 {
                    res = Some((c1, self.ag_id_str(ox[0]) != self.ag_id_str(ox[1])));
                    break;
                }
            }
            found = Some(res?);
        }
        let (ac, linking) = found?;
        let same = |x: u32| self.ag_id_str(x) == self.ag_id_str(ac) || self.rg_id_str(x) == self.rg_id_str(ac);
        let ring_with = |e: &str| nb(ac).iter().copied().find(|&x| el(x) == e && same(x));
        let ring_oxygen = ring_with("O").or_else(|| ring_with("C"));
        let ring_carbon = ring_with("C");
        let anomeric_h = ring_with("H");
        let other_group = |x: u32| self.ag_id_str(x) != self.ag_id_str(ac) || self.rg_id_str(x) != self.rg_id_str(ac);
        let mut link_oxygen = nb(ac).iter().copied().find(|&x| el(x) == "O" && other_group(x));
        let on_distance = || -> Option<u32> {
            let ga = self.ag_of[ac as usize];
            let link_group = if self.ags[g2 as usize].atoms.contains(&ac) || ga == g2 {
                g1
            } else {
                g2
            };
            self.ags[link_group as usize].atoms.iter().copied().find(|&x| el(x) == "O" && self.d2(x, ac) < 5.0)
        };
        if link_oxygen.is_none() {
            link_oxygen = on_distance();
        }
        let lo = link_oxygen?;
        let mut link_carbon =
            nb(lo).iter().copied().find(|&x| el(x) == "C" && x != ac && self.ag_id_str(x) != self.ag_id_str(ac));
        if link_carbon.is_none() {
            // get_link_carbon_on_distance returns an oxygen, as in cctbx
            link_carbon = on_distance();
        }
        Some(GlycoAtoms { anomeric_carbon: ac, ring_oxygen, ring_carbon, link_oxygen: lo, link_carbon, anomeric_h, linking })
    }

    /// `glyco_utils.apply_glyco_link_using_proxies_and_atoms`.
    fn apply_glyco_link(&self, it: &mut Interp, g1: u32, g2: u32) -> Option<(u32, u32)> {
        let mut gla = self.glyco_link_atoms(g1, g2, CARBOHYDRATE_CUTOFF);
        if let Some(g) = &gla {
            if !g.is_correct() {
                gla = self.glyco_link_atoms(g2, g1, CARBOHYDRATE_CUTOFF);
            }
        }
        let g = gla?;
        if !g.is_correct() || !g.linking {
            return None; // Sorry in cctbx
        }
        let isomer = self.isomer(&g);
        let origin = link_origin(&isomer).unwrap_or(ORIGIN_GLYCO_CUSTOM);
        it.upsert_bond(g.anomeric_carbon, g.link_oxygen, 1.439, origin);
        let lc = g.link_carbon;
        for (atoms, v) in [
            ([lc, Some(g.link_oxygen), Some(g.anomeric_carbon)], 108.7),
            ([Some(g.link_oxygen), Some(g.anomeric_carbon), g.ring_oxygen], 112.3),
            ([Some(g.link_oxygen), Some(g.anomeric_carbon), g.ring_carbon], 109.47),
            ([Some(g.link_oxygen), Some(g.anomeric_carbon), g.anomeric_h], 109.47),
        ] {
            if let [Some(a), Some(b), Some(c)] = atoms {
                it.add_angle(a, b, c, v, 3.0);
            }
        }
        Some((g.anomeric_carbon, g.link_oxygen))
    }

    /// `glyco_link_class.get_isomer`.
    fn isomer(&self, g: &GlycoAtoms) -> String {
        let rn = self.ag(g.anomeric_carbon).resname.clone();
        let mut s = match linkdata::GLYCO_VOLUMES.iter().find(|(k, _)| *k == rn) {
            Some(&(_, v)) if v < 0.0 => "ALPHA".to_string(),
            Some(_) => "BETA".to_string(),
            None => "?".to_string(),
        };
        let cn = self.name(g.anomeric_carbon).trim();
        match cn.chars().last() {
            Some(ch) if ch.is_ascii_digit() => s.push(ch),
            _ => s.push_str(&format!(" {} ", cn)),
        }
        let on = self.name(g.link_oxygen).trim();
        match on.chars().last() {
            Some(ch) if ch.is_ascii_digit() => {
                s.push('-');
                s.push(ch);
            }
            _ => s.push_str(&format!("- {} ", on)),
        }
        s
    }

    /// `process_atom_groups_for_linking_single_link`: the link key, or None.
    fn single_link_key(&self, a1: u32, a2: u32) -> Option<String> {
        let (mut a1, mut a2) = (a1, a2);
        let t1 = self.ccd_type(&self.ag(a1).resname);
        let t2 = self.ccd_type(&self.ag(a2).resname);
        let sugar = |t: &Option<String>| t.as_deref().map_or(false, |t| SUGAR_TYPES.contains(&t));
        let amino = |t: &Option<String>| t.as_deref().map_or(false, |t| AMINO_TYPES.contains(&t));
        let both = t1.is_some() && t2.is_some();
        let glyco = both && sugar(&t1) && sugar(&t2);
        let glyco_amino = both && {
            let s = sugar(&t1) as u8 + sugar(&t2) as u8;
            let am = (!sugar(&t1) && amino(&t1)) as u8 + (!sugar(&t2) && amino(&t2)) as u8;
            s == 1 && am == 1
        };
        if glyco {
            if self.name(a1).contains('C') {
                std::mem::swap(&mut a1, &mut a2);
            }
        } else if glyco_amino && self.name(a2).contains('C') {
            std::mem::swap(&mut a1, &mut a2);
        }
        let (r1, r2) = (self.ag(a1).resname.trim().to_string(), self.ag(a2).resname.trim().to_string());
        let long_key = format!("{}:{}-{}:{}", r1, self.name(a1).trim(), r2, self.name(a2).trim());
        let tmp_key = format!("{}-{}", r1, r2);
        let (rn1, rn2) = (&self.ag(a1).resname, &self.ag(a2).resname);
        let (t1, t2) = (self.ccd_type(rn1), self.ccd_type(rn2));
        let both = t1.is_some() && t2.is_some();
        let count = |names: &[&str]| -> bool {
            let s = sugar(&t1) as u8 + sugar(&t2) as u8;
            let l = (!sugar(&t1) && names.contains(&rn1.as_str())) as u8 + (!sugar(&t2) && names.contains(&rn2.as_str())) as u8;
            s == 1 && l == 1
        };
        if both && count(&["ASN"]) {
            return Some(tmp_key);
        }
        if both && count(&["SER", "THR"]) {
            return Some(tmp_key);
        }
        if both && sugar(&t1) && sugar(&t2) {
            let c_atom = if self.name(a1).contains('C') {
                Some(a1)
            } else if self.name(a2).contains('C') {
                Some(a2)
            } else {
                None
            };
            let o_atom = if self.name(a1).contains('O') {
                Some(a1)
            } else if self.name(a2).contains('O') {
                Some(a2)
            } else {
                None
            };
            if let (Some(c), Some(o)) = (c_atom, o_atom) {
                // the hand (ALPHA/BETA) does not change what is applied
                let cd = self.name(c).trim().chars().last().unwrap_or(' ');
                let od = self.name(o).trim().chars().last().unwrap_or(' ');
                return Some(format!("{}{}-{}", self.glyco_hand(c, o), cd, od));
            } else if self.d2(a1, a2) <= OTHER_CUTOFF * OTHER_CUTOFF {
                return Some(long_key);
            }
            return None;
        }
        Some(long_key)
    }

    /// `get_hand` via `get_angles_from_included_bonds` and `get_chiral_volume`.
    fn glyco_hand(&self, c: u32, o: u32) -> &'static str {
        let (a1, a2) = if self.name(c).contains('C') && c != o { (o, c) } else { (c, o) };
        let mut angles: Vec<[u32; 3]> = Vec::new();
        for (i, &x) in [a1, a2].iter().enumerate() {
            let other = if i == 1 { a1 } else { a2 };
            for r in self.get_bonded(x, 1.75) {
                angles.push([other, x, r]);
            }
        }
        let c_resseq = &self.st.residue_group(self.flat.path[c as usize]).resseq;
        let mut others: Vec<u32> = Vec::new();
        for ang in &angles {
            for &x in ang {
                if self.flat.is_h[x as usize] && matches!(self.elem(x), "H" | "D" | "T") {
                    continue;
                }
                if matches!(self.elem(x), "H" | "D" | "T") {
                    continue;
                }
                if &self.st.residue_group(self.flat.path[x as usize]).resseq != c_resseq {
                    continue;
                }
                if self.name(x) == self.name(c) || self.name(x) == self.name(o) {
                    continue;
                }
                others.push(x);
            }
        }
        others.sort_by(|&x, &y| self.name(x).cmp(self.name(y)));
        if others.len() != 2 {
            return "ALPHA";
        }
        let pc = self.pos(c);
        let v = (self.pos(o) - pc).dot((self.pos(others[0]) - pc).cross(self.pos(others[1]) - pc));
        if v < 0.0 {
            "BETA"
        } else {
            "ALPHA"
        }
    }

    /// `check_for_peptide_links`: Some((key, swap)), or None/false.
    fn peptide_link(&self, a1: u32, a2: u32, cl1: &Classes, cl2: &Classes) -> Option<(&'static str, bool)> {
        if !(cl1.has("common_amino_acid") || cl2.has("common_amino_acid")) {
            return None;
        }
        let other = if cl1.has("common_amino_acid") { self.ag_of[a2 as usize] } else { self.ag_of[a1 as usize] };
        let og = &self.ags[other as usize];
        if resclass::get_class(&og.resname).name() != "other" {
            return None;
        }
        let (n1, n2) = (self.name(a1).trim(), self.name(a2).trim());
        if n1 == "SG" && n2 == "SG" {
            return Some(("SS", false));
        }
        let count = og.atoms.iter().filter(|&&x| matches!(self.name(x).trim(), "C" | "N" | "O")).count();
        if count != 3 {
            return None;
        }
        if n2 == "C" && n1 == "N" {
            return Some(("TRANS", false));
        }
        if n1 == "C" && n2 == "N" {
            return Some(("TRANS", true));
        }
        None
    }
}

struct GlycoAtoms {
    anomeric_carbon: u32,
    ring_oxygen: Option<u32>,
    ring_carbon: Option<u32>,
    link_oxygen: u32,
    link_carbon: Option<u32>,
    anomeric_h: Option<u32>,
    linking: bool,
}

impl GlycoAtoms {
    fn is_correct(&self) -> bool {
        self.ring_oxygen.is_some() && self.link_carbon.is_some()
    }
}

/// Heavy-atom pairs closer than `MAX_BONDED_CUTOFF` within each model, in
/// order of distance (then atom indices). Pairs the main loop would skip
/// without side effects whatever happened before them are left out: atoms
/// bonded before linking starts, atoms of one residue group, and different
/// non-blank altlocs.
fn candidate_pairs(cx: &Ctx, initial: &[Vec<u32>]) -> Vec<(f64, u32, u32)> {
    let flat = cx.flat;
    let n = flat.pos.len();
    let cell = MAX_BONDED_CUTOFF;
    let ckey = |p: Vec3| ((p.x / cell).floor() as i64, (p.y / cell).floor() as i64, (p.z / cell).floor() as i64);
    let pack = |m: u32, x: i64, y: i64, z: i64| -> u64 {
        // 16 bits model, 16 bits per coordinate (offset), plenty for any model
        ((m as u64 & 0xffff) << 48) | (((x + 32768) as u64 & 0xffff) << 32) | (((y + 32768) as u64 & 0xffff) << 16) | ((z + 32768) as u64 & 0xffff)
    };
    let mut cells: Vec<(u64, u32)> = (0..n)
        .filter(|&a| !flat.is_h[a])
        .map(|a| {
            let (x, y, z) = ckey(flat.pos[a]);
            (pack(flat.path[a].model, x, y, z), a as u32)
        })
        .collect();
    crate::par::sort_unstable_by_key(&mut cells, |&c| c);
    let mut ranges: FxHashMap<u64, (usize, usize)> = FxHashMap::default();
    let mut k = 0;
    while k < cells.len() {
        let mut e = k;
        while e < cells.len() && cells[e].0 == cells[k].0 {
            e += 1;
        }
        ranges.insert(cells[k].0, (k, e));
        k = e;
    }
    let rg_of = |a: u32| {
        let p = flat.path[a as usize];
        (p.model, p.chain, p.rg)
    };
    let mut out: Vec<(f64, u32, u32)> = crate::par::flat_map_collect(&cells, |&(_, a)| {
            let pa = flat.pos[a as usize];
            let (x, y, z) = ckey(pa);
            let m = flat.path[a as usize].model;
            let ra = rg_of(a);
            let alt_a = flat.altloc[a as usize].trim();
            let mut v: Vec<(f64, u32, u32)> = Vec::new();
            for dx in -1..=1 {
                for dy in -1..=1 {
                    for dz in -1..=1 {
                        let Some(&(s, e)) = ranges.get(&pack(m, x + dx, y + dy, z + dz)) else { continue };
                        for &(_, b) in &cells[s..e] {
                            if b <= a {
                                continue;
                            }
                            let d = (flat.pos[b as usize] - pa).length();
                            if d >= MAX_BONDED_CUTOFF || rg_of(b) == ra {
                                continue;
                            }
                            let alt_b = flat.altloc[b as usize].trim();
                            if !(alt_a == alt_b || alt_a.is_empty() || alt_b.is_empty()) {
                                continue;
                            }
                            if initial[a as usize].contains(&b) {
                                continue;
                            }
                            v.push((d, a, b));
                        }
                    }
                }
            }
            v.into_iter()
        });
    crate::par::sort_unstable_by_key(&mut out, |&(d, i, j)| (d.to_bits(), i, j));
    out
}

/// Make the automatic links. Bonds (with their angles, dihedrals and planes)
/// go into `it`; a line per link goes to `log`.
pub fn auto_link(it: &mut Interp, st: &Structure, flat: &FlatAtoms, ml: &MonLib, p: &AutoLinkParams, log: &mut String) {
    let cx = Ctx::new(st, flat, ml, p);
    let n = flat.pos.len();
    let mut initial: Vec<Vec<u32>> = vec![Vec::new(); n];
    for b in &it.bonds {
        initial[b.i as usize].push(b.j);
        initial[b.j as usize].push(b.i);
    }
    // `done`: atoms that took their link, and per residue pair the atom-name
    // pairs considered
    let mut atom_done: Vec<bool> = vec![false; n];
    let mut res_done: FxHashMap<(u32, u32), smallvec::SmallVec<[(u32, u32); 2]>> = FxHashMap::default();
    let mut n_links: Vec<u32> = vec![0; n];
    let mut custom: Vec<(u32, u32, u16)> = Vec::new();
    let mut custom_partners: FxHashMap<u32, Vec<u32>> = FxHashMap::default();

    for (dist, i, j) in candidate_pairs(&cx, &initial) {
        if it.has_bond(i, j) {
            continue;
        }
        let (gi, gj) = (cx.ag(i), cx.ag(j));
        let moved = ["SF4", "F3S", "FES"];
        if moved.contains(&gi.resname.as_str()) || moved.contains(&gj.resname.as_str()) {
            continue;
        }
        let (ri, rj) = (gi.resname.trim(), gj.resname.trim());
        if matches!(ri, "ZN" | "CYS") && matches!(rj, "ZN" | "CYS") {
            continue;
        }
        if matches!(ri, "ZN" | "HIS") && matches!(rj, "ZN" | "HIS") {
            continue;
        }
        if (gi.model, gi.chain, gi.rg) == (gj.model, gj.chain, gj.rg) {
            continue;
        }
        let (alt_i, alt_j) = (gi.altloc.trim(), gj.altloc.trim());
        if !(alt_i == alt_j || alt_i.is_empty() || alt_j.is_empty()) {
            continue;
        }
        let cl1 = cx.classes(i);
        let cl2 = cx.classes(j);
        // link_small_molecules = False
        if cl1.has("common_small_molecule") || cl2.has("common_small_molecule") {
            continue;
        }
        let peptides = ["common_amino_acid", "d_amino_acid", "uncommon_amino_acid"];
        let mut only_cutoff = false;
        if peptides.contains(&cl1.important) && peptides.contains(&cl2.important) && cx.possible_cyclic_peptide(i, j) {
            only_cutoff = true;
        }
        // bonded atoms can't link to the same atom
        if !custom.is_empty() {
            let mut bonded = false;
            if let Some(v) = custom_partners.get(&i) {
                bonded |= v.iter().any(|&t| it.has_bond(j, t));
            }
            if let Some(v) = custom_partners.get(&j) {
                bonded |= v.iter().any(|&t| it.has_bond(i, t));
            }
            if bonded {
                continue;
            }
        }
        let key = cx.ordered(cx.res_key_id(i), cx.res_key_id(j));
        if !cx.is_atom_pair_linked(i, j, cl1.important, cl2.important, dist, only_cutoff) {
            // recorded unsorted, as the original does
            res_done.entry(key).or_default().push((cx.name_id(i), cx.name_id(j)));
            continue;
        }
        if !(cl1.has("common_element") || cl2.has("common_element")) {
            if !cx.check_valence(i) || !cx.check_valence(j) {
                log.push_str(&format!("  Atom rejected from bonding due to valence issues: {}\n", cx.id_str(if cx.check_valence(i) { j } else { i })));
                continue;
            }
        }
        let (c1, c2) = (cl1.important, cl2.important);
        let class_key = if c1 <= c2 { (c1, c2) } else { (c2, c1) };
        // link_metals = Auto
        if class_key.0 == "metal" || class_key.1 == "metal" {
            continue;
        }
        if class_key.0 == "common_amino_acid" && class_key.1 == "common_rna_dna" {
            continue;
        }
        let names = cx.ordered(cx.name_id(i), cx.name_id(j));
        if res_done.get(&key).map_or(false, |v| v.contains(&names)) {
            continue;
        }
        if gi.altloc == gj.altloc {
            let may_link_again = |a: u32, n_links: &Vec<u32>| -> bool {
                if dist > OTHER_CUTOFF {
                    return false;
                }
                if cx.elem(a) != "N" {
                    return false;
                }
                let heavy = initial[a as usize].iter().filter(|&&x| !flat.is_h[x as usize]).count();
                heavy + (n_links[a as usize] as usize) < 3
            };
            if has_maximum_per_atom_links(c1) {
                if atom_done[i as usize] && !may_link_again(i, &n_links) {
                    continue;
                }
                atom_done[i as usize] = true;
                n_links[i as usize] += 1;
            }
            if has_maximum_per_atom_links(c2) {
                if atom_done[j as usize] && !may_link_again(j, &n_links) {
                    continue;
                }
                atom_done[j as usize] = true;
                n_links[j as usize] += 1;
            }
        }
        let entry = res_done.entry(key).or_default();
        if entry.len() >= maximum_inter_residue_links(class_key.0, class_key.1) {
            continue;
        }
        entry.push(names);
        // a bond between the two atom groups already (before auto-linking)?
        let (agi, agj) = (cx.ag_of[i as usize], cx.ag_of[j as usize]);
        let link_found = cx.ags[agi as usize]
            .atoms
            .iter()
            .any(|&x| initial[x as usize].iter().any(|&y| cx.ag_of[y as usize] == agj));
        if link_found {
            continue;
        }
        let (mut g1, mut g2) = (agi, agj);
        // predefined residue-pair link, e.g. NAG-ASN
        let mut found: Option<(&ChemLink, String)> = None;
        let k1 = format!("{}-{}", gi.resname, gj.resname);
        let k2 = format!("{}-{}", gj.resname, gi.resname);
        if let Some(l) = ml.link(&k1) {
            found = Some((l, k1));
        } else if let Some(l) = ml.link(&k2) {
            found = Some((l, k2));
            std::mem::swap(&mut g1, &mut g2);
        } else if peptides.contains(&c1) && peptides.contains(&c2) && ((cx.name(i) == " N  " && cx.name(j) == " C  ") || (cx.name(i) == " C  " && cx.name(j) == " N  ")) {
            if let Some(l) = ml.link("TRANS") {
                if cx.name(i) == " C  " {
                    std::mem::swap(&mut g1, &mut g2);
                }
                found = Some((l, "TRANS".to_string()));
            }
        } else {
            let k1 = format!("{}_{}-ANY_{}", gi.resname, cx.name(i).trim(), cx.name(j).trim());
            let k2 = format!("{}_{}-ANY_{}", gj.resname, cx.name(j).trim(), cx.name(i).trim());
            if let Some(l) = ml.link(&k1) {
                found = Some((l, k1));
            } else if let Some(l) = ml.link(&k2) {
                found = Some((l, k2));
                std::mem::swap(&mut g1, &mut g2);
            }
        }
        if let Some((link, lkey)) = found {
            let Some(origin) = link_origin(&lkey) else { continue };
            let made = cx.apply_link(it, link, g1, g2, origin);
            if !made.is_empty() {
                log.push_str(&format!("  Link {}: {} - {}\n", lkey, cx.id_str(made[0].0), cx.id_str(made[0].1)));
            }
            continue;
        }
        let Some(lkey) = cx.single_link_key(i, j) else {
            if let Some(v) = res_done.get_mut(&key) {
                if let Some(pos) = v.iter().position(|x| *x == names) {
                    v.remove(pos);
                }
            }
            continue;
        };
        let is_glyco = ["ALPHA1", "BETA1", "ALPHA2", "BETA2"].iter().any(|s| lkey.contains(s));
        if is_glyco {
            if let Some((a, b)) = cx.apply_glyco_link(it, agj, agi) {
                log.push_str(&format!("  Glycosidic link: {} - {}\n", cx.id_str(a), cx.id_str(b)));
            }
            continue;
        }
        if let Some(link) = ml.link(&lkey) {
            if let Some(origin) = link_origin(&lkey) {
                cx.apply_link(it, link, agi, agj, origin);
            }
            continue;
        }
        if let Some((pkey, swap)) = cx.peptide_link(i, j, &cl1, &cl2) {
            if let Some(link) = ml.link(pkey) {
                let (g1, g2) = if swap { (agj, agi) } else { (agi, agj) };
                cx.apply_link(it, link, g1, g2, link_origin(pkey).unwrap_or(ORIGIN_MISC));
            }
            continue;
        }
        let origin = if c1 == "metal" || c2 == "metal" { ORIGIN_METAL } else { ORIGIN_MISC };
        custom.push((i, j, origin));
        custom_partners.entry(i).or_default().push(j);
        custom_partners.entry(j).or_default().push(i);
    }
    // custom bonds go in after the loop, in atom order
    custom.sort_by_key(|&(i, j, _)| (i, j));
    for (i, j, origin) in custom {
        let ideal = default_bond_length(cx.elem(i), cx.elem(j));
        if it.add_bond(i, j, ideal, origin) {
            log.push_str(&format!("  Custom bond: {} - {}\n", cx.id_str(i), cx.id_str(j)));
        }
    }
}

/// `bondlength_defaults.get_default_bondlength` for a single bond (2.3 A when
/// the pair is not tabulated, as `process_nonbonded_for_links` does).
fn default_bond_length(e1: &str, e2: &str) -> f64 {
    let fix = |e: &str| -> String { e.to_string() };
    let (s1, s2) = (fix(e1), fix(e2));
    for (a, b, o) in [(&s1, &s2, 1u8), (&s2, &s1, 1), (&s1, &s2, 0), (&s2, &s1, 0)] {
        if let Some(&(_, _, _, v)) = linkdata::QM_DEFAULTS.iter().find(|(x, y, oo, _)| *x == a.as_str() && *y == b.as_str() && *oo == o) {
            return v;
        }
    }
    2.3
}
