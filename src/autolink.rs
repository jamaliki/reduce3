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
use rustc_hash::{FxHashMap, FxHashSet};
use std::cell::RefCell;
use std::sync::Arc;

#[cfg(test)]
mod tests;

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
struct Ag<'a> {
    model: u32,
    chain: u32,
    rg: u32,
    /// The flat indices of the atoms, which are consecutive.
    atoms: std::ops::Range<u32>,
    /// All atoms of the residue group (every altloc), also consecutive.
    rg_atoms: std::ops::Range<u32>,
    resname: &'a str,
    altloc: &'a str,
}

pub struct AutoLinkParams {
    pub compat: bool,
}

/// Maps a residue name and a flag for groups of more than one atom to a
/// class and a saccharide flag.
type GroupClassCache<'a> = FxHashMap<(&'a str, bool), (&'static str, bool)>;

struct Ctx<'a> {
    st: &'a Structure,
    flat: &'a FlatAtoms,
    ml: &'a MonLib,
    p: &'a AutoLinkParams,
    ag_of: Vec<u32>,
    ags: Vec<Ag<'a>>,
    /// 0 or 1 for atoms in the first or last residue group of a chain.
    first_last: Vec<i8>,
    type_cache: RefCell<FxHashMap<String, Option<String>>>,
    class_cache: RefCell<FxHashMap<u32, (&'static str, bool, bool)>>,
    /// The class and saccharide flag of `classes` by residue name and
    /// whether the atom group has more than one atom, on which they depend.
    group_class_cache: RefCell<GroupClassCache<'a>>,
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
        let mut first_last = vec![-1i8; n];
        let mut ags: Vec<Ag<'a>> = Vec::new();
        let mut last: Option<(u32, u32, u32, u32)> = None;
        // `first_last` of the atoms of the current atom group
        let mut fl = -1i8;
        for (k, pth) in flat.path.iter().enumerate() {
            let key = (pth.model, pth.chain, pth.rg, pth.ag);
            if last != Some(key) {
                let g = st.atom_group(*pth);
                let nrg = st.chain(*pth).residue_groups.len() as u32;
                fl = if nrg == 0 { -1 } else if pth.rg == nrg - 1 { 1 } else if pth.rg == 0 { 0 } else { -1 };
                ags.push(Ag {
                    model: pth.model,
                    chain: pth.chain,
                    rg: pth.rg,
                    atoms: k as u32..k as u32,
                    rg_atoms: 0..0,
                    resname: &g.resname,
                    altloc: &g.altloc,
                });
                last = Some(key);
            }
            let gi = ags.len() - 1;
            ags[gi].atoms.end = k as u32 + 1;
            ag_of[k] = gi as u32;
            first_last[k] = fl;
        }
        // residue-group atom ranges: consecutive atom groups of one residue
        // group hold consecutive atoms
        let mut s = 0;
        while s < ags.len() {
            let mut e = s;
            while e < ags.len() && (ags[e].model, ags[e].chain, ags[e].rg) == (ags[s].model, ags[s].chain, ags[s].rg) {
                e += 1;
            }
            let all = ags[s].atoms.start..ags[e - 1].atoms.end;
            for g in &mut ags[s..e] {
                g.rg_atoms = all.clone();
            }
            s = e;
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
            group_class_cache: RefCell::new(FxHashMap::default()),
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
    fn ag(&self, a: u32) -> &Ag<'a> {
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

    /// The class of an atom group's residue for `classes`, before metals,
    /// and whether `classes` counts it as a saccharide.
    fn group_class(&self, gi: u32) -> (&'static str, bool) {
        let g = &self.ags[gi as usize];
        let consider = g.atoms.len() != 1;
        if let Some(&v) = self.group_class_cache.borrow().get(&(g.resname, consider)) {
            return v;
        }
        let mut gc = resclass::get_class_ext(g.resname, consider).name();
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
            self.is_sugar_type(g.resname)
        };
        self.group_class_cache.borrow_mut().insert((g.resname, consider), (gc, sacch));
        (gc, sacch)
    }

    /// `linking_utils.get_classes(atom)`.
    fn classes(&self, a: u32) -> Classes {
        let gi = self.ag_of[a as usize];
        let cached = self.class_cache.borrow().get(&gi).copied();
        let (gc, sacch, uncommon) = match cached {
            Some(v) => v,
            None => {
                let g = &self.ags[gi as usize];
                let (gc, sacch) = self.group_class(gi);
                let bb = g
                    .rg_atoms
                    .clone()
                    .filter(|&x| matches!(self.name(x).trim(), "C" | "CA" | "N" | "O" | "OXT"))
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
            if &rg.resseq != resseq || g.resname != resname {
                continue;
            }
            let mut out = Vec::new();
            for x in g.atoms.clone() {
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

    /// Sum the present, compatible dictionary bonds of a carbon in half units.
    fn carbon_valence(&self, it: &Interp, initial: &Initial, a: u32, partner: u32) -> Option<u8> {
        if self.elem(a) != "C" { return None; }
        let info = &it.dict[a as usize];
        let (Some(comp), Some(id)) = (&info.comp, &info.dict_id) else { return None };
        if comp.atom(id).is_none_or(|atom| atom.type_symbol != "C") { return None; }
        let alt = if self.ag(a).altloc.is_empty() { &self.ag(partner).altloc } else { &self.ag(a).altloc };
        let mut seen = FxHashSet::default();
        let mut names = FxHashSet::default();
        let mut value = 0u8;
        for &other in initial.of(a) {
            let group = self.ag(other);
            let parent = self.ag(a);
            if (group.model, group.chain, group.rg) != (parent.model, parent.chain, parent.rg)
                || (!alt.is_empty() && !group.altloc.is_empty() && alt != &group.altloc) { continue; }
            let info = &it.dict[other as usize];
            let (Some(source), Some(other_id)) = (&info.comp, &info.dict_id) else { return None };
            if !Arc::ptr_eq(source, comp) || it.dict[other as usize].instance != it.dict[a as usize].instance { return None; }
            if !seen.insert((other_id.as_str(), group.altloc)) { return None; }
            if !names.insert(other_id.as_str()) { continue; }
            let bonds: Vec<_> = comp.bonds.iter().filter(|b|
                (&b.a1 == id && &b.a2 == other_id) || (&b.a2 == id && &b.a1 == other_id)).collect();
            let [bond] = bonds[..] else { return None };
            value = value.checked_add(crate::interp::source_bond_order(self.ml, comp, bond, parent.resname)?)?;
        }
        Some(value)
    }

    fn named_link(&self, i: u32, j: u32, c1: &Classes, c2: &Classes) -> Option<(&'a ChemLink, String, u32, u32)> {
        let (gi, gj) = (self.ag(i), self.ag(j));
        let (g1, g2) = (self.ag_of[i as usize], self.ag_of[j as usize]);
        let keys = [format!("{}-{}", gi.resname, gj.resname), format!("{}-{}", gj.resname, gi.resname)];
        for (key, a, b) in [(&keys[0], g1, g2), (&keys[1], g2, g1)] {
            if let Some(link) = self.ml.link(key) { return Some((link, key.clone(), a, b)); }
        }
        let peptides = ["common_amino_acid", "d_amino_acid", "uncommon_amino_acid"];
        if peptides.contains(&c1.important) && peptides.contains(&c2.important)
            && ((self.name(i) == " N  " && self.name(j) == " C  ") || (self.name(i) == " C  " && self.name(j) == " N  ")) {
            return self.ml.link("TRANS").map(|link| {
                let (a, b) = if self.name(i) == " C  " { (g2, g1) } else { (g1, g2) };
                (link, "TRANS".into(), a, b)
            });
        }
        let keys = [format!("{}_{}-ANY_{}", gi.resname, self.name(i).trim(), self.name(j).trim()),
            format!("{}_{}-ANY_{}", gj.resname, self.name(j).trim(), self.name(i).trim())];
        for (key, a, b) in [(&keys[0], g1, g2), (&keys[1], g2, g1)] {
            if let Some(link) = self.ml.link(key) { return Some((link, key.clone(), a, b)); }
        }
        None
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
        self.ags[g as usize].atoms.clone().find(|&x| self.name(x).trim() == n)
    }

    fn peptide_pair(&self, it: &Interp, a: u32, b: u32) -> bool {
        let (ga, gb) = (self.ag(a), self.ag(b));
        ga.model == gb.model && (ga.chain, ga.rg) != (gb.chain, gb.rg)
            && matches!((self.name(a).trim(), self.elem(a), self.name(b).trim(), self.elem(b)),
                ("C", "C", "N", "N") | ("N", "N", "C", "C"))
            && [a, b].iter().all(|&index|
                resclass::get_class(self.ag(index).resname.trim()).is_amino_acid()
                    || it.dict[index as usize].comp.as_ref().is_some_and(|comp| comp.test_for_peptide()))
    }

    fn has_peptide_partner(&self, it: &Interp, initial: &Initial, initial_bonds: usize,
                           atom: u32, candidate: u32) -> bool {
        let compatible = |other: u32| {
            let (a, b) = (&self.ag(other).altloc, &self.ag(candidate).altloc);
            a.is_empty() || b.is_empty() || a == b
        };
        // peptide bonds made before linking
        initial.of(atom).iter().zip(initial.bonds_of(atom)).any(|(&other, &k)|
            peptide_origin(it.bonds[k as usize].origin) && self.peptide_pair(it, atom, other) && compatible(other))
            || it.bonds[initial_bonds..].iter().any(|bond| {
                let other = if bond.i == atom { bond.j } else if bond.j == atom { bond.i } else { return false };
                peptide_origin(bond.origin) && self.peptide_pair(it, atom, other) && compatible(other)
            })
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
        if !bonds.is_empty() {
            let anchors = [self.ags[g1 as usize].atoms.start, self.ags[g2 as usize].atoms.start];
            it.automatic_links.push(crate::interp::AutomaticLink::Dictionary {
                id: link.id.clone(), anchors, origin,
            });
        }
        bonds
    }

    /// `glyco_utils.get_glyco_link_atoms(atom_group1, atom_group2)`.
    fn glyco_link_atoms(&self, g1: u32, g2: u32, cutoff: f64) -> Option<GlycoAtoms> {
        // bonds by distance over the atoms of both residue groups
        let all: Vec<u32> =
            self.ags[g1 as usize].rg_atoms.clone().chain(self.ags[g2 as usize].rg_atoms.clone()).collect();
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
            self.ags[g1 as usize].atoms.clone().chain(self.ags[g2 as usize].atoms.clone()).collect();
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
            self.ags[link_group as usize].atoms.clone().find(|&x| el(x) == "O" && self.d2(x, ac) < 5.0)
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
        let mut angles = Vec::new();
        for (atoms, v) in [
            ([lc, Some(g.link_oxygen), Some(g.anomeric_carbon)], 108.7),
            ([Some(g.link_oxygen), Some(g.anomeric_carbon), g.ring_oxygen], 112.3),
            ([Some(g.link_oxygen), Some(g.anomeric_carbon), g.ring_carbon], 109.47),
            ([Some(g.link_oxygen), Some(g.anomeric_carbon), g.anomeric_h], 109.47),
        ] {
            if let [Some(a), Some(b), Some(c)] = atoms {
                it.add_angle(a, b, c, v, 3.0);
                angles.push(crate::interp::AngleProxy { i: [a, b, c], ideal: v, esd: 3.0 });
            }
        }
        let chirality = g.ring_oxygen.zip(g.ring_carbon).and_then(|(oxygen, carbon)| {
            let rn = self.ag(g.anomeric_carbon).resname.trim();
            linkdata::GLYCO_VOLUMES.iter().find(|(name, _)| *name == rn)
                .map(|(_, value)| ([g.anomeric_carbon, g.link_oxygen, oxygen, carbon], *value))
        });
        it.automatic_links.push(crate::interp::AutomaticLink::Glycosidic {
            id: isomer, bond: [g.anomeric_carbon, g.link_oxygen], angles, chirality, origin,
        });
        Some((g.anomeric_carbon, g.link_oxygen))
    }

    /// `glyco_link_class.get_isomer`.
    fn isomer(&self, g: &GlycoAtoms) -> String {
        let rn = self.ag(g.anomeric_carbon).resname;
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
        let t1 = self.ccd_type(self.ag(a1).resname);
        let t2 = self.ccd_type(self.ag(a2).resname);
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
        let (rn1, rn2) = (self.ag(a1).resname, self.ag(a2).resname);
        let (t1, t2) = (self.ccd_type(rn1), self.ccd_type(rn2));
        let both = t1.is_some() && t2.is_some();
        let count = |names: &[&str]| -> bool {
            let s = sugar(&t1) as u8 + sugar(&t2) as u8;
            let l = (!sugar(&t1) && names.contains(&rn1)) as u8 + (!sugar(&t2) && names.contains(&rn2)) as u8;
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
        if resclass::get_class(og.resname).name() != "other" {
            return None;
        }
        let (n1, n2) = (self.name(a1).trim(), self.name(a2).trim());
        if n1 == "SG" && n2 == "SG" {
            return Some(("SS", false));
        }
        let count = og.atoms.clone().filter(|&x| matches!(self.name(x).trim(), "C" | "N" | "O")).count();
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

/// The bonds made before linking: each atom's partners and bond indices, in
/// bond order.
struct Initial {
    start: Vec<u32>,
    partner: Vec<u32>,
    bond: Vec<u32>,
}

impl Initial {
    fn new(it: &Interp, n: usize) -> Initial {
        let mut start = vec![0u32; n + 1];
        for b in &it.bonds {
            start[b.i as usize + 1] += 1;
            start[b.j as usize + 1] += 1;
        }
        for k in 0..n {
            start[k + 1] += start[k];
        }
        let mut fill = start.clone();
        let mut partner = vec![0u32; start[n] as usize];
        let mut bond = vec![0u32; start[n] as usize];
        for (k, b) in it.bonds.iter().enumerate() {
            for (x, y) in [(b.i, b.j), (b.j, b.i)] {
                let slot = &mut fill[x as usize];
                partner[*slot as usize] = y;
                bond[*slot as usize] = k as u32;
                *slot += 1;
            }
        }
        Initial { start, partner, bond }
    }
    /// Adjacency lists given directly, without bond indices.
    #[cfg(test)]
    fn from_lists(lists: &[Vec<u32>]) -> Initial {
        let mut start = vec![0u32];
        for l in lists {
            start.push(start[start.len() - 1] + l.len() as u32);
        }
        let partner: Vec<u32> = lists.iter().flatten().copied().collect();
        let bond = vec![u32::MAX; partner.len()];
        Initial { start, partner, bond }
    }
    fn of(&self, a: u32) -> &[u32] {
        &self.partner[self.start[a as usize] as usize..self.start[a as usize + 1] as usize]
    }
    fn bonds_of(&self, a: u32) -> &[u32] {
        &self.bond[self.start[a as usize] as usize..self.start[a as usize + 1] as usize]
    }
}

/// Heavy-atom pairs closer than `MAX_BONDED_CUTOFF` within each model, in
/// order of distance (then atom indices). Pairs the main loop would skip
/// without side effects whatever happened before them are left out: atoms
/// bonded before linking starts, atoms of one residue group, and different
/// non-blank altlocs.
fn candidate_pairs(cx: &Ctx, initial: &Initial) -> Vec<(f64, u32, u32)> {
    let flat = cx.flat;
    let n = flat.pos.len();
    let rg_of = |a: u32| {
        let p = flat.path[a as usize];
        (p.model, p.chain, p.rg)
    };
    // altloc code per heavy atom: 0 when blank, else 1 + index of the
    // trimmed altloc among those seen
    let mut alt_names: Vec<&str> = Vec::new();
    let mut alt = vec![0u8; n];
    // heavy atoms per model
    let mut by_model: Vec<Vec<u32>> = Vec::new();
    for a in (0..n).filter(|&a| !flat.is_h[a]) {
        let t = flat.altloc[a].trim();
        if !t.is_empty() {
            let k = match alt_names.iter().position(|&x| x == t) {
                Some(k) => k,
                None => {
                    alt_names.push(t);
                    alt_names.len() - 1
                }
            };
            alt[a] = u8::try_from(k + 1).unwrap_or(u8::MAX);
        }
        let m = flat.path[a].model as usize;
        if by_model.len() <= m {
            by_model.resize(m + 1, Vec::new());
        }
        by_model[m].push(a as u32);
    }
    // more than 254 altlocs share a code; compare their names then
    let same_alt = |a: u32, b: u32| {
        let (x, y) = (alt[a as usize], alt[b as usize]);
        x == 0 || y == 0 || (x == y && (x != u8::MAX || flat.altloc[a as usize].trim() == flat.altloc[b as usize].trim()))
    };
    // cells a little wider than the cutoff, so that atoms closer than it lie
    // in the same or adjacent cells; each pair of cells is visited once
    let cell = MAX_BONDED_CUTOFF * (1.0 + 1e-6);
    let limit2 = cell * cell;
    let mut half: Vec<[isize; 3]> = Vec::new();
    for dz in -1..=1isize {
        for dy in -1..=1isize {
            for dx in -1..=1isize {
                if (dz, dy, dx) > (0, 0, 0) {
                    half.push([dx, dy, dz]);
                }
            }
        }
    }
    let mut out: Vec<(f64, u32, u32)> = Vec::new();
    for atoms in by_model.iter().filter(|v| !v.is_empty()) {
        let mut lo = v3(f64::INFINITY, f64::INFINITY, f64::INFINITY);
        let mut hi = v3(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
        for &a in atoms {
            let p = flat.pos[a as usize];
            lo = v3(lo.x.min(p.x), lo.y.min(p.y), lo.z.min(p.z));
            hi = v3(hi.x.max(p.x), hi.y.max(p.y), hi.z.max(p.z));
        }
        let dims = [
            ((hi.x - lo.x) / cell).floor() as usize + 1,
            ((hi.y - lo.y) / cell).floor() as usize + 1,
            ((hi.z - lo.z) / cell).floor() as usize + 1,
        ];
        let axis = |v: f64, l: f64, k: usize| (((v - l) / cell).floor() as usize).min(dims[k] - 1);
        let cell_of = |p: Vec3| axis(p.x, lo.x, 0) + dims[0] * (axis(p.y, lo.y, 1) + dims[1] * axis(p.z, lo.z, 2));
        let ncell = dims[0] * dims[1] * dims[2];
        let mut starts = vec![0u32; ncell + 1];
        let cells: Vec<usize> = atoms.iter().map(|&a| cell_of(flat.pos[a as usize])).collect();
        for &c in &cells {
            starts[c + 1] += 1;
        }
        for c in 0..ncell {
            starts[c + 1] += starts[c];
        }
        let mut fill = starts.clone();
        let mut items = vec![0u32; atoms.len()];
        for (k, &a) in atoms.iter().enumerate() {
            items[fill[cells[k]] as usize] = a;
            fill[cells[k]] += 1;
        }
        let range = |c: usize| starts[c] as usize..starts[c + 1] as usize;
        let planes: Vec<usize> = (0..dims[2]).collect();
        let found = crate::par::flat_map_collect(&planes, |&z| {
            let mut v: Vec<(f64, u32, u32)> = Vec::new();
            let mut test = |p: u32, q: u32| {
                let (a, b) = if p < q { (p, q) } else { (q, p) };
                let dv = flat.pos[b as usize] - flat.pos[a as usize];
                if dv.x * dv.x + dv.y * dv.y + dv.z * dv.z > limit2 {
                    return;
                }
                let d = dv.length();
                if d >= MAX_BONDED_CUTOFF || rg_of(b) == rg_of(a) || !same_alt(a, b) || initial.of(a).contains(&b) {
                    return;
                }
                v.push((d, a, b));
            };
            for y in 0..dims[1] {
                for x in 0..dims[0] {
                    let c = x + dims[0] * (y + dims[1] * z);
                    let own = range(c);
                    if own.is_empty() {
                        continue;
                    }
                    for k in own.clone() {
                        for l in k + 1..own.end {
                            test(items[k], items[l]);
                        }
                    }
                    for d in &half {
                        let (nx, ny, nz) = (x as isize + d[0], y as isize + d[1], z as isize + d[2]);
                        if nx < 0 || ny < 0 || nx >= dims[0] as isize || ny >= dims[1] as isize || nz >= dims[2] as isize {
                            continue;
                        }
                        let other = range(nx as usize + dims[0] * (ny as usize + dims[1] * nz as usize));
                        for k in own.clone() {
                            for l in other.clone() {
                                test(items[k], items[l]);
                            }
                        }
                    }
                }
            }
            v.into_iter()
        });
        out.extend(found);
    }
    crate::par::sort_unstable_by_key(&mut out, |&(d, i, j)| (d.to_bits(), i, j));
    out
}

/// Atoms of amino acids, nucleotides and waters other than sulfur link to
/// each other only when closer than this: `is_atom_pair_linked` passes such
/// a pair only for amino acid with amino acid or nucleotide, within
/// `AMINO_ACID_CUTOFF` plus the second-row buffer (1.965 A).
const NARROW_CUTOFF: f64 = 2.0;

/// Fixed mode: the pairs of `candidate_pairs` whose processing can change the
/// result. The main loop changes its state only for a pair that passes
/// `is_atom_pair_linked`. A pair that fails it only records its atom names
/// for its pair of atom groups, and only later pairs of the same two groups
/// that pass read the record. So the pairs kept are all those of each pair
/// of atom groups with at least one passing pair, in the same order.
fn linking_candidates(cx: &Ctx, initial: &Initial) -> Vec<(f64, u32, u32)> {
    let flat = cx.flat;
    let n = flat.pos.len();
    // Groups whose atoms other than sulfur link only within NARROW_CUTOFF:
    // those whose `classes(..).important` is an amino acid, a nucleotide or
    // water, which it is for every atom of the group. `group_class` gives
    // the same class without the backbone count.
    let narrow_kind: Vec<Narrow> = (0..cx.ags.len() as u32).map(|g| match cx.group_class(g) {
        (_, true) => Narrow::No,
        ("common_amino_acid", _) => Narrow::AminoAcid,
        ("common_rna_dna", _) => Narrow::Nucleotide,
        ("common_water", _) => Narrow::Water,
        _ => Narrow::No,
    }).collect();
    // What the narrow pass reads of an atom. A pair with water never passes
    // `is_atom_pair_linked`, nor does a pair of two nucleotides. An amino acid
    // passes with a nucleotide, and with another amino acid only as N with C
    // by element or as atoms named N and C (`possible_cyclic_peptide`).
    const NARROW: u8 = 1;
    const AMINO: u8 = 2;
    const NUCLEOTIDE: u8 = 4;
    const AMIDE_N: u8 = 8;
    const AMIDE_C: u8 = 16;
    let kind = |a: u32| -> u8 {
        if cx.elem(a) == "S" {
            return 0;
        }
        match narrow_kind[cx.ag_of[a as usize] as usize] {
            Narrow::No => 0,
            Narrow::Water => NARROW,
            Narrow::Nucleotide => NARROW | NUCLEOTIDE,
            Narrow::AminoAcid => {
                let (e, name) = (cx.elem(a), cx.name(a));
                let n_bit = if e == "N" || name == " N  " { AMIDE_N } else { 0 };
                let c_bit = if e == "C" || name == " C  " { AMIDE_C } else { 0 };
                NARROW | AMINO | n_bit | c_bit
            }
        }
    };
    let same_alt = |a: u32, b: u32| {
        let (x, y) = (flat.altloc[a as usize].trim(), flat.altloc[b as usize].trim());
        x.is_empty() || y.is_empty() || x == y
    };
    let rg_of = |a: u32| {
        let p = flat.path[a as usize];
        (p.model, p.chain, p.rg)
    };
    let full = MAX_BONDED_CUTOFF * (1.0 + 1e-6);
    // the pair as `candidate_pairs` records it, if it does
    let candidate = |p: u32, q: u32| -> Option<(f64, u32, u32)> {
        let (a, b) = if p < q { (p, q) } else { (q, p) };
        let dv = flat.pos[b as usize] - flat.pos[a as usize];
        if dv.x * dv.x + dv.y * dv.y + dv.z * dv.z > full * full {
            return None;
        }
        let d = dv.length();
        if d >= MAX_BONDED_CUTOFF || rg_of(b) == rg_of(a) || !same_alt(a, b) || initial.of(a).contains(&b) {
            return None;
        }
        Some((d, a, b))
    };
    let peptides = ["common_amino_acid", "d_amino_acid", "uncommon_amino_acid"];
    let passes = |d: f64, i: u32, j: u32| {
        let (cl1, cl2) = (cx.classes(i), cx.classes(j));
        let only_cutoff = peptides.contains(&cl1.important) && peptides.contains(&cl2.important) && cx.possible_cyclic_peptide(i, j);
        cx.is_atom_pair_linked(i, j, cl1.important, cl2.important, d, only_cutoff)
    };
    let mut by_model: Vec<Vec<u32>> = Vec::new();
    for a in (0..n).filter(|&a| !flat.is_h[a]) {
        let m = flat.path[a].model as usize;
        if by_model.len() <= m {
            by_model.resize(m + 1, Vec::new());
        }
        by_model[m].push(a as u32);
    }
    // residue number of each atom; atoms of one residue are contiguous
    let mut rg_id = vec![0u32; n];
    let (mut id, mut prev) = (0u32, None);
    for g in &cx.ags {
        let key = (g.model, g.chain, g.rg);
        if prev.is_some_and(|p| p != key) {
            id += 1;
        }
        prev = Some(key);
        rg_id[g.atoms.start as usize..g.atoms.end as usize].fill(id);
    }
    let cell = NARROW_CUTOFF * (1.0 + 1e-6);
    let mut keys: FxHashSet<(u32, u32)> = FxHashSet::default();
    let mut found = |d: f64, a: u32, b: u32| {
        if passes(d, a, b) {
            let (x, y) = (cx.ag_of[a as usize], cx.ag_of[b as usize]);
            keys.insert((x.min(y), x.max(y)));
        }
    };
    for atoms in by_model.iter().filter(|v| !v.is_empty()) {
        let mut lo = v3(f64::INFINITY, f64::INFINITY, f64::INFINITY);
        let mut hi = v3(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
        for &a in atoms {
            let p = flat.pos[a as usize];
            lo = v3(lo.x.min(p.x), lo.y.min(p.y), lo.z.min(p.z));
            hi = v3(hi.x.max(p.x), hi.y.max(p.y), hi.z.max(p.z));
        }
        let dims = [
            ((hi.x - lo.x) / cell).floor() as usize + 1,
            ((hi.y - lo.y) / cell).floor() as usize + 1,
            ((hi.z - lo.z) / cell).floor() as usize + 1,
        ];
        let axis = |v: f64, l: f64, k: usize| (((v - l) / cell).floor() as usize).min(dims[k] - 1);
        let coords = |p: Vec3| [axis(p.x, lo.x, 0), axis(p.y, lo.y, 1), axis(p.z, lo.z, 2)];
        // The atoms are sorted by row of cells (y, z), then by x cell, then
        // by index: a counting sort by x cell, then a stable one by row. Each
        // row is a run of the list, found from `row_start`.
        let nrow = dims[1] * dims[2];
        // x cell, row and kind of each atom, in index order
        let mut cell_x = Vec::with_capacity(atoms.len());
        let mut cell_row = Vec::with_capacity(atoms.len());
        let mut flags = Vec::with_capacity(atoms.len());
        let lists_of = [AMINO, NUCLEOTIDE, AMIDE_N, AMIDE_C];
        let mut sizes = [0usize; 4];
        for &a in atoms {
            let c = coords(flat.pos[a as usize]);
            cell_x.push(c[0] as u32);
            cell_row.push((c[1] + dims[1] * c[2]) as u32);
            let f = kind(a);
            for (size, bit) in sizes.iter_mut().zip(lists_of) {
                *size += usize::from(f & bit != 0);
            }
            flags.push(f);
        }
        let mut by_x = vec![0u32; dims[0] + 1];
        for &x in &cell_x {
            by_x[x as usize + 1] += 1;
        }
        for x in 0..dims[0] {
            by_x[x + 1] += by_x[x];
        }
        let mut order = vec![0u32; atoms.len()];
        for (k, &x) in cell_x.iter().enumerate() {
            let slot = &mut by_x[x as usize];
            order[*slot as usize] = k as u32;
            *slot += 1;
        }
        let mut row_start = vec![0u32; nrow + 1];
        for &r in &cell_row {
            row_start[r as usize + 1] += 1;
        }
        for r in 0..nrow {
            row_start[r + 1] += row_start[r];
        }
        let mut fill = row_start.clone();
        let mut sorted = vec![0u32; atoms.len()];
        for &k in &order {
            let slot = &mut fill[cell_row[k as usize] as usize];
            sorted[*slot as usize] = k;
            *slot += 1;
        }
        let items: Vec<u32> = sorted.iter().map(|&k| atoms[k as usize]).collect();
        let xs: Vec<u32> = sorted.iter().map(|&k| cell_x[k as usize]).collect();
        let is_narrow: Vec<bool> = sorted.iter().map(|&k| flags[k as usize] & NARROW != 0).collect();
        // the atoms of each narrow list in the same order
        let mut lists = sizes.map(|size| RowAtoms::with_capacity(nrow, size));
        for r in 0..nrow {
            for &k in &sorted[row_start[r] as usize..row_start[r + 1] as usize] {
                let f = flags[k as usize];
                if f & (AMINO | NUCLEOTIDE) == 0 {
                    continue;
                }
                let a = atoms[k as usize];
                let (x, p, rg) = (cell_x[k as usize], flat.pos[a as usize], rg_id[a as usize]);
                for (list, bit) in lists.iter_mut().zip(lists_of) {
                    if f & bit != 0 {
                        list.push(a, x, p, rg);
                    }
                }
            }
            for list in &mut lists {
                list.end_row();
            }
        }
        let [amino, nucleotide, amide_n, amide_c] = lists;
        let near2 = cell * cell;
        // narrow pairs from different residues in the same or adjacent cells
        for (own, other) in [(&amino, &nucleotide), (&amide_n, &amide_c)] {
            adjacent_pairs(own, other, dims, |k, l| {
                if other.rgs[l] == own.rgs[k] || (other.pos[l] - own.pos[k]).length_sq() > near2 {
                    return;
                }
                if let Some((d, a, b)) = candidate(own.items[k], other.items[l]) {
                    found(d, a, b);
                }
            });
        }
        // pairs with an atom that is not narrow, up to the full cutoff
        let reach = (MAX_BONDED_CUTOFF / cell).ceil() as usize;
        for k in (0..items.len()).filter(|&k| !is_narrow[k]) {
            let c = coords(flat.pos[items[k] as usize]);
            let (x0, x1) = (c[0].saturating_sub(reach) as u32, (c[0] + reach).min(dims[0] - 1) as u32);
            for z in c[2].saturating_sub(reach)..=(c[2] + reach).min(dims[2] - 1) {
                for y in c[1].saturating_sub(reach)..=(c[1] + reach).min(dims[1] - 1) {
                    let r = y + dims[1] * z;
                    let (first, last) = (row_start[r] as usize, row_start[r + 1] as usize);
                    // the atoms of cells x0 to x1 of the row
                    let from = first + xs[first..last].partition_point(|&x| x < x0);
                    let to = first + xs[first..last].partition_point(|&x| x <= x1);
                    for l in from..to {
                        // a pair of two such atoms is taken from its first atom
                        if l == k || (!is_narrow[l] && items[l] < items[k]) {
                            continue;
                        }
                        if let Some((d, a, b)) = candidate(items[k], items[l]) {
                            found(d, a, b);
                        }
                    }
                }
            }
        }
    }
    // every candidate pair of the atom groups with a passing pair
    let mut out: Vec<(f64, u32, u32)> = Vec::new();
    for &(x, y) in &keys {
        for a in cx.ags[x as usize].atoms.clone() {
            if flat.is_h[a as usize] {
                continue;
            }
            for b in cx.ags[y as usize].atoms.clone() {
                if flat.is_h[b as usize] {
                    continue;
                }
                out.extend(candidate(a, b));
            }
        }
    }
    out.sort_unstable_by_key(|&(d, i, j)| (d.to_bits(), i, j));
    out
}

/// What an atom is for the narrow pass of `linking_candidates`.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Narrow {
    No,
    AminoAcid,
    Nucleotide,
    Water,
}

/// Atoms sorted by row of grid cells and then by x cell, with the data the
/// pair test reads.
struct RowAtoms {
    /// The first atom of each row, and the end of the list.
    start: Vec<u32>,
    items: Vec<u32>,
    xs: Vec<u32>,
    pos: Vec<Vec3>,
    rgs: Vec<u32>,
}

impl RowAtoms {
    fn with_capacity(nrow: usize, atoms: usize) -> RowAtoms {
        let mut start = Vec::with_capacity(nrow + 1);
        start.push(0);
        RowAtoms {
            start,
            items: Vec::with_capacity(atoms),
            xs: Vec::with_capacity(atoms),
            pos: Vec::with_capacity(atoms),
            rgs: Vec::with_capacity(atoms),
        }
    }
    fn push(&mut self, a: u32, x: u32, p: Vec3, rg: u32) {
        self.items.push(a);
        self.xs.push(x);
        self.pos.push(p);
        self.rgs.push(rg);
    }
    /// End the current row.
    fn end_row(&mut self) {
        self.start.push(self.items.len() as u32);
    }
    fn row(&self, r: usize) -> std::ops::Range<usize> {
        self.start[r] as usize..self.start[r + 1] as usize
    }
}

/// Calls `f(k, l)` for each atom `k` of `own` and atom `l` of `other` in the
/// same or adjacent cells of a grid with `dims` cells. For each row of `own`
/// and each of the nine rows of `other` next to it, a pointer that only
/// moves forward finds the atoms of cells x - 1 to x + 1.
fn adjacent_pairs(own: &RowAtoms, other: &RowAtoms, dims: [usize; 3], mut f: impl FnMut(usize, usize)) {
    for z in 0..dims[2] {
        for y in 0..dims[1] {
            let mine = own.row(y + dims[1] * z);
            if mine.is_empty() {
                continue;
            }
            for nz in z.saturating_sub(1)..=(z + 1).min(dims[2] - 1) {
                for ny in y.saturating_sub(1)..=(y + 1).min(dims[1] - 1) {
                    let theirs = other.row(ny + dims[1] * nz);
                    if theirs.is_empty() {
                        continue;
                    }
                    let mut lo = theirs.start;
                    for k in mine.clone() {
                        let x = own.xs[k];
                        while lo < theirs.end && other.xs[lo] + 1 < x {
                            lo += 1;
                        }
                        for l in lo..theirs.end {
                            if other.xs[l] > x + 1 {
                                break;
                            }
                            f(k, l);
                        }
                    }
                }
            }
        }
    }
}

/// Keep unique records for custom-link guards and resolved records for peptide links.
fn declared_covalent_pairs(st: &Structure, flat: &FlatAtoms) -> (FxHashSet<(u32, u32)>, FxHashSet<(u32, u32)>) {
    if !st.links.iter().any(|link| link.kind == "covale") { return (FxHashSet::default(), FxHashSet::default()); }
    let labels = st.links.iter().filter(|link| link.kind == "covale").flat_map(|link| [&link.atom1, &link.atom2]);
    let resolver = EndpointResolver::for_labels(st, &flat.path, labels);
    let mut pairs = FxHashMap::default();
    for link in &st.links {
        if link.kind != "covale" { continue; }
        let matches = |label: &AtomLabel| resolver.matching_indices(label).collect::<Vec<_>>();
        let (first, second) = (matches(&link.atom1), matches(&link.atom2));
        let ([i], [j]) = (&first[..], &second[..]) else { continue };
        let (a, b) = (flat.path[*i as usize], flat.path[*j as usize]);
        if a.model != b.model || (a.chain == b.chain && a.rg == b.rg)
            || (!flat.altloc[*i as usize].is_empty() && !flat.altloc[*j as usize].is_empty()
                && flat.altloc[*i as usize] != flat.altloc[*j as usize]) { continue; }
        let entry = pairs.entry(((*i).min(*j), (*i).max(*j))).or_insert((0usize, true));
        entry.0 += 1;
        entry.1 &= link.order.as_deref().is_none_or(|order|
            order.eq_ignore_ascii_case("sing") || order.eq_ignore_ascii_case("single"));
    }
    let mut unique = FxHashSet::default();
    let mut resolved = FxHashSet::default();
    for (pair, (count, valid)) in pairs {
        if valid {
            resolved.insert(pair);
            if count == 1 { unique.insert(pair); }
        }
    }
    (unique, resolved)
}

fn peptide_origin(origin: u16) -> bool {
    // Polymer bonds and the CIS, NMCIS, NMTRANS, PCIS, PTRANS and TRANS dictionaries.
    matches!(origin, crate::interp::ORIGIN_COVALENT | 35 | 63 | 65 | 66 | 70 | 79)
}

/// Make the automatic links. Bonds (with their angles, dihedrals and planes)
/// go into `it`; a line per link goes to `log`.
pub fn auto_link(it: &mut Interp, st: &Structure, flat: &FlatAtoms, ml: &MonLib, p: &AutoLinkParams, log: &mut String) {
    let cx = Ctx::new(st, flat, ml, p);
    let n = flat.pos.len();
    let initial = Initial::new(it, n);
    let initial_bonds = it.bonds.len();
    // `done`: atoms that took their link, and per residue pair the atom-name
    // pairs considered
    let mut atom_done: Vec<bool> = vec![false; n];
    let mut res_done: FxHashMap<(u32, u32), smallvec::SmallVec<[(u32, u32); 2]>> = FxHashMap::default();
    let mut n_links: Vec<u32> = vec![0; n];
    let mut custom: Vec<(u32, u32, u16)> = Vec::new();
    let mut custom_partners: FxHashMap<u32, Vec<u32>> = FxHashMap::default();
    let (declared, declared_peptides) = if p.compat { Default::default() } else { declared_covalent_pairs(st, flat) };

    let candidates = if p.compat { candidate_pairs(&cx, &initial) } else { linking_candidates(&cx, &initial) };
    for (dist, i, j) in candidates {
        if it.has_bond(i, j) {
            continue;
        }
        let (gi, gj) = (cx.ag(i), cx.ag(j));
        let moved = ["SF4", "F3S", "FES"];
        if moved.contains(&gi.resname) || moved.contains(&gj.resname) {
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
        let found = cx.named_link(i, j, &cl1, &cl2);
        let lkey = if found.is_none() { cx.single_link_key(i, j) } else { None };
        let is_glyco = lkey.as_ref().is_some_and(|key|
            ["ALPHA1", "BETA1", "ALPHA2", "BETA2"].iter().any(|s| key.contains(s)));
        let peptide = if found.is_none() && !is_glyco && lkey.as_ref().is_some_and(|key| ml.link(key).is_none()) {
            cx.peptide_link(i, j, &cl1, &cl2)
        } else { None };
        let annotated = declared.contains(&(i.min(j), i.max(j)));
        let selected_key = found.as_ref().map(|(_, key, _, _)| key.as_str())
            .or(peptide.map(|(key, _)| key)).or(lkey.as_deref());
        if !p.compat && selected_key.is_some_and(|key| matches!(key, "CIS" | "NMCIS" | "NMTRANS" | "PCIS" | "PTRANS" | "TRANS"))
            && cx.peptide_pair(it, i, j) && !declared_peptides.contains(&(i.min(j), i.max(j)))
            && (cx.has_peptide_partner(it, &initial, initial_bonds, i, j)
                || cx.has_peptide_partner(it, &initial, initial_bonds, j, i)) {
            continue;
        }
        let custom_fallback = found.is_none() && !is_glyco && peptide.is_none()
            && lkey.as_ref().is_some_and(|key| ml.link(key).is_none());
        let his_coordination = !p.compat && custom_fallback && !annotated && [(i, j), (j, i)].iter().any(|&(donor, metal)|
            cx.ag(donor).resname.trim() == "HIS" && cx.elem(donor) == "N"
                && matches!(cx.name(donor).trim(), "ND1" | "NE2")
                && METALS.contains(&cx.elem(metal)) && cx.ag(metal).rg_atoms.len() > 1);
        if !p.compat && custom_fallback && !annotated && !his_coordination
            && (METALS.contains(&cx.elem(i)) || METALS.contains(&cx.elem(j))
                || cx.carbon_valence(it, &initial, i, j).is_some_and(|value| value >= 8)
                || cx.carbon_valence(it, &initial, j, i).is_some_and(|value| value >= 8)) {
            continue;
        }
        // link_metals = Auto
        if (class_key.0 == "metal" || class_key.1 == "metal") && (p.compat || (found.is_none() && !annotated)) {
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
                let heavy = initial.of(a).iter().filter(|&&x| !flat.is_h[x as usize]).count();
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
            .clone()
            .any(|x| initial.of(x).iter().any(|&y| cx.ag_of[y as usize] == agj));
        if link_found {
            continue;
        }
        if let Some((link, lkey, g1, g2)) = found {
            let Some(origin) = link_origin(&lkey) else { continue };
            let made = cx.apply_link(it, link, g1, g2, origin);
            if !made.is_empty() {
                log.push_str(&format!("  Link {}: {} - {}\n", lkey, cx.id_str(made[0].0), cx.id_str(made[0].1)));
            }
            continue;
        }
        let Some(lkey) = lkey else {
            if let Some(v) = res_done.get_mut(&key) {
                if let Some(pos) = v.iter().position(|x| *x == names) {
                    v.remove(pos);
                }
            }
            continue;
        };
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
        if let Some((pkey, swap)) = peptide {
            if let Some(link) = ml.link(pkey) {
                let (g1, g2) = if swap { (agj, agi) } else { (agi, agj) };
                cx.apply_link(it, link, g1, g2, link_origin(pkey).unwrap_or(ORIGIN_MISC));
            }
            continue;
        }
        let origin = if his_coordination || c1 == "metal" || c2 == "metal" { ORIGIN_METAL } else { ORIGIN_MISC };
        custom.push((i, j, origin));
        custom_partners.entry(i).or_default().push(j);
        custom_partners.entry(j).or_default().push(i);
    }
    // custom bonds go in after the loop, in atom order
    custom.sort_by_key(|&(i, j, _)| (i, j));
    for (i, j, origin) in custom {
        let ideal = default_bond_length(cx.elem(i), cx.elem(j));
        if it.add_bond(i, j, ideal, origin) {
            it.automatic_links.push(crate::interp::AutomaticLink::Custom { indices: [i, j], origin });
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
