//! Restraint interpretation (the parts of cctbx `pdb_interpretation` that
//! Reduce2's hydrogen placement and optimizer depend on): per-residue
//! dictionaries and modifications, polymer links, disulfides, metal
//! coordination and automatic links, producing bond/angle/dihedral/plane
//! proxies in cctbx's registry order, plus energy types.

use crate::geom::*;
use crate::model::*;
use crate::monlib::{Comp, MonLib};
use crate::names;
use crate::resclass::{self, ResClass};
use rustc_hash::{FxHashMap, FxHashSet};
use std::sync::Arc;

pub const ORIGIN_COVALENT: u16 = 0;
pub const ORIGIN_SS: u16 = 1;
pub const ORIGIN_METAL: u16 = 3;
pub const ORIGIN_GLYCO_CUSTOM: u16 = 5;
pub const ORIGIN_MISC: u16 = 10;

#[derive(Clone, Copy, Debug)]
pub struct BondProxy {
    pub i: u32,
    pub j: u32,
    pub ideal: f64,
    pub origin: u16,
}

#[derive(Clone, Copy, Debug)]
pub struct AngleProxy {
    pub i: [u32; 3],
    pub ideal: f64,
    pub esd: f64,
}

#[derive(Clone, Debug)]
pub struct DihedralProxy {
    pub i: [u32; 4],
    pub ideal: f64,
    pub period: i32,
}

/// Energy type of an atom as cctbx records it.
#[derive(Clone, Debug, PartialEq)]
pub enum EType {
    Typed(String),
    /// Dictionary atom without a type (`'None'`).
    NoneType,
    /// Unexpected or duplicate atom (`'False'`).
    Unexpected,
    /// Atom of an unknown residue (`''`).
    Unknown,
}

#[derive(Clone, Debug, Default)]
pub struct AtomDictInfo {
    /// The residue dictionary used for the atom and its dictionary atom id.
    pub comp: Option<Arc<Comp>>,
    pub dict_id: Option<String>,
}

pub struct Interp {
    pub bonds: Vec<BondProxy>,
    bond_keys: FxHashMap<(u32, u32), usize>,
    pub angles: Vec<AngleProxy>,
    angle_keys: FxHashMap<(u32, u32, u32), usize>,
    pub dihedrals: Vec<DihedralProxy>,
    dih_keys: FxHashMap<[u32; 4], usize>,
    pub const_dihedrals: Vec<DihedralProxy>,
    const_keys: FxHashMap<[u32; 4], usize>,
    pub planes: Vec<Vec<u32>>,
    plane_keys: FxHashSet<Vec<u32>>,
    /// Bonds to a symmetry copy (cctbx `bond_asu_proxy`): disulfides across
    /// crystal symmetry, with `j` the atom whose copy is bonded.
    pub sym_bonds: Vec<BondProxy>,
    pub etype: Vec<EType>,
    pub dict: Vec<AtomDictInfo>,
    /// Residue dictionaries that could not be found (residue names).
    pub missing_residues: Vec<String>,
    pub log: String,
}

impl Interp {
    /// An interpretation with no restraints (to be filled by the caller).
    pub fn empty(n: usize) -> Interp {
        Interp::new(n)
    }
    fn new(n: usize) -> Interp {
        Interp {
            bonds: Vec::new(),
            bond_keys: FxHashMap::default(),
            angles: Vec::new(),
            angle_keys: FxHashMap::default(),
            dihedrals: Vec::new(),
            dih_keys: FxHashMap::default(),
            const_dihedrals: Vec::new(),
            const_keys: FxHashMap::default(),
            planes: Vec::new(),
            plane_keys: FxHashSet::default(),
            sym_bonds: Vec::new(),
            etype: vec![EType::Unknown; n],
            dict: vec![AtomDictInfo::default(); n],
            missing_residues: Vec::new(),
            log: String::new(),
        }
    }

    pub fn add_bond(&mut self, a: u32, b: u32, ideal: f64, origin: u16) -> bool {
        let (i, j) = if a < b { (a, b) } else { (b, a) };
        if i == j || self.bond_keys.contains_key(&(i, j)) {
            return false;
        }
        self.bond_keys.insert((i, j), self.bonds.len());
        self.bonds.push(BondProxy { i, j, ideal, origin });
        true
    }
    /// `bond_params_table.update`: add the bond, or replace the parameters of
    /// an existing one.
    pub fn upsert_bond(&mut self, a: u32, b: u32, ideal: f64, origin: u16) {
        let (i, j) = if a < b { (a, b) } else { (b, a) };
        if i == j {
            return;
        }
        match self.bond_keys.get(&(i, j)) {
            Some(&k) => {
                self.bonds[k].ideal = ideal;
                self.bonds[k].origin = origin;
            }
            None => {
                self.bond_keys.insert((i, j), self.bonds.len());
                self.bonds.push(BondProxy { i, j, ideal, origin });
            }
        }
    }
    /// Ideal length of a bond, if there is one.
    pub fn bond_ideal(&self, a: u32, b: u32) -> Option<f64> {
        let (i, j) = if a < b { (a, b) } else { (b, a) };
        self.bond_keys.get(&(i, j)).map(|&k| self.bonds[k].ideal)
    }
    pub fn has_bond(&self, a: u32, b: u32) -> bool {
        let (i, j) = if a < b { (a, b) } else { (b, a) };
        self.bond_keys.contains_key(&(i, j))
    }

    pub(crate) fn add_angle(&mut self, a: u32, b: u32, c: u32, ideal: f64, esd: f64) {
        let (x, z) = if a < c { (a, c) } else { (c, a) };
        let key = (b, x, z);
        if self.angle_keys.contains_key(&key) {
            return;
        }
        self.angle_keys.insert(key, self.angles.len());
        self.angles.push(AngleProxy { i: [x, b, z], ideal, esd });
    }
    pub fn has_angle(&self, a: u32, b: u32, c: u32) -> bool {
        let (x, z) = if a < c { (a, c) } else { (c, a) };
        self.angle_keys.contains_key(&(b, x, z))
    }
    pub fn push_angle_unchecked(&mut self, a: u32, b: u32, c: u32, ideal: f64, esd: f64) {
        self.add_angle(a, b, c, ideal, esd);
    }
    pub fn remove_angle(&mut self, a: u32, b: u32, c: u32) {
        let (x, z) = if a < c { (a, c) } else { (c, a) };
        if let Some(idx) = self.angle_keys.remove(&(b, x, z)) {
            self.angles.remove(idx);
            self.angle_keys.clear();
            for (k, p) in self.angles.iter().enumerate() {
                self.angle_keys.insert((p.i[1], p.i[0], p.i[2]), k);
            }
        }
    }

    /// Canonicalize like cctbx's dihedral `sort_i_seqs`: swap the ends if
    /// i0 > i3 and, independently, the middle pair if i1 > i2, negating the
    /// ideal angle for each swap.
    fn canon_dihedral(idx: [u32; 4], ideal: f64) -> ([u32; 4], f64) {
        let mut i = idx;
        let mut v = ideal;
        if i[0] > i[3] {
            i.swap(0, 3);
            v = -v;
        }
        if i[1] > i[2] {
            i.swap(1, 2);
            v = -v;
        }
        (i, v)
    }

    pub(crate) fn add_dihedral(&mut self, idx: [u32; 4], ideal: f64, period: i32, is_const: bool) {
        let (i, v) = Self::canon_dihedral(idx, ideal);
        if is_const {
            if self.const_keys.contains_key(&i) {
                return;
            }
            self.const_keys.insert(i, self.const_dihedrals.len());
            self.const_dihedrals.push(DihedralProxy { i, ideal: v, period: 0 });
        } else {
            if self.dih_keys.contains_key(&i) {
                return;
            }
            self.dih_keys.insert(i, self.dihedrals.len());
            self.dihedrals.push(DihedralProxy { i, ideal: v, period });
        }
    }

    pub(crate) fn add_plane(&mut self, atoms: Vec<u32>) {
        if atoms.len() < 4 {
            return;
        }
        let mut a = atoms;
        a.sort_unstable();
        a.dedup();
        if a.len() < 4 || self.plane_keys.contains(&a) {
            return;
        }
        self.plane_keys.insert(a.clone());
        self.planes.push(a);
    }
}

/// One residue of one conformer: atom indices (flat) and their model names.
#[derive(Clone)]
struct ConfResidue {
    atoms: Vec<u32>,
    resname: String,
    #[allow(dead_code)]
    rg: (usize, usize, usize), // model, chain, residue group
    link_to_previous: bool,
}

/// Interpreted residue: dictionary and id -> atom map.
struct ResInterp {
    comp: Option<Arc<Comp>>,
    id_to_atom: FxHashMap<String, u32>,
    is_peptide: bool,
    is_rna_dna: bool,
    is_rna2p: Option<bool>,
    is_water: bool,
}

pub struct InterpParams {
    pub neutron: bool,
    pub link_distance_cutoff: f64,
    /// Reproduce Reduce2's quirks (see `autolink`).
    pub compat: bool,
    /// Dictionaries built from the CCD during H placement, by residue name
    /// (Reduce2's `auto_<resname>` restraint objects).
    pub auto_comps: FxHashMap<String, Arc<Comp>>,
}

/// Flat per-atom view used during interpretation.
pub struct FlatAtoms {
    pub pos: Vec<Vec3>,
    pub element: Vec<String>,
    pub name: Vec<String>,
    pub altloc: Vec<String>,
    pub is_h: Vec<bool>,
    pub path: Vec<AtomPath>,
}

impl FlatAtoms {
    pub fn from_structure(st: &Structure) -> FlatAtoms {
        let paths = st.atom_paths();
        let mut f = FlatAtoms {
            pos: Vec::with_capacity(paths.len()),
            element: Vec::with_capacity(paths.len()),
            name: Vec::with_capacity(paths.len()),
            altloc: Vec::with_capacity(paths.len()),
            is_h: Vec::with_capacity(paths.len()),
            path: paths.clone(),
        };
        for p in &paths {
            let a = st.atom(*p);
            f.pos.push(a.xyz);
            f.element.push(a.elem().to_ascii_uppercase());
            f.name.push(a.name.clone());
            f.altloc.push(st.atom_group(*p).altloc.clone());
            f.is_h.push(a.is_hydrogen());
        }
        f
    }
}

fn dihedral_model(pos: &[Vec3], i: [u32; 4]) -> Option<f64> {
    dihedral_deg(pos[i[0] as usize], pos[i[1] as usize], pos[i[2] as usize], pos[i[3] as usize])
}

/// `angle_delta_deg(a1, a2, periodicity)` from cctbx geometry_restraints.
pub fn angle_delta_deg(a1: f64, a2: f64, periodicity: i32) -> f64 {
    let half = 180.0 / (periodicity.abs().max(1) as f64);
    let mut d = (a2 - a1) % (2.0 * half);
    if d < -half {
        d += 2.0 * half;
    } else if d > half {
        d -= 2.0 * half;
    }
    d
}

/// Build the conformer residues of every chain, in cctbx iteration order.
fn conformer_residues(st: &Structure, flat_index: &FxHashMap<AtomPath, u32>) -> Vec<Vec<Vec<ConfResidue>>> {
    // [model][chain-conformer][residue]
    let mut out = Vec::new();
    for (mi, m) in st.models.iter().enumerate() {
        let mut per_model = Vec::new();
        for (ci, c) in m.chains.iter().enumerate() {
            let alts = c.conformer_altlocs();
            for alt in &alts {
                let mut residues = Vec::new();
                for (ri, rg) in c.residue_groups.iter().enumerate() {
                    // atom groups with blank altloc or this altloc, grouped by resname
                    let mut by_res: Vec<(String, Vec<u32>)> = Vec::new();
                    for (gi, ag) in rg.atom_groups.iter().enumerate() {
                        if !(ag.altloc.is_empty() || ag.altloc == *alt) {
                            continue;
                        }
                        let entry = match by_res.iter_mut().position(|e| e.0 == ag.resname) {
                            Some(k) => k,
                            None => {
                                by_res.push((ag.resname.clone(), Vec::new()));
                                by_res.len() - 1
                            }
                        };
                        for k in 0..ag.atoms.len() {
                            let p = AtomPath { model: mi as u32, chain: ci as u32, rg: ri as u32, ag: gi as u32, atom: k as u32 };
                            by_res[entry].1.push(flat_index[&p]);
                        }
                    }
                    for (resname, atoms) in by_res {
                        residues.push(ConfResidue { atoms, resname, rg: (mi, ci, ri), link_to_previous: rg.link_to_previous });
                    }
                }
                per_model.push(residues);
            }
        }
        out.push(per_model);
    }
    out
}

/// Look up an atom id in a residue, retrying with '<->* swaps (bonds,
/// dihedrals, planes) like cctbx.
fn lookup(r: &ResInterp, id: &str) -> Option<u32> {
    if let Some(&a) = r.id_to_atom.get(id) {
        return Some(a);
    }
    if id.contains('\'') {
        if let Some(&a) = r.id_to_atom.get(&id.replace('\'', "*")) {
            return Some(a);
        }
    }
    if id.contains('*') {
        if let Some(&a) = r.id_to_atom.get(&id.replace('*', "'")) {
            return Some(a);
        }
    }
    None
}

fn v3_to_v2(id: &str) -> String {
    match id {
        "OP1" => "O1P".into(),
        "OP2" => "O2P".into(),
        _ => id.replace('\'', "*"),
    }
}
fn v2_to_v3(id: &str) -> String {
    match id {
        "O1P" => "OP1".into(),
        "O2P" => "OP2".into(),
        _ => id.replace('*', "'"),
    }
}

fn lookup_angle(r: &ResInterp, id: &str) -> Option<u32> {
    r.id_to_atom
        .get(id)
        .or_else(|| r.id_to_atom.get(&v3_to_v2(id)))
        .or_else(|| r.id_to_atom.get(&v2_to_v3(id)))
        .copied()
}

/// RNA sugar pucker analysis (`rna_sugar_pucker_analysis`): Some(true) for
/// 2'-endo (rna2p), Some(false) for 3'-endo, None when undecided.
fn rna_pucker(flat: &FlatAtoms, cur: &ResInterp, next_p: Option<u32>) -> Option<bool> {
    let g = |n: &str| lookup(cur, n).map(|a| flat.pos[a as usize]);
    let (c1, c2, c3, c4, o4, o3) = (g("C1'")?, g("C2'")?, g("C3'")?, g("C4'")?, g("O4'")?, g("O3'")?);
    let c5 = g("C5'");
    for (a, b) in [(c1, c2), (c2, c3), (c3, c4), (c4, o4), (o4, c1), (c3, o3)] {
        let d = a.dist(b);
        if !(1.2..=1.8).contains(&d) {
            return None;
        }
    }
    if let Some(c5) = c5 {
        let d = c4.dist(c5);
        if !(1.2..=1.8).contains(&d) {
            return None;
        }
    }
    // outbound atom from C1': closest N or C (no prime) within 1.463 + 0.5
    let mut out: Option<(f64, Vec3)> = None;
    for &a in cur.id_to_atom.values() {
        let nm = flat.name[a as usize].trim();
        let el = flat.element[a as usize].as_str();
        if nm.contains('\'') || nm.contains('*') || (el != "N" && el != "C") {
            continue;
        }
        let d = flat.pos[a as usize].dist(c1);
        if d <= 1.463 + 0.5 && out.map(|o| d < o.0).unwrap_or(true) {
            out = Some((d, flat.pos[a as usize]));
        }
    }
    let (_, outb) = out?;
    let line_dist = |p: Vec3| -> f64 {
        let u = (outb - c1).normalize();
        let v = p - c1;
        (v - u * v.dot(u)).length()
    };
    if let Some(p) = next_p {
        let pp = flat.pos[p as usize];
        if pp.dist(o3) <= 1.8 {
            return Some(line_dist(pp) < 2.9);
        }
    }
    Some(line_dist(o3) < 2.4)
}

/// Interpret one conformer residue: pick the dictionary, map names, apply
/// terminal and other modifications.
fn interpret_residue(
    ml: &MonLib,
    auto_comps: &FxHashMap<String, Arc<Comp>>,
    flat: &FlatAtoms,
    res: &ConfResidue,
    neutron_unused: bool,
    log: &mut String,
) -> ResInterp {
    let _ = neutron_unused;
    let resname = res.resname.trim().to_ascii_uppercase();
    let names: Vec<String> = res.atoms.iter().map(|&a| flat.name[a as usize].clone()).collect();
    let mut is_water = resclass::get_class(&res.resname) == ResClass::CommonWater;
    let mut work = resname.clone();
    let mut d_aa = false;
    if let Some(l) = names::l_given_d(&resname) {
        work = l.to_string();
        d_aa = true;
    }
    let ani = names::protein_mon_lib_names(&work, &names);
    let mut na_interp = false;
    if ani.is_none() {
        let has_o2 = names.iter().any(|n| {
            let t = n.trim().to_ascii_uppercase();
            t == "O2'" || t == "O2*" || t == "HO2'" || t == "2HO*" || t == "HO2*"
        });
        if let Some(w) = names::rna_dna_mon_lib_name(&resname, has_o2) {
            work = w.to_string();
            na_interp = true;
        }
    }
    let comp0 = ml.comp(&work).or_else(|| auto_comps.get(&resname).cloned());
    let Some(comp0) = comp0 else {
        return ResInterp {
            comp: None,
            id_to_atom: FxHashMap::default(),
            is_peptide: false,
            is_rna_dna: false,
            is_rna2p: None,
            is_water,
        };
    };
    let _ = d_aa;
    // shared dictionary; a modified copy is made only when a mod applies
    let mut comp: Arc<Comp> = comp0;
    let is_peptide_dict = comp.test_for_peptide();
    let is_na = na_interp || comp.test_for_rna_dna().is_some() || resclass::get_class(&res.resname).is_rna_dna();
    if comp.atoms.len() == 3 && comp.has_atom("O") && comp.has_atom("H1") && comp.has_atom("H2") {
        is_water = true;
    }
    let syn = ml.atom_synonyms.get(&comp.id);
    let mut mapping = names::map_atoms(&comp, &names, ani.as_deref(), is_na, syn);
    // ---------------- resolve_unexpected (terminal modifications)
    // `_get_incomplete_info`: no unexpected atoms, some heavy atoms missing,
    // and the atoms present are exactly a backbone fragment (or a lone P)
    let skip_resolve = {
        let mut present: Vec<&str> = mapping.ids.iter().flatten().map(|s| s.as_str()).collect();
        present.sort_unstable();
        present.dedup();
        let missing_heavy = comp.atoms.iter().any(|a| !comp.is_h(&a.id) && !present.contains(&a.id.as_str()));
        if !mapping.unexpected.is_empty() || !missing_heavy {
            false
        } else if is_peptide_dict {
            let ids = present.join(" ");
            matches!(ids.as_str(), "CA" | "C CA N" | "C CA N O" | "C CA CB N O")
        } else if is_na {
            present == ["P"]
        } else {
            false
        }
    };
    let apply = |comp: &mut Arc<Comp>, mod_id: &str, mapping: &mut names::Mapping, log: &mut String| {
        if let Some(c) = ml.apply_mod(comp, mod_id) {
            *comp = Arc::new(c);
            *mapping = names::map_atoms(comp, &names, ani.as_deref(), is_na, syn);
            log.push_str(&format!("  mod {} on {}\n", mod_id, resname));
        }
    };
    if !skip_resolve {
        let u: FxHashSet<String> = mapping.unexpected.iter().cloned().collect();
        if is_peptide_dict {
            let n_term = ["H1", "H2", "H3"].iter().filter(|n| u.contains(**n)).count();
            if n_term == 3 && comp.has_atom("H") {
                apply(&mut comp, "NH3", &mut mapping, log);
            } else if n_term == 2 {
                if comp.id == "PRO" {
                    apply(&mut comp, "NH2", &mut mapping, log);
                } else {
                    apply(&mut comp, "NH2NOTPRO", &mut mapping, log);
                }
                rename_terminal(&mut mapping, &names, ani.as_deref(), &["HN1", "HN2"]);
            } else if n_term == 1 {
                if comp.id == "PRO" {
                    apply(&mut comp, "NH1", &mut mapping, log);
                } else {
                    apply(&mut comp, "NH1NOTPRO", &mut mapping, log);
                }
                rename_terminal(&mut mapping, &names, ani.as_deref(), &["HN"]);
            }
            let u: FxHashSet<String> = mapping.unexpected.iter().cloned().collect();
            if u.contains("HXT") {
                apply(&mut comp, "COOH", &mut mapping, log);
            } else if u.contains("OXT") {
                apply(&mut comp, "COO", &mut mapping, log);
            }
            if comp.id == "GLU" && u.contains("HE2") {
                apply(&mut comp, "ACID-GLU", &mut mapping, log);
            }
            if comp.id == "ASP" && u.contains("HD2") {
                apply(&mut comp, "ACID-ASP", &mut mapping, log);
            }
        } else if is_na {
            let has = |n: &str| names.iter().any(|x| names::na_reference_name(x) == n);
            if has("OP3") || has("HOP3") {
                apply(&mut comp, "p5*END", &mut mapping, log);
            } else if !(has("P") && (has("OP1") || has("OP2"))) && !has("P") {
                apply(&mut comp, "5*END", &mut mapping, log);
            }
            if has("HO3'") {
                apply(&mut comp, "3*END", &mut mapping, log);
            }
        }
    }
    let mut id_to_atom: FxHashMap<String, u32> = FxHashMap::default();
    for (k, id) in mapping.ids.iter().enumerate() {
        if let Some(id) = id {
            id_to_atom.insert(id.clone(), res.atoms[k]);
        }
    }
    let is_rna_dict = comp.test_for_rna_dna();
    ResInterp {
        comp: Some(comp),
        id_to_atom,
        is_peptide: is_peptide_dict && !is_na,
        is_rna_dna: is_na || is_rna_dict.is_some(),
        is_rna2p: None,
        is_water,
    }
}

/// After NH1/NH2 mods, the model's unexpected H1/H2/H3 map to HN/HN1/HN2 in order.
fn rename_terminal(mapping: &mut names::Mapping, names: &[String], ani: Option<&[Option<String>]>, targets: &[&str]) {
    let mut t = targets.iter();
    for (k, n) in names.iter().enumerate() {
        let ml = ani.and_then(|a| a[k].clone()).unwrap_or_else(|| n.trim().to_string());
        if ["H1", "H2", "H3"].contains(&ml.as_str()) && mapping.ids[k].is_none() {
            if let Some(&tn) = t.next() {
                mapping.ids[k] = Some(tn.to_string());
                mapping.unexpected.retain(|x| x != &ml);
            }
        }
    }
}

/// `get_lib_link` for two consecutive residues.
fn lib_link_id(ml: &MonLib, prev: &ResInterp, cur: &ResInterp) -> Option<String> {
    if prev.is_water || cur.is_water {
        return None;
    }
    let (pc, cc) = (prev.comp.as_ref()?, cur.comp.as_ref()?);
    if prev.is_peptide && cur.is_peptide {
        if cur.id_to_atom.contains_key("CN") {
            return Some("NMTRANS".into());
        }
        if cc.id == "PRO" {
            return Some("PTRANS".into());
        }
        return Some("TRANS".into());
    }
    if prev.is_rna_dna && cur.is_rna_dna {
        return Some(if prev.is_rna2p == Some(true) { "rna2p".into() } else { "rna3p".into() });
    }
    // generic link search over the library list (simplified best match)
    let norm = |g: &str| -> String {
        match g.to_ascii_lowercase().as_str() {
            "l-peptide" | "d-peptide" | "peptide" => "peptide".into(),
            "dna" | "rna" => "DNA/RNA".into(),
            other => other.into(),
        }
    };
    let pg = norm(&pc.group);
    let cg = norm(&cc.group);
    let mut best: Option<(usize, usize, usize, String)> = None;
    for l in &ml.links {
        if l.name.contains("SS-bridge") || l.bonds.is_empty() {
            continue;
        }
        if l.comp_id_1.is_empty() && l.comp_id_2.is_empty() && l.group_comp_1.is_empty() && l.group_comp_2.is_empty() {
            continue;
        }
        let side_ok = |cid: &str, grp: &str, comp: &Comp, g: &str| -> bool {
            let comp_ok = cid.is_empty() || cid.eq_ignore_ascii_case(&comp.id);
            let grp_ok = grp.is_empty() || norm(grp) == g;
            comp_ok && grp_ok
        };
        if !side_ok(&l.comp_id_1, &l.group_comp_1, pc, &pg) || !side_ok(&l.comp_id_2, &l.group_comp_2, cc, &cg) {
            continue;
        }
        let specificity = l.comp_id_1.len() + l.comp_id_2.len();
        let unresolved = l
            .bonds
            .iter()
            .filter(|b| {
                let r1 = if b.c1 == 1 { pc } else { cc };
                let r2 = if b.c2 == 1 { pc } else { cc };
                !r1.has_atom(&b.a1) || !r2.has_atom(&b.a2)
            })
            .count();
        let key = (unresolved, usize::MAX - specificity, 0, l.id.clone());
        if best.as_ref().map(|b| (key.0, key.1) < (b.0, b.1)).unwrap_or(true) {
            best = Some(key);
        }
    }
    let b = best?;
    let l = ml.link(&b.3)?;
    if l.comp_id_1.is_empty() && l.comp_id_2.is_empty() && l.group_comp_1.is_empty() && l.group_comp_2.is_empty() {
        return None;
    }
    Some(b.3)
}

/// Add the proxies of a chain link between `prev` and `cur`. Returns false
/// when a link bond is broken (chain break).
fn add_link(it: &mut Interp, ml: &MonLib, flat: &FlatAtoms, link_id: &str, prev: &ResInterp, cur: &ResInterp, cutoff: f64) -> bool {
    let Some(l) = ml.link(link_id) else { return false };
    let side = |c: u8| if c == 1 { prev } else { cur };
    let mut broken: Vec<(u32, u32)> = Vec::new();
    let mut any = false;
    for b in &l.bonds {
        let (Some(i), Some(j)) = (lookup(side(b.c1), &b.a1), lookup(side(b.c2), &b.a2)) else { continue };
        let Some(v) = b.value_dist else { continue };
        if b.esd.unwrap_or(0.0) <= 0.0 {
            continue;
        }
        if flat.pos[i as usize].dist(flat.pos[j as usize]) > cutoff {
            broken.push((i.min(j), i.max(j)));
            continue;
        }
        it.add_bond(i, j, v, ORIGIN_COVALENT);
        any = true;
    }
    let involves_broken = |ids: &[u32]| broken.iter().any(|(a, b)| ids.contains(a) && ids.contains(b));
    for a in &l.angles {
        let ids: Vec<Option<u32>> = (0..3).map(|k| lookup_angle(side(a.c[k]), &a.a[k])).collect();
        if ids.iter().any(|x| x.is_none()) {
            continue;
        }
        let ids: Vec<u32> = ids.into_iter().flatten().collect();
        let (Some(v), Some(e)) = (a.value, a.esd) else { continue };
        if e <= 0.0 || involves_broken(&ids) {
            continue;
        }
        it.add_angle(ids[0], ids[1], ids[2], v, e);
    }
    for t in &l.tors {
        if t.id == "psi" || t.id == "phi" {
            continue;
        }
        let ids: Vec<Option<u32>> = (0..4).map(|k| lookup(side(t.c[k]), &t.a[k])).collect();
        if ids.iter().any(|x| x.is_none()) {
            continue;
        }
        let ids: [u32; 4] = [ids[0].unwrap(), ids[1].unwrap(), ids[2].unwrap(), ids[3].unwrap()];
        if involves_broken(&ids) {
            continue;
        }
        let Some(mut v) = t.value else { continue };
        if t.id == "omega" && matches!(link_id, "TRANS" | "PTRANS" | "NMTRANS") {
            if let Some(m) = dihedral_model(&flat.pos, ids) {
                if angle_delta_deg(m, 180.0, 1).abs() > 180.0 - 45.0 {
                    v = 0.0;
                }
            }
        }
        match t.esd {
            Some(e) if e > 0.0 => it.add_dihedral(ids, v, t.period, false),
            _ => it.add_dihedral(ids, v, 0, true),
        }
    }
    let mut planes: Vec<(String, Vec<u32>)> = Vec::new();
    for p in &l.planes {
        let Some(a) = lookup(side(p.c), &p.atom) else { continue };
        match planes.iter_mut().find(|x| x.0 == p.plane_id) {
            Some(x) => x.1.push(a),
            None => planes.push((p.plane_id.clone(), vec![a])),
        }
    }
    for (_, atoms) in planes {
        if !involves_broken(&atoms) {
            it.add_plane(atoms);
        }
    }
    any || broken.is_empty()
}

fn add_residue_restraints(it: &mut Interp, flat: &FlatAtoms, r: &ResInterp, neutron: bool) {
    let Some(comp) = r.comp.as_ref() else { return };
    let _ = flat;
    for b in &comp.bonds {
        let (Some(i), Some(j)) = (lookup(r, &b.a1), lookup(r, &b.a2)) else { continue };
        let v = if neutron { b.value_dist_neutron.or(b.value_dist) } else { b.value_dist };
        let Some(v) = v else { continue };
        match b.esd {
            Some(e) if e > 0.0 => {}
            _ => continue,
        }
        it.add_bond(i, j, v, ORIGIN_COVALENT);
    }
    for a in &comp.angles {
        let (Some(i), Some(j), Some(k)) = (lookup_angle(r, &a.a1), lookup_angle(r, &a.a2), lookup_angle(r, &a.a3)) else { continue };
        let (Some(v), Some(e)) = (a.value, a.esd) else { continue };
        if e <= 0.0 {
            continue;
        }
        it.add_angle(i, j, k, v, e);
    }
    for t in &comp.tors {
        let ids: Vec<Option<u32>> = t.a.iter().map(|x| lookup(r, x)).collect();
        if ids.iter().any(|x| x.is_none()) {
            continue;
        }
        let ids = [ids[0].unwrap(), ids[1].unwrap(), ids[2].unwrap(), ids[3].unwrap()];
        let Some(v) = t.value else { continue };
        match t.esd {
            Some(e) if e > 0.0 => it.add_dihedral(ids, v, t.period, false),
            _ => it.add_dihedral(ids, v, 0, true),
        }
    }
    let mut planes: Vec<(String, Vec<u32>)> = Vec::new();
    for p in &comp.planes {
        let Some(a) = lookup(r, &p.atom) else { continue };
        match planes.iter_mut().find(|x| x.0 == p.plane_id) {
            Some(x) => x.1.push(a),
            None => planes.push((p.plane_id.clone(), vec![a])),
        }
    }
    for (_, atoms) in planes {
        it.add_plane(atoms);
    }
}

/// Run the interpretation over the whole structure.
pub fn interpret(st: &Structure, flat: &FlatAtoms, ml: &MonLib, p: &InterpParams) -> Interp {
    let n = flat.pos.len();
    let mut it = Interp::new(n);
    let mut flat_index: FxHashMap<AtomPath, u32> = FxHashMap::default();
    for (k, path) in flat.path.iter().enumerate() {
        flat_index.insert(*path, k as u32);
    }
    let confs = conformer_residues(st, &flat_index);
    let mut log = String::new();
    let mut cys_sg: Vec<u32> = Vec::new();
    let mut cys_sg_seen: FxHashSet<u32> = FxHashSet::default();
    for model_confs in &confs {
        for residues in model_confs {
            // interpret all residues of this chain conformer first (pucker needs next P)
            let mut interps: Vec<ResInterp> =
                residues.iter().map(|r| interpret_residue(ml, &p.auto_comps, flat, r, p.neutron, &mut log)).collect();
            for k in 0..interps.len() {
                if interps[k].is_rna_dna && interps[k].comp.is_some() && interps[k].id_to_atom.contains_key("O2'") {
                    let next_p = interps.get(k + 1).and_then(|nx| nx.id_to_atom.get("P").copied());
                    let pk = rna_pucker(flat, &interps[k], next_p);
                    interps[k].is_rna2p = pk;
                    if let Some(two) = pk {
                        let comp = interps[k].comp.as_ref().unwrap();
                        let base = residues[k].resname.trim().to_ascii_uppercase();
                        let kind = match base.as_str() {
                            "A" | "G" => "_pur",
                            "C" | "U" => "_pyr",
                            _ => "",
                        };
                        let mod_id = format!("{}{}", if two { "rna2p" } else { "rna3p" }, kind);
                        if let Some(c) = ml.apply_mod(comp, &mod_id) {
                            interps[k].comp = Some(Arc::new(c));
                        }
                    }
                }
            }
            let mut prev: Option<usize> = None;
            for (k, r) in residues.iter().enumerate() {
                let ri = &interps[k];
                if ri.comp.is_none() {
                    let rn = r.resname.trim().to_string();
                    if !it.missing_residues.contains(&rn) {
                        it.missing_residues.push(rn);
                    }
                    for &a in &r.atoms {
                        it.etype[a as usize] = EType::Unknown;
                    }
                    prev = None;
                    continue;
                }
                if let Some(pk) = prev {
                    if r.link_to_previous {
                        if let Some(lid) = lib_link_id(ml, &interps[pk], ri) {
                            add_link(&mut it, ml, flat, &lid, &interps[pk], ri, p.link_distance_cutoff);
                        }
                    }
                }
                // energy types
                let comp = ri.comp.as_ref().unwrap();
                let mapped: FxHashMap<u32, &String> = ri.id_to_atom.iter().map(|(id, &a)| (a, id)).collect();
                for &a in &r.atoms {
                    match mapped.get(&a) {
                        Some(id) => {
                            let ca = comp.atom(id).unwrap();
                            it.etype[a as usize] = match &ca.type_energy {
                                Some(t) => EType::Typed(t.clone()),
                                // Reduce2's CCD-built dictionaries leave the type empty
                                None if comp.from_ccd => EType::Unknown,
                                None => EType::NoneType,
                            };
                            it.dict[a as usize] = AtomDictInfo { comp: Some(comp.clone()), dict_id: Some((*id).clone()) };
                        }
                        None => {
                            it.etype[a as usize] = EType::Unexpected;
                            it.dict[a as usize] = AtomDictInfo { comp: Some(comp.clone()), dict_id: None };
                        }
                    }
                }
                add_residue_restraints(&mut it, flat, ri, p.neutron);
                if comp.id == "CYS" {
                    if let Some(&sg) = ri.id_to_atom.get("SG") {
                        if cys_sg_seen.insert(sg) {
                            cys_sg.push(sg);
                        }
                    }
                }
                prev = Some(k);
            }
        }
    }
    it.log = log;
    add_disulfides(&mut it, st, flat, &cys_sg);
    let mut link_log = String::new();
    crate::autolink::auto_link(
        &mut it,
        st,
        flat,
        ml,
        &crate::autolink::AutoLinkParams { compat: p.compat },
        &mut link_log,
    );
    it.log.push_str(&link_log);
    add_iron_sulfur_coordination(&mut it, st, flat, p.compat);
    add_zinc_coordination(&mut it, st, flat);
    adjust_for_ph_variants(&mut it, st, flat, &flat_index, ml, p);
    it
}

/// cctbx `pH_dependent_restraints.adjust_geometry_proxies_registeries`, which
/// pdb_interpretation runs when hydrogens have no energy type (residues with
/// CCD-built dictionaries): bonds and angles of such a residue take the values
/// of the GeoStd neutron or low-pH variant of its dictionary, where there is
/// one, and missing ones are added. Reduce2 adjusts only the first residue of
/// each name; fixed mode adjusts every residue, and in neutron mode uses the
/// neutron distance for the bonds it changes as well as for those it adds.
fn adjust_for_ph_variants(
    it: &mut Interp,
    st: &Structure,
    flat: &FlatAtoms,
    flat_index: &FxHashMap<AtomPath, u32>,
    ml: &MonLib,
    p: &InterpParams,
) {
    let unknown: Vec<u32> = (0..flat.pos.len() as u32).filter(|&k| it.etype[k as usize] == EType::Unknown).collect();
    if !unknown.iter().any(|&k| flat.is_h[k as usize]) {
        return;
    }
    let mut checked: Vec<String> = Vec::new();
    let mut done: FxHashSet<(u32, u32, u32)> = FxHashSet::default();
    // atoms of added bonds -> their dictionary energy type (None: not in the dictionary)
    let mut added_energy: Vec<(u32, Option<Option<String>>)> = Vec::new();
    for &k in &unknown {
        let pth = flat.path[k as usize];
        let resname = st.atom_group(pth).resname.clone();
        if p.compat {
            if checked.contains(&resname) {
                continue;
            }
            checked.push(resname.clone());
        } else if !done.insert((pth.model, pth.chain, pth.rg)) {
            continue;
        }
        let Some(v) = ml.ph_variant(&resname) else { continue };
        // the residue group's atoms by atom group: (stripped name, flat index)
        let rg = st.residue_group(pth);
        let groups: Vec<Vec<(String, u32)>> = rg
            .atom_groups
            .iter()
            .enumerate()
            .map(|(g, ag)| {
                ag.atoms
                    .iter()
                    .enumerate()
                    .filter_map(|(i, a)| {
                        let path = AtomPath { model: pth.model, chain: pth.chain, rg: pth.rg, ag: g as u32, atom: i as u32 };
                        flat_index.get(&path).map(|&f| (a.name.trim().to_string(), f))
                    })
                    .collect()
            })
            .collect();
        // `_get_atom_neutron`: the atom by name, or its D form for bonds under 1 A
        let atom = |g: usize, name: &str, bondlength: Option<f64>| -> Option<(u32, bool)> {
            let find = |n: &str| groups[g].iter().find(|(x, _)| x == n).map(|x| x.1);
            if let Some(a) = find(name) {
                return Some((a, false));
            }
            if bondlength.is_some_and(|b| b < 1.0) {
                return find(&name.replace('H', "D")).map(|a| (a, true));
            }
            None
        };
        let ng = groups.len();
        for b in &v.bonds {
            let Some(dist) = b.value_dist else { continue };
            for g1 in 0..ng {
                for g2 in 0..=g1 {
                    let Some((a1, n1)) = atom(g1, &b.a1, b.value_dist) else { continue };
                    let Some((a2, n2)) = atom(g2, &b.a2, b.value_dist) else { continue };
                    let neutron_value = b.value_dist_neutron.unwrap_or(dist);
                    if it.has_bond(a1, a2) {
                        let ideal = if !p.compat && p.neutron { neutron_value } else { dist };
                        it.upsert_bond(a1, a2, ideal, ORIGIN_COVALENT);
                    } else if a1 != a2 {
                        it.add_bond(a1, a2, if n1 || n2 { neutron_value } else { dist }, ORIGIN_COVALENT);
                        for (a, name) in [(a1, &b.a1), (a2, &b.a2)] {
                            added_energy.retain(|x| x.0 != a);
                            added_energy.push((a, v.atom(name).map(|c| c.type_energy.clone())));
                        }
                    }
                }
            }
        }
        // angles: atoms -> (value, esd), both directions, in first-insertion order
        let mut lookup: Vec<([u32; 3], (f64, f64))> = Vec::new();
        let mut index: FxHashMap<[u32; 3], usize> = FxHashMap::default();
        for a in &v.angles {
            let (Some(value), Some(esd)) = (a.value, a.esd) else { continue };
            for g1 in 0..ng {
                for g2 in 0..=g1 {
                    let Some((i1, _)) = atom(g1, &a.a1, Some(0.5)) else { continue };
                    let Some((i3, _)) = atom(g2, &a.a3, Some(0.5)) else { continue };
                    let Some((i2, _)) = [g1, g2].iter().find_map(|&g| atom(g, &a.a2, None)) else { continue };
                    for key in [[i1, i2, i3], [i3, i2, i1]] {
                        match index.get(&key) {
                            Some(&x) => lookup[x].1 = (value, esd),
                            None => {
                                index.insert(key, lookup.len());
                                lookup.push((key, (value, esd)));
                            }
                        }
                    }
                }
            }
        }
        let mut used: FxHashSet<[u32; 3]> = FxHashSet::default();
        for proxy in it.angles.iter_mut() {
            if used.contains(&proxy.i) {
                continue;
            }
            if let Some(&x) = index.get(&proxy.i) {
                proxy.ideal = lookup[x].1 .0;
                proxy.esd = lookup[x].1 .1;
                used.insert(proxy.i);
                used.insert([proxy.i[2], proxy.i[1], proxy.i[0]]);
            }
        }
        let mut reversed_done: FxHashSet<[u32; 3]> = FxHashSet::default();
        for (key, (value, esd)) in lookup {
            if used.contains(&key) || reversed_done.contains(&key) {
                continue;
            }
            it.add_angle(key[0], key[1], key[2], value, esd);
            reversed_done.insert([key[2], key[1], key[0]]);
        }
    }
    let unknown_set: FxHashSet<u32> = unknown.into_iter().collect();
    for (a, te) in added_energy {
        // an atom missing from the variant dictionary stops cctbx (AttributeError)
        let Some(te) = te else { continue };
        if unknown_set.contains(&a) {
            it.etype[a as usize] = match te {
                Some(t) => EType::Typed(t),
                None => EType::NoneType,
            };
        }
    }
}


fn add_disulfides(it: &mut Interp, st: &Structure, flat: &FlatAtoms, sgs: &[u32]) {
    // exclude SG near atoms that are not H D T S O P N C SE (typically metals)
    let ok_el = ["H", "D", "T", "S", "O", "P", "N", "C", "SE"];
    let n = flat.pos.len();
    let others: Vec<usize> = (0..n).filter(|&a| !ok_el.contains(&flat.element[a].as_str())).collect();
    let excluded: FxHashSet<u32> = sgs
        .iter()
        .copied()
        .filter(|&sg| {
            let p = flat.pos[sg as usize];
            others.iter().any(|&a| flat.pos[a].dist(p) <= 3.0 && flat.path[a].model == flat.path[sg as usize].model)
        })
        .collect();
    // candidate pairs within 3 A through a coarse grid (pairs in input order)
    let cell = |p: Vec3| ((p.x / 3.0).floor() as i64, (p.y / 3.0).floor() as i64, (p.z / 3.0).floor() as i64);
    let mut grid: FxHashMap<(i64, i64, i64), Vec<usize>> = FxHashMap::default();
    for (k, &a) in sgs.iter().enumerate() {
        grid.entry(cell(flat.pos[a as usize])).or_default().push(k);
    }
    for (x, &a) in sgs.iter().enumerate() {
        if excluded.contains(&a) {
            continue;
        }
        let pa = flat.pos[a as usize];
        let (cx, cy, cz) = cell(pa);
        let mut partners: Vec<usize> = Vec::new();
        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    if let Some(v) = grid.get(&(cx + dx, cy + dy, cz + dz)) {
                        partners.extend(v.iter().copied().filter(|&y| y > x));
                    }
                }
            }
        }
        partners.sort_unstable();
        for y in partners {
            let b = sgs[y];
            if excluded.contains(&b) {
                continue;
            }
            if flat.path[a as usize].model != flat.path[b as usize].model {
                continue;
            }
            let (aa, ab) = (&flat.altloc[a as usize], &flat.altloc[b as usize]);
            if !aa.is_empty() && !ab.is_empty() && aa != ab {
                continue;
            }
            if pa.dist(flat.pos[b as usize]) <= 3.0 {
                it.add_bond(a, b, 2.031, ORIGIN_SS);
            }
        }
    }    // disulfides to symmetry copies (`create_disulfides` over the asu mappings)
    let Some((uc, ops)) = crate::cell::crystal_operators(st) else { return };
    let frac: Vec<Vec3> = sgs.iter().map(|&a| uc.fractionalize(flat.pos[a as usize])).collect();
    for x in 0..sgs.len() {
        let a = sgs[x];
        if excluded.contains(&a) {
            continue;
        }
        for y in x..sgs.len() {
            let b = sgs[y];
            if excluded.contains(&b) || flat.path[a as usize].model != flat.path[b as usize].model {
                continue;
            }
            let (aa, ab) = (&flat.altloc[a as usize], &flat.altloc[b as usize]);
            if !aa.is_empty() && !ab.is_empty() && aa != ab {
                continue;
            }
            for op in &ops {
                let img = op.apply(frac[y]);
                let d = frac[x] - img;
                let base = v3(d.x.round(), d.y.round(), d.z.round());
                let mut hit = false;
                for dx in -1..=1 {
                    for dy in -1..=1 {
                        for dz in -1..=1 {
                            let shift = base + v3(dx as f64, dy as f64, dz as f64);
                            if op.is_identity() && shift.x == 0.0 && shift.y == 0.0 && shift.z == 0.0 {
                                continue; // the same copy: a simple disulfide
                            }
                            if uc.orthogonalize(img + shift).dist(flat.pos[a as usize]) <= 3.0 {
                                hit = true;
                            }
                        }
                    }
                }
                if hit {
                    it.sym_bonds.push(BondProxy { i: a, j: b, ideal: 2.031, origin: ORIGIN_SS });
                }
            }
        }
    }
}


/// Metal Coordination Library, iron-sulfur cluster part
/// (`mcl_sf4_coordination.get_sulfur_iron_cluster_coordination`): over the
/// non-bonded pairs within 3.5 A between an SF4/F3S/FES residue and another
/// residue, closest first, each other residue's first S or N bonds to the
/// cluster Fe. cctbx keys the residues by atom-group id, without the model;
/// fixed mode tells models apart.
fn add_iron_sulfur_coordination(it: &mut Interp, st: &Structure, flat: &FlatAtoms, compat: bool) {
    let cluster = |a: usize| matches!(st.atom_group(flat.path[a]).resname.trim(), "SF4" | "F3S" | "FES");
    let n = flat.pos.len();
    let in_cluster: Vec<usize> = (0..n).filter(|&a| cluster(a)).collect();
    if in_cluster.is_empty() {
        return;
    }
    let mut nbrs: Vec<Vec<u32>> = vec![Vec::new(); n];
    for b in &it.bonds {
        nbrs[b.i as usize].push(b.j);
        nbrs[b.j as usize].push(b.i);
    }
    // 1-2 and 1-3 pairs are not non-bonded pairs
    let excluded = |a: usize, b: usize| -> bool {
        nbrs[a].contains(&(b as u32)) || nbrs[a].iter().any(|&m| nbrs[m as usize].contains(&(b as u32)))
    };
    let mut pairs: Vec<(f64, usize, usize)> = Vec::new();
    for &c in &in_cluster {
        for o in 0..n {
            if cluster(o) || flat.path[o].model != flat.path[c].model {
                continue;
            }
            let (ac, ao) = (&flat.altloc[c], &flat.altloc[o]);
            if !ac.is_empty() && !ao.is_empty() && ac != ao {
                continue;
            }
            let d = flat.pos[c].dist(flat.pos[o]);
            if d < 3.5 && !excluded(c, o) {
                pairs.push((d, c.min(o), c.max(o)));
            }
        }
    }
    pairs.sort_by(|x, y| x.0.partial_cmp(&y.0).unwrap().then(x.1.cmp(&y.1)).then(x.2.cmp(&y.2)));
    let mut done: Vec<(u32, String, String, String, String, String)> = Vec::new();
    for (_, i, j) in pairs {
        let (fe, aa) = if cluster(i) { (i, j) } else { (j, i) };
        if !matches!(flat.element[aa].as_str(), "S" | "N") || flat.element[fe] != "FE" {
            continue;
        }
        let pth = flat.path[aa];
        let ag = st.atom_group(pth);
        let rg = st.residue_group(pth);
        let key = (
            if compat { 0 } else { pth.model },
            st.chain(pth).id.clone(),
            rg.resseq.clone(),
            rg.icode.clone(),
            ag.altloc.clone(),
            ag.resname.clone(),
        );
        if done.contains(&key) {
            continue;
        }
        done.push(key);
        // the S-Fe-S angles cctbx adds as well do not involve hydrogens
        it.add_bond(fe as u32, aa as u32, 2.268, ORIGIN_METAL);
    }
}

/// Metal Coordination Library, zinc part: ZN with SG/ND1/NE2 partners within 3 A.
fn add_zinc_coordination(it: &mut Interp, st: &Structure, flat: &FlatAtoms) {
    let n = flat.pos.len();
    let cands: Vec<usize> = (0..n).filter(|&a| matches!(flat.name[a].trim(), "SG" | "ND1" | "NE2")).collect();
    for z in 0..n {
        if flat.name[z].trim() != "ZN" {
            continue;
        }
        let pz = flat.pos[z];
        let mut per_res: FxHashSet<(u32, u32, u32)> = FxHashSet::default();
        let mut partners: Vec<(f64, u32)> = Vec::new();
        for &a in &cands {
            if a == z || flat.path[a].model != flat.path[z].model {
                continue;
            }
            let d = flat.pos[a].dist(pz);
            if d <= 3.0 {
                partners.push((d, a as u32));
            }
        }
        partners.sort_by(|x, y| x.0.partial_cmp(&y.0).unwrap());
        for (_, a) in partners {
            let pth = flat.path[a as usize];
            let key = (pth.model, pth.chain, pth.rg);
            if per_res.insert(key) {
                it.add_bond(z as u32, a, 2.30, ORIGIN_METAL);
            }
        }
    }
    let _ = st;
}
