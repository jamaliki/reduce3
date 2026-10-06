//! Hydrogen placement (port of mmtbx.hydrogens.reduce_hydrogen.place_hydrogens).

use crate::geom::*;
use crate::interp::{self, EType, FlatAtoms, Interp, InterpParams};
use crate::model::*;
use crate::monlib::{Comp, MonLib};
use crate::names;
use crate::resclass::{self, ResClass};
use crate::riding::{self, RidingAtoms, RidingCoef};
use rustc_hash::{FxHashMap, FxHashSet};
use std::sync::Arc;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum NTermCharge {
    ResidueOne,
    FirstInChain,
    NoCharge,
}

#[derive(Clone, Debug)]
pub struct HPlaceParams {
    pub neutron: bool,
    pub n_terminal_charge: NTermCharge,
    pub exclude_water: bool,
    pub keep_existing_h: bool,
    pub adp_scale: f64,
    pub compat: bool,
    /// The cell cctbx processes the model with (compat mode repeats its
    /// coordinate round-off; see `cell`).
    pub cell: Option<crate::cell::UnitCell>,
}

/// Result of hydrogen placement: the updated structure plus everything the
/// optimizer needs, indexed by the final flat atom order.
pub struct Placed {
    pub bonds: Vec<(u32, u32, u16)>,
    pub etype: Vec<EType>,
    pub riding: Vec<Option<RidingCoef>>,
    pub dict: Vec<interp::AtomDictInfo>,
    pub no_h_placed: Vec<String>,
    pub site_labels_no_para: Vec<String>,
    /// H whose third neighbor lies on the parent bond axis (compat mode,
    /// where Reduce2 divides by zero on them).
    pub axial_reference: Vec<String>,
    pub removed_on_links: Vec<(String, String)>,
    pub n_h_initial: usize,
    pub n_h_final: usize,
    pub log: String,
}

fn label(st: &Structure, p: AtomPath) -> String {
    let a = st.atom(p);
    let ag = st.atom_group(p);
    let rg = st.residue_group(p);
    let ch = st.chain(p);
    format!("{}{}{:>3}{:>2}{:>4}{}", a.name, if ag.altloc.is_empty() { " " } else { &ag.altloc }, ag.resname, ch.id, rg.resseq, rg.icode)
}

/// `mon_lib_query`: the residue dictionary (unmodified) for an atom group.
pub fn residue_dictionary(ml: &MonLib, resname: &str, names: &[String]) -> Option<Arc<Comp>> {
    let r = resname.trim().to_ascii_uppercase();
    if r == "UNL" || r.is_empty() {
        return None;
    }
    let mut work = r.clone();
    if let Some(l) = names::l_given_d(&r) {
        work = l.to_string();
    }
    if !names::has_protein_interpreter(&work) {
        let has_o2 = names.iter().any(|n| {
            let t = n.trim().to_ascii_uppercase();
            t == "O2'" || t == "O2*" || t == "HO2'" || t == "2HO*"
        });
        if let Some(w) = names::rna_dna_mon_lib_name(&r, has_o2) {
            work = w.to_string();
        }
    }
    ml.comp(&work)
}

/// `_terminal_h`: free-form-only H of a peptide-like dictionary.
fn terminal_h(comp: &Comp) -> Vec<String> {
    let mut on_n: Vec<String> = Vec::new();
    let mut on_oxt: Vec<String> = Vec::new();
    for b in &comp.bonds {
        for (x, y) in [(&b.a1, &b.a2), (&b.a2, &b.a1)] {
            let Some(a) = comp.atom(y) else { continue };
            if a.type_symbol != "H" && a.type_symbol != "D" {
                continue;
            }
            if x == "N" {
                on_n.push(y.clone());
            } else if x == "OXT" {
                on_oxt.push(y.clone());
            }
        }
    }
    let mut remove = on_oxt;
    if on_n.len() > 1 {
        let order: Vec<&str> = comp.atoms.iter().map(|a| a.id.as_str()).collect();
        let mut sorted = on_n.clone();
        sorted.sort_by_key(|h| (!(h == "H" || h == "D"), order.iter().position(|o| o == h).unwrap_or(usize::MAX)));
        let keep = sorted[0].clone();
        remove.extend(on_n.into_iter().filter(|h| *h != keep));
    }
    remove
}

const ALT_NAMES: &[[&str; 3]] = &[
    ["HA1", "HA2", "HA3"],
    ["HB1", "HB2", "HB3"],
    ["HG1", "HG2", "HG3"],
    ["HD1", "HD2", "HD3"],
    ["HE1", "HE2", "HE3"],
    ["HG11", "HG12", "HG13"],
];

fn pad_name(n: &str) -> String {
    if n.len() < 4 { format!("{:<4}", format!(" {}", n)) } else { n.to_string() }
}

/// Add every missing dictionary H at a bogus position (mean of the atom
/// group + 0.5 on each axis).
fn add_missing_h(
    st: &mut Structure,
    ml: &MonLib,
    p: &HPlaceParams,
    no_h_placed: &mut Vec<String>,
    auto_comps: &mut FxHashMap<String, Arc<Comp>>,
) {
    for m in &mut st.models {
        for c in &mut m.chains {
            for rg in &mut c.residue_groups {
                let n_ag = rg.atom_groups.len();
                for ag in &mut rg.atom_groups {
                    if n_ag > 2 && ag.altloc.is_empty() {
                        continue;
                    }
                    let cls = resclass::get_class(&ag.resname);
                    if p.exclude_water && cls == ResClass::CommonWater {
                        continue;
                    }
                    let names: Vec<String> = ag.atoms.iter().map(|a| a.name.clone()).collect();
                    let actual: FxHashSet<String> = ag.atoms.iter().map(|a| a.name.trim().to_ascii_uppercase()).collect();
                    // `mon_lib_query`: the libraries, else restraints built from the CCD
                    let comp = match residue_dictionary(ml, &ag.resname, &names) {
                        Some(c) => c,
                        None => match ml.ccd_comp(&ag.resname, p.compat) {
                            Some(c) => {
                                auto_comps.entry(ag.resname.trim().to_ascii_uppercase()).or_insert_with(|| c.clone());
                                c
                            }
                            None => {
                                no_h_placed.push(ag.resname.clone());
                                continue;
                            }
                        },
                    };
                    let mut removed: Vec<String> = Vec::new();
                    if comp.test_for_peptide() {
                        removed = terminal_h(&comp);
                    }
                    let mut expected: Vec<String> = comp
                        .atoms
                        .iter()
                        .filter(|a| a.type_symbol == "H" && !removed.contains(&a.id))
                        .map(|a| a.id.clone())
                        .collect();
                    let rename_classes: &[ResClass] = if p.compat {
                        &[ResClass::CommonAminoAcid, ResClass::ModifiedAminoAcid, ResClass::DAminoAcid]
                    } else {
                        // Modified amino acids follow CCD names; renaming them
                        // made the H unexpected and it was silently dropped.
                        &[ResClass::CommonAminoAcid, ResClass::DAminoAcid]
                    };
                    if rename_classes.contains(&cls) {
                        for alt in ALT_NAMES {
                            if expected.iter().any(|e| e == alt[0]) && expected.iter().any(|e| e == alt[1]) {
                                let te = |n: &str| comp.atom(n).and_then(|a| a.type_energy.clone()).unwrap_or_default();
                                if te(alt[0]) == "HCH2" && te(alt[1]) == "HCH2" {
                                    expected.push(alt[2].to_string());
                                    expected.retain(|e| e != alt[0]);
                                }
                            }
                        }
                    }
                    let missing: Vec<String> = expected.into_iter().filter(|e| !actual.contains(e)).collect();
                    if missing.is_empty() {
                        continue;
                    }
                    let mut mean = Vec3::ZERO;
                    for a in &ag.atoms {
                        mean += a.xyz;
                    }
                    let nat = ag.atoms.len().max(1) as f64;
                    let bogus = mean / nat + v3(0.5, 0.5, 0.5);
                    let hetero = ag.atoms.first().map(|a| a.hetero).unwrap_or(false);
                    let segid = ag.atoms.first().map(|a| a.segid.clone()).unwrap_or_default();
                    for mh in missing {
                        let mut a = Atom::new(&pad_name(&mh), "H", bogus);
                        a.element = "H".into();
                        a.hetero = hetero;
                        a.segid = segid.clone();
                        a.occ = 0.0;
                        a.b = 0.0;
                        ag.atoms.push(a);
                    }
                }
            }
        }
    }
}

/// CCD heavy-atom bond partners of an atom (`bonds_in_restraints`).
fn ccd_heavy_partners(ml: &MonLib, resname: &str, atom: &str) -> Option<usize> {
    let e = ml.ccd(resname)?;
    let elem: FxHashMap<&str, &str> = e.atoms.iter().map(|a| (a.0.as_str(), a.1.as_str())).collect();
    let mut n = 0;
    for (a1, a2, _) in &e.bonds {
        let other = if a1 == atom {
            a2
        } else if a2 == atom {
            a1
        } else {
            continue;
        };
        let el = elem.get(other.as_str()).copied().unwrap_or("");
        if el != "H" && el != "D" {
            n += 1;
        }
    }
    Some(n)
}

fn construct_xyz(n: Vec3, bv: f64, ca: Vec3, av: f64, c: Vec3, dv: f64) -> [Vec3; 3] {
    let rcca = c - ca;
    let e0 = (n - ca).normalize();
    let e1 = (rcca - e0 * rcca.dot(e0)).normalize();
    let e2 = e0.cross(e1);
    let alpha = av.to_radians();
    let phi = dv.to_radians();
    let mut out = [Vec3::ZERO; 3];
    for k in 0..3 {
        let ang = phi + k as f64 * 2.0 * std::f64::consts::PI / 3.0;
        out[k] = n + ((e1 * ang.cos() + e2 * ang.sin()) * alpha.sin() - e0 * alpha.cos()) * bv;
    }
    out
}

fn place_n_terminal_propeller(st: &mut Structure, ml: &MonLib, p: &HPlaceParams) {
    if p.n_terminal_charge == NTermCharge::NoCharge {
        return;
    }
    for m in &mut st.models {
        for c in &mut m.chains {
            let Some(rg) = c.residue_groups.first_mut() else { continue };
            if p.n_terminal_charge == NTermCharge::ResidueOne && rg.resseq_as_int() != 1 {
                continue;
            }
            let Some(nag) = rg.atom_groups.iter().position(|ag| ag.get_atom("N").is_some()) else { continue };
            let resname = rg.atom_groups[nag].resname.trim().to_string();
            let heavies = if resname == "PRO" { 3 } else { 2 };
            let Some(nb) = ccd_heavy_partners(ml, &resname, "N") else { continue };
            if nb >= heavies {
                continue;
            }
            let cls = resclass::get_class(&rg.atom_groups[nag].resname);
            if cls.is_amino_acid() {
                for ag in &mut rg.atom_groups {
                    ag.atoms.retain(|a| a.name.trim() != "H");
                }
            }
            // altloc combinations with N, CA and C
            let mut by_alt: Vec<(String, Vec<(usize, usize)>)> = Vec::new();
            for (gi, ag) in rg.atom_groups.iter().enumerate() {
                for k in 0..ag.atoms.len() {
                    match by_alt.iter_mut().find(|x| x.0 == ag.altloc) {
                        Some(x) => x.1.push((gi, k)),
                        None => by_alt.push((ag.altloc.clone(), vec![(gi, k)])),
                    }
                }
            }
            if by_alt.len() > 1 {
                if let Some(bi) = by_alt.iter().position(|x| x.0.is_empty()) {
                    let blank = by_alt.remove(bi).1;
                    for x in by_alt.iter_mut() {
                        x.1.extend(blank.iter().copied());
                    }
                }
            }
            for (_, members) in by_alt {
                let find = |nm: &str| members.iter().copied().find(|&(g, k)| rg.atom_groups[g].atoms[k].name.trim() == nm);
                let (Some(np), Some(cap), Some(cp)) = (find("N"), find("CA"), find("C")) else { continue };
                let gi = np.0;
                let ag = &rg.atom_groups[gi];
                let n_atom = ag.atoms[np.1].clone();
                let ca = rg.atom_groups[cap.0].atoms[cap.1].xyz;
                let cc = rg.atom_groups[cp.0].atoms[cp.1].xyz;
                // proton element: D if the group carries only D
                let hs: Vec<&str> = ag.atoms.iter().filter(|a| a.is_hydrogen()).map(|a| a.elem()).collect();
                let only_d = !hs.is_empty() && hs.iter().all(|e| *e == "D");
                let pe = if only_d { "D" } else { "H" };
                let mut dihedral = 120.0;
                if let Some(h) = ag.get_atom(pe) {
                    if let Some(d) = dihedral_deg(h.xyz, n_atom.xyz, ca, cc) {
                        dihedral = d;
                    }
                }
                let rh3 = construct_xyz(n_atom.xyz, 1.0, ca, 109.5, cc, dihedral);
                let possible = if pe == "H" { vec!["H", "H1", "H2", "H3", "HT1", "HT2"] } else { vec!["D", "D1", "D2", "D3"] };
                let mut count = 0;
                for h in &possible {
                    if ag.get_atom(h).is_some() {
                        count += 1;
                    }
                    if ag.get_atom(&h.replace('H', "D")).is_some() {
                        count += 1;
                    }
                }
                if count >= 3 {
                    continue;
                }
                let hetero = ag.atoms.first().map(|a| a.hetero).unwrap_or(false);
                let is_pro = ag.resname.trim() == "PRO";
                let has_proton = ag.get_atom(pe).is_some();
                let mut j = 0usize;
                let mut new_atoms = Vec::new();
                for i in 0..3 {
                    let name = format!(" {}{} ", pe, i + 1);
                    if i == 0 && has_proton {
                        continue;
                    }
                    if i == 1 && has_proton {
                        let r = ag.get_atom(pe).unwrap();
                        if r.xyz.dist_sq(rh3[j]) < 0.5 {
                            j += 1;
                        }
                    }
                    if ag.get_atom(name.trim()).is_some() {
                        continue;
                    }
                    if is_pro && i == 0 {
                        continue;
                    }
                    let mut a = Atom::new(&name, pe, rh3[j]);
                    a.occ = n_atom.occ;
                    a.b = n_atom.b;
                    a.segid = "    ".into();
                    a.hetero = hetero;
                    new_atoms.push(a);
                    j += 1;
                    if j == 3 {
                        j = 0;
                    }
                }
                rg.atom_groups[gi].atoms.extend(new_atoms);
            }
        }
    }
}

/// `hierarchy.flip_symmetric_amino_acids()`. With `per_atom_group`, each atom
/// group is tested on its own geometry (the original flips every atom group
/// of a residue once any one of them needs it).
fn flip_symmetric_amino_acids(st: &mut Structure, per_atom_group: bool) {
    struct FlipData {
        dihedral: Option<[&'static str; 4]>,
        chiral: Option<[&'static str; 4]>,
        pairs: &'static [[&'static str; 2]],
    }
    const ARG: &[[&str; 2]] = &[["NH1", "NH2"], ["HH11", "HH21"], ["HH12", "HH22"], ["DH11", "DH21"], ["DH12", "DH22"]];
    const ASP: &[[&str; 2]] = &[["OD1", "OD2"]];
    const GLU: &[[&str; 2]] = &[["OE1", "OE2"]];
    const PHE: &[[&str; 2]] = &[["CD1", "CD2"], ["CE1", "CE2"], ["HD1", "HD2"], ["HE1", "HE2"], ["DD1", "DD2"], ["DE1", "DE2"]];
    const VAL: &[[&str; 2]] =
        &[["CG1", "CG2"], ["HG11", "HG21"], ["HG12", "HG22"], ["HG13", "HG23"], ["DG11", "DG21"], ["DG12", "DG22"], ["DG13", "DG23"]];
    const LEU: &[[&str; 2]] =
        &[["CD1", "CD2"], ["HD11", "HD21"], ["HD12", "HD22"], ["HD13", "HD23"], ["DD11", "DD21"], ["DD12", "DD22"], ["DD13", "DD23"]];
    let data = |rn: &str| -> Option<FlipData> {
        Some(match rn {
            "ARG" => FlipData { dihedral: Some(["CD", "NE", "CZ", "NH1"]), chiral: None, pairs: ARG },
            "ASP" => FlipData { dihedral: Some(["CA", "CB", "CG", "OD1"]), chiral: None, pairs: ASP },
            "GLU" => FlipData { dihedral: Some(["CB", "CG", "CD", "OE1"]), chiral: None, pairs: GLU },
            "PHE" | "TYR" => FlipData { dihedral: Some(["CA", "CB", "CG", "CD1"]), chiral: None, pairs: PHE },
            "VAL" => FlipData { dihedral: None, chiral: Some(["CB", "CA", "CG1", "CG2"]), pairs: VAL },
            "LEU" => FlipData { dihedral: None, chiral: Some(["CG", "CB", "CD1", "CD2"]), pairs: LEU },
            _ => return None,
        })
    };
    for m in &mut st.models {
        for c in &mut m.chains {
            for rg in &mut c.residue_groups {
                let mut flip_it = false;
                for ag in &mut rg.atom_groups {
                    let Some(fd) = data(ag.resname.as_str()) else { continue };
                    if per_atom_group {
                        flip_it = false;
                    }
                    if !flip_it {
                        let names = fd.dihedral.or(fd.chiral).unwrap();
                        let sites: Vec<Vec3> = names.iter().filter_map(|n| ag.get_atom(n).map(|a| a.xyz)).collect();
                        if sites.len() != 4 {
                            continue;
                        }
                        if fd.dihedral.is_some() {
                            if let Some(d) = dihedral_deg(sites[0], sites[1], sites[2], sites[3]) {
                                if d.abs() > 90.0 {
                                    flip_it = true;
                                }
                            }
                        } else {
                            let d01 = sites[1] - sites[0];
                            let d02 = sites[2] - sites[0];
                            let d03 = sites[3] - sites[0];
                            let vol = d01.dot(d02.cross(d03));
                            let delta = -2.5 - vol;
                            if delta.abs() > 2.0 {
                                flip_it = true;
                            }
                        }
                    }
                    if flip_it {
                        let mut swaps: Vec<(usize, usize)> = Vec::new();
                        let mut incomplete = false;
                        for pr in fd.pairs {
                            let i1 = ag.atoms.iter().position(|a| a.name.trim() == pr[0]);
                            let i2 = ag.atoms.iter().position(|a| a.name.trim() == pr[1]);
                            match (i1, i2) {
                                (None, None) => continue,
                                (Some(x), Some(y)) => swaps.push((x, y)),
                                _ => {
                                    incomplete = true;
                                    break;
                                }
                            }
                        }
                        if incomplete {
                            swaps.clear();
                        }
                        for (x, y) in swaps {
                            let (px, bx) = (ag.atoms[x].xyz, ag.atoms[x].b);
                            ag.atoms[x].xyz = ag.atoms[y].xyz;
                            ag.atoms[x].b = ag.atoms[y].b;
                            ag.atoms[y].xyz = px;
                            ag.atoms[y].b = bx;
                        }
                    }
                }
            }
        }
    }
}

/// Valences by element (mmtbx.ligands.chemistry.get_valences).
fn valences(el: &str) -> Vec<i32> {
    let v = match el {
        "H" | "LI" | "NA" | "K" | "RB" | "CS" | "F" | "CL" | "BR" | "I" | "SC" | "Y" | "CU" | "AG" | "AU" => 1,
        "BE" | "MG" | "CA" | "SR" | "BA" | "O" | "S" | "SE" | "TE" | "NI" | "PD" | "PT" | "TI" | "ZR" => 2,
        "B" | "N" | "P" | "AS" | "SB" | "AL" | "GA" | "IN" | "V" | "NB" | "CO" | "RH" | "IR" => 3,
        "C" | "SI" | "GE" | "SN" | "CR" | "MO" | "FE" | "RU" | "OS" | "HF" => 4,
        "MN" | "TC" | "RE" | "TA" | "W" => 5,
        "ZN" | "CD" | "HG" | "HE" | "NE" | "AR" | "KR" | "XE" | "RN" => 0,
        _ => return vec![],
    };
    vec![v]
}

/// Bond orders by atom-name pair for a residue (dictionary, then CCD overrides).
fn bond_orders(ml: &MonLib, resname: &str, cache: &mut FxHashMap<String, FxHashMap<(String, String), i32>>) -> FxHashMap<(String, String), i32> {
    if let Some(c) = cache.get(resname) {
        return c.clone();
    }
    let mut orders: FxHashMap<(String, String), i32> = FxHashMap::default();
    let key = |a: &str, b: &str| if a < b { (a.to_string(), b.to_string()) } else { (b.to_string(), a.to_string()) };
    if let Some(cc) = ml.comp(resname) {
        for b in &cc.bonds {
            let o = match b.type_.trim().to_ascii_lowercase().as_str() {
                "double" => 2,
                "triple" => 3,
                _ => 0,
            };
            if o > 0 {
                orders.insert(key(&b.a1, &b.a2), o);
            }
        }
    }
    if let Some(e) = ml.ccd(resname) {
        for (a1, a2, o) in &e.bonds {
            let k = key(a1, a2);
            match o.trim().to_ascii_uppercase().as_str() {
                "DOUB" => {
                    orders.insert(k, 2);
                }
                "TRIP" => {
                    orders.insert(k, 3);
                }
                _ => {
                    orders.remove(&k);
                }
            }
        }
    }
    cache.insert(resname.to_string(), orders.clone());
    orders
}

/// Reference CH2 centres (n_h=2, n_heavy=2) or propellers (3, 1) from the CCD.
fn h_references(ml: &MonLib, resname: &str, n_h: usize, n_heavy: usize) -> Vec<(String, Vec<String>, Vec<String>, FxHashMap<String, Vec3>)> {
    let mut result = Vec::new();
    let groups_from = |elements: &FxHashMap<String, String>, pairs: &[(String, String)], sites: &FxHashMap<String, Vec3>| {
        let mut nb: Vec<(String, Vec<String>)> = Vec::new();
        for (a, b) in pairs {
            for (x, y) in [(a, b), (b, a)] {
                match nb.iter_mut().find(|e| e.0 == *x) {
                    Some(e) => e.1.push(y.clone()),
                    None => nb.push((x.clone(), vec![y.clone()])),
                }
            }
        }
        let is_h = |n: &str| matches!(elements.get(n).map(|s| s.as_str()), Some("H") | Some("D"));
        let mut out = Vec::new();
        for (p, ns) in nb {
            if is_h(&p) {
                continue;
            }
            let mut hs: Vec<String> = ns.iter().filter(|n| is_h(n)).cloned().collect();
            let mut hv: Vec<String> = ns.iter().filter(|n| !is_h(n)).cloned().collect();
            hs.sort();
            hv.sort();
            if hs.len() == n_h && hv.len() == n_heavy && std::iter::once(&p).chain(hv.iter()).chain(hs.iter()).all(|n| sites.contains_key(n)) {
                out.push((p, hv, hs, sites.clone()));
            }
        }
        out
    };
    let comp = ml.comp(resname);
    let describes = comp.as_ref().map(|c| c.source.to_string_lossy().contains("chem_data")).unwrap_or(true);
    if describes {
        if let Some(e) = ml.ccd(resname) {
            let elements: FxHashMap<String, String> = e.atoms.iter().map(|a| (a.0.clone(), a.1.clone())).collect();
            let model: FxHashMap<String, Vec3> =
                e.atoms.iter().filter_map(|a| a.2.map(|x| (a.0.clone(), Vec3::from_array(x)))).collect();
            let ideal: FxHashMap<String, Vec3> =
                e.atoms.iter().filter_map(|a| a.3.map(|x| (a.0.clone(), Vec3::from_array(x)))).collect();
            let pairs: Vec<(String, String)> = e.bonds.iter().map(|b| (b.0.clone(), b.1.clone())).collect();
            let fix: &[&str] = if n_h == 2 {
                match resname {
                    "ARG" => &["CB", "CG"],
                    "ILE" => &["CG1"],
                    "LEU" => &["CB"],
                    "MET" | "MSE" => &["CB", "CG"],
                    _ => &[],
                }
            } else {
                &[]
            };
            if !fix.is_empty() {
                result.extend(groups_from(&elements, &pairs, &model).into_iter().filter(|g| fix.contains(&g.0.as_str())));
            }
            result.extend(groups_from(&elements, &pairs, &ideal));
        }
    }
    if let Some(c) = comp {
        if let Some(sites) = dictionary_sites(&c, resname) {
            let elements: FxHashMap<String, String> = c.atoms.iter().map(|a| (a.id.clone(), a.type_symbol.clone())).collect();
            let pairs: Vec<(String, String)> = c.bonds.iter().map(|b| (b.a1.clone(), b.a2.clone())).collect();
            result.extend(groups_from(&elements, &pairs, &sites));
        }
    }
    result
}

/// Ideal sites from a restraint dictionary that carries coordinates.
fn dictionary_sites(c: &Comp, resname: &str) -> Option<FxHashMap<String, Vec3>> {
    let text = std::fs::read_to_string(&c.source).ok()?;
    let doc = crate::cif::parse(&text);
    for b in &doc.blocks {
        let Some(cat) = b.category("_chem_comp_atom") else { continue };
        let (Some(ai), Some(xi), Some(yi), Some(zi)) = (cat.col("atom_id"), cat.col("x"), cat.col("y"), cat.col("z")) else { continue };
        let ci = cat.col("comp_id");
        if ci.is_none() && b.name != format!("comp_{}", resname) {
            continue;
        }
        let mut sites = FxHashMap::default();
        for r in 0..cat.nrows() {
            if let Some(ci) = ci {
                if cat.get(r, ci).trim() != resname {
                    continue;
                }
            }
            let (x, y, z) = (crate::cif::parse_f64(cat.get(r, xi)), crate::cif::parse_f64(cat.get(r, yi)), crate::cif::parse_f64(cat.get(r, zi)));
            if let (Some(x), Some(y), Some(z)) = (x, y, z) {
                sites.insert(cat.get(r, ai).trim_matches('"').to_string(), v3(x, y, z));
            }
        }
        if !sites.is_empty() {
            return Some(sites);
        }
    }
    None
}

fn chiral_volume(c: Vec3, a: Vec3, b: Vec3, h: Vec3) -> f64 {
    (a - c).dot((b - c).cross(h - c))
}

/// A residue's dictionary as the monomer server has it after interpretation:
/// the libraries, or the CCD-built one registered during placement.
fn server_comp(ml: &MonLib, auto_comps: &FxHashMap<String, Arc<Comp>>, resname: &str) -> Option<Arc<Comp>> {
    ml.comp(resname).or_else(|| auto_comps.get(&resname.trim().to_ascii_uppercase()).cloned())
}

/// Per residue name: the number of heavy-atom bonds of each dictionary atom
/// (what riding uses to tell a missing heavy neighbor from a complete parent).
pub fn expected_heavy_table(
    ml: &MonLib,
    auto_comps: &FxHashMap<String, Arc<Comp>>,
    resnames: &[String],
) -> FxHashMap<String, Option<FxHashMap<String, usize>>> {
    let mut heavy_cache: FxHashMap<String, Option<FxHashMap<String, usize>>> = FxHashMap::default();
    for r in resnames.iter() {
        if heavy_cache.contains_key(r) {
            continue;
        }
        let v = server_comp(ml, auto_comps, r).map(|c| {
            let mut counts: FxHashMap<String, usize> = c.atoms.iter().map(|a| (a.id.trim().to_string(), 0)).collect();
            for b in &c.bonds {
                let (x, y) = (c.atom(&b.a1), c.atom(&b.a2));
                let (Some(x), Some(y)) = (x, y) else { continue };
                if x.type_symbol != "H" {
                    *counts.entry(b.a2.trim().to_string()).or_insert(0) += 1;
                }
                if y.type_symbol != "H" {
                    *counts.entry(b.a1.trim().to_string()).or_insert(0) += 1;
                }
            }
            counts
        });
        heavy_cache.insert(r.clone(), v);
    }
    heavy_cache
}

/// `model.process()` moves every site through fractional space and back
/// twice (`apply_symmetry_sites`); repeat that round-off.
pub fn round_off_like_cctbx(st: &mut Structure, uc: &crate::cell::UnitCell) {
    for m in &mut st.models {
        for c in &mut m.chains {
            for rg in &mut c.residue_groups {
                for ag in &mut rg.atom_groups {
                    for a in &mut ag.atoms {
                        let mut x = [a.xyz];
                        crate::cell::roundtrip_sites(&mut x, uc, 2);
                        a.xyz = x[0];
                    }
                }
            }
        }
    }
}

/// Working state over the flat atom list during placement.
struct Work {
    flat: FlatAtoms,
    removed: Vec<bool>,
}

/// Run hydrogen placement on a structure in place.
pub fn place_hydrogens(st: &mut Structure, ml: &MonLib, p: &HPlaceParams) -> Placed {
    let mut log = String::new();
    // element X atoms are removed by the caller (reduce2 does it too)
    let n_h_initial = st.models.iter().flat_map(|m| &m.chains).flat_map(|c| &c.residue_groups).flat_map(|r| &r.atom_groups).flat_map(|g| &g.atoms).filter(|a| a.is_hydrogen()).count();
    if !p.keep_existing_h {
        st.retain_atoms(|a| !a.is_hydrogen());
    }
    let mut no_h_placed = Vec::new();
    let mut auto_comps: FxHashMap<String, Arc<Comp>> = FxHashMap::default();
    add_missing_h(st, ml, p, &mut no_h_placed, &mut auto_comps);
    place_n_terminal_propeller(st, ml, p);
    st.sort_atoms_in_place();
    st.reset_serial();
    st.reset_i_seq();
    flip_symmetric_amino_acids(st, !p.compat);
    if p.compat {
        if let Some(uc) = &p.cell {
            round_off_like_cctbx(st, uc);
        }
    }

    crate::model::mem_checkpoint("placement: H added");
    let flat = FlatAtoms::from_structure(st);
    let n = flat.pos.len();
    let ip = InterpParams { neutron: p.neutron, link_distance_cutoff: 3.0, compat: p.compat, auto_comps: auto_comps.clone() };
    let mut it = interp::interpret(st, &flat, ml, &ip);
    log += &it.log;
    let mut w = Work { flat, removed: vec![false; n] };

    add_link_h_restraints(&mut it, ml, &w.flat);

    // ---------------- riding
    let occ: Vec<f64> = w.flat.path.iter().map(|&pp| st.atom(pp).occ).collect();
    let rg_of: Vec<u64> = w.flat.path.iter().map(|pp| ((pp.model as u64) << 40) | ((pp.chain as u64) << 20) | pp.rg as u64).collect();
    let resnames: Vec<String> = w.flat.path.iter().map(|&pp| st.atom_group(pp).resname.trim().to_string()).collect();
    let heavy_cache = expected_heavy_table(ml, &auto_comps, &resnames);
    let expected_heavy = |a: u32| -> Option<usize> {
        let m = heavy_cache.get(&resnames[a as usize])?.as_ref()?;
        m.get(w.flat.name[a as usize].trim()).copied()
    };
    let mut sites = w.flat.pos.clone();
    let ra = RidingAtoms {
        is_h: &w.flat.is_h,
        altloc: &w.flat.altloc,
        name: &w.flat.name,
        occ: &occ,
        rg_of: &rg_of,
        expected_heavy: &expected_heavy,
        dictionary_nh2_torsion: !p.compat,
        reroute_axial_reference: !p.compat,
    };
    crate::model::mem_checkpoint("placement: interpreted");
    let rr = riding::riding(&it, &mut sites, &ra, true);
    crate::model::mem_checkpoint("placement: riding done");
    for msg in &rr.warnings {
        log.push_str(msg);
        log.push('\n');
    }
    let axial_reference: Vec<String> = rr.axial_reference.iter().map(|&h| label(st, w.flat.path[h as usize])).collect();
    let mut coefs = rr.coef;
    w.flat.pos = sites;

    // ---------------- remove H that could not be parameterized (waters keep theirs)
    let bonds_all: Vec<(u32, u32)> = it.bonds.iter().map(|b| (b.i, b.j)).collect();
    let mut nbrs: Vec<Vec<u32>> = vec![Vec::new(); n];
    for &(i, j) in &bonds_all {
        nbrs[i as usize].push(j);
        nbrs[j as usize].push(i);
    }
    let heavy_neighbors = |nbrs: &Vec<Vec<u32>>, flat: &FlatAtoms, i: u32| -> Vec<u32> {
        let mut v: Vec<u32> = nbrs[i as usize].iter().copied().filter(|&m| !flat.is_h[m as usize]).collect();
        v.sort_unstable();
        v.dedup();
        v
    };
    let mut site_labels_no_para = Vec::new();
    let unpara: Vec<u32> = rr
        .unparameterized
        .iter()
        .copied()
        .filter(|&h| resclass::get_class(&st.atom_group(w.flat.path[h as usize]).resname) != ResClass::CommonWater)
        .collect();
    let placed_anchorless = place_anchorless_h(&it, ml, &auto_comps, st, &mut w, &unpara, &nbrs);
    for &h in &unpara {
        if placed_anchorless.contains(&h) {
            continue;
        }
        let isolated = nbrs[h as usize].is_empty();
        if !isolated {
            let parents = heavy_neighbors(&nbrs, &w.flat, h);
            let missing_neighbor = parents.len() == 1 && heavy_neighbors(&nbrs, &w.flat, parents[0]).len() < 2;
            if !missing_neighbor {
                site_labels_no_para.push(label(st, w.flat.path[h as usize]));
            }
        }
        w.removed[h as usize] = true;
        coefs[h as usize] = None;
    }

    // ---------------- write positions; reset ADP and occupancy of H
    for k in 0..n {
        let pth = w.flat.path[k];
        st.atom_mut(pth).xyz = w.flat.pos[k];
    }
    reset_h_adp_occ(st, &w, &it, p.adp_scale, if p.compat { p.cell.as_ref() } else { None });

    // ---------------- prochiral naming, links, esterified O
    let mut name_done: FxHashSet<u32> = FxHashSet::default();
    name_prochiral_h(st, ml, &w, &[2], &mut name_done);
    let removed_on_links = exclude_h_on_links(st, ml, &mut w, &it);
    exclude_h_on_esterified_o(&mut w, &it);
    name_prochiral_h(st, ml, &w, &[3], &mut name_done);
    // write back any positions changed by exclude_h_on_links
    for k in 0..n {
        let pth = w.flat.path[k];
        st.atom_mut(pth).xyz = w.flat.pos[k];
    }

    // ---------------- compact: drop removed atoms and remap everything
    let mut new_index = vec![u32::MAX; n];
    let mut next = 0u32;
    for k in 0..n {
        if !w.removed[k] {
            new_index[k] = next;
            next += 1;
        }
    }
    let removed_paths: FxHashSet<AtomPath> = (0..n).filter(|&k| w.removed[k]).map(|k| w.flat.path[k]).collect();
    // remove by path (mark via a sentinel serial)
    {
        let mut idx = 0usize;
        let paths = st.atom_paths();
        let mut keep = vec![true; paths.len()];
        for (k, pth) in paths.iter().enumerate() {
            if removed_paths.contains(pth) {
                keep[k] = false;
            }
        }
        st.retain_atoms(|_| {
            let r = keep[idx];
            idx += 1;
            r
        });
    }
    st.reset_i_seq();
    let mut bonds = Vec::new();
    for b in &it.bonds {
        let (i, j) = (new_index[b.i as usize], new_index[b.j as usize]);
        if i != u32::MAX && j != u32::MAX {
            bonds.push((i, j, b.origin));
        }
    }
    let mut etype = Vec::with_capacity(next as usize);
    let mut dict = Vec::with_capacity(next as usize);
    let mut riding_out = Vec::with_capacity(next as usize);
    for k in 0..n {
        if w.removed[k] {
            continue;
        }
        etype.push(it.etype[k].clone());
        dict.push(it.dict[k].clone());
        riding_out.push(coefs[k].map(|mut c| {
            let map = |x: i64| if x < 0 { x } else { new_index[x as usize] as i64 };
            c.ih = new_index[c.ih as usize];
            c.a0 = new_index[c.a0 as usize];
            c.a1 = new_index[c.a1 as usize];
            c.a2 = map(c.a2);
            c.a3 = map(c.a3);
            c
        }));
    }
    let n_h_final = st.models.iter().flat_map(|m| &m.chains).flat_map(|c| &c.residue_groups).flat_map(|r| &r.atom_groups).flat_map(|g| &g.atoms).filter(|a| a.is_hydrogen()).count();
    let mut uniq_missing = Vec::new();
    for r in no_h_placed {
        if !uniq_missing.contains(&r) {
            uniq_missing.push(r);
        }
    }
    Placed {
        bonds,
        etype,
        riding: riding_out,
        dict,
        no_h_placed: uniq_missing,
        site_labels_no_para,
        axial_reference,
        removed_on_links,
        n_h_initial,
        n_h_final,
        log,
    }
}

/// `add_link_h_restraints`: amide H named other than H/D get the TRANS
/// C-N-H angle (and plane), and their amine CA-N-H angle is replaced.
fn add_link_h_restraints(it: &mut Interp, ml: &MonLib, flat: &FlatAtoms) {
    let Some(link) = ml.link("TRANS") else { return };
    let Some(def) = link.angles.iter().find(|a| a.a[1] == "N" && (a.a[2] == "H" || a.a[2] == "D") && a.value.is_some() && a.esd.unwrap_or(0.0) > 0.0)
    else {
        return;
    };
    let (ideal, esd) = (def.value.unwrap(), def.esd.unwrap());
    let ca_val = link.angles.iter().find(|a| a.a[1] == "N" && a.a[2] == "CA" && a.value.is_some()).map(|a| a.value.unwrap());
    let ca_ideal = ca_val.map(|v| 360.0 - ideal - v);
    let n = flat.pos.len();
    let mut nbrs: Vec<Vec<u32>> = vec![Vec::new(); n];
    for b in &it.bonds {
        nbrs[b.i as usize].push(b.j);
        nbrs[b.j as usize].push(b.i);
    }
    let rg_key = |a: u32| {
        let p = flat.path[a as usize];
        (p.model, p.chain, p.rg)
    };
    let mut new_angles: Vec<(u32, u32, u32, f64, f64)> = Vec::new();
    let mut stale: Vec<(u32, u32, u32)> = Vec::new();
    for h in 0..n as u32 {
        if !flat.is_h[h as usize] {
            continue;
        }
        let heavy: Vec<u32> = nbrs[h as usize].iter().copied().filter(|&x| !flat.is_h[x as usize]).collect();
        if heavy.len() != 1 {
            continue;
        }
        let nn = heavy[0];
        if flat.name[nn as usize].trim() != "N" {
            continue;
        }
        for &c in &nbrs[nn as usize] {
            if flat.name[c as usize].trim() != "C" || rg_key(c) == rg_key(nn) {
                continue;
            }
            if it.has_angle(c, nn, h) {
                continue;
            }
            new_angles.push((c, nn, h, ideal, esd));
            let ca: Vec<u32> = nbrs[nn as usize].iter().copied().filter(|&x| flat.name[x as usize].trim() == "CA" && rg_key(x) == rg_key(nn)).collect();
            if ca_ideal.is_some() && !ca.is_empty() {
                stale.push((ca[0], nn, h));
            }
        }
    }
    for (ca, nn, h) in stale {
        if let Some(ap) = it.angles.iter().find(|a| a.i[1] == nn && ((a.i[0] == ca && a.i[2] == h) || (a.i[0] == h && a.i[2] == ca))) {
            let esd = ap.esd;
            it.remove_angle(ca, nn, h);
            new_angles.push((ca, nn, h, ca_ideal.unwrap(), esd));
        }
    }
    for (a, b, c, v, e) in new_angles {
        it.push_angle_unchecked(a, b, c, v, e);
    }
}

fn place_anchorless_h(
    it: &Interp,
    ml: &MonLib,
    auto_comps: &FxHashMap<String, Arc<Comp>>,
    st: &Structure,
    w: &mut Work,
    iseqs: &[u32],
    nbrs: &[Vec<u32>],
) -> Vec<u32> {
    let heavy = |i: u32| -> Vec<u32> {
        let mut v: Vec<u32> = nbrs[i as usize].iter().copied().filter(|&m| !w.flat.is_h[m as usize]).collect();
        v.sort_unstable();
        v.dedup();
        v
    };
    let dict_heavy_degree = |a: u32| -> Option<usize> {
        let ag = st.atom_group(w.flat.path[a as usize]);
        let c = server_comp(ml, auto_comps, ag.resname.trim())?;
        let name = w.flat.name[a as usize].trim();
        c.atom(name)?;
        let mut n = 0;
        for b in &c.bonds {
            for (x, y) in [(&b.a1, &b.a2), (&b.a2, &b.a1)] {
                if x == name {
                    if let Some(at) = c.atom(y) {
                        if at.type_symbol != "H" && at.type_symbol != "D" {
                            n += 1;
                        }
                    }
                }
            }
        }
        Some(n)
    };
    let mut by_parent: Vec<((u32, u32), Vec<u32>)> = Vec::new();
    for &ih in iseqs {
        let parents = heavy(ih);
        if parents.len() != 1 {
            continue;
        }
        let i0 = parents[0];
        let partners = heavy(i0);
        if partners.len() != 1 {
            continue;
        }
        let i1 = partners[0];
        if heavy(i1) != vec![i0] {
            continue;
        }
        if dict_heavy_degree(i0) != Some(1) || dict_heavy_degree(i1) != Some(1) {
            continue;
        }
        match by_parent.iter_mut().find(|x| x.0 == (i0, i1)) {
            Some(x) => x.1.push(ih),
            None => by_parent.push(((i0, i1), vec![ih])),
        }
    }
    if by_parent.is_empty() {
        return vec![];
    }
    let targets: FxHashSet<u32> = by_parent.iter().flat_map(|x| x.1.iter().copied()).collect();
    let mut dist: FxHashMap<u32, f64> = FxHashMap::default();
    let mut ang: FxHashMap<u32, f64> = FxHashMap::default();
    for b in &it.bonds {
        for x in [b.i, b.j] {
            if targets.contains(&x) {
                dist.insert(x, b.ideal);
            }
        }
    }
    for a in &it.angles {
        if targets.contains(&a.i[0]) {
            ang.insert(a.i[0], a.ideal);
        }
        if targets.contains(&a.i[2]) {
            ang.insert(a.i[2], a.ideal);
        }
    }
    let mut placed = Vec::new();
    for ((i0, i1), mut hs) in by_parent {
        let r0 = w.flat.pos[i0 as usize];
        let r1 = w.flat.pos[i1 as usize];
        let axis = if i0 < i1 { r1 - r0 } else { r0 - r1 };
        let pv = axis.ortho().normalize();
        let q = axis.normalize().cross(pv);
        let u = (r1 - r0).normalize();
        let offset = if i0 < i1 { 0.0 } else { 180.0 };
        hs.sort_unstable();
        for (k, &ih) in hs.iter().enumerate() {
            let (Some(&d), Some(&a)) = (dist.get(&ih), ang.get(&ih)) else { continue };
            let a = a.to_radians();
            let phi = (offset + 120.0 * k as f64).to_radians();
            let dir = u * a.cos() + (pv * phi.cos() + q * phi.sin()) * a.sin();
            w.flat.pos[ih as usize] = r0 + dir * d;
            placed.push(ih);
        }
    }
    placed
}

fn reset_h_adp_occ(st: &mut Structure, w: &Work, it: &Interp, scale: f64, cell: Option<&crate::cell::UnitCell>) {
    let n = w.flat.pos.len();
    // xh connectivity: (heavy, H) pairs from bonds
    let mut parent: Vec<Option<u32>> = vec![None; n];
    for b in &it.bonds {
        let (i, j) = (b.i as usize, b.j as usize);
        if w.removed[i] || w.removed[j] {
            continue;
        }
        if w.flat.is_h[i] && !w.flat.is_h[j] {
            parent[i] = Some(b.j);
        } else if w.flat.is_h[j] && !w.flat.is_h[i] {
            parent[j] = Some(b.i);
        }
    }
    // conformer occupancy of altloc atom groups (mean over non-H atoms)
    let mut conf_occ: FxHashMap<(u32, u32, u32, u32), f64> = FxHashMap::default();
    for (mi, m) in st.models.iter().enumerate() {
        for (ci, c) in m.chains.iter().enumerate() {
            for (ri, rg) in c.residue_groups.iter().enumerate() {
                for (gi, ag) in rg.atom_groups.iter().enumerate() {
                    if ag.altloc.trim().is_empty() {
                        continue;
                    }
                    let v: Vec<f64> = ag.atoms.iter().filter(|a| !a.is_hydrogen()).map(|a| a.occ).collect();
                    if !v.is_empty() {
                        conf_occ.insert((mi as u32, ci as u32, ri as u32, gi as u32), v.iter().sum::<f64>() / v.len() as f64);
                    }
                }
            }
        }
    }
    let eight_pi_sq = 8.0 * std::f64::consts::PI * std::f64::consts::PI;
    // The model's X-ray structure holds U* for anisotropic atoms; the
    // hierarchy gets U back from it (compat: with cctbx's round-off).
    if let Some(uc) = cell {
        for k in 0..n {
            if w.removed[k] {
                continue;
            }
            let a = st.atom_mut(w.flat.path[k]);
            if let Some(u) = a.uij {
                a.uij = Some(uc.u_cart_roundtrip(u));
            }
        }
    }
    for h in 0..n {
        let Some(px) = parent[h] else { continue };
        if w.removed[h] {
            continue;
        }
        let pa = st.atom(w.flat.path[px as usize]).clone();
        let u_eq = match pa.uij {
            Some(u) => (u[0] + u[1] + u[2]) / 3.0,
            None => pa.b / eight_pi_sq,
        };
        let hp = w.flat.path[h];
        let hag = (hp.model, hp.chain, hp.rg, hp.ag);
        let h_alt = !st.atom_group(hp).altloc.trim().is_empty();
        let p_blank = st.atom_group(w.flat.path[px as usize]).altloc.trim().is_empty();
        let ha = st.atom_mut(hp);
        ha.b = u_eq * eight_pi_sq * scale;
        ha.uij = None;
        ha.occ = if h_alt && p_blank && conf_occ.contains_key(&hag) { conf_occ[&hag] } else { pa.occ };
    }
    // Syncing the hierarchy from the X-ray structure rewrites the B of every
    // anisotropic atom as its equivalent isotropic B.
    for k in 0..n {
        if w.removed[k] {
            continue;
        }
        let a = st.atom_mut(w.flat.path[k]);
        if let Some(u) = a.uij {
            a.b = (u[0] + u[1] + u[2]) / 3.0 * eight_pi_sq;
        }
    }
}

fn name_prochiral_h(st: &mut Structure, ml: &MonLib, w: &Work, kinds: &[usize], done: &mut FxHashSet<u32>) {
    let mut cache: FxHashMap<(String, usize), Vec<(String, Vec<String>, Vec<String>, FxHashMap<String, Vec3>)>> = FxHashMap::default();
    // flat index by path, to honor removals
    let mut idx: FxHashMap<AtomPath, u32> = FxHashMap::default();
    for (k, pth) in w.flat.path.iter().enumerate() {
        idx.insert(*pth, k as u32);
    }
    for mi in 0..st.models.len() {
        for ci in 0..st.models[mi].chains.len() {
            let alts = st.models[mi].chains[ci].conformer_altlocs();
            for alt in &alts {
                for ri in 0..st.models[mi].chains[ci].residue_groups.len() {
                    // conformer residues: blank + this altloc, grouped by resname
                    let rg = &st.models[mi].chains[ci].residue_groups[ri];
                    let mut groups: Vec<(String, Vec<AtomPath>)> = Vec::new();
                    for (gi, ag) in rg.atom_groups.iter().enumerate() {
                        if !(ag.altloc.is_empty() || ag.altloc == *alt) {
                            continue;
                        }
                        let e = match groups.iter_mut().position(|g| g.0 == ag.resname) {
                            Some(k) => k,
                            None => {
                                groups.push((ag.resname.clone(), Vec::new()));
                                groups.len() - 1
                            }
                        };
                        for k in 0..ag.atoms.len() {
                            let pth = AtomPath { model: mi as u32, chain: ci as u32, rg: ri as u32, ag: gi as u32, atom: k as u32 };
                            if let Some(&fi) = idx.get(&pth) {
                                if w.removed[fi as usize] {
                                    continue;
                                }
                            }
                            groups[e].1.push(pth);
                        }
                    }
                    for (resname, paths) in groups {
                        let rn = resname.trim().to_string();
                        let mut refs = Vec::new();
                        for &k in kinds {
                            let n_heavy = if k == 2 { 2 } else { 1 };
                            let key = (rn.clone(), k);
                            if !cache.contains_key(&key) {
                                cache.insert(key.clone(), h_references(ml, &rn, k, n_heavy));
                            }
                            for g in &cache[&key] {
                                refs.push((k, g.clone()));
                            }
                        }
                        if refs.is_empty() {
                            continue;
                        }
                        let mut atoms: FxHashMap<String, AtomPath> = FxHashMap::default();
                        for &pth in &paths {
                            atoms.insert(st.atom(pth).name.trim().to_string(), pth);
                        }
                        for (k, (p, hv, hs, sites)) in refs {
                            if !std::iter::once(&p).chain(hv.iter()).chain(hs.iter()).all(|n| atoms.contains_key(n)) {
                                continue;
                            }
                            let pp = atoms[&p];
                            let pidx = idx[&pp];
                            if done.contains(&pidx) {
                                continue;
                            }
                            let ppos = st.atom(pp).xyz;
                            if hv.iter().chain(hs.iter()).any(|n| st.atom(atoms[n]).xyz.dist(ppos) > 2.4) {
                                continue;
                            }
                            done.insert(pidx);
                            let (refn, swap) = if k == 2 {
                                ([hv[0].clone(), hv[1].clone(), hs[0].clone()], (hs[0].clone(), hs[1].clone()))
                            } else {
                                ([hs[0].clone(), hs[1].clone(), hv[0].clone()], (hs[1].clone(), hs[2].clone()))
                            };
                            let v_ideal = chiral_volume(sites[&p], sites[&refn[0]], sites[&refn[1]], sites[&refn[2]]);
                            let v_model = chiral_volume(ppos, st.atom(atoms[&refn[0]]).xyz, st.atom(atoms[&refn[1]]).xyz, st.atom(atoms[&refn[2]]).xyz);
                            if v_ideal.abs() < 0.5 || v_model.abs() < 0.5 {
                                continue;
                            }
                            if (v_ideal > 0.0) != (v_model > 0.0) {
                                let (a1, a2) = (atoms[&swap.0], atoms[&swap.1]);
                                let n1 = st.atom(a1).name.clone();
                                let n2 = st.atom(a2).name.clone();
                                st.atom_mut(a1).name = n2;
                                st.atom_mut(a2).name = n1;
                            }
                        }
                    }
                }
            }
        }
    }
}

/// `exclude_H_on_links`: drop H on atoms that take part in a link (bond with
/// origin != 0) when the parent's valence is exceeded.
fn exclude_h_on_links(st: &mut Structure, ml: &MonLib, w: &mut Work, it: &Interp) -> Vec<(String, String)> {
    let n = w.flat.pos.len();
    // live bond proxies in cctbx order (sorted by i, then j)
    let mut live: Vec<(u32, u32, f64, u16)> = it
        .bonds
        .iter()
        .filter(|b| !w.removed[b.i as usize] && !w.removed[b.j as usize])
        .map(|b| (b.i, b.j, b.ideal, b.origin))
        .collect();
    live.sort_unstable_by_key(|x| (x.0, x.1));
    // then the bonds to symmetry copies (cctbx's asu proxies)
    let n_simple = live.len();
    live.extend(
        it.sym_bonds
            .iter()
            .filter(|b| !w.removed[b.i as usize] && !w.removed[b.j as usize])
            .map(|b| (b.i, b.j, b.ideal, b.origin)),
    );
    // origin of the link an atom takes part in (0: none; the last proxy wins)
    let mut exclusion = vec![0u16; n];
    let mut link_partners: FxHashMap<u32, Vec<u32>> = FxHashMap::default();
    for (k, &(i, j, _, o)) in live.iter().enumerate() {
        if o != 0 {
            exclusion[i as usize] = o;
            exclusion[j as usize] = o;
            if k < n_simple {
                link_partners.entry(i).or_default().push(j);
                link_partners.entry(j).or_default().push(i);
            }
        }
    }
    // adjacency in proxy order, and the atoms in order of first appearance
    // (the key order of the original's `bonds` dict)
    let mut start_of = vec![0u32; n + 1];
    for &(i, j, _, _) in &live {
        start_of[i as usize + 1] += 1;
        start_of[j as usize + 1] += 1;
    }
    for k in 0..n {
        start_of[k + 1] += start_of[k];
    }
    let mut fill = start_of.clone();
    let mut adj = vec![0u32; start_of[n] as usize];
    let mut key_order: Vec<u32> = Vec::new();
    let mut seen = vec![false; n];
    for &(i, j, _, _) in &live {
        for (x, y) in [(i, j), (j, i)] {
            adj[fill[x as usize] as usize] = y;
            fill[x as usize] += 1;
        }
        for x in [i, j] {
            if !seen[x as usize] {
                seen[x as usize] = true;
                key_order.push(x);
            }
        }
    }
    let nbrs = |x: u32| &adj[start_of[x as usize] as usize..start_of[x as usize + 1] as usize];
    fn name_of<'a>(st: &'a Structure, path: &[AtomPath], k: u32) -> &'a str {
        st.atom(path[k as usize]).name.trim()
    }
    fn resname_of<'a>(st: &'a Structure, path: &[AtomPath], k: u32) -> &'a str {
        st.atom_group(path[k as usize]).resname.trim()
    }
    let name = |st: &Structure, k: u32| name_of(st, &w.flat.path, k).to_string();
    let resname = |st: &Structure, k: u32| resname_of(st, &w.flat.path, k).to_string();
    let his_exception = |st: &Structure, h: u32, parent: u32| -> bool {
        if !w.flat.is_h[h as usize] || resname_of(st, &w.flat.path, h) != "HIS" {
            return false;
        }
        if !["HD1", "DD1", "HE2", "DE2"].contains(&name_of(st, &w.flat.path, h)) {
            return false;
        }
        exclusion[parent as usize] == 0 || exclusion[parent as usize] == interp::ORIGIN_METAL
    };
    let mut sel_remove: Vec<u32> = Vec::new();
    let mut in_sel = vec![false; n];
    let mut parent_of: FxHashMap<u32, u32> = FxHashMap::default();
    let mut bond_len: FxHashMap<u32, f64> = FxHashMap::default();
    let mut removed_why: FxHashMap<u32, String> = FxHashMap::default();
    for &(i, j, ideal, _) in &live {
        if his_exception(st, i, j) || his_exception(st, j, i) {
            continue;
        }
        for (h, x) in [(i, j), (j, i)] {
            if w.flat.is_h[h as usize] && exclusion[x as usize] != 0 && !in_sel[h as usize] {
                in_sel[h as usize] = true;
                sel_remove.push(h);
                bond_len.insert(h, ideal);
                removed_why.insert(h, origin_name(exclusion[x as usize]));
                parent_of.insert(h, x);
            }
        }
    }
    // tertiary amide H
    for &nn in &key_order {
        if name_of(st, &w.flat.path, nn) != "N" {
            continue;
        }
        let cls = resclass::get_class(resname_of(st, &w.flat.path, nn));
        if !cls.is_amino_acid() {
            continue;
        }
        let mut uniq: Vec<u32> = nbrs(nn).to_vec();
        uniq.sort_unstable();
        uniq.dedup();
        for &k in &uniq {
            if !w.flat.is_h[k as usize] {
                continue;
            }
            let h_alt = &w.flat.altloc[k as usize];
            let n_heavy = uniq.iter().filter(|&&m| !w.flat.is_h[m as usize] && (w.flat.altloc[m as usize].is_empty() || w.flat.altloc[m as usize] == *h_alt)).count();
            if n_heavy >= 3 && !in_sel[k as usize] {
                in_sel[k as usize] = true;
                sel_remove.push(k);
                removed_why.insert(k, "a tertiary amide".into());
                parent_of.insert(k, nn);
            }
        }
    }
    let drop_order = |h: u32| -> f64 {
        let parent = parent_of[&h];
        let partners = link_partners.get(&parent).cloned().unwrap_or_default();
        let mut heavy: Vec<u32> = nbrs(parent).iter().copied().filter(|&x| !w.flat.is_h[x as usize] && !partners.contains(&x)).collect();
        heavy.sort_unstable();
        heavy.dedup();
        if partners.is_empty() || heavy.len() < 2 {
            return f64::INFINITY;
        }
        partners.iter().map(|&p| w.flat.pos[h as usize].dist(w.flat.pos[p as usize])).fold(f64::INFINITY, f64::min)
    };
    let mut order: Vec<u32> = sel_remove.iter().rev().copied().collect();
    order.sort_by(|a, b| drop_order(*a).partial_cmp(&drop_order(*b)).unwrap());
    let ideal_of = |i: u32, j: u32| -> f64 { it.bond_ideal(i, j).unwrap_or(0.0) };
    let mut removed_edges: FxHashSet<(u32, u32)> = FxHashSet::default();
    let mut order_cache: FxHashMap<String, FxHashMap<(String, String), i32>> = FxHashMap::default();
    let mut keep: Vec<u32> = Vec::new();
    for h in order {
        let j = parent_of[&h];
        let vals = valences(&w.flat.element[j as usize]);
        let mut nb_count = 0;
        for &x in nbrs(j) {
            if removed_edges.contains(&(j.min(x), j.max(x))) {
                continue;
            }
            nb_count += bond_order(st, ml, w, j, x, &ideal_of, &mut order_cache);
        }
        if vals.contains(&nb_count) {
            keep.push(h);
            removed_why.remove(&h);
        } else {
            removed_edges.insert((j.min(h), j.max(h)));
        }
    }
    // geometry fix-ups for H kept on linked parents (uses the original
    // connectivity, like the restraint tables in the original)
    let orig_nb = |x: u32| -> Vec<u32> {
        let mut v = nbrs(x).to_vec();
        v.sort_unstable();
        v.dedup();
        v
    };
    let parent_fsc = |h: u32| -> Option<u32> { orig_nb(h).first().copied() };
    let mut renames: Vec<u32> = Vec::new();
    for &h in &keep {
        let Some(parent) = parent_fsc(h) else { continue };
        let first_neighbors: Vec<u32> = orig_nb(parent).into_iter().filter(|&x| x != h).collect();
        let fn_filtered: Vec<u32> = first_neighbors.into_iter().filter(|&x| !in_sel[x as usize]).collect();
        let siblings: Vec<u32> = {
            let mut s: Vec<u32> = keep.iter().copied().filter(|&k| parent_fsc(k) == Some(parent)).collect();
            s.sort_unstable();
            s
        };
        let n_kept = siblings.len();
        let coordp = w.flat.pos[parent as usize];
        if n_kept == 2 && fn_filtered.len() == 2 {
            let u1 = (w.flat.pos[fn_filtered[0] as usize] - coordp).normalize();
            let u2 = (w.flat.pos[fn_filtered[1] as usize] - coordp).normalize();
            let anti = -(u1 + u2).normalize();
            let perp = u1.cross(u2).normalize();
            let beta = 54.735f64.to_radians();
            let s = if siblings.iter().position(|&x| x == h) == Some(0) { 1.0 } else { -1.0 };
            let d = anti * beta.cos() + perp * (s * beta.sin());
            w.flat.pos[h as usize] = coordp + d.normalize() * bond_len.get(&h).copied().unwrap_or(1.0);
            continue;
        }
        if n_kept > 1 {
            continue;
        }
        renames.push(h);
        let bl = bond_len.get(&h).copied().unwrap_or(1.0);
        if fn_filtered.len() == 3 {
            let c1 = w.flat.pos[fn_filtered[0] as usize];
            let c2 = w.flat.pos[fn_filtered[1] as usize];
            let c3 = w.flat.pos[fn_filtered[2] as usize];
            let mut orth = (c2 - c1).cross(c3 - c1).normalize();
            if orth.dot(coordp - c1) > 0.0 {
                orth = -orth;
            }
            w.flat.pos[h as usize] = coordp - orth * bl;
        } else if fn_filtered.len() == 2 {
            let c1 = w.flat.pos[fn_filtered[0] as usize];
            let c2 = w.flat.pos[fn_filtered[1] as usize];
            let half = (c1 - coordp).normalize() + (c2 - coordp).normalize();
            w.flat.pos[h as usize] = coordp - half.normalize() * bl;
        }
    }
    let keep_set: FxHashSet<u32> = keep.iter().copied().collect();
    let mut out = Vec::new();
    for &h in &sel_remove {
        if keep_set.contains(&h) {
            continue;
        }
        w.removed[h as usize] = true;
        out.push((label(st, w.flat.path[h as usize]), removed_why.get(&h).cloned().unwrap_or_default()));
    }
    // rename a kept H1/H2/H3 on a linked backbone N to H
    for &h in &keep {
        let parent = parent_of[&h];
        let kept_on_parent = keep.iter().filter(|&&k| parent_of[&k] == parent).count();
        if kept_on_parent != 1 {
            continue;
        }
        let pn = name(st, parent);
        if pn == "N" && resclass::get_class(&resname(st, parent)).is_amino_acid() {
            let hn = name(st, h);
            if ["H1", "H2", "H3", "D1", "D2", "D3"].contains(&hn.as_str()) {
                let el = w.flat.element[h as usize].clone();
                st.atom_mut(w.flat.path[h as usize]).name = format!(" {:<3}", el);
            }
        }
    }
    let _ = n;
    out
}

fn origin_name(o: u16) -> String {
    match o {
        1 => "SS BOND".into(),
        3 => "metal coordination".into(),
        5 => "glycosidic custom".into(),
        10 => "Misc. bond".into(),
        x => format!("link {}", x),
    }
}

fn bond_order(
    st: &Structure,
    ml: &MonLib,
    w: &Work,
    i: u32,
    j: u32,
    ideal: &dyn Fn(u32, u32) -> f64,
    cache: &mut FxHashMap<String, FxHashMap<(String, String), i32>>,
) -> i32 {
    let (pi, pj) = (w.flat.path[i as usize], w.flat.path[j as usize]);
    let ri = (pi.model, pi.chain, pi.rg, st.atom_group(pi).resname.trim().to_string());
    let rj = (pj.model, pj.chain, pj.rg, st.atom_group(pj).resname.trim().to_string());
    if ri != rj {
        return 1;
    }
    if !cache.contains_key(ri.3.as_str()) {
        bond_orders(ml, &ri.3, cache);
    }
    let orders = &cache[ri.3.as_str()];
    let (a, b) = (st.atom(pi).name.trim().to_string(), st.atom(pj).name.trim().to_string());
    let key = if a < b { (a, b) } else { (b, a) };
    let o = *orders.get(&key).unwrap_or(&1);
    let o = if o == 0 { 1 } else { o };
    if o > 1 {
        let d = w.flat.pos[i as usize].dist(w.flat.pos[j as usize]);
        if d > ideal(i, j) + 0.1 {
            return 1;
        }
    }
    o
}

/// `exclude_H_on_esterified_O`: an O with two heavy neighbours loses its H.
fn exclude_h_on_esterified_o(w: &mut Work, it: &Interp) {
    let mut h_on: FxHashMap<u32, Vec<u32>> = FxHashMap::default();
    let mut heavy_on: FxHashMap<u32, FxHashSet<u32>> = FxHashMap::default();
    for b in &it.bonds {
        if w.removed[b.i as usize] || w.removed[b.j as usize] {
            continue;
        }
        for (o, other) in [(b.i, b.j), (b.j, b.i)] {
            if w.flat.element[o as usize] != "O" {
                continue;
            }
            if w.flat.is_h[other as usize] {
                h_on.entry(o).or_default().push(other);
            } else {
                heavy_on.entry(o).or_default().insert(other);
            }
        }
    }
    for (o, hs) in h_on {
        if heavy_on.get(&o).map(|s| s.len()).unwrap_or(0) > 1 {
            for h in hs {
                w.removed[h as usize] = true;
            }
        }
    }
}
