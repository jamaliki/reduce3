//! Report deposited bonds whose targets lack a supported chemical environment.

use super::{FlatAtoms, ORIGIN_COVALENT};
use crate::model::{Atom, AtomLabel, Structure};
use rustc_hash::{FxHashMap, FxHashSet};

/// Why a deposited bond cannot supply a full chemical-state target.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BondUncertaintyReason {
    /// The deposited endpoint does not resolve to one atom.
    AmbiguousEndpoint,
    /// Hydrogen placement did not accept this deposited bond.
    UninterpretedBond,
    /// A cysteine sulfur has more than one deposited covalent partner.
    MultipleSulfurPartners,
    /// The modified residue environment is absent from the state catalogue.
    UnsupportedResidueEnvironment,
    /// An accepted heavy bond has identical deposited endpoint coordinates.
    DegenerateHeavyBond,
}

/// One affected original heavy atom; other targets remain unchanged.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BondUncertainty {
    /// Original atom-site row of the affected heavy atom.
    pub parent_source: u32,
    /// Whether generated H targets on this parent also lack support.
    pub hydrogens_unknown: bool,
    /// The source or representation limit that prevents a full target.
    pub reason: BondUncertaintyReason,
}

fn canonical(component: &str) -> bool {
    ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS",
     "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"].contains(&component)
}

fn cysteine_sulfur(st: &Structure, flat: &FlatAtoms, index: u32) -> bool {
    flat.element[index as usize] == "S" && flat.name[index as usize].trim() == "SG"
        && st.atom_group(flat.path[index as usize]).resname.trim() == "CYS"
}

fn peptide(st: &Structure, flat: &FlatAtoms, first: u32, second: u32) -> bool {
    let names = (flat.name[first as usize].trim(), flat.name[second as usize].trim());
    !same_residue(flat, first, second) && matches!(names, ("N", "C") | ("C", "N")) && [first, second].iter().all(|&index|
        crate::resclass::get_class(st.atom_group(flat.path[index as usize]).resname.trim()).is_amino_acid())
}

fn peptide_origin(origin: Option<u16>) -> bool {
    // Named CIS, NMCIS, NMTRANS, PCIS, PTRANS and TRANS links use dictionary origins.
    matches!(origin, Some(ORIGIN_COVALENT | 35 | 63 | 65 | 66 | 70 | 79))
}

fn pair(a: u32, b: u32) -> (u32, u32) { (a.min(b), a.max(b)) }

fn same_residue(flat: &FlatAtoms, a: u32, b: u32) -> bool {
    let (a, b) = (flat.path[a as usize], flat.path[b as usize]);
    (a.model, a.chain, a.rg) == (b.model, b.chain, b.rg)
}

fn priority(reason: BondUncertaintyReason) -> u8 {
    match reason {
        BondUncertaintyReason::UnsupportedResidueEnvironment => 0,
        BondUncertaintyReason::UninterpretedBond => 1,
        BondUncertaintyReason::MultipleSulfurPartners => 2,
        BondUncertaintyReason::AmbiguousEndpoint => 3,
        BondUncertaintyReason::DegenerateHeavyBond => 4,
    }
}

fn add(result: &mut FxHashMap<u32, BondUncertainty>, st: &Structure, flat: &FlatAtoms,
       index: u32, hydrogens_unknown: bool, reason: BondUncertaintyReason) {
    let parent_source = st.atom(flat.path[index as usize]).src;
    if parent_source == Atom::NEW || flat.is_h[index as usize] { return; }
    let entry = result.entry(parent_source).or_insert(BondUncertainty {
        parent_source, hydrogens_unknown, reason,
    });
    entry.hydrogens_unknown |= hydrogens_unknown;
    if priority(reason) > priority(entry.reason) { entry.reason = reason; }
}

fn deposited_pairs(st: &Structure, flat: &FlatAtoms,
    result: &mut FxHashMap<u32, BondUncertainty>) -> FxHashSet<(u32, u32)> {
    let mut pairs = FxHashSet::default();
    let matches = |label: &AtomLabel| flat.path.iter().enumerate()
        .filter_map(|(index, &p)| label.matches(st, p).then_some(index as u32)).collect::<Vec<_>>();
    for link in st.links.iter().filter(|link| matches!(link.kind.as_str(), "covale" | "disulf")) {
        let (first, second) = (matches(&link.atom1), matches(&link.atom2));
        let ([a], [b]) = (&first[..], &second[..]) else {
            for &index in first.iter().chain(&second) {
                add(result, st, flat, index, true, BondUncertaintyReason::AmbiguousEndpoint);
            }
            continue;
        };
        if a == b || flat.path[*a as usize].model != flat.path[*b as usize].model {
            for index in [*a, *b] {
                add(result, st, flat, index, true, BondUncertaintyReason::AmbiguousEndpoint);
            }
            continue;
        }
        pairs.insert(pair(*a, *b));
    }
    pairs
}

fn multiple_sulfur_partners(st: &Structure, flat: &FlatAtoms,
    pairs: &FxHashSet<(u32, u32)>) -> FxHashSet<u32> {
    let mut counts: FxHashMap<u32, usize> = FxHashMap::default();
    for &(a, b) in pairs {
        for index in [a, b] {
            if cysteine_sulfur(st, flat, index) { *counts.entry(index).or_default() += 1; }
        }
    }
    counts.into_iter().filter_map(|(index, count)| (count > 1).then_some(index)).collect()
}

fn hydrogen_parents(st: &Structure, flat: &FlatAtoms, bonds: &[(u32, u32, u16)]) -> FxHashMap<u32, bool> {
    let mut parents = FxHashMap::default();
    for &(a, b, _) in bonds {
        for (parent, h) in [(a, b), (b, a)] {
            if flat.is_h[h as usize] {
                let generated = st.atom(flat.path[h as usize]).src == Atom::NEW;
                *parents.entry(parent).or_default() |= generated;
            }
        }
    }
    parents
}

fn report_pair(st: &Structure, flat: &FlatAtoms, edge: (u32, u32), origin: Option<u16>,
    multiple: &FxHashSet<u32>, h_parents: &FxHashMap<u32, bool>, result: &mut FxHashMap<u32, BondUncertainty>) {
    let (a, b) = edge;
    let reason = if multiple.contains(&a) || multiple.contains(&b) {
        Some(BondUncertaintyReason::MultipleSulfurPartners)
    } else if origin.is_none() { Some(BondUncertaintyReason::UninterpretedBond) } else { None };
    if let Some(reason) = reason {
        for index in [a, b] { add(result, st, flat, index, true, reason); }
    } else if !peptide_origin(origin) || !peptide(st, flat, a, b) {
        for index in [a, b] {
            if canonical(st.atom_group(flat.path[index as usize]).resname.trim()) {
                let sulfur = cysteine_sulfur(st, flat, index);
                if !sulfur || h_parents.contains_key(&index) {
                    add(result, st, flat, index, sulfur && h_parents[&index], BondUncertaintyReason::UnsupportedResidueEnvironment);
                }
            }
        }
    }
}

fn terminal_oxygen(st: &Structure, flat: &FlatAtoms, intrinsic: &FxHashSet<(u32, u32)>,
    result: &mut FxHashMap<u32, BondUncertainty>) {
    for &(a, b) in intrinsic {
        for (carbon, oxygen) in [(a, b), (b, a)] {
            if flat.name[carbon as usize].trim() != "C" || flat.element[carbon as usize] != "C"
                || flat.name[oxygen as usize].trim() != "OXT" || flat.element[oxygen as usize] != "O"
                || !canonical(st.atom_group(flat.path[carbon as usize]).resname.trim()) { continue; }
            let source = st.atom(flat.path[carbon as usize]).src;
            if let Some(uncertainty) = result.get(&source).cloned() {
                add(result, st, flat, oxygen, uncertainty.hydrogens_unknown, uncertainty.reason);
            }
        }
    }
}

/// Compare deposited records with the exact graph used for H placement.
/// This reports uncertainty; it never adds a bond or changes an atom.
pub fn inspect(st: &Structure, bonds: &[(u32, u32, u16)]) -> Vec<BondUncertainty> {
    if !st.links.iter().any(|link| matches!(link.kind.as_str(), "covale" | "disulf")) { return Vec::new(); }
    let flat = FlatAtoms::from_structure(st);
    let graph: FxHashMap<_, _> = bonds.iter().map(|&(a, b, origin)| (pair(a, b), origin)).collect();
    let intrinsic: FxHashSet<_> = graph.iter().filter_map(|(&(a, b), &origin)|
        (origin == ORIGIN_COVALENT && same_residue(&flat, a, b)).then_some(pair(a, b))).collect();
    let mut result = FxHashMap::default();
    let mut pairs = deposited_pairs(st, &flat, &mut result);
    pairs.retain(|edge| !intrinsic.contains(edge));
    let multiple = multiple_sulfur_partners(st, &flat, &pairs);
    let h_parents = hydrogen_parents(st, &flat, bonds);
    for edge in pairs {
        report_pair(st, &flat, edge, graph.get(&edge).copied(), &multiple, &h_parents, &mut result);
    }
    terminal_oxygen(st, &flat, &intrinsic, &mut result);
    let mut result: Vec<_> = result.into_values().collect();
    result.sort_by_key(|entry| entry.parent_source);
    result
}

/// Add exact coincident heavy endpoints from the accepted placement graph.
/// Call this only after an explicit uncertainty policy encounters unplaced H.
pub fn add_degenerate_heavy_bonds(st: &Structure, bonds: &[(u32, u32, u16)],
    reports: &mut Vec<BondUncertainty>) {
    let flat = FlatAtoms::from_structure(st);
    let mut result = reports.drain(..).map(|report| (report.parent_source, report)).collect();
    for &(a, b, _) in bonds {
        if a == b || flat.is_h[a as usize] || flat.is_h[b as usize]
            || flat.path[a as usize].model != flat.path[b as usize].model { continue; }
        let (first, second) = (st.atom(flat.path[a as usize]), st.atom(flat.path[b as usize]));
        if first.src == Atom::NEW || second.src == Atom::NEW || first.src == second.src
            || first.xyz != second.xyz || first.xyz.as_array().iter().any(|x| !x.is_finite()) { continue; }
        for index in [a, b] {
            add(&mut result, st, &flat, index, true, BondUncertaintyReason::DegenerateHeavyBond);
        }
    }
    *reports = result.into_values().collect();
    reports.sort_by_key(|report| report.parent_source);
}

/// Retain observed CYS HG on an exact accepted, deposited external SG bond.
pub(crate) fn add_observed_sulfur_targets(st: &Structure, bonds: &[(u32, u32, u16)],
    observed: &FxHashSet<u32>, reports: &mut Vec<BondUncertainty>) {
    if observed.is_empty() { return; }
    let flat = FlatAtoms::from_structure(st);
    let graph: FxHashSet<_> = bonds.iter().map(|&(a, b, _)| pair(a, b)).collect();
    let mut result = reports.drain(..).map(|report| (report.parent_source, report)).collect();
    for (a, b) in deposited_pairs(st, &flat, &mut FxHashMap::default()) {
        if same_residue(&flat, a, b) || !graph.contains(&pair(a, b))
            || [a, b].iter().any(|&index| flat.is_h[index as usize]
                || st.atom(flat.path[index as usize]).src == Atom::NEW) { continue; }
        for index in [a, b] {
            if observed.contains(&st.atom(flat.path[index as usize]).src) && cysteine_sulfur(st, &flat, index) {
                add(&mut result, st, &flat, index, false, BondUncertaintyReason::UnsupportedResidueEnvironment);
            }
        }
    }
    *reports = result.into_values().collect();
    reports.sort_by_key(|report| report.parent_source);
}
