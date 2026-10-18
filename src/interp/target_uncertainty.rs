//! Report deposited bonds whose targets lack a supported chemical environment.

use super::ORIGIN_COVALENT;
use crate::model::{Atom, AtomLabel, AtomPath, EndpointResolver, Structure};
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

/// False for every atom whose trimmed name is not OXT, without trimming.
fn maybe_oxt(atom: &Atom) -> bool { atom.name.as_bytes().contains(&b'X') }

// Atoms are addressed by their index in `paths`, the structure's atom order.
fn atom<'a>(st: &'a Structure, paths: &[AtomPath], index: u32) -> &'a Atom {
    st.atom(paths[index as usize])
}

/// The trimmed atom name, and whether the upper-case element is `element`.
fn name_and_element<'a>(st: &'a Structure, paths: &[AtomPath], index: u32, element: &str) -> (&'a str, bool) {
    let atom = atom(st, paths, index);
    (atom.name.trim(), atom.elem().eq_ignore_ascii_case(element))
}

fn cysteine_sulfur(st: &Structure, paths: &[AtomPath], index: u32) -> bool {
    let (name, sulfur) = name_and_element(st, paths, index, "S");
    sulfur && name == "SG" && st.atom_group(paths[index as usize]).resname.trim() == "CYS"
}

fn peptide(st: &Structure, paths: &[AtomPath], first: u32, second: u32) -> bool {
    let names = (atom(st, paths, first).name.trim(), atom(st, paths, second).name.trim());
    !same_residue(paths, first, second) && matches!(names, ("N", "C") | ("C", "N")) && [first, second].iter().all(|&index|
        crate::resclass::get_class(st.atom_group(paths[index as usize]).resname.trim()).is_amino_acid())
}

fn peptide_origin(origin: Option<u16>) -> bool {
    // Named CIS, NMCIS, NMTRANS, PCIS, PTRANS and TRANS links use dictionary origins.
    matches!(origin, Some(ORIGIN_COVALENT | 35 | 63 | 65 | 66 | 70 | 79))
}

fn pair(a: u32, b: u32) -> (u32, u32) { (a.min(b), a.max(b)) }

fn same_residue(paths: &[AtomPath], a: u32, b: u32) -> bool {
    let (a, b) = (paths[a as usize], paths[b as usize]);
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

fn add(result: &mut FxHashMap<u32, BondUncertainty>, st: &Structure, paths: &[AtomPath],
       index: u32, hydrogens_unknown: bool, reason: BondUncertaintyReason) {
    let parent = atom(st, paths, index);
    let parent_source = parent.src;
    if parent_source == Atom::NEW || parent.is_hydrogen() { return; }
    let entry = result.entry(parent_source).or_insert(BondUncertainty {
        parent_source, hydrogens_unknown, reason,
    });
    entry.hydrogens_unknown |= hydrogens_unknown;
    if priority(reason) > priority(entry.reason) { entry.reason = reason; }
}

fn deposited_pairs(st: &Structure, paths: &[AtomPath],
    result: &mut FxHashMap<u32, BondUncertainty>) -> FxHashSet<(u32, u32)> {
    let mut pairs = FxHashSet::default();
    if !st.links.iter().any(|link| matches!(link.kind.as_str(), "covale" | "disulf")) { return pairs; }
    let labels = st.links.iter().filter(|link| matches!(link.kind.as_str(), "covale" | "disulf"))
        .flat_map(|link| [&link.atom1, &link.atom2]);
    let resolver = EndpointResolver::for_labels(st, paths, labels);
    let matches = |label: &AtomLabel| resolver.matching_indices(label).collect::<Vec<_>>();
    for link in st.links.iter().filter(|link| matches!(link.kind.as_str(), "covale" | "disulf")) {
        let (first, second) = (matches(&link.atom1), matches(&link.atom2));
        let ([a], [b]) = (&first[..], &second[..]) else {
            for &index in first.iter().chain(&second) {
                add(result, st, paths, index, true, BondUncertaintyReason::AmbiguousEndpoint);
            }
            continue;
        };
        if a == b || paths[*a as usize].model != paths[*b as usize].model {
            for index in [*a, *b] {
                add(result, st, paths, index, true, BondUncertaintyReason::AmbiguousEndpoint);
            }
            continue;
        }
        pairs.insert(pair(*a, *b));
    }
    pairs
}

fn multiple_sulfur_partners(st: &Structure, paths: &[AtomPath],
    pairs: &FxHashSet<(u32, u32)>) -> FxHashSet<u32> {
    let mut counts: FxHashMap<u32, usize> = FxHashMap::default();
    for &(a, b) in pairs {
        for index in [a, b] {
            if cysteine_sulfur(st, paths, index) { *counts.entry(index).or_default() += 1; }
        }
    }
    counts.into_iter().filter_map(|(index, count)| (count > 1).then_some(index)).collect()
}


fn report_pair(st: &Structure, paths: &[AtomPath], edge: (u32, u32), origin: Option<u16>,
    multiple: &FxHashSet<u32>, h_parents: &FxHashMap<u32, bool>, result: &mut FxHashMap<u32, BondUncertainty>) {
    let (a, b) = edge;
    let reason = if multiple.contains(&a) || multiple.contains(&b) {
        Some(BondUncertaintyReason::MultipleSulfurPartners)
    } else if origin.is_none() { Some(BondUncertaintyReason::UninterpretedBond) } else { None };
    if let Some(reason) = reason {
        for index in [a, b] { add(result, st, paths, index, true, reason); }
    } else if !peptide_origin(origin) || !peptide(st, paths, a, b) {
        for index in [a, b] {
            if canonical(st.atom_group(paths[index as usize]).resname.trim()) {
                let sulfur = cysteine_sulfur(st, paths, index);
                if !sulfur || h_parents.contains_key(&index) {
                    add(result, st, paths, index, sulfur && h_parents[&index], BondUncertaintyReason::UnsupportedResidueEnvironment);
                }
            }
        }
    }
}

/// The (carbon, oxygen) orders of a bond that joins C to OXT in a canonical residue.
fn terminal_oxygen_bond<'a>(st: &'a Structure, paths: &'a [AtomPath], a: u32, b: u32) -> impl Iterator<Item = (u32, u32)> + 'a {
    [(a, b), (b, a)].into_iter().filter(move |&(carbon, oxygen)|
        name_and_element(st, paths, carbon, "C") == ("C", true)
            && name_and_element(st, paths, oxygen, "O") == ("OXT", true)
            && canonical(st.atom_group(paths[carbon as usize]).resname.trim()))
}

fn terminal_oxygen(st: &Structure, paths: &[AtomPath], intrinsic: &[(u32, u32)],
    result: &mut FxHashMap<u32, BondUncertainty>) {
    for &(a, b) in intrinsic {
        for (carbon, oxygen) in terminal_oxygen_bond(st, paths, a, b) {
            let source = atom(st, paths, carbon).src;
            if let Some(uncertainty) = result.get(&source).cloned() {
                add(result, st, paths, oxygen, uncertainty.hydrogens_unknown, uncertainty.reason);
            }
        }
    }
}

/// Compare deposited records with the exact graph used for H placement.
/// This reports uncertainty; it never adds a bond or changes an atom.
pub fn inspect(st: &Structure, bonds: &[(u32, u32, u16)]) -> Vec<BondUncertainty> {
    if !st.links.iter().any(|link| matches!(link.kind.as_str(), "covale" | "disulf")) { return Vec::new(); }
    let paths = st.atom_paths();
    let mut result = FxHashMap::default();
    let mut pairs = deposited_pairs(st, &paths, &mut result);
    // One pass over the bonds gathers what the deposited pairs and the C-OXT
    // bonds need: the last origin given for each such pair, and whether each
    // pair endpoint has an H, and a generated one.
    let mut endpoint = vec![false; paths.len()];
    for &(a, b) in &pairs { endpoint[a as usize] = true; endpoint[b as usize] = true; }
    let mut origins: FxHashMap<(u32, u32), u16> = FxHashMap::default();
    let mut terminal: FxHashMap<(u32, u32), u16> = FxHashMap::default();
    let mut h_parents: FxHashMap<u32, bool> = FxHashMap::default();
    let oxt: Vec<bool> = st.atom_refs().into_iter().map(maybe_oxt).collect();
    for &(a, b, origin) in bonds {
        let edge = pair(a, b);
        if endpoint[a as usize] && endpoint[b as usize] && pairs.contains(&edge) { origins.insert(edge, origin); }
        if (oxt[a as usize] || oxt[b as usize])
            && terminal_oxygen_bond(st, &paths, a, b).next().is_some() { terminal.insert(edge, origin); }
        for (parent, h) in [(a, b), (b, a)] {
            if !endpoint[parent as usize] { continue; }
            let h = atom(st, &paths, h);
            if h.is_hydrogen() { *h_parents.entry(parent).or_default() |= h.src == Atom::NEW; }
        }
    }
    let intrinsic = |&(a, b): &(u32, u32), origin: u16| origin == ORIGIN_COVALENT && same_residue(&paths, a, b);
    pairs.retain(|edge| !origins.get(edge).is_some_and(|&origin| intrinsic(edge, origin)));
    let multiple = multiple_sulfur_partners(st, &paths, &pairs);
    for edge in pairs {
        report_pair(st, &paths, edge, origins.get(&edge).copied(), &multiple, &h_parents, &mut result);
    }
    let terminal: Vec<_> = terminal.into_iter().filter(|(edge, origin)| intrinsic(edge, *origin))
        .map(|(edge, _)| edge).collect();
    terminal_oxygen(st, &paths, &terminal, &mut result);
    let mut result: Vec<_> = result.into_values().collect();
    result.sort_by_key(|entry| entry.parent_source);
    result
}

/// Add exact coincident heavy endpoints from the accepted placement graph.
/// Call this only after an explicit uncertainty policy encounters unplaced H.
pub fn add_degenerate_heavy_bonds(st: &Structure, bonds: &[(u32, u32, u16)],
    reports: &mut Vec<BondUncertainty>) {
    let paths = st.atom_paths();
    let atoms = st.atom_refs();
    let mut result = reports.drain(..).map(|report| (report.parent_source, report)).collect();
    for &(a, b, _) in bonds {
        let (first, second) = (atoms[a as usize], atoms[b as usize]);
        if a == b || first.is_hydrogen() || second.is_hydrogen()
            || paths[a as usize].model != paths[b as usize].model { continue; }
        if first.src == Atom::NEW || second.src == Atom::NEW || first.src == second.src
            || first.xyz != second.xyz || first.xyz.as_array().iter().any(|x| !x.is_finite()) { continue; }
        for index in [a, b] {
            add(&mut result, st, &paths, index, true, BondUncertaintyReason::DegenerateHeavyBond);
        }
    }
    *reports = result.into_values().collect();
    reports.sort_by_key(|report| report.parent_source);
}

#[cfg(test)]
mod endpoint_tests {
    use super::*;
    use crate::model::LinkRecord;

    #[test]
    fn absent_or_unrelated_links_leave_existing_uncertainties_unchanged() {
        let mut st = crate::mmcif::read_mmcif(concat!("data_no_links\nloop_\n",
            "_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_alt_id\n_atom_site.auth_comp_id\n",
            "_atom_site.auth_asym_id\n_atom_site.auth_seq_id\n",
            "_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n",
            "_atom_site.occupancy\n_atom_site.B_iso_or_equiv\nS SG . CYS A 1 0 0 0 1 10\n")).unwrap();
        let paths = st.atom_paths();
        let expected: FxHashMap<_, _> = [(7, BondUncertainty { parent_source: 7,
            hydrogens_unknown: false, reason: BondUncertaintyReason::UninterpretedBond })].into_iter().collect();
        let label = AtomLabel { name: "SG".into(), resname: "CYS".into(), chain: "A".into(),
            resseq: "1".into(), icode: "".into(), altloc: "".into() };
        for kind in [None, Some("metalc"), Some("hydrog")] {
            st.links = kind.into_iter().map(|kind| LinkRecord { atom1: label.clone(), atom2: label.clone(),
                distance: None, order: None, kind: kind.into() }).collect();
            let mut result = expected.clone();
            assert!(deposited_pairs(&st, &paths, &mut result).is_empty());
            assert_eq!(result, expected);
        }
    }
}

/// Retain observed CYS HG on an exact accepted, deposited external SG bond.
pub(crate) fn add_observed_sulfur_targets(st: &Structure, bonds: &[(u32, u32, u16)],
    observed: &FxHashSet<u32>, reports: &mut Vec<BondUncertainty>) {
    if observed.is_empty() { return; }
    let paths = st.atom_paths();
    let graph: FxHashSet<_> = bonds.iter().map(|&(a, b, _)| pair(a, b)).collect();
    let mut result = reports.drain(..).map(|report| (report.parent_source, report)).collect();
    for (a, b) in deposited_pairs(st, &paths, &mut FxHashMap::default()) {
        if same_residue(&paths, a, b) || !graph.contains(&pair(a, b))
            || [a, b].iter().any(|&index| atom(st, &paths, index).is_hydrogen()
                || atom(st, &paths, index).src == Atom::NEW) { continue; }
        for index in [a, b] {
            if observed.contains(&atom(st, &paths, index).src) && cysteine_sulfur(st, &paths, index) {
                add(&mut result, st, &paths, index, false, BondUncertaintyReason::UnsupportedResidueEnvironment);
            }
        }
    }
    *reports = result.into_values().collect();
    reports.sort_by_key(|report| report.parent_source);
}
