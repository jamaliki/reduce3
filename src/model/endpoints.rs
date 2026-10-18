//! Exact deposited-endpoint lookup for one unchanged atom hierarchy.

use super::{Atom, AtomGroup, AtomLabel, AtomPath, Structure};
use rustc_hash::{FxHashMap, FxHashSet};

/// Residue part of an endpoint label: the chain id as stored, and the trimmed
/// residue number and insertion code.
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
struct ResKey<'a> {
    chain: &'a str,
    resseq: &'a str,
    icode: &'a str,
}

impl<'a> ResKey<'a> {
    fn of_label(label: &'a AtomLabel) -> Self {
        ResKey { chain: &label.chain, resseq: label.resseq.trim(), icode: label.icode.trim() }
    }
}

/// Keep every candidate in the caller's atom-path order, including duplicates.
/// Atoms are filed by residue; a lookup compares the atom fields of that
/// residue's candidates.
pub(crate) struct EndpointResolver<'a> {
    by_residue: FxHashMap<ResKey<'a>, u32>,
    /// Per residue key: (index in `paths`, atom group, atom) in path order.
    candidates: Vec<Vec<(u32, &'a AtomGroup, &'a Atom)>>,
}

impl<'a> EndpointResolver<'a> {
    #[cfg(test)]
    pub(crate) fn new(st: &'a Structure, paths: &[AtomPath]) -> Self {
        Self::build(st, paths, None)
    }

    /// Like `new`, but files only the residues that `labels` name. A lookup
    /// gives the same result as with `new` for any of these labels.
    pub(crate) fn for_labels<'l>(st: &'a Structure, paths: &[AtomPath],
        labels: impl IntoIterator<Item = &'l AtomLabel>) -> Self {
        let wanted: FxHashSet<ResKey<'l>> = labels.into_iter().map(ResKey::of_label).collect();
        Self::build(st, paths, Some(&wanted))
    }

    fn build(st: &'a Structure, paths: &[AtomPath], wanted: Option<&FxHashSet<ResKey<'_>>>) -> Self {
        let mut by_residue: FxHashMap<ResKey<'a>, u32> = FxHashMap::default();
        let mut candidates: Vec<Vec<(u32, &'a AtomGroup, &'a Atom)>> = Vec::new();
        // consecutive paths usually share a residue group
        let mut last: Option<((u32, u32, u32), Option<u32>)> = None;
        for (index, &p) in paths.iter().enumerate() {
            let rg_id = (p.model, p.chain, p.rg);
            let slot = match last {
                Some((id, slot)) if id == rg_id => slot,
                _ => {
                    let rg = st.residue_group(p);
                    let key = ResKey { chain: &st.chain(p).id, resseq: rg.resseq.trim(), icode: rg.icode.trim() };
                    let slot = wanted.is_none_or(|w| w.contains(&key)).then(|| *by_residue.entry(key).or_insert_with(|| {
                        candidates.push(Vec::new());
                        (candidates.len() - 1) as u32
                    }));
                    last = Some((rg_id, slot));
                    slot
                }
            };
            if let Some(slot) = slot {
                candidates[slot as usize].push((index as u32, st.atom_group(p), st.atom(p)));
            }
        }
        Self { by_residue, candidates }
    }

    pub(crate) fn matching_indices<'b>(&'b self, label: &'b AtomLabel) -> impl Iterator<Item = u32> + 'b {
        let key = ResKey::of_label(label);
        let (name, resname) = (label.name.trim(), label.resname.trim());
        self.by_residue.get(&key).into_iter()
            .flat_map(move |&slot| self.candidates[slot as usize].iter())
            .filter(move |(_, ag, atom)| atom.name.trim() == name && ag.resname.trim() == resname
                && (label.altloc.is_empty() || ag.altloc == label.altloc))
            .map(|&(index, _, _)| index)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::Vec3;
    use crate::model::{Atom, AtomGroup, Chain, Model, ResidueGroup};

    fn structure() -> Structure {
        Structure { models: (0..3).map(|model| Model { id: model.to_string(), chains:
            ["A", " A", "?", ".", ""].into_iter().map(|chain| Chain { id: chain.into(), residue_groups:
                [" 12 ", "12A", "?", ".", ""].into_iter().map(|resseq| ResidueGroup {
                    resseq: resseq.into(), icode: if model == 1 { " B " } else { " " }.into(),
                    link_to_previous: false, atom_groups:
                        ["", "A", "B", " A", "?", "."].into_iter().map(|altloc| AtomGroup {
                            altloc: altloc.into(), resname: " LIG ".into(),
                            atoms: [" CA ", "CA", " D1 ", "?", "."].into_iter()
                                .map(|name| Atom::new(name, "C", Vec3::ZERO)).collect(),
                        }).collect(),
                }).collect(),
            }).collect(),
        }).collect(), ..Structure::default() }
    }

    fn compare(st: &Structure, paths: &[AtomPath], label: &AtomLabel) {
        let resolver = EndpointResolver::new(st, paths);
        let expected: Vec<_> = paths.iter().enumerate().filter_map(|(index, &p)|
            label.matches(st, p).then_some(index as u32)).collect();
        assert_eq!(resolver.matching_indices(label).collect::<Vec<_>>(), expected);
        for model in 0..st.models.len() as u32 {
            let expected: Vec<_> = expected.iter().copied().filter(|&i| paths[i as usize].model == model).collect();
            assert_eq!(resolver.matching_indices(label).filter(|&i| paths[i as usize].model == model)
                .collect::<Vec<_>>(), expected);
        }
    }

    #[test]
    fn indexed_matches_preserve_scan_fields_candidates_and_order() {
        let st = structure();
        let mut paths = st.atom_paths();
        // Flat order may interleave models, repeat paths, or omit atoms.
        paths.reverse();
        paths.truncate(paths.len() - 3);
        paths.extend_from_within(..7);
        let resolver = EndpointResolver::new(&st, &paths);
        for name in ["CA", " CA ", "D1", "H1", "?", ".", "", "ca"] {
            for chain in ["A", " A", "A ", "?", ".", "", "absent"] {
                for resseq in ["12", " 12 ", "12A", "?", ".", ""] {
                    for icode in ["", " ", "B", " B ", "?"] {
                        for altloc in ["", "A", "B", " A", "A ", "?", "."] {
                            let label = AtomLabel { name: name.into(), resname: "LIG".into(),
                                chain: chain.into(), resseq: resseq.into(), icode: icode.into(), altloc: altloc.into() };
                            let expected: Vec<_> = paths.iter().enumerate().filter_map(|(index, &p)|
                                label.matches(&st, p).then_some(index as u32)).collect();
                            assert_eq!(resolver.matching_indices(&label).collect::<Vec<_>>(), expected, "{label:?}");
                        }
                    }
                }
            }
        }
        for resname in ["LIG", " LIG ", "lig", "?", ".", ""] {
            compare(&st, &paths, &AtomLabel { name: "CA".into(), resname: resname.into(),
                chain: "A".into(), resseq: "12".into(), icode: "".into(), altloc: "".into() });
        }
    }

    #[test]
    fn label_subset_matches_full_index() {
        let st = structure();
        let mut paths = st.atom_paths();
        paths.reverse();
        paths.extend_from_within(..9);
        let labels: Vec<AtomLabel> = [("CA", "A", "12", ""), ("D1", " A", "12A", "B"), ("CA", "?", "?", " B ")]
            .into_iter().map(|(name, chain, resseq, icode)| AtomLabel { name: name.into(), resname: "LIG".into(),
                chain: chain.into(), resseq: resseq.into(), icode: icode.into(), altloc: "".into() }).collect();
        let full = EndpointResolver::new(&st, &paths);
        let subset = EndpointResolver::for_labels(&st, &paths, &labels[..2]);
        for label in &labels[..2] {
            assert_eq!(subset.matching_indices(label).collect::<Vec<_>>(), full.matching_indices(label).collect::<Vec<_>>());
            assert!(full.matching_indices(label).next().is_some());
        }
        assert_eq!(subset.matching_indices(&labels[2]).count(), 0);
    }

    #[test]
    fn rebuild_tracks_changed_hierarchy_and_empty_input() {
        let mut st = structure();
        let paths = st.atom_paths();
        let label = AtomLabel { name: "CA".into(), resname: "LIG".into(), chain: "A".into(),
            resseq: "12".into(), icode: "".into(), altloc: "".into() };
        compare(&st, &paths, &label);
        st.atom_mut(paths[0]).name = "renamed".into();
        compare(&st, &paths, &label);
        compare(&st, &[], &label);
    }

    #[test]
    fn whitespace_altloc_and_numeric_ids_keep_literal_meaning() {
        let mut st = structure();
        let p = st.atom_paths()[0];
        st.atom_mut(p).name = "\u{2003}CA\u{2003}".into();
        let chain = &mut st.models[p.model as usize].chains[p.chain as usize];
        let residue = &mut chain.residue_groups[p.rg as usize];
        residue.resseq = "\u{2003}01\u{2003}".into();
        residue.icode = "\u{2003}B\u{2003}".into();
        let group = &mut residue.atom_groups[p.ag as usize];
        group.resname = "\u{2003}LIG\u{2003}".into();
        group.altloc = " ".into();
        let paths = vec![p, AtomPath { ag: 1, ..p }];
        let resolver = EndpointResolver::new(&st, &paths);
        let mut label = AtomLabel { name: "CA".into(), resname: "LIG".into(), chain: "A".into(),
            resseq: "01".into(), icode: "B".into(), altloc: " ".into() };
        compare(&st, &paths, &label);
        assert_eq!(resolver.matching_indices(&label).collect::<Vec<_>>(), vec![0]);
        label.altloc = "A".into();
        assert_eq!(resolver.matching_indices(&label).collect::<Vec<_>>(), vec![1]);
        label.altloc.clear();
        assert_eq!(resolver.matching_indices(&label).collect::<Vec<_>>(), vec![0, 1]);
        label.resseq = "1".into();
        assert_eq!(resolver.matching_indices(&label).count(), 0);
    }

}
