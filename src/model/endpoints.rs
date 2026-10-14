//! Exact deposited-endpoint lookup for one unchanged atom hierarchy.

use super::{AtomLabel, AtomPath, Structure};
use rustc_hash::FxHashMap;
use smallvec::SmallVec;

#[derive(Clone, Copy, PartialEq, Eq, Hash)]
struct Key<'a> {
    name: &'a str,
    resname: &'a str,
    chain: &'a str,
    resseq: &'a str,
    icode: &'a str,
}

impl<'a> Key<'a> {
    fn from_label(label: &'a AtomLabel) -> Self {
        Self { name: label.name.trim(), resname: label.resname.trim(), chain: &label.chain,
            resseq: label.resseq.trim(), icode: label.icode.trim() }
    }

    fn from_path(st: &'a Structure, p: AtomPath) -> Self {
        Self { name: st.atom(p).name.trim(), resname: st.atom_group(p).resname.trim(),
            chain: &st.chain(p).id, resseq: st.residue_group(p).resseq.trim(),
            icode: st.residue_group(p).icode.trim() }
    }
}

/// Keep every candidate in the caller's atom-path order, including duplicates.
pub(crate) struct EndpointResolver<'a> {
    by_label: FxHashMap<Key<'a>, SmallVec<[(u32, &'a str); 1]>>,
}

impl<'a> EndpointResolver<'a> {
    pub(crate) fn new(st: &'a Structure, paths: &[AtomPath]) -> Self {
        let mut by_label: FxHashMap<_, SmallVec<[_; 1]>> = FxHashMap::default();
        for (index, &p) in paths.iter().enumerate() {
            by_label.entry(Key::from_path(st, p)).or_default()
                .push((index as u32, st.atom_group(p).altloc.as_str()));
        }
        Self { by_label }
    }

    pub(crate) fn matching_indices<'b>(&'b self, label: &'b AtomLabel) -> impl Iterator<Item = u32> + 'b {
        self.by_label.get(&Key::from_label(label)).into_iter().flatten()
            .filter(move |(_, altloc)| label.altloc.is_empty() || *altloc == label.altloc)
            .map(|&(index, _)| index)
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
