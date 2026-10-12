use super::*;

/// Source rows whose coordinates and B factors exchange during interpretation.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SiteSwap {
    pub source_rows: [u32; 2],
}

/// Apply reference label conventions to a copy, without changing source atoms.
pub fn normalized_copy(source: &Structure) -> (Structure, Vec<SiteSwap>) {
    let mut result = source.clone();
    let mut changes = Vec::new();
    for model in &mut result.models {
        for chain in &mut model.chains {
            for residue in &mut chain.residue_groups {
                let mut flip = false;
                for group in &mut residue.atom_groups {
                    let Some((names, pairs, chiral)) = convention(&group.resname) else {
                        continue;
                    };
                    if !flip {
                        flip = needs_flip(group, names, chiral);
                    }
                    if !flip {
                        continue;
                    }
                    let Some(mapped) = paired_indices(group, pairs) else {
                        continue;
                    };
                    for (a, b) in mapped {
                        let (xyz, bf) = (group.atoms[a].xyz, group.atoms[a].b);
                        group.atoms[a].xyz = group.atoms[b].xyz;
                        group.atoms[a].b = group.atoms[b].b;
                        group.atoms[b].xyz = xyz;
                        group.atoms[b].b = bf;
                        changes.push(SiteSwap {
                            source_rows: [group.atoms[a].src, group.atoms[b].src],
                        });
                    }
                }
            }
        }
    }
    (result, changes)
}

fn needs_flip(group: &AtomGroup, names: [&str; 4], chiral: bool) -> bool {
    let sites = names.map(|name| {
        group
            .atoms
            .iter()
            .find(|a| a.name.trim() == name)
            .map(|a| a.xyz)
    });
    if sites.iter().any(Option::is_none) {
        return false;
    }
    let sites = sites.map(Option::unwrap);
    if chiral {
        let volume = (sites[1] - sites[0]).dot((sites[2] - sites[0]).cross(sites[3] - sites[0]));
        (-2.5 - volume).abs() > 2.0
    } else {
        crate::geom::dihedral_deg(sites[0], sites[1], sites[2], sites[3])
            .is_some_and(|d| d.abs() > 90.)
    }
}

fn paired_indices(group: &AtomGroup, pairs: &[(&str, &str)]) -> Option<Vec<(usize, usize)>> {
    let mut mapped = Vec::new();
    for &(left, right) in pairs {
        for isotope in [false, true] {
            if isotope && !left.starts_with('H') {
                continue;
            }
            let left = if isotope {
                format!("D{}", &left[1..])
            } else {
                left.into()
            };
            let right = if isotope {
                format!("D{}", &right[1..])
            } else {
                right.into()
            };
            let a = group.atoms.iter().position(|a| a.name.trim() == left);
            let b = group.atoms.iter().position(|a| a.name.trim() == right);
            match (a, b) {
                (Some(a), Some(b)) => mapped.push((a, b)),
                (None, None) => {}
                _ => return None,
            }
        }
    }
    Some(mapped)
}

type Convention = (
    [&'static str; 4],
    &'static [(&'static str, &'static str)],
    bool,
);

fn convention(name: &str) -> Option<Convention> {
    Some(match name {
        "ARG" => (
            ["CD", "NE", "CZ", "NH1"],
            &[("NH1", "NH2"), ("HH11", "HH21"), ("HH12", "HH22")],
            false,
        ),
        "ASP" => (["CA", "CB", "CG", "OD1"], &[("OD1", "OD2")], false),
        "GLU" => (["CB", "CG", "CD", "OE1"], &[("OE1", "OE2")], false),
        "PHE" | "TYR" => (
            ["CA", "CB", "CG", "CD1"],
            &[
                ("CD1", "CD2"),
                ("CE1", "CE2"),
                ("HD1", "HD2"),
                ("HE1", "HE2"),
            ],
            false,
        ),
        "VAL" => (
            ["CB", "CA", "CG1", "CG2"],
            &[
                ("CG1", "CG2"),
                ("HG11", "HG21"),
                ("HG12", "HG22"),
                ("HG13", "HG23"),
            ],
            true,
        ),
        "LEU" => (
            ["CG", "CB", "CD1", "CD2"],
            &[
                ("CD1", "CD2"),
                ("HD11", "HD21"),
                ("HD12", "HD22"),
                ("HD13", "HD23"),
            ],
            true,
        ),
        _ => return None,
    })
}

/// CCTBX swaps coordinates under these labels before building restraints.
pub(super) fn needs_reference_swaps(st: &Structure) -> bool {
    st.models.iter().any(|model| {
        model.chains.iter().any(|chain| {
            chain.residue_groups.iter().any(|residue| {
                let mut flip = false;
                residue.atom_groups.iter().any(|group| {
                    let names = match group.resname.as_str() {
                        "ARG" => ["CD", "NE", "CZ", "NH1"],
                        "ASP" => ["CA", "CB", "CG", "OD1"],
                        "GLU" => ["CB", "CG", "CD", "OE1"],
                        "PHE" | "TYR" => ["CA", "CB", "CG", "CD1"],
                        "VAL" => ["CB", "CA", "CG1", "CG2"],
                        "LEU" => ["CG", "CB", "CD1", "CD2"],
                        _ => return false,
                    };
                    let sites = names.map(|name| {
                        group
                            .atoms
                            .iter()
                            .find(|a| a.name.trim() == name)
                            .map(|a| a.xyz)
                    });
                    if !flip && sites.iter().all(Option::is_some) {
                        let sites = sites.map(Option::unwrap);
                        flip = if matches!(group.resname.as_str(), "VAL" | "LEU") {
                            let volume = (sites[1] - sites[0])
                                .dot((sites[2] - sites[0]).cross(sites[3] - sites[0]));
                            (-2.5 - volume).abs() > 2.0
                        } else {
                            crate::geom::dihedral_deg(sites[0], sites[1], sites[2], sites[3])
                                .is_some_and(|d| d.abs() > 90.)
                        };
                    }
                    flip
                })
            })
        })
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::v3;
    use crate::model::{Atom, AtomGroup, Chain, Model, ResidueGroup};

    #[test]
    fn detects_reference_coordinate_swaps_without_changing_sites() {
        let mut st = Structure {
            models: vec![Model {
                id: "1".into(),
                chains: vec![Chain {
                    id: "A".into(),
                    residue_groups: vec![ResidueGroup {
                        resseq: "1".into(),
                        icode: "".into(),
                        link_to_previous: false,
                        atom_groups: vec![AtomGroup {
                            resname: "VAL".into(),
                            altloc: "".into(),
                            atoms: vec![
                                Atom::new("CB", "C", v3(0., 0., 0.)),
                                Atom::new("CA", "C", v3(1., 0., 0.)),
                                Atom::new("CG1", "C", v3(0., 1., 0.)),
                                Atom::new("CG2", "C", v3(0., 0., 2.5)),
                            ],
                        }],
                    }],
                }],
            }],
            ..Default::default()
        };
        let before = FlatAtoms::from_structure(&st).pos;
        assert!(needs_reference_swaps(&st));
        assert_eq!(before, FlatAtoms::from_structure(&st).pos);
        let (normalized, changes) = normalized_copy(&st);
        assert!(!needs_reference_swaps(&normalized));
        assert_eq!(changes.len(), 1);
        assert_eq!(FlatAtoms::from_structure(&normalized).pos[2], before[3]);
        assert_eq!(before, FlatAtoms::from_structure(&st).pos);
        let mut incomplete = st.clone();
        incomplete.models[0].chains[0].residue_groups[0].atom_groups[0]
            .atoms
            .push(Atom::new("HG11", "H", v3(2., 2., 2.)));
        let (unchanged, changes) = normalized_copy(&incomplete);
        assert!(changes.is_empty());
        assert_eq!(
            FlatAtoms::from_structure(&incomplete).pos,
            FlatAtoms::from_structure(&unchanged).pos
        );
        st.models[0].chains[0].residue_groups[0].atom_groups[0].atoms[3]
            .xyz
            .z = -2.5;
        assert!(!needs_reference_swaps(&st));
    }
}
