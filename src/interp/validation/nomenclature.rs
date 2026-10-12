use super::*;

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
        st.models[0].chains[0].residue_groups[0].atom_groups[0].atoms[3]
            .xyz
            .z = -2.5;
        assert!(!needs_reference_swaps(&st));
    }
}
