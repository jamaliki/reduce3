use super::*;

mod metal_defaults;

/// Remove only absent leaving atoms named by the selected link's mods.
pub(super) fn prepare(graph: &mut ValidationGraph, residues: &mut [ResInterp], ml: &MonLib) {
    for application in graph.topology.automatic_links.clone() {
        let (id, anchors) = match application {
            AutomaticLink::Dictionary { id, anchors, .. } => (id, anchors),
            AutomaticLink::Glycosidic { id, bond, .. } => (id, [bond[1], bond[0]]),
            _ => continue,
        };
        let Some(link) = ml.link(&id) else {
            continue;
        };
        for (atom, modification) in anchors.into_iter().zip([&link.mod_id_1, &link.mod_id_2]) {
            if modification.is_empty() {
                continue;
            }
            let Some(residue) = residues
                .iter_mut()
                .find(|r| r.id_to_atom.values().any(|&i| i == atom))
            else {
                continue;
            };
            let Some(comp) = residue.comp.as_ref() else {
                continue;
            };
            let Some(modified) = ml.apply_mod(comp, modification) else {
                graph.unsupported(
                    ALL,
                    format!("Missing automatic link modification {modification}"),
                );
                continue;
            };
            if comp.atoms.iter().all(|a| modified.atom(&a.id).is_some()) {
                continue;
            }
            // Automatic link discovery must not delete an observed atom.
            if comp
                .atoms
                .iter()
                .any(|a| modified.atom(&a.id).is_none() && lookup(residue, &a.id).is_some())
            {
                graph.unsupported(ALL, "Automatic link leaving atom is still observed");
                continue;
            }
            let deleted: FxHashSet<_> = comp
                .atoms
                .iter()
                .filter(|a| modified.atom(&a.id).is_none())
                .map(|a| a.id.clone())
                .collect();
            // CCTBX emits monomer targets before automatic linking. Preserve
            // those targets; use the mod only to identify absent leaving atoms.
            let mut retained = (**comp).clone();
            retained.atoms.retain(|a| !deleted.contains(&a.id));
            retained
                .bonds
                .retain(|b| !deleted.contains(&b.a1) && !deleted.contains(&b.a2));
            retained.angles.retain(|a| {
                !deleted.contains(&a.a1) && !deleted.contains(&a.a2) && !deleted.contains(&a.a3)
            });
            retained.chirs.retain(|c| {
                !deleted.contains(&c.centre) && c.a.iter().all(|a| !deleted.contains(a))
            });
            retained.planes.retain(|p| !deleted.contains(&p.atom));
            retained
                .applied_mods
                .push(format!("{modification}:absent-leaving-atoms"));
            residue.comp = Some(Arc::new(retained));
        }
    }
}

/// Automatic links use the selected dictionary; they do not apply its mods.
pub(super) fn emit(
    graph: &mut ValidationGraph,
    residues: &[ResInterp],
    ml: &MonLib,
    st: &Structure,
    flat: &FlatAtoms,
    options: &ValidationOptions,
) {
    for application in graph.topology.automatic_links.clone() {
        match application {
            AutomaticLink::Custom { indices, origin } => {
                if graph
                    .bond_keys
                    .contains_key(&[indices[0].min(indices[1]), indices[0].max(indices[1])])
                {
                    continue;
                }
                let (target, sigma) = default_target(st, flat, indices);
                graph.bond(indices, target, sigma, origin);
            }
            AutomaticLink::Dictionary {
                id,
                anchors,
                origin,
            } => {
                let find = |atom| {
                    residues
                        .iter()
                        .find(|r| r.id_to_atom.values().any(|&i| i == atom))
                };
                let Some((first, second)) = find(anchors[0]).zip(find(anchors[1])) else {
                    graph.unsupported(ALL, "Automatic dictionary link atom mapping is unresolved");
                    continue;
                };
                let Some(definition) = ml.link(&id) else {
                    graph.unsupported(ALL, format!("Missing automatic link dictionary {id}"));
                    continue;
                };
                let (angles, chiralities) = (graph.angles.len(), graph.chiralities.len());
                link::emit(graph, first, second, definition, options);
                for angle in &mut graph.angles[angles..] {
                    angle.origin = origin;
                }
                for chirality in &mut graph.chiralities[chiralities..] {
                    chirality.origin = origin;
                }
                for bond in &definition.bonds {
                    let side = |c| if c == 1 { first } else { second };
                    if let Some((i, j)) =
                        lookup(side(bond.c1), &bond.a1).zip(lookup(side(bond.c2), &bond.a2))
                    {
                        if let Some(&k) = graph.bond_keys.get(&[i.min(j), i.max(j)]) {
                            graph.bonds[k].origin = origin;
                            graph
                                .topology
                                .upsert_bond(i, j, graph.bonds[k].ideal, origin);
                        }
                    }
                }
            }
            AutomaticLink::Glycosidic {
                id: _,
                bond,
                angles,
                chirality,
                origin,
            } => {
                graph.bond(bond, 1.439, 0.02, origin);
                for angle in angles {
                    graph.angle_with_origin(angle.i, angle.ideal, angle.esd, origin);
                }
                if let Some((indices, target)) = chirality {
                    let before = graph.chiralities.len();
                    graph.chirality(indices, target, 1. / 0.02_f64.powi(2), false, false);
                    for chirality in &mut graph.chiralities[before..] {
                        chirality.origin = origin;
                    }
                }
            }
        }
    }
}

fn default_target(st: &Structure, flat: &FlatAtoms, [first, second]: [u32; 2]) -> (f64, f64) {
    for (metal, other) in [
        (first as usize, second as usize),
        (second as usize, first as usize),
    ] {
        let (m, e) = (flat.element[metal].as_str(), flat.element[other].as_str());
        let residue = st.atom_group(flat.path[other]).resname.trim();
        if e == "O"
            && flat.name[other].trim() == "O"
            && resclass::get_class(residue) == ResClass::CommonAminoAcid
        {
            if let Some(&(_, _, target, sigma)) = metal_defaults::CARBONYL
                .iter()
                .find(|&&(a, b, _, _)| a == m && b == e)
            {
                return (target, sigma);
            }
        }
        if let Some(&(_, _, _, target, sigma)) = metal_defaults::RESIDUE
            .iter()
            .find(|&&(a, b, r, _, _)| a == m && b == e && r == residue)
        {
            return (target, sigma);
        }
        if let Some(&(_, _, target, sigma)) = metal_defaults::OTHER
            .iter()
            .find(|&&(a, b, _, _)| a == m && b == e)
        {
            return (target, sigma);
        }
    }
    let element = |i: u32| {
        let value = &flat.element[i as usize];
        if value.len() < 2 {
            value.to_string()
        } else {
            format!("{}{}", &value[..1], value[1..].to_ascii_lowercase())
        }
    };
    let (a, b) = (element(first), element(second));
    for (a, b, order) in [(&a, &b, 1), (&b, &a, 1), (&a, &b, 0), (&b, &a, 0)] {
        if let Some(&(_, _, _, target)) = crate::linkdata::QM_DEFAULTS
            .iter()
            .find(|&&(x, y, o, _)| x == a && y == b && o == order)
        {
            return (target, 0.01);
        }
    }
    (2.3, 0.02)
}
