use super::*;

pub(super) fn prepare_declared(
    graph: &mut ValidationGraph,
    residues: &mut [ResInterp],
    ml: &MonLib,
    options: &ValidationOptions,
) -> Result<Vec<(usize, usize, String)>, String> {
    let mut result = Vec::new();
    for request in &options.dictionary_links {
        if request
            .atom_indices
            .iter()
            .any(|&i| i as usize >= graph.sites.len())
        {
            return Err("A dictionary link atom index is outside the model".into());
        }
        let find = |atom| {
            residues
                .iter()
                .enumerate()
                .filter_map(|(i, r)| r.id_to_atom.values().any(|&a| a == atom).then_some(i))
                .collect::<Vec<_>>()
        };
        let (first, second) = (find(request.atom_indices[0]), find(request.atom_indices[1]));
        if first.len() != 1 || second.len() != 1 || first[0] == second[0] {
            graph.unsupported(
                ALL,
                "A declared dictionary link needs two unique residue conformers",
            );
            continue;
        }
        let Some(link) = ml.link(&request.id) else {
            return Err(format!("Unknown dictionary link {}", request.id));
        };
        if !link_matches(
            link,
            &residues[first[0]],
            &residues[second[0]],
            request.atom_indices[0],
            request.atom_indices[1],
        ) {
            return Err(format!(
                "Dictionary link {} does not match the selected residues and atoms",
                request.id
            ));
        }
        if result
            .iter()
            .any(|(a, b, id)| *a == first[0] && *b == second[0] && *id == request.id)
        {
            return Err("Duplicate dictionary link application".into());
        }
        for (index, id) in [(first[0], &link.mod_id_1), (second[0], &link.mod_id_2)] {
            if id.is_empty() {
                continue;
            }
            apply_modification(graph, &mut residues[index], ml, id)?;
        }
        result.push((first[0], second[0], request.id.clone()));
    }
    Ok(result)
}

fn apply_modification(
    graph: &mut ValidationGraph,
    residue: &mut ResInterp,
    ml: &MonLib,
    id: &str,
) -> Result<(), String> {
    let base = residue.comp.as_ref().unwrap();
    let Some(comp) = ml.apply_mod(base, id) else {
        return Err(format!("Unknown dictionary link modification {id}"));
    };
    let comp = Arc::new(comp);
    for (name, &atom) in &residue.id_to_atom {
        if let Some(definition) = comp.atom(name) {
            graph.topology.etype[atom as usize] = definition
                .type_energy
                .as_ref()
                .map(|value| EType::Typed(value.into()))
                .unwrap_or(EType::NoneType);
            graph.topology.dict[atom as usize] = AtomDictInfo {
                comp: Some(Arc::clone(&comp)),
                dict_id: Some(name.into()),
                instance: 0,
            };
        }
    }
    residue.comp = Some(comp);

    Ok(())
}

pub(super) fn emit_declared(
    graph: &mut ValidationGraph,
    residues: &[ResInterp],
    links: &[(usize, usize, String)],
    ml: &MonLib,
    options: &ValidationOptions,
) {
    for (first, second, id) in links {
        let link = ml.link(id).unwrap();
        let (angles, chiralities) = (graph.angles.len(), graph.chiralities.len());
        super::link::emit(graph, &residues[*first], &residues[*second], link, options);
        for angle in &mut graph.angles[angles..] {
            angle.origin = ORIGIN_MISC;
        }
        for chirality in &mut graph.chiralities[chiralities..] {
            chirality.origin = ORIGIN_MISC;
        }
        let side = |c| {
            if c == 1 {
                &residues[*first]
            } else {
                &residues[*second]
            }
        };
        for bond in &link.bonds {
            let (Some(i), Some(j)) = (
                lookup(side(bond.c1), &bond.a1),
                lookup(side(bond.c2), &bond.a2),
            ) else {
                continue;
            };
            if let Some(&index) = graph.bond_keys.get(&[i.min(j), i.max(j)]) {
                graph.bonds[index].origin = ORIGIN_MISC;
                graph
                    .topology
                    .upsert_bond(i, j, graph.bonds[index].ideal, ORIGIN_MISC);
            }
        }
    }
}

pub(super) fn disulfides(
    graph: &mut ValidationGraph,
    residues: &[ResInterp],
    ml: &MonLib,
    options: &ValidationOptions,
) {
    let pairs: Vec<(u32, u32)> = graph
        .topology
        .bonds
        .iter()
        .filter(|b| b.origin == ORIGIN_SS)
        .map(|b| (b.i, b.j))
        .collect();
    if pairs.is_empty() {
        return;
    }
    let Some(link) = ml.link("SS") else {
        graph.unsupported(ALL, "Disulfide dictionary link is missing");
        return;
    };
    for (a, b) in pairs {
        let first = residues.iter().find(|r| lookup(r, "SG") == Some(a));
        let second = residues.iter().find(|r| lookup(r, "SG") == Some(b));
        if let (Some(first), Some(second)) = (first, second) {
            let before = graph.angles.len();
            super::link::emit(graph, first, second, link, options);
            for angle in &mut graph.angles[before..] {
                angle.origin = ORIGIN_SS;
            }
            let key = [a.min(b), a.max(b)];
            if let Some(&k) = graph.bond_keys.get(&key) {
                graph.bonds[k].origin = ORIGIN_SS;
                graph
                    .topology
                    .upsert_bond(a, b, graph.bonds[k].ideal, ORIGIN_SS);
            }
        } else {
            graph.unsupported(ALL, "Disulfide atom mapping is unresolved");
        }
    }
}

pub(super) fn explicit_links(graph: &mut ValidationGraph, st: &Structure, flat: &FlatAtoms) {
    for link in &st.links {
        let find = |label: &AtomLabel| {
            flat.path
                .iter()
                .enumerate()
                .filter_map(|(i, &p)| {
                    let ag = st.atom_group(p);
                    let chain = &st.models[p.model as usize].chains[p.chain as usize];
                    let residue = &chain.residue_groups[p.rg as usize];
                    (flat.name[i].trim() == label.name.trim()
                        && ag.resname.trim() == label.resname.trim()
                        && ag.altloc.trim() == label.altloc.trim()
                        && chain.id == label.chain
                        && residue.resseq.trim() == label.resseq.trim()
                        && residue.icode.trim() == label.icode.trim())
                    .then_some(i as u32)
                })
                .collect::<Vec<_>>()
        };
        let (left, right) = (find(&link.atom1), find(&link.atom2));
        if left.len() != 1 || right.len() != 1 {
            graph.unsupported(ALL, "Explicit link atom identity is missing or ambiguous");
            continue;
        }
        let key = [left[0].min(right[0]), left[0].max(right[0])];
        if !graph.bond_keys.contains_key(&key) {
            graph.unsupported(
                ALL,
                "Explicit links need a declared chemical dictionary application",
            );
        }
    }
}

fn link_matches(
    link: &crate::monlib::ChemLink,
    first: &ResInterp,
    second: &ResInterp,
    a: u32,
    b: u32,
) -> bool {
    let Some((left, right)) = first.comp.as_ref().zip(second.comp.as_ref()) else {
        return false;
    };
    if link_match_one(&link.comp_id_1, &link.group_comp_1, &left.id, &left.group).is_none()
        || link_match_one(&link.comp_id_2, &link.group_comp_2, &right.id, &right.group).is_none()
    {
        return false;
    }
    let side = |c: u8| if c == 1 { first } else { second };
    link.bonds.iter().any(|bond| {
        match (
            lookup(side(bond.c1), &bond.a1),
            lookup(side(bond.c2), &bond.a2),
        ) {
            (Some(i), Some(j)) => (i == a && j == b) || (i == b && j == a),
            _ => false,
        }
    })
}
