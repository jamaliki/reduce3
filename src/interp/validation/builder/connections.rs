use super::*;

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
            super::restraints::link(graph, first, second, link, options);
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
                "Explicit covalent link lacks a validated chemical dictionary link",
            );
        }
    }
}
