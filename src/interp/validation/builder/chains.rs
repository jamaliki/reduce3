use super::*;

pub(super) fn prepare(
    graph: &mut ValidationGraph,
    chain: &[ConfResidue],
    flat: &FlatAtoms,
    ml: &MonLib,
    options: &ValidationOptions,
) -> Vec<ResInterp> {
    let mut log = String::new();
    let mut interpreted: Vec<ResInterp> = chain
        .iter()
        .map(|r| {
            interpret_residue(
                ml,
                &FxHashMap::default(),
                flat,
                r,
                options.nuclear,
                &mut log,
            )
            .0
        })
        .collect();
    nucleotide_targets(graph, &mut interpreted, chain, flat, ml);
    for (residue, source) in interpreted.iter_mut().zip(chain) {
        if names::l_given_d(source.resname.trim()).is_some() {
            let mut mods = vec!["PEPT-D"];
            if source.resname.trim() == "DIL" {
                mods.push("DIL_chir_02_both");
            }
            if source.resname.trim() == "DTH" {
                mods.push("DTH_chir_02_both");
            }
            for id in mods {
                if let Some(comp) = &residue.comp {
                    if let Some(modified) = ml.apply_mod(comp, id) {
                        residue.comp = Some(Arc::new(modified));
                    } else {
                        graph.unsupported(ALL, format!("Missing D-amino acid modification {id}"));
                    }
                }
            }
        }
    }
    super::enol::prepare(graph, &mut interpreted, chain, flat);
    interpreted
}

pub(super) fn emit(
    graph: &mut ValidationGraph,
    residues: &mut [ResInterp],
    sources: &[ConfResidue],
    flat: &FlatAtoms,
    ml: &MonLib,
    options: &ValidationOptions,
) {
    if let Some(first) = residues.first() {
        super::restraints::residue(graph, first, &sources[0], ml, options);
    }
    for current in 1..residues.len() {
        if let Some(id) =
            apply_polymer_modification(graph, residues, sources, flat, ml, options, current)
        {
            super::link::emit(
                graph,
                &residues[current - 1],
                &residues[current],
                ml.link(&id).unwrap(),
                options,
            );
        }
        restraints::residue(graph, &residues[current], &sources[current], ml, options);
    }
    if options.use_cdl {
        super::cdl::apply(graph, residues, sources);
    }
    super::enol::links(graph, residues, sources);
}

fn nucleotide_targets(
    graph: &mut ValidationGraph,
    residues: &mut [ResInterp],
    chain: &[ConfResidue],
    flat: &FlatAtoms,
    ml: &MonLib,
) {
    for k in 0..residues.len() {
        if !residues[k].is_rna_dna || lookup(&residues[k], "O2'").is_none() {
            continue;
        }
        let next_p = residues.get(k + 1).and_then(|r| lookup(r, "P"));
        let Some(two) = rna_pucker(flat, &residues[k], next_p) else {
            continue;
        };
        residues[k].is_rna2p = Some(two);
        let base = names::rna_dna_mon_lib_name(chain[k].resname.trim(), true);
        let family = match base {
            Some("A" | "G") => "_pur",
            Some("C" | "U") => "_pyr",
            _ => "",
        };
        let id = format!("{}{}", if two { "rna2p" } else { "rna3p" }, family);
        if let Some(comp) = &residues[k].comp {
            if let Some(modified) = ml.apply_mod(comp, &id) {
                residues[k].comp = Some(Arc::new(modified));
            } else {
                graph.unsupported(ALL, format!("Missing nucleotide modification {id}"));
            }
        }
    }
}

fn apply_polymer_modification(
    graph: &mut ValidationGraph,
    residues: &mut [ResInterp],
    chain: &[ConfResidue],
    flat: &FlatAtoms,
    ml: &MonLib,
    options: &ValidationOptions,
    k: usize,
) -> Option<String> {
    if !chain[k].link_to_previous {
        return None;
    }
    let mut id = lib_link_id(ml, &residues[k - 1], &residues[k])?;
    if matches!(id.as_str(), "TRANS" | "PTRANS") {
        let omega = [
            lookup(&residues[k - 1], "CA"),
            lookup(&residues[k - 1], "C"),
            lookup(&residues[k], "N"),
            lookup(&residues[k], "CA"),
        ];
        if omega.iter().all(Option::is_some)
            && dihedral_model(&flat.pos, omega.map(Option::unwrap)).is_some_and(|v| v.abs() < 45.)
        {
            id = if id == "PTRANS" { "PCIS" } else { "CIS" }.into();
        }
    }
    let Some(link) = ml.link(&id) else {
        graph.unsupported(ALL, format!("Missing polymer link {id}"));
        return None;
    };
    let side = |c: u8| {
        if c == 1 {
            &residues[k - 1]
        } else {
            &residues[k]
        }
    };
    if link.bonds.iter().any(|bond| {
        match (
            lookup(side(bond.c1), &bond.a1),
            lookup(side(bond.c2), &bond.a2),
        ) {
            (Some(i), Some(j)) => {
                flat.pos[i as usize].dist(flat.pos[j as usize]) > options.link_distance_cutoff
            }
            _ => true,
        }
    }) {
        return None;
    }
    for (slot, id) in [(k - 1, &link.mod_id_1), (k, &link.mod_id_2)] {
        if id.is_empty() || id == "." || id == "?" {
            continue;
        }
        if let Some(comp) = &residues[slot].comp {
            if let Some(modified) = ml.apply_mod(comp, id) {
                residues[slot].comp = Some(Arc::new(modified));
            } else {
                graph.unsupported(ALL, format!("Missing link modification {id}"));
            }
        }
    }
    Some(id)
}
