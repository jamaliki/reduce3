use super::*;

mod automatic;
mod chains;
mod connections;
mod enol;
mod link;
mod metals;
mod restraints;

pub(super) const ALL: &[&str] = &[
    "bonds",
    "angles",
    "chiralities",
    "planarities",
    "hydrogens",
    "atom_types",
];

pub(super) fn build(
    st: &Structure,
    flat: &FlatAtoms,
    ml: &MonLib,
    options: &ValidationOptions,
) -> Result<ValidationGraph, String> {
    let input = FlatAtoms::from_structure(st);
    if flat.path != input.path
        || flat.pos != input.pos
        || flat.name != input.name
        || flat.element != input.element
        || flat.altloc != input.altloc
        || flat.is_h != input.is_h
        || flat
            .pos
            .iter()
            .any(|p| !p.x.is_finite() || !p.y.is_finite() || !p.z.is_finite())
    {
        return Err("Validation sites and atom order must match the input structure".into());
    }
    if !options.link_distance_cutoff.is_finite()
        || options.link_distance_cutoff <= 0.
        || !options.chirality_sigma.is_finite()
        || options.chirality_sigma <= 0.
    {
        return Err("Validation target parameters must be finite and positive".into());
    }
    let params = InterpParams {
        neutron: options.nuclear,
        link_distance_cutoff: options.link_distance_cutoff,
        compat: false,
        auto_comps: FxHashMap::default(),
    };
    let topology = interpret_validation(st, flat, ml, &params);
    let mut graph = ValidationGraph {
        sites: flat.pos.clone(),
        atom_paths: flat.path.clone(),
        hydrogen_flags: flat.is_h.clone(),
        typed_charges: flat
            .path
            .iter()
            .map(|&p| charges::normalized(st.atom(p).elem(), st.atom(p).charge_value()))
            .collect(),
        bonds: Vec::new(),
        symmetry_bonds: Vec::new(),
        angles: Vec::new(),
        chiralities: Vec::new(),
        planes: Vec::new(),
        dictionaries: Vec::new(),
        topology,
        bond_keys: FxHashMap::default(),
        angle_keys: FxHashMap::default(),
        chiral_keys: FxHashMap::default(),
        plane_keys: FxHashMap::default(),
        coverage: ALL
            .iter()
            .map(|&check| Coverage {
                check,
                unsupported: Vec::new(),
                unresolved: 0,
            })
            .collect(),
    };
    if st.models.len() != 1 {
        graph.unsupported(ALL, "Multiple models require separate validation views");
    }
    symmetry_bonds(&mut graph, ml);
    let index: FxHashMap<AtomPath, u32> = flat
        .path
        .iter()
        .enumerate()
        .map(|(i, &p)| (p, i as u32))
        .collect();
    let chains = conformer_residues(st, Some(&index));
    let mut all_residues = Vec::new();
    let mut prepared_chains = Vec::new();
    for model in &chains {
        for chain in model {
            let prepared = chains::prepare(&mut graph, chain, flat, ml, options);
            let start = all_residues.len();
            all_residues.extend(prepared);
            prepared_chains.push((start..all_residues.len(), chain));
        }
    }
    let links = connections::prepare_declared(&mut graph, &mut all_residues, ml, options)?;
    connections::emit_declared(&mut graph, &all_residues, &links, ml, options);
    automatic::prepare(&mut graph, &mut all_residues, ml);
    metals::emit(&mut graph, st, flat, ml, options);
    for (range, sources) in prepared_chains {
        chains::emit(
            &mut graph,
            &mut all_residues[range],
            sources,
            flat,
            ml,
            options,
        );
    }
    automatic::emit(&mut graph, &all_residues, ml, st, flat, options);
    connections::disulfides(&mut graph, &all_residues, ml, options);
    connections::explicit_links(&mut graph, st, flat);
    let missing = graph.topology.bonds.iter().any(|bond| {
        !graph
            .bond_keys
            .contains_key(&[bond.i.min(bond.j), bond.i.max(bond.j)])
    });
    if missing {
        graph.unsupported(
            ALL,
            "Automatic or metal/disulfide links need reference target weights",
        );
    }
    Ok(graph)
}

fn symmetry_bonds(graph: &mut ValidationGraph, ml: &MonLib) {
    let images = graph.topology.symmetry_images.clone();
    for image in images {
        let Some(definition) = ml.link("SS").and_then(|link| link.bonds.first()) else {
            graph.unsupported(ALL, "Missing symmetry disulfide dictionary target");
            continue;
        };
        let Some((ideal, sigma)) = definition.value_dist.zip(definition.esd) else {
            graph.unsupported(ALL, "Undefined symmetry disulfide target uncertainty");
            continue;
        };
        let Some(weight) = bond_weight(ideal,sigma) else {
            graph.unsupported(&["bonds", "hydrogens"], "Invalid symmetry disulfide target uncertainty");
            continue;
        };
        if image.bond.i != image.bond.j {
            let Some(reverse) = image.image.inverse() else {
                graph.unsupported(ALL, "Symmetry bond operator is singular");
                continue;
            };
            graph.symmetry_bonds.push(WeightedSymmetryBond {
                indices: [image.bond.j, image.bond.i],
                image: reverse,
                ideal,
                weight,
                origin: ORIGIN_SS,
            });
        }
        graph.symmetry_bonds.push(WeightedSymmetryBond {
            indices: [image.bond.i, image.bond.j],
            image: image.image,
            ideal,
            weight,
            origin: ORIGIN_SS,
        });
    }
    if graph
        .topology
        .sym_bonds
        .iter()
        .any(|bond| bond.origin != ORIGIN_SS)
    {
        graph.unsupported(
            ALL,
            "Symmetry metal links need weighted image interpretation",
        );
    }
}
