use super::*;

pub(super) fn residue(
    graph: &mut ValidationGraph,
    r: &ResInterp,
    source: &ConfResidue,
    _ml: &MonLib,
    options: &ValidationOptions,
) {
    let Some(comp) = &r.comp else {
        graph.unsupported(
            ALL,
            format!("Unknown chemical component {}", source.resname),
        );
        return;
    };
    if comp.from_ccd {
        graph.unsupported(
            ALL,
            "CCD-derived restraints have not passed reference interpretation parity",
        );
    }
    let disulfide = comp.id == "CYS"
        && lookup(r, "SG").is_some_and(|sg| {
            graph
                .topology
                .bonds
                .iter()
                .chain(&graph.topology.sym_bonds)
                .any(|b| b.origin == ORIGIN_SS && (b.i == sg || b.j == sg))
        });
    let substituted_thiol = comp.id == "CYS"
        && lookup(r, "SG").is_some_and(|sg| {
            let cb = lookup(r, "CB");
            graph
                .topology
                .bonds
                .iter()
                .any(|b| (b.i == sg && Some(b.j) != cb) || (b.j == sg && Some(b.i) != cb))
        });
    let optional = |name: &str| {
        if lookup(r, name).is_some() {
            return false;
        }
        if (disulfide || substituted_thiol) && name == "HG" {
            return true;
        }
        let is_h = comp
            .atom(name)
            .is_some_and(|a| matches!(a.type_symbol.as_str(), "H" | "D"));
        is_h && (r.is_water
            || (comp.id == "HIS"
                && matches!(name, "HD1" | "HE2")
                && (lookup(r, "HD1").is_some() ^ lookup(r, "HE2").is_some())))
    };
    let known_missing = |names: &[&str]| {
        names
            .iter()
            .filter(|name| lookup(r, name).is_none())
            .all(|name| optional(name))
    };
    graph.dictionaries.push(DictionaryUse {
        residue: source.resname.clone(),
        source: comp.source.clone(),
        modifications: comp.applied_mods.clone(),
    });
    for atom in &comp.atoms {
        if matches!(atom.type_symbol.as_str(), "H" | "D") && lookup(r, &atom.id).is_none() {
            if optional(&atom.id) {
                continue;
            }
            graph.unresolved("hydrogens");
        }
    }
    for &atom in &source.atoms {
        let mapped = r
            .id_to_atom
            .iter()
            .find(|(_, i)| **i == atom)
            .and_then(|(id, _)| comp.atom(id));
        if mapped.is_some() {
            // CCTBX fixes energy types before applying polymer link mods.
            // Keep that map; the modified dictionary below supplies targets.
            if !matches!(graph.topology.etype[atom as usize], EType::Typed(_)) {
                graph.unsupported(&["atom_types"], "Dictionary atom lacks an energy type");
            }
        } else {
            graph.unsupported(
                ALL,
                format!("Unmapped or deleted observed atom {}", graph.topology.log),
            );
        }
    }
    for b in &comp.bonds {
        let (Some(i), Some(j)) = (lookup(r, &b.a1), lookup(r, &b.a2)) else {
            if known_missing(&[&b.a1, &b.a2]) {
                continue;
            }
            graph.unresolved("bonds");
            continue;
        };
        let target = if options.nuclear {
            b.value_dist_neutron.or(b.value_dist)
        } else {
            b.value_dist
        };
        if let (Some(target), Some(sigma)) = (target, b.esd) {
            graph.bond([i, j], target, sigma, ORIGIN_COVALENT);
        } else {
            graph.unsupported(&["bonds"], "Undefined dictionary bond targets");
        }
    }
    for a in &comp.angles {
        let indices = [
            lookup_angle(r, &a.a1),
            lookup_angle(r, &a.a2),
            lookup_angle(r, &a.a3),
        ];
        if indices.iter().any(Option::is_none) {
            if known_missing(&[&a.a1, &a.a2, &a.a3]) {
                continue;
            }
            graph.unresolved("angles");
            continue;
        }
        if let (Some(target), Some(sigma)) = (a.value, a.esd) {
            graph.angle(indices.map(Option::unwrap), target, sigma);
        } else {
            graph.unsupported(&["angles"], "Undefined dictionary angle targets");
        }
    }
    for chir in &comp.chirs {
        let ids = [
            lookup(r, &chir.centre),
            lookup(r, &chir.a[0]),
            lookup(r, &chir.a[1]),
            lookup(r, &chir.a[2]),
        ];
        if ids.iter().any(Option::is_none) {
            graph.unresolved("chiralities");
            continue;
        }
        let Some(ideal) = chirality::target(comp, chir) else {
            graph.unsupported(&["chiralities"], "Undefined chiral target geometry");
            continue;
        };
        let sign = chir.volume_sign.to_ascii_lowercase();
        if !["pos", "neg", "both"].iter().any(|s| sign.starts_with(s)) {
            graph.unsupported(&["chiralities"], "Unknown chiral volume sign");
            continue;
        }
        graph.chirality(
            ids.map(Option::unwrap),
            ideal,
            1. / options.chirality_sigma.powi(2),
            sign.starts_with("both"),
            chirality::pseudo(&comp.id, &chir.centre),
        );
    }
    let mut planes: Vec<(String, Vec<(u32, f64)>)> = Vec::new();
    for p in &comp.planes {
        let Some(atom) = lookup(r, &p.atom) else {
            if optional(&p.atom) {
                continue;
            }
            graph.unresolved("planarities");
            continue;
        };
        let Some(sigma) = p.esd.filter(|v| v.is_finite() && *v > 0.) else {
            graph.unsupported(&["planarities"], "Undefined plane target uncertainty");
            continue;
        };
        match planes.iter_mut().find(|row| row.0 == p.plane_id) {
            Some(row) => row.1.push((atom, sigma)),
            None => planes.push((p.plane_id.clone(), vec![(atom, sigma)])),
        }
    }
    for (_, atoms) in planes {
        graph.plane(atoms);
    }
}
