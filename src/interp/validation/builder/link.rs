use super::*;
use crate::monlib::ChemLink;

pub(super) fn emit(
    graph: &mut ValidationGraph,
    previous: &ResInterp,
    current: &ResInterp,
    link: &ChemLink,
    options: &ValidationOptions,
) {
    bonds(graph, previous, current, link);
    angles(graph, previous, current, link);
    chiralities(graph, previous, current, link, options);
    planes(graph, previous, current, link);
}

fn bonds(graph: &mut ValidationGraph, previous: &ResInterp, current: &ResInterp, link: &ChemLink) {
    let side = |c| if c == 1 { previous } else { current };
    for b in &link.bonds {
        if let (Some(i), Some(j), Some(ideal), Some(sigma)) = (
            lookup(side(b.c1), &b.a1),
            lookup(side(b.c2), &b.a2),
            b.value_dist,
            b.esd,
        ) {
            graph.bond([i, j], ideal, sigma, ORIGIN_COVALENT);
        } else {
            graph.unresolved("bonds");
        }
    }
}

fn angles(graph: &mut ValidationGraph, previous: &ResInterp, current: &ResInterp, link: &ChemLink) {
    let side = |c| if c == 1 { previous } else { current };
    for a in &link.angles {
        let ids = std::array::from_fn::<_, 3, _>(|i| lookup_angle(side(a.c[i]), &a.a[i]));
        if let (Some(ideal), Some(sigma)) = (a.value, a.esd) {
            if ids.iter().all(Option::is_some) {
                graph.angle(ids.map(Option::unwrap), ideal, sigma);
            } else {
                let missing_enol = lookup(previous, "HNO").is_some()
                    && a.c
                        .iter()
                        .zip(&a.a)
                        .zip(&ids)
                        .filter(|(_, id)| id.is_none())
                        .all(|((&side, name), _)| side == 2 && name == "H");
                if !missing_enol {
                    graph.unresolved("angles");
                }
            }
        } else {
            graph.unsupported(&["angles"], "Undefined link angle uncertainty");
        }
    }
}

fn chiralities(
    graph: &mut ValidationGraph,
    previous: &ResInterp,
    current: &ResInterp,
    link: &ChemLink,
    options: &ValidationOptions,
) {
    let side = |c| if c == 1 { previous } else { current };
    for chir in &link.chirs {
        let sign = chir.volume_sign.to_ascii_lowercase();
        if !["pos", "neg", "both"]
            .iter()
            .any(|prefix| sign.starts_with(prefix))
        {
            graph.unsupported(&["chiralities"], "Unknown cross-residue chiral volume sign");
            continue;
        }
        let ids = std::array::from_fn::<_, 4, _>(|i| lookup(side(chir.c[i]), &chir.a[i]));
        if ids.iter().any(Option::is_none) {
            graph.unresolved("chiralities");
            continue;
        }
        let Some((first, second)) = previous.comp.as_ref().zip(current.comp.as_ref()) else {
            graph.unresolved("chiralities");
            continue;
        };
        let Some(ideal) = chirality::link_target(first, second, link, chir) else {
            graph.unsupported(&["chiralities"], "Undefined cross-residue chiral geometry");
            continue;
        };
        graph.chirality(
            ids.map(Option::unwrap),
            ideal,
            1. / options.chirality_sigma.powi(2),
            sign.starts_with("both"),
            false,
        );
    }
}

fn planes(graph: &mut ValidationGraph, previous: &ResInterp, current: &ResInterp, link: &ChemLink) {
    let side = |c| if c == 1 { previous } else { current };
    let mut planes: Vec<(String, Vec<(u32, f64)>)> = Vec::new();
    for p in &link.planes {
        let Some(atom) = lookup(side(p.c), &p.atom) else {
            if p.c == 2 && p.atom == "H" && lookup(previous, "HNO").is_some() {
                continue;
            }
            graph.unresolved("planarities");
            continue;
        };
        let Some(sigma) = p.esd.filter(|v| v.is_finite() && *v > 0.) else {
            graph.unsupported(&["planarities"], "Undefined link plane uncertainty");
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
