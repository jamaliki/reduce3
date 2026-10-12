use super::*;
mod candidates;
mod images;
mod zinc_targets;
use candidates::{Candidate, group};

/// Weighted targets from CCTBX's Metal Coordination Library.
pub(super) fn emit(
    graph: &mut ValidationGraph,
    st: &Structure,
    flat: &FlatAtoms,
    ml: &MonLib,
    options: &ValidationOptions,
) {
    let candidates = match candidates::collect(graph, st, flat, ml, options) {
        Ok(pairs) => pairs,
        Err(reason) => {
            graph.unsupported(ALL, reason);
            return;
        }
    };
    let zincs = zinc_partners(&candidates, flat);
    let clusters = cluster_partners(&candidates, st, flat);
    iron_sulfur(graph, st, flat, &clusters);
    for (metal, others) in zincs {
        zinc(graph, st, flat, metal, &others);
    }
}

fn zinc_partners(candidates: &[Candidate], flat: &FlatAtoms) -> FxHashMap<u32, Vec<Candidate>> {
    let mut zincs: FxHashMap<u32, Vec<Candidate>> = FxHashMap::default();
    for &pair in candidates {
        let [i, j] = pair.indices.map(|index| index as usize);
        if pair.distance > 3. {
            continue;
        }
        for (metal, other) in [(i, j), (j, i)] {
            if flat.name[metal].trim() != "ZN"
                || flat.element[metal] != "ZN"
                || !matches!(flat.name[other].trim(), "SG" | "ND1" | "NE2")
            {
                continue;
            }
            let donors = zincs.entry(metal as u32).or_default();
            if !donors
                .iter()
                .any(|existing| group(flat, existing.indices[1]) == group(flat, other as u32))
            {
                donors.push(directed(pair, metal as u32));
            }
        }
    }
    zincs
}

fn cluster_partners(candidates: &[Candidate], st: &Structure, flat: &FlatAtoms) -> Vec<(u32, u32)> {
    let mut clusters = Vec::new();
    let mut cluster_residues = FxHashSet::default();
    for pair in candidates {
        let [i, j] = pair.indices.map(|index| index as usize);
        if cluster_name(st, flat, i).is_some() && cluster_name(st, flat, j).is_some() {
            // Different cluster names also do not define an external partner.
            if cluster_name(st, flat, i) != cluster_name(st, flat, j) {
                continue;
            }
        }
        let cluster = if cluster_name(st, flat, i).is_some() {
            i
        } else if cluster_name(st, flat, j).is_some() {
            j
        } else {
            continue;
        };
        // CCTBX clears all cluster partners on any symmetry contact, before
        // it tests whether the candidate is a coordinating S or N atom.
        if !pair.image.is_identity() {
            clusters.clear();
            break;
        }
        let other = if cluster == i { j } else { i };
        if flat.element[cluster] == "FE"
            && matches!(flat.element[other].as_str(), "S" | "N")
            && cluster_name(st, flat, other).is_none()
            && cluster_residues.insert(group(flat, other as u32))
        {
            clusters.push((cluster as u32, other as u32));
        }
    }
    clusters
}

fn directed(mut pair: Candidate, metal: u32) -> Candidate {
    if pair.indices[0] != metal {
        pair.indices.swap(0, 1);
        pair.image = pair.image.inverse().unwrap();
    }
    pair
}

fn cluster_name(st: &Structure, flat: &FlatAtoms, i: usize) -> Option<&'static str> {
    match st.atom_group(flat.path[i]).resname.trim() {
        "SF4" => Some("SF4"),
        "FES" => Some("FES"),
        "F3S" => Some("F3S"),
        _ => None,
    }
}

fn zinc(
    graph: &mut ValidationGraph,
    st: &Structure,
    flat: &FlatAtoms,
    metal: u32,
    others: &[Candidate],
) {
    let count = |name| {
        others
            .iter()
            .filter(|pair| {
                st.atom_group(flat.path[pair.indices[1] as usize])
                    .resname
                    .trim()
                    == name
            })
            .count()
    };
    let counts = (count("CYS"), count("HIS"));
    let target = |names: &[&str]| {
        zinc_targets::TARGETS
            .iter()
            .find(|(key, atoms, _, _)| *key == counts && *atoms == names)
            .map(|(_, _, target, sigma)| (*target, *sigma))
    };
    for (k, pair) in others.iter().enumerate() {
        let other = pair.indices[1];
        let name = flat.name[other as usize].trim();
        if let Some((target, sigma)) = if others.len() < 4 {
            Some((2.3, 0.03))
        } else {
            target(&["ZN", name])
        } {
            if pair.image.is_identity() {
                graph.bond([metal, other], target, sigma, ORIGIN_METAL);
            } else {
                graph.image_bond([metal, other], pair.image, target, sigma, ORIGIN_METAL);
            }
        }
        for previous_pair in &others[..k] {
            let previous = previous_pair.indices[1];
            let previous_name = flat.name[previous as usize].trim();
            if others.len() < 4 {
                continue;
            }
            if let Some((target, sigma)) = target(&[name, "ZN", previous_name]) {
                graph.angle_with_origin([previous, metal, other], target, sigma, ORIGIN_METAL);
            }
        }
    }
}

type ClusterTarget = (f64, f64, Option<(f64, f64)>);

fn cluster_targets(cluster: &str, ligand: &str, element: &str) -> Option<ClusterTarget> {
    match (cluster, ligand, element) {
        ("SF4", "HIS", "N") => Some((2.04, 0.05, None)),
        ("FES", "HIS", "N") => Some((2.14, 0.05, None)),
        ("SF4", "MET", "S") => Some((2.311, 0.012, Some((113.97, 17.528)))),
        ("SF4", _, "S") => Some((2.268, 0.034, Some((114.24, 11.5)))),
        ("FES", _, "S") => Some((2.305, 0.044, Some((111.20, 8.1)))),
        ("F3S", _, "S") => Some((2.318, 0.016, Some((112.23, 12.06)))),
        _ => None,
    }
}

fn iron_sulfur(
    graph: &mut ValidationGraph,
    st: &Structure,
    flat: &FlatAtoms,
    pairs: &[(u32, u32)],
) {
    for &(metal, other) in pairs {
        let m = metal as usize;
        let o = other as usize;
        let cluster = cluster_name(st, flat, m).unwrap();
        let ligand = st.atom_group(flat.path[o]).resname.trim();
        let Some((target, sigma, angle)) = cluster_targets(cluster, ligand, &flat.element[o])
        else {
            continue;
        };
        graph.bond([metal, other], target, sigma, ORIGIN_METAL);
        let Some((target, sigma)) = angle else {
            continue;
        };
        let Some(number) = flat.name[m]
            .trim()
            .chars()
            .last()
            .and_then(|c| c.to_digit(10))
        else {
            graph.unsupported(ALL, "Iron-sulfur cluster metal name lacks its site number");
            continue;
        };
        let excluded = if cluster == "F3S" {
            match number {
                1 => 4,
                3 => 2,
                4 => 1,
                _ => 0,
            }
        } else {
            number
        };
        let group = st.atom_group(flat.path[m]);
        for (k, atom) in group.atoms.iter().enumerate() {
            if atom.elem() != "S" || atom.name.trim() == format!("S{excluded}") {
                continue;
            }
            let mut atom_path = flat.path[m];
            atom_path.atom = k as u32;
            if let Some(index) = flat.path.iter().position(|&p| p == atom_path) {
                graph.angle_with_origin([index as u32, metal, other], target, sigma, ORIGIN_METAL);
            }
        }
        if cluster == "FES" && ligand == "CYS" {
            for &(_, partner) in pairs {
                if partner == other
                    || flat.name[partner as usize].trim() != "SG"
                    || flat.pos[partner as usize].dist(flat.pos[m]) >= 3.5
                {
                    continue;
                }
                graph.angle_with_origin([partner, metal, other], 107.77, 8.16, ORIGIN_METAL);
            }
        }
    }
}
