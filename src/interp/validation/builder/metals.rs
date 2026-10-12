use super::*;

/// Weighted targets from CCTBX's Metal Coordination Library.
pub(super) fn emit(graph: &mut ValidationGraph, st: &Structure, flat: &FlatAtoms) {
    let pairs: Vec<_> = graph
        .topology
        .bonds
        .iter()
        .copied()
        .filter(|b| b.origin == ORIGIN_METAL)
        .collect();
    let mut zincs: FxHashMap<u32, Vec<u32>> = FxHashMap::default();
    let mut clusters = Vec::new();
    for bond in pairs {
        let (i, j) = (bond.i as usize, bond.j as usize);
        if !flat.altloc[i].is_empty()
            && !flat.altloc[j].is_empty()
            && flat.altloc[i] != flat.altloc[j]
        {
            graph.unsupported(
                ALL,
                "Metal candidate joins incompatible alternate conformers",
            );
            continue;
        }
        if flat.name[i].trim() == "ZN" && flat.element[i] == "ZN" {
            zincs.entry(bond.i).or_default().push(bond.j);
        } else if flat.name[j].trim() == "ZN" && flat.element[j] == "ZN" {
            zincs.entry(bond.j).or_default().push(bond.i);
        } else if cluster_name(st, flat, i).is_some() && flat.element[i] == "FE" {
            clusters.push((bond.i, bond.j));
        } else if cluster_name(st, flat, j).is_some() && flat.element[j] == "FE" {
            clusters.push((bond.j, bond.i));
        }
    }
    for (metal, others) in zincs {
        zinc(graph, st, flat, metal, &others);
    }
    iron_sulfur(graph, st, flat, &clusters);
}

fn cluster_name(st: &Structure, flat: &FlatAtoms, i: usize) -> Option<&'static str> {
    match st.atom_group(flat.path[i]).resname.trim() {
        "SF4" => Some("SF4"),
        "FES" => Some("FES"),
        "F3S" => Some("F3S"),
        _ => None,
    }
}

fn zinc(graph: &mut ValidationGraph, st: &Structure, flat: &FlatAtoms, metal: u32, others: &[u32]) {
    if others
        .iter()
        .any(|&i| !matches!(flat.name[i as usize].trim(), "SG" | "ND1" | "NE2"))
    {
        return;
    }
    let count = |name| {
        others
            .iter()
            .filter(|&&i| st.atom_group(flat.path[i as usize]).resname.trim() == name)
            .count()
    };
    let (sulfur, nitrogen, sulfur_angle, nitrogen_angle) = if others.len() < 4 {
        (Some((2.3, 0.03)), Some((2.3, 0.03)), None, None)
    } else {
        match (count("CYS"), count("HIS")) {
            (4, 0) => (Some((2.330, 0.029)), None, Some((109.45, 5.46)), None),
            (3, 1) => (
                Some((2.318, 0.027)),
                Some((2.074, 0.056)),
                Some((112.15, 3.96)),
                None,
            ),
            (2, 2) => (
                Some((2.306, 0.029)),
                Some((2.040, 0.050)),
                Some((116.23, 4.58)),
                Some((102.38, 5.44)),
            ),
            (1, 3) => (
                Some((2.298, 0.017)),
                Some((2.002, 0.045)),
                None,
                Some((107.23, 4.78)),
            ),
            _ => (None, None, None, None),
        }
    };
    for (k, &other) in others.iter().enumerate() {
        let name = flat.name[other as usize].trim();
        if let Some((target, sigma)) = if name == "SG" { sulfur } else { nitrogen } {
            graph.bond([metal, other], target, sigma, ORIGIN_METAL);
        }
        for &previous in &others[..k] {
            let previous_name = flat.name[previous as usize].trim();
            // The reference has separate ND1 and NE2 keys, no mixed-N key.
            if name != previous_name {
                continue;
            }
            if let Some((target, sigma)) = if name == "SG" {
                sulfur_angle
            } else {
                nitrogen_angle
            } {
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
