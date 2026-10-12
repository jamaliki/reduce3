//! Protein backbone CDL v1.2, nearest-ten-degree lookup without interpolation.
//! Table values come from CCTBX 2026.9; see the provenance and license files.

use super::*;

const DATA: &[u8] = include_bytes!("cdl_targets.bin");

fn group(current: &str, next: &str) -> Option<usize> {
    let family = match current {
        "GLY" => 0,
        "ILE" | "VAL" => 2,
        "PRO" => 6,
        "ALA" | "ARG" | "ASN" | "ASP" | "CYS" | "GLN" | "GLU" | "HIS" | "LEU" | "LYS" | "MET"
        | "PHE" | "SER" | "THR" | "TRP" | "TYR" => 4,
        _ => return None,
    };
    Some(family + usize::from(next == "PRO"))
}

fn row(group: usize, phi: f64, psi: f64) -> [f64; 24] {
    let bin = |v: f64| ((v / 10.).round_ties_even() as i32 + 18).rem_euclid(36) as usize;
    let offset = ((group * 36 + bin(phi)) * 36 + bin(psi)) * 24 * 8;
    std::array::from_fn(|i| {
        f64::from_le_bytes(
            DATA[offset + i * 8..offset + (i + 1) * 8]
                .try_into()
                .unwrap(),
        )
    })
}

pub(super) fn apply(graph: &mut ValidationGraph, residues: &[ResInterp], chain: &[ConfResidue]) {
    for centre in 1..residues.len().saturating_sub(1) {
        let previous = &residues[centre - 1];
        let current = &residues[centre];
        let next = &residues[centre + 1];
        if !chain[centre].link_to_previous
            || !chain[centre + 1].link_to_previous
            || !previous.is_peptide
            || !current.is_peptide
            || !next.is_peptide
        {
            continue;
        }
        let Some(family) = group(
            chain[centre].resname.trim(),
            chain[centre + 1].resname.trim(),
        ) else {
            continue;
        };
        let ids = [
            lookup(previous, "C"),
            lookup(previous, "CA"),
            lookup(current, "N"),
            lookup(current, "CA"),
            lookup(current, "CB"),
            lookup(current, "C"),
            lookup(current, "O"),
            lookup(next, "N"),
        ];
        let (Some(c0), Some(ca0), Some(n), Some(ca), Some(c), Some(n1)) =
            (ids[0], ids[1], ids[2], ids[3], ids[5], ids[7])
        else {
            continue;
        };
        if !graph.topology.has_bond(c0, n) || !graph.topology.has_bond(c, n1) {
            continue;
        }
        // CCTBX skips CDL for cis groups and handles enol targets separately.
        let Some(omega) = dihedral_model(&graph.sites, [ca0, c0, n, ca]) else {
            continue;
        };
        if omega.abs() < 45. {
            continue;
        }
        if lookup(current, "HNO").is_some() {
            graph.unsupported(
                &["bonds", "angles", "planarities"],
                "Enol peptide targets need reference interpretation",
            );
            continue;
        }
        let (Some(phi), Some(psi)) = (
            dihedral_model(&graph.sites, [c0, n, ca, c]),
            dihedral_model(&graph.sites, [n, ca, c, n1]),
        ) else {
            continue;
        };
        let targets = row(family, phi, psi);
        let angles = [
            [ids[0], ids[2], ids[3]],
            [ids[2], ids[3], ids[4]],
            [ids[2], ids[3], ids[5]],
            [ids[4], ids[3], ids[5]],
            [ids[3], ids[5], ids[6]],
            [ids[3], ids[5], ids[7]],
            [ids[6], ids[5], ids[7]],
        ];
        for (i, atoms) in angles.iter().enumerate() {
            if atoms.iter().any(Option::is_none) {
                continue;
            }
            let mut key = atoms.map(Option::unwrap);
            if key[0] > key[2] {
                key.swap(0, 2);
            }
            if let Some(&k) = graph.angle_keys.get(&key) {
                graph.angles[k].ideal = targets[2 * i];
                graph.angles[k].weight = 1. / targets[2 * i + 1].powi(2);
            }
        }
        let bonds = [
            [ids[0], ids[2]],
            [ids[2], ids[3]],
            [ids[3], ids[4]],
            [ids[3], ids[5]],
            [ids[5], ids[6]],
        ];
        for (i, atoms) in bonds.iter().enumerate() {
            if atoms.iter().any(Option::is_none) {
                continue;
            }
            let mut key = atoms.map(Option::unwrap);
            key.sort_unstable();
            if let Some(&k) = graph.bond_keys.get(&key) {
                graph.bonds[k].ideal = targets[14 + 2 * i];
                graph.bonds[k].weight = 1. / targets[15 + 2 * i].powi(2);
            }
        }
        if let Some(h) = lookup(current, "H") {
            let canonical = |mut key: [u32; 3]| {
                if key[0] > key[2] {
                    key.swap(0, 2);
                }
                key
            };
            let keys = [
                canonical([c0, n, ca]),
                canonical([c0, n, h]),
                canonical([ca, n, h]),
            ];
            if let (Some(&a), Some(&b), Some(&c)) = (
                graph.angle_keys.get(&keys[0]),
                graph.angle_keys.get(&keys[1]),
                graph.angle_keys.get(&keys[2]),
            ) {
                let difference =
                    (graph.angles[a].ideal + graph.angles[b].ideal + graph.angles[c].ideal - 360.)
                        / 2.;
                graph.angles[b].ideal -= difference;
                graph.angles[c].ideal -= difference;
            } else {
                graph.unsupported(
                    &["angles"],
                    "CDL requires both C-N-H and CA-N-H angle proxies",
                );
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn grid_wraps_and_uses_ties_even() {
        assert_eq!(row(4, 180., 180.), row(4, -180., -180.));
        assert_eq!(row(4, 5., 5.), row(4, 0., 0.));
        assert_eq!(row(4, 15., 15.), row(4, 20., 20.));
        assert_eq!(DATA.len(), 8 * 36 * 36 * 24 * 8);
    }
}
