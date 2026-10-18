use super::*;
use crate::cell::SymOp;
use crate::clashcheck::neighbors::Query;

#[derive(Clone, Copy, Debug)]
pub(super) struct Candidate {
    pub indices: [u32; 2],
    pub distance: f64,
    pub delta: f64,
    /// Transform the second atom into the first atom's crystal image.
    pub image: SymOp,
}

fn radii(graph: &ValidationGraph, ml: &MonLib, nuclear: bool) -> Vec<Option<f64>> {
    let missing_h = !graph.hydrogen_flags.iter().any(|&h| h);
    graph
        .topology
        .etype
        .iter()
        .enumerate()
        .map(|(i, kind)| {
            let EType::Typed(name) = kind else {
                return None;
            };
            let atom = ml.ener.get(name.as_str())?;
            let vdw = if nuclear {
                atom.vdw_radius_neutron.or(atom.vdw_radius)
            } else {
                atom.vdw_radius
            };
            let normal = if missing_h {
                atom.vdwh_radius.or(vdw)
            } else {
                vdw.or(atom.vdwh_radius)
            };
            if graph.typed_charges[i] != 0 {
                atom.ion_radius.or(normal)
            } else {
                normal
            }
        })
        .collect()
}

fn compatible(flat: &FlatAtoms, first: usize, second: usize) -> bool {
    flat.path[first].model == flat.path[second].model
        && (flat.altloc[first].is_empty()
            || flat.altloc[second].is_empty()
            || flat.altloc[first] == flat.altloc[second])
}

fn image_key(image: SymOp) -> [u64; 12] {
    let mut key = [0; 12];
    for (out, value) in key
        .iter_mut()
        .zip(image.r.into_iter().flatten().chain(image.t))
    {
        // Space-group operators and cell shifts have exact rational values.
        *out = (if value == 0. { 0. } else { value }).to_bits();
    }
    key
}

fn exclusions(graph: &ValidationGraph) -> (FxHashSet<[u32; 2]>, FxHashSet<[u32; 2]>) {
    let mut neighbors: Vec<Vec<u32>> = vec![Vec::new(); graph.sites.len()];
    let mut result = FxHashSet::default();
    for bond in &graph.topology.bonds {
        neighbors[bond.i as usize].push(bond.j);
        neighbors[bond.j as usize].push(bond.i);
        result.insert([bond.i.min(bond.j), bond.i.max(bond.j)]);
    }
    for group in &neighbors {
        for (k, &first) in group.iter().enumerate() {
            for &second in &group[..k] {
                result.insert([first.min(second), first.max(second)]);
            }
        }
    }
    let mut shell3 = FxHashSet::default();
    for bond in &graph.topology.bonds {
        for &first in &neighbors[bond.i as usize] {
            for &second in &neighbors[bond.j as usize] {
                let pair = [first.min(second), first.max(second)];
                if first != second && !result.contains(&pair) {
                    shell3.insert(pair);
                }
            }
        }
    }
    (result, shell3)
}

struct Nonbonded<'a> {
    graph: &'a ValidationGraph,
    flat: &'a FlatAtoms,
    radius: Vec<Option<f64>>,
    excluded: FxHashSet<[u32; 2]>,
    shell3: FxHashSet<[u32; 2]>,
}

impl Nonbonded<'_> {
    fn candidate(
        &self,
        anchor: usize,
        index: u32,
        distance_sq: f64,
        mut image: SymOp,
    ) -> Option<Candidate> {
        let other = index as usize;
        let simple = image.is_identity();
        if !compatible(self.flat, anchor, other) || (simple && anchor == other) {
            return None;
        }
        let indices = [index.min(anchor as u32), index.max(anchor as u32)];
        if indices[0] == anchor as u32 {
            image = image.inverse()?;
        }
        if simple && self.excluded.contains(&indices) {
            return None;
        }
        if self
            .graph
            .symmetry_bonds
            .iter()
            .any(|bond| bond.indices == indices && image_key(bond.image) == image_key(image))
        {
            return None;
        }
        let distance = distance_sq.sqrt();
        let mut ideal = self.radius[anchor]
            .zip(self.radius[other])
            .map_or(1., |(a, b)| a + b);
        if simple && self.shell3.contains(&indices) {
            ideal *= 0.8;
        }
        ideal = ideal.max(1.);
        Some(Candidate {
            indices,
            distance,
            delta: distance - ideal,
            image,
        })
    }
}

/// Initial nonbonded candidates shared by both MCL hooks, before either adds bonds.
pub(super) fn collect(
    graph: &ValidationGraph,
    st: &Structure,
    flat: &FlatAtoms,
    ml: &MonLib,
    options: &ValidationOptions,
) -> Result<Vec<Candidate>, String> {
    let anchors: Vec<_> = (0..flat.pos.len())
        .filter(|&i| {
            (flat.element[i] == "ZN" && flat.name[i].trim() == "ZN")
                || cluster_name(st, flat, i).is_some()
        })
        .collect();
    if anchors.is_empty() {
        return Ok(Vec::new());
    }
    let query = Query::new(&flat.pos)?;
    let (excluded, shell3) = exclusions(graph);
    let nonbonded = Nonbonded {
        graph,
        flat,
        radius: radii(graph, ml, options.nuclear),
        excluded,
        shell3,
    };
    let mut seen = FxHashSet::default();
    let mut result = Vec::new();
    let identity = SymOp {
        r: [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
        t: [0.; 3],
    };
    let mut visit = |anchor: usize, position: Vec3, image: SymOp| {
        query.for_each_within(position, 0., 3.5, |index, _, distance_sq| {
            if let Some(pair) = nonbonded.candidate(anchor, index, distance_sq, image) {
                if seen.insert((pair.indices, image_key(pair.image))) {
                    result.push(pair);
                }
            }
        });
    };
    for &anchor in &anchors {
        visit(anchor, flat.pos[anchor], identity);
    }
    images::for_each(st, flat, &anchors, &mut visit)?;
    result.sort_by(|a, b| {
        a.delta
            .total_cmp(&b.delta)
            .then(a.indices.cmp(&b.indices))
            .then(image_key(a.image).cmp(&image_key(b.image)))
    });
    Ok(result)
}

pub(super) fn group(flat: &FlatAtoms, index: u32) -> AtomPath {
    let mut result = flat.path[index as usize];
    result.atom = 0;
    result
}
