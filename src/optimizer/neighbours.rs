//! Neighbor search used while placing Movers.

use super::*;

/// Neighbor search used while placing Movers.
///
/// Fixed mode finds atoms at their current positions. Compat mode emulates
/// Reduce2's `SpatialQuery`, which files each atom in the grid cell of the
/// position it had when last added: atoms that Movers or methyl staggering
/// move without removing and re-adding them stay filed in their old cell, so
/// queries can miss them.
pub(super) struct NeighborIndex {
    loose: Option<SpatialGrid>,
    /// Fixed mode: the atoms added after `loose` was built (phantoms).
    late: Option<SpatialGrid>,
    margin: f64,
    reg: Option<RegGrid>,
}

struct RegGrid {
    geom: SpatialGrid,
    cells: Vec<SmallVec<[u32; 8]>>,
}

impl NeighborIndex {
    pub(super) fn new(w: &World, atoms: &[u32], compat: bool) -> Self {
        if compat {
            let (lo, bin, dims) = reduce2_grid_geometry(atoms.iter().map(|&a| w.pos[a as usize]));
            let geom = SpatialGrid::with_geometry(lo, bin, dims, &[]);
            let mut cells = vec![SmallVec::new(); dims[0] * dims[1] * dims[2]];
            for &a in atoms {
                let c = geom.cell_of(w.pos[a as usize]);
                if !cells[c].contains(&a) {
                    cells[c].push(a);
                }
            }
            NeighborIndex { loose: None, late: None, margin: 0.0, reg: Some(RegGrid { geom, cells }) }
        } else {
            let pts: Vec<(u32, Vec3)> = atoms.iter().map(|&a| (a, w.pos[a as usize])).collect();
            NeighborIndex { loose: Some(SpatialGrid::new(&pts, 3.0)), late: None, margin: 2.5, reg: None }
        }
    }

    /// Fixed mode: add atoms `ids` to the index, as a rebuild over all atoms
    /// would. No atom may have moved since the index was built.
    pub(super) fn add_late(&mut self, w: &World, ids: &[u32]) {
        if self.loose.is_some() && !ids.is_empty() {
            let pts: Vec<(u32, Vec3)> = ids.iter().map(|&a| (a, w.pos[a as usize])).collect();
            self.late = Some(SpatialGrid::new(&pts, 3.0));
        }
    }

    /// `SpatialQuery::remove(a)` at the atom's current position.
    pub(super) fn remove(&mut self, w: &World, a: u32) {
        if let Some(r) = &mut self.reg {
            let c = r.geom.cell_of(w.pos[a as usize]);
            r.cells[c].retain(|x| *x != a);
        }
    }
    /// `SpatialQuery::add(a)` at the atom's current position.
    pub(super) fn add(&mut self, w: &World, a: u32) {
        if let Some(r) = &mut self.reg {
            let c = r.geom.cell_of(w.pos[a as usize]);
            if !r.cells[c].contains(&a) {
                r.cells[c].push(a);
            }
        } else if let Some(g) = &self.loose {
            // atoms added later (phantoms) are not in the loose grid; callers
            // keep them out of placement queries, as the original does not
            // search for them there either except for ion checks.
            let _ = g;
        }
    }

    pub(super) fn neighbors(&self, w: &World, p: Vec3, min_d: f64, max_d: f64, removed: &FxHashSet<u32>) -> Vec<u32> {
        let mut out = Vec::new();
        let (min2, max2) = (min_d * min_d, max_d * max_d);
        if let Some(r) = &self.reg {
            let [(x0, x1), (y0, y1), (z0, z1)] = r.geom.cell_ranges(p, max_d);
            let dims = r.geom.dims();
            for z in z0..=z1 {
                for y in y0..=y1 {
                    for x in x0..=x1 {
                        for &id in &r.cells[x + dims[0] * (y + dims[1] * z)] {
                            let d2 = w.pos[id as usize].dist_sq(p);
                            if d2 >= min2 && d2 <= max2 {
                                out.push(id);
                            }
                        }
                    }
                }
            }
        } else {
            // the distance test first: it is cheaper than the set lookup and
            // rejects most of the grid's candidates
            for g in self.loose.iter().chain(&self.late) {
                g.for_each_within(p, 0.0, max_d + self.margin, |id, _, _| {
                    let d2 = w.pos[id as usize].dist_sq(p);
                    if d2 >= min2 && d2 <= max2 && !removed.contains(&id) {
                        out.push(id);
                    }
                });
            }
        }
        out.sort_unstable();
        out.dedup();
        out
    }

    /// Grid of the static scoring targets. Compat mode files them in their
    /// registered cells with Reduce2's geometry.
    pub(super) fn static_grid(&self, w: &World, include: impl Fn(u32) -> bool) -> Option<SpatialGrid> {
        let r = self.reg.as_ref()?;
        let mut pts: Vec<(u32, Vec3, Vec3)> = Vec::new();
        for (ci, cell) in r.cells.iter().enumerate() {
            if cell.is_empty() {
                continue;
            }
            let center = r.geom.cell_center(ci);
            for &a in cell {
                if include(a) {
                    pts.push((a, center, w.pos[a as usize]));
                }
            }
        }
        let g = &r.geom;
        Some(SpatialGrid::with_geometry(g.lower_bound(), g.bins(), g.dims(), &pts))
    }
}
