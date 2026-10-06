//! Probe-style contact-dot scoring (port of mmtbx/probe Scoring.cpp,
//! DotSpheres.cpp and SpatialQuery.cpp), organized for speed.

use crate::geom::{v3, Vec3};
use rustc_hash::FxHashMap;
use std::f64::consts::PI;

/// Per-atom information Probe needs beyond the hierarchy (`ExtraAtomInfo`).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AtomInfo {
    pub vdw_radius: f64,
    pub is_acceptor: bool,
    pub is_donor: bool,
    pub is_dummy_hydrogen: bool,
    pub is_ion: bool,
    pub charge: i8,
    /// Alternate location character, or 0 for none.
    pub alt: u8,
}

impl Default for AtomInfo {
    fn default() -> Self {
        AtomInfo {
            vdw_radius: 0.0,
            is_acceptor: false,
            is_donor: false,
            is_dummy_hydrogen: false,
            is_ion: false,
            charge: 0,
            alt: 0,
        }
    }
}

/// `DotScorer::compatible_conformations`.
#[inline(always)]
pub fn compatible_alts(a: u8, b: u8) -> bool {
    a == 0 || a == b' ' || b == 0 || b == b' ' || a == b
}

/// Probe scoring parameters (the `probe` PHIL scope, as used by Reduce2).
#[derive(Clone, Debug)]
pub struct ProbeParams {
    pub probe_radius: f64,
    pub density: f64,
    pub worse_clash_cutoff: f64,
    pub clash_cutoff: f64,
    pub contact_cutoff: f64,
    pub uncharged_hydrogen_cutoff: f64,
    pub charged_hydrogen_cutoff: f64,
    pub bump_weight: f64,
    pub hydrogen_bond_weight: f64,
    pub gap_weight: f64,
    pub allow_weak_hydrogen_bonds: bool,
    pub ignore_ion_interactions: bool,
    pub set_polar_hydrogen_radius: bool,
}

impl ProbeParams {
    /// Reduce2 defaults: probe2 defaults with bump and H-bond weights raised 10x.
    pub fn reduce2_defaults() -> ProbeParams {
        ProbeParams {
            probe_radius: 0.25,
            density: 16.0,
            worse_clash_cutoff: 0.5,
            clash_cutoff: 0.4,
            contact_cutoff: 0.25,
            uncharged_hydrogen_cutoff: 0.6,
            charged_hydrogen_cutoff: 0.8,
            bump_weight: 100.0,
            hydrogen_bond_weight: 40.0,
            gap_weight: 0.25,
            allow_weak_hydrogen_bonds: false,
            ignore_ion_interactions: false,
            set_polar_hydrogen_radius: true,
        }
    }
}

// ----------------------------------------------------------------------------
// Dot spheres

/// Dots on a sphere, generated exactly as `molprobity::probe::DotSphere`.
pub fn dot_sphere(radius: f64, density: f64) -> Vec<Vec3> {
    let rad = radius.max(0.0);
    let dens = density.max(0.0);
    let mut out = Vec::new();
    if rad == 0.0 || dens == 0.0 {
        return out;
    }
    let num_dots = (4.0 * PI * dens * (rad * rad)).floor() as usize;
    let offset = 0.2;
    let nequator = ((num_dots as f64) * PI).sqrt().floor() as i32;
    let ang = 5.0 * PI / 360.0;
    let cosang = ang.cos();
    let sinang = ang.sin();
    let nvert = nequator / 2;
    let mut odd = true;
    for j in 0..=nvert {
        let phi = (PI * j as f64) / nvert as f64;
        let z0 = phi.cos() * rad;
        let xy0 = phi.sin() * rad;
        let mut nhoriz = (nequator as f64 * phi.sin()).floor() as i32;
        if nhoriz < 1 {
            nhoriz = 1;
        }
        for k in 0..nhoriz {
            let theta = if odd {
                (2.0 * PI * k as f64 + offset) / nhoriz as f64
            } else {
                (2.0 * PI * k as f64) / nhoriz as f64
            };
            let x0 = theta.cos() * xy0;
            let y0 = theta.sin() * xy0;
            out.push(v3(x0, y0 * cosang - z0 * sinang, y0 * sinang + z0 * cosang));
        }
        odd = !odd;
    }
    out
}

/// Cache of dot spheres keyed by exact radius (`DotSphereCache`).
#[derive(Default)]
pub struct DotSphereCache {
    density: f64,
    spheres: FxHashMap<u64, std::sync::Arc<Vec<Vec3>>>,
}

impl DotSphereCache {
    pub fn new(density: f64) -> Self {
        DotSphereCache { density, spheres: FxHashMap::default() }
    }
    pub fn get(&mut self, radius: f64) -> std::sync::Arc<Vec<Vec3>> {
        let d = self.density;
        self.spheres
            .entry(radius.to_bits())
            .or_insert_with(|| std::sync::Arc::new(dot_sphere(radius, d)))
            .clone()
    }
}

// ----------------------------------------------------------------------------
// Static spatial grid (CSR layout)

pub struct SpatialGrid {
    lower: Vec3,
    bin: [f64; 3],
    dims: [usize; 3],
    starts: Vec<u32>,
    items: Vec<u32>,
    pos: Vec<Vec3>,
}

impl SpatialGrid {
    /// Build a grid over `(id, position)` pairs with the given cell size.
    pub fn new(points: &[(u32, Vec3)], cell: f64) -> SpatialGrid {
        let mut lo = v3(f64::INFINITY, f64::INFINITY, f64::INFINITY);
        let mut hi = v3(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
        for &(_, p) in points {
            lo = v3(lo.x.min(p.x), lo.y.min(p.y), lo.z.min(p.z));
            hi = v3(hi.x.max(p.x), hi.y.max(p.y), hi.z.max(p.z));
        }
        if points.is_empty() {
            lo = Vec3::ZERO;
            hi = Vec3::ZERO;
        }
        let dims = [
            (((hi.x - lo.x) / cell).floor() as usize + 1).max(1),
            (((hi.y - lo.y) / cell).floor() as usize + 1).max(1),
            (((hi.z - lo.z) / cell).floor() as usize + 1).max(1),
        ];
        let pts: Vec<(u32, Vec3, Vec3)> = points.iter().map(|&(i, p)| (i, p, p)).collect();
        Self::with_geometry(lo, [cell; 3], dims, &pts)
    }

    /// Build a grid with explicit geometry over `(id, bucket position, actual
    /// position)`: each point is filed in the cell of its bucket position but
    /// distances use its actual position.
    pub fn with_geometry(lower: Vec3, bin: [f64; 3], dims: [usize; 3], points: &[(u32, Vec3, Vec3)]) -> SpatialGrid {
        let mut g = SpatialGrid { lower, bin, dims, starts: Vec::new(), items: Vec::new(), pos: Vec::new() };
        let ncell = dims[0] * dims[1] * dims[2];
        let mut counts = vec![0u32; ncell + 1];
        let cells: Vec<usize> = points.iter().map(|&(_, b, _)| g.cell_of(b)).collect();
        for &c in &cells {
            counts[c + 1] += 1;
        }
        for i in 0..ncell {
            counts[i + 1] += counts[i];
        }
        let mut fill = counts.clone();
        let mut items = vec![0u32; points.len()];
        let mut pos = vec![Vec3::ZERO; points.len()];
        for (k, &(id, _, p)) in points.iter().enumerate() {
            let c = cells[k];
            let slot = fill[c] as usize;
            items[slot] = id;
            pos[slot] = p;
            fill[c] += 1;
        }
        g.starts = counts;
        g.items = items;
        g.pos = pos;
        g
    }

    #[inline]
    fn axis_index(&self, v: f64, k: usize) -> usize {
        // As molprobity::probe::SpatialQuery::grid_index: below the lower bound -> 0,
        // past the end -> last cell.
        let lo = self.lower[k];
        let i = if v < lo { 0 } else { ((v - lo) / self.bin[k]).floor() as usize };
        i.min(self.dims[k] - 1)
    }

    #[inline]
    pub fn cell_of(&self, p: Vec3) -> usize {
        let ix = self.axis_index(p.x, 0);
        let iy = self.axis_index(p.y, 1);
        let iz = self.axis_index(p.z, 2);
        ix + self.dims[0] * (iy + self.dims[1] * iz)
    }

    /// Inclusive cell index ranges along x, y, z covering `p +/- max_d`.
    pub fn cell_ranges(&self, p: Vec3, max_d: f64) -> [(usize, usize); 3] {
        [
            (self.axis_index(p.x - max_d, 0), self.axis_index(p.x + max_d, 0)),
            (self.axis_index(p.y - max_d, 1), self.axis_index(p.y + max_d, 1)),
            (self.axis_index(p.z - max_d, 2), self.axis_index(p.z + max_d, 2)),
        ]
    }
    #[inline]
    pub fn dims(&self) -> [usize; 3] {
        self.dims
    }
    #[inline]
    pub fn lower_bound(&self) -> Vec3 {
        self.lower
    }
    #[inline]
    pub fn bins(&self) -> [f64; 3] {
        self.bin
    }
    /// Center of a cell (inside the grid, so `cell_of(center) == index`).
    pub fn cell_center(&self, index: usize) -> Vec3 {
        let ix = index % self.dims[0];
        let iy = (index / self.dims[0]) % self.dims[1];
        let iz = index / (self.dims[0] * self.dims[1]);
        v3(
            self.lower.x + (ix as f64 + 0.5) * self.bin[0],
            self.lower.y + (iy as f64 + 0.5) * self.bin[1],
            self.lower.z + (iz as f64 + 0.5) * self.bin[2],
        )
    }

    /// Call `f(id, pos, dist_sq)` for every point with `min_d <= dist <= max_d` from `p`.
    #[inline]
    pub fn for_each_within<F: FnMut(u32, Vec3, f64)>(&self, p: Vec3, min_d: f64, max_d: f64, mut f: F) {
        let (x0, x1) = (self.axis_index(p.x - max_d, 0), self.axis_index(p.x + max_d, 0));
        let (y0, y1) = (self.axis_index(p.y - max_d, 1), self.axis_index(p.y + max_d, 1));
        let (z0, z1) = (self.axis_index(p.z - max_d, 2), self.axis_index(p.z + max_d, 2));
        let min2 = min_d * min_d;
        let max2 = max_d * max_d;
        for z in z0..=z1 {
            for y in y0..=y1 {
                let row = self.dims[0] * (y + self.dims[1] * z);
                let s = self.starts[row + x0] as usize;
                let e = self.starts[row + x1 + 1] as usize;
                for k in s..e {
                    let d2 = self.pos[k].dist_sq(p);
                    if d2 >= min2 && d2 <= max2 {
                        f(self.items[k], self.pos[k], d2);
                    }
                }
            }
        }
    }
}

/// Geometry of `molprobity::probe::SpatialQuery(atoms)`: bins of at least 3 A
/// (or 1/50 of the extent), covering the bounding box of the atoms.
pub fn reduce2_grid_geometry(points: impl Iterator<Item = Vec3>) -> (Vec3, [f64; 3], [usize; 3]) {
    let mut lo = v3(1e10, 1e10, 1e10);
    let mut hi = v3(-1e10, -1e10, -1e10);
    for p in points {
        lo = v3(lo.x.min(p.x), lo.y.min(p.y), lo.z.min(p.z));
        hi = v3(hi.x.max(p.x), hi.y.max(p.y), hi.z.max(p.z));
    }
    let mut bin = [3.0f64; 3];
    let mut dims = [1usize; 3];
    for k in 0..3 {
        let (l, mut u) = (lo[k], hi[k]);
        if u < l {
            u = l;
        }
        let min_size = (u - l) / 50.0;
        if bin[k] < min_size {
            bin[k] = min_size;
        }
        dims[k] = (((u - l) / bin[k]).ceil() as usize).max(1);
    }
    (lo, bin, dims)
}

// ----------------------------------------------------------------------------
// Dot scoring

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum OverlapType {
    Ignore,
    Clash,
    NoOverlap,
    HydrogenBond,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum InteractionType {
    WideContact,
    CloseContact,
    WeakHydrogenBond,
    SmallOverlap,
    Bump,
    BadBump,
    StandardHydrogenBond,
    Invalid,
}

#[derive(Clone, Copy, Debug, Default)]
pub struct ScoreDotsResult {
    pub bump: f64,
    pub hbond: f64,
    pub attract: f64,
    pub has_bad_bump: bool,
}

impl ScoreDotsResult {
    #[inline]
    pub fn total(&self) -> f64 {
        self.bump + self.hbond + self.attract
    }
}

/// An atom that may interact with the dots of a source atom.
#[derive(Clone, Copy, Debug)]
pub struct Target {
    pub pos: Vec3,
    pub info: AtomInfo,
}

/// Precomputed scorer constants.
#[derive(Clone, Debug)]
pub struct DotScorer {
    pub gap_scale: f64,
    pub bump_weight: f64,
    pub hbond_weight: f64,
    pub max_regular_h_overlap: f64,
    pub max_charged_h_overlap: f64,
    pub bump_overlap: f64,
    pub bad_bump_overlap: f64,
    pub contact_cutoff: f64,
    pub weak_hbonds: bool,
    pub ignore_ions: bool,
}

impl DotScorer {
    pub fn new(p: &ProbeParams) -> DotScorer {
        DotScorer {
            gap_scale: p.gap_weight.max(p.probe_radius),
            bump_weight: p.bump_weight,
            hbond_weight: p.hydrogen_bond_weight,
            max_regular_h_overlap: p.uncharged_hydrogen_cutoff,
            max_charged_h_overlap: p.charged_hydrogen_cutoff,
            bump_overlap: p.clash_cutoff,
            bad_bump_overlap: p.worse_clash_cutoff,
            contact_cutoff: p.contact_cutoff,
            weak_hbonds: p.allow_weak_hydrogen_bonds,
            ignore_ions: p.ignore_ion_interactions,
        }
    }

    #[inline]
    pub fn interaction_type(&self, ot: OverlapType, gap: f64, separate_bad_bumps: bool) -> InteractionType {
        match ot {
            OverlapType::NoOverlap => {
                if gap > self.contact_cutoff {
                    InteractionType::WideContact
                } else {
                    InteractionType::CloseContact
                }
            }
            OverlapType::Clash => {
                if gap > -self.bump_overlap {
                    InteractionType::SmallOverlap
                } else if separate_bad_bumps {
                    if gap > -self.bad_bump_overlap { InteractionType::Bump } else { InteractionType::BadBump }
                } else {
                    InteractionType::Bump
                }
            }
            OverlapType::HydrogenBond => {
                if self.weak_hbonds && gap > 0.0 {
                    InteractionType::WeakHydrogenBond
                } else {
                    InteractionType::StandardHydrogenBond
                }
            }
            OverlapType::Ignore => InteractionType::Invalid,
        }
    }

    /// Score the dots of a source atom against already-selected interacting
    /// targets. `dots` must already exclude dots inside excluded atoms
    /// (`preTrimmedDots` path of `DotScorer::score_dots`).
    ///
    /// Mirrors `score_dots`/`check_dot` exactly, including the order in which
    /// targets are tested (ties on gap keep the first target).
    pub fn score_dots(
        &self,
        src_pos: Vec3,
        src: &AtomInfo,
        dots: &[Vec3],
        targets: &[Target],
        probe_radius: f64,
        density: f64,
        only_bumps: bool,
    ) -> ScoreDotsResult {
        let mut ret = ScoreDotsResult::default();
        if density <= 0.0 || probe_radius < 0.0 {
            return ret;
        }
        if src.is_ion && self.ignore_ions {
            return ret;
        }
        let src_r = src.vdw_radius;
        let probe_extra = if src_r > 0.0 { (src_r + probe_radius) / src_r } else { 1.0 };
        for d in dots {
            let dot_abs = src_pos + *d;
            // probe location: same direction, further out by the probe radius
            let probe_loc = if src_r > 0.0 {
                src_pos + *d * probe_extra
            } else {
                let l = d.length();
                if l > 0.0 { src_pos + *d * ((l + probe_radius) / l) } else { src_pos }
            };
            let r = self.check_dot_core(src, src_pos, dot_abs, probe_loc, targets, probe_radius);
            let Some((ot, gap, overlap, annular)) = r else { continue };
            let it = self.interaction_type(ot, gap, true);
            match it {
                InteractionType::Invalid => {}
                InteractionType::WideContact | InteractionType::CloseContact | InteractionType::WeakHydrogenBond => {
                    if !only_bumps && !annular {
                        let sg = gap / self.gap_scale;
                        ret.attract += (-sg * sg).exp();
                    }
                }
                InteractionType::SmallOverlap | InteractionType::Bump | InteractionType::BadBump => {
                    ret.bump += -self.bump_weight * overlap;
                    if it == InteractionType::BadBump {
                        ret.has_bad_bump = true;
                    }
                }
                InteractionType::StandardHydrogenBond => {
                    if !only_bumps {
                        ret.hbond += self.hbond_weight * overlap;
                    } else {
                        ret.bump += -self.bump_weight * overlap;
                    }
                }
            }
        }
        ret.bump /= density;
        ret.hbond /= density;
        ret.attract /= density;
        ret
    }

    /// Core of `DotScorer::check_dot`, returning (overlap type, gap, overlap,
    /// annular) or None when no target is in range of the probe.
    #[inline(always)]
    fn check_dot_core(
        &self,
        src: &AtomInfo,
        src_pos: Vec3,
        dot_abs: Vec3,
        probe_loc: Vec3,
        targets: &[Target],
        probe_radius: f64,
    ) -> Option<(OverlapType, f64, f64, bool)> {
        let mut best_gap = 1e100;
        let mut best: usize = usize::MAX;
        let mut is_hbond = false;
        let mut too_close_hbond = false;
        let mut hbond_min_dist = 0.0;
        let mut keep = false;
        let mut cause_is_dummy = false;
        for (k, b) in targets.iter().enumerate() {
            let bi = &b.info;
            if bi.is_ion && self.ignore_ions {
                continue;
            }
            let vdwb = bi.vdw_radius;
            let pr = vdwb + probe_radius;
            if probe_loc.dist_sq(b.pos) > pr * pr {
                continue;
            }
            if !compatible_alts(src.alt, bi.alt) {
                continue;
            }
            let gap = dot_abs.dist(b.pos) - vdwb;
            if gap < best_gap {
                let charge_s = src.charge as i32;
                let charge_b = bi.charge as i32;
                let both_charged = charge_s != 0 && charge_b != 0;
                let complement = both_charged && charge_s * charge_b < 0;
                let could_hbond = (src.is_donor && bi.is_acceptor) || (src.is_acceptor && bi.is_donor);
                if could_hbond && (!both_charged || complement) {
                    is_hbond = true;
                    hbond_min_dist = if both_charged { self.max_charged_h_overlap } else { self.max_regular_h_overlap };
                    too_close_hbond = gap < -hbond_min_dist;
                } else {
                    if src.is_dummy_hydrogen || bi.is_dummy_hydrogen {
                        continue;
                    }
                    is_hbond = false;
                    too_close_hbond = false;
                }
                keep = true;
                cause_is_dummy = bi.is_dummy_hydrogen;
                best_gap = gap;
                best = k;
            }
        }
        let mut ot = OverlapType::Ignore;
        let mut overlap = 0.0;
        let mut gap = best_gap;
        let mut annular = false;
        if keep {
            if gap > 0.0 {
                overlap = 0.0;
                ot = if self.weak_hbonds && is_hbond { OverlapType::HydrogenBond } else { OverlapType::NoOverlap };
            } else if is_hbond {
                overlap = -0.5 * gap;
                if too_close_hbond {
                    gap += hbond_min_dist;
                    overlap = -0.5 * gap;
                    ot = OverlapType::Clash;
                } else {
                    ot = OverlapType::HydrogenBond;
                }
            } else {
                overlap = -0.5 * gap;
                ot = OverlapType::Clash;
            }
            let b = &targets[best];
            annular = annular_dots(dot_abs, src_pos, src.vdw_radius, b.pos, b.info.vdw_radius, probe_radius);
        }
        if (src.is_dummy_hydrogen || cause_is_dummy) && (too_close_hbond || ot != OverlapType::HydrogenBond) {
            ot = OverlapType::Ignore;
        }
        if ot == OverlapType::Ignore {
            return if keep { Some((OverlapType::Ignore, gap, overlap, annular)) } else { None };
        }
        Some((ot, gap, overlap, annular))
    }
}

#[inline(always)]
fn annular_dots(dot: Vec3, src: Vec3, src_r: f64, targ: Vec3, targ_r: f64, probe_r: f64) -> bool {
    // dot2srcCenter > kissEdge2bullsEye
    let src2targ = (targ - src).normalize() * src_r;
    let surface = src2targ + src;
    let d = (surface - dot).length();
    let kiss = 2.0 * src_r * (targ_r * probe_r / ((src_r + targ_r) * (src_r + probe_r))).sqrt();
    d > kiss
}

/// Remove dots that fall inside any excluded atom (`DotScorer::trim_dots`).
pub fn trim_dots(src_pos: Vec3, dots: &[Vec3], excluded: &[(Vec3, f64)], out: &mut Vec<Vec3>) {
    out.clear();
    'dot: for d in dots {
        let p = src_pos + *d;
        for &(ep, er) in excluded {
            if p.dist_sq(ep) < er * er {
                continue 'dot;
            }
        }
        out.push(*d);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn dot_counts_close_to_expected() {
        for &(r, den) in &[(1.0, 16.0), (1.7, 16.0), (1.1, 16.0)] {
            let d = dot_sphere(r, den);
            let expected = 4.0 * PI * den * r * r;
            assert!((d.len() as f64 - expected).abs() < expected * 0.2, "{} {}", d.len(), expected);
            for p in &d {
                assert!((p.length() - r).abs() < 1e-9);
            }
        }
    }
}
