//! Scoring context

use super::*;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Cfg {
    Coarse(u16),
    Fine(u16, u16),
}
impl Cfg {
    pub(super) fn coarse(self) -> usize {
        match self {
            Cfg::Coarse(c) | Cfg::Fine(c, _) => c as usize,
        }
    }
}

#[derive(Clone, Copy)]
pub(super) struct DotPair {
    pub(super) d: Vec3,
    pub(super) probe: Vec3,
}

pub(super) const NONE: u32 = u32::MAX;

pub(super) struct Ctx<'a> {
    pub(super) w: &'a World,
    pub(super) movers: &'a [Mover],
    pub(super) p: &'a OptParams,
    pub(super) scorer: DotScorer,
    /// Scorer view of atom info in compat mode (Reduce2's DotScorer keeps a copy).
    pub(super) snapshot: Option<Vec<AtomInfo>>,
    pub(super) grid: SpatialGrid,
    /// (mover, slot) for atoms that move with a Mover (slot < n_moved).
    pub(super) dyn_of: DynOf,
    pub(super) atom_movers: FxHashMap<u32, SmallVec<[u32; 4]>>,
    pub(super) exclude: FxHashMap<u32, Vec<u32>>,
    pub(super) dots: FxHashMap<u32, Arc<Vec<DotPair>>>,
    /// The same dots as coordinate arrays.
    pub(super) dots_soa: FxHashMap<u32, Arc<DotsSoa>>,
    /// Largest distance of an atom's dots from its center.
    pub(super) dot_reach: FxHashMap<u32, f64>,
    /// Per Mover atom: the static atoms near any of its positions that can
    /// be its targets (not excluded, occupied), in grid order with their
    /// info (center, radius covered, range of `near_entries`), filtered per
    /// state instead of querying the grid for each.
    pub(super) near: FxHashMap<u32, (Vec3, f64, u32, u32)>,
    pub(super) near_entries: Vec<(Vec3, AtomInfo)>,
    pub(super) max_vdw: f64,
    /// Local index of each mover inside its clique.
    pub(super) local: Vec<u32>,
}

/// (mover, slot) of each atom that moves with a Mover, by atom index.
pub(super) struct DynOf(pub(super) Vec<(u32, u16)>);

impl DynOf {
    #[inline]
    pub(super) fn get(&self, a: &u32) -> Option<&(u32, u16)> {
        self.0.get(*a as usize).filter(|e| e.0 != NONE)
    }
    #[inline]
    pub(super) fn contains_key(&self, a: &u32) -> bool {
        self.get(a).is_some()
    }
    pub(super) fn insert(&mut self, a: u32, v: (u32, u16)) {
        if self.0.len() <= a as usize {
            self.0.resize(a as usize + 1, (NONE, 0));
        }
        self.0[a as usize] = v;
    }
}

/// The configuration of the Movers of one clique.
pub(super) struct CliqueState<'s> {
    pub(super) movers: &'s [u32],
    pub(super) cfg: Vec<Cfg>,
}

impl<'a> Ctx<'a> {
    #[inline]
    fn cfg_of(&self, st: &CliqueState, m: u32) -> Cfg {
        st.cfg[self.local[m as usize] as usize]
    }

    #[inline]
    pub(super) fn mover_atom_pos(&self, m: u32, slot: usize, cfg: Cfg) -> Vec3 {
        let mv = &self.movers[m as usize];
        match cfg {
            Cfg::Coarse(c) => mv.coarse_pos[c as usize][slot],
            Cfg::Fine(c, f) => mv.fine.pos(c as usize, f as usize)[slot],
        }
    }

    #[inline]
    pub(super) fn mover_atom_deleted(&self, m: u32, slot: usize, cfg: Cfg) -> bool {
        let d = in_state(&self.movers[m as usize].coarse_del, cfg.coarse());
        slot < d.len() && d[slot]
    }

    #[inline]
    pub(super) fn mover_atom_info(&self, m: u32, slot: usize, cfg: Cfg, a: u32) -> AtomInfo {
        if let Some(s) = &self.snapshot {
            return s[a as usize];
        }
        let inf = in_state(&self.movers[m as usize].coarse_info, cfg.coarse());
        if slot < inf.len() { inf[slot] } else { self.w.info[a as usize] }
    }

    #[inline]
    fn static_info(&self, a: u32) -> AtomInfo {
        match &self.snapshot {
            Some(s) => s[a as usize],
            None => self.w.info[a as usize],
        }
    }

    /// Position, deletion and info of any atom under the clique state.
    #[inline]
    pub(super) fn atom_state(&self, st: &CliqueState, a: u32) -> (Vec3, bool, AtomInfo) {
        match self.dyn_of.get(&a) {
            Some(&(m, slot)) if self.local[m as usize] != NONE && st.movers.get(self.local[m as usize] as usize) == Some(&m) => {
                let cfg = self.cfg_of(st, m);
                (
                    self.mover_atom_pos(m, slot as usize, cfg),
                    self.mover_atom_deleted(m, slot as usize, cfg),
                    self.mover_atom_info(m, slot as usize, cfg, a),
                )
            }
            _ => (self.w.pos[a as usize], false, self.static_info(a)),
        }
    }

    /// The static targets of atom `a` at `pa` (`for_each_within` over the grid
    /// with the target filters), from its cached neighborhood when that
    /// covers the position.
    #[inline]
    pub(super) fn static_targets(&self, a: u32, pa: Vec3, ia: &AtomInfo, nearby: f64, excl: &[u32], out: &mut Vec<Target>) {
        let pr = self.p.probe.probe_radius;
        if let Some(&(c, r, lo, hi)) = self.near.get(&a) {
            if pa.dist(c) + nearby + NEAR_MARGIN <= r {
                let (min2, max2) = (1e-5 * 1e-5, nearby * nearby);
                for &(pb, ib) in &self.near_entries[lo as usize..hi as usize] {
                    let d2 = pb.dist_sq(pa);
                    if d2 >= min2 && d2 <= max2 && d2.sqrt() <= ia.vdw_radius + ib.vdw_radius + 2.0 * pr {
                        out.push(Target { pos: pb, info: ib });
                    }
                }
                return;
            }
        }
        let (w, min_occ) = (self.w, self.p.min_occupancy);
        self.grid.for_each_within(pa, 1e-5, nearby, |b, pb, d2| {
            if excl.contains(&b) || w.occ[b as usize].abs() < min_occ {
                return;
            }
            let ib = self.static_info(b);
            if d2.sqrt() <= ia.vdw_radius + ib.vdw_radius + 2.0 * pr {
                out.push(Target { pos: pb, info: ib });
            }
        });
    }

    /// One evaluation of the dots of Mover atom `slot` of the Mover at local
    /// index `i` for all of its fine poses in coarse state `c`, with the rest
    /// of the clique as in `st`. None unless the atom keeps its place and info
    /// in every pose and its own Mover offers it no targets.
    pub(super) fn fine_shared(&self, st: &CliqueState, i: usize, c: usize, slot: usize, buf: &mut ScoreBuf) -> Option<FineShared> {
        let (w, p) = (self.w, self.p);
        let m = st.movers[i];
        let mv = &self.movers[m as usize];
        let a = mv.atoms[slot];
        let nf = mv.fine.count();
        let coarse = Cfg::Coarse(c as u16);
        if nf < 2 || p.compat || w.occ[a as usize] < p.min_occupancy || self.mover_atom_deleted(m, slot, coarse) {
            return None;
        }
        let ia = self.mover_atom_info(m, slot, coarse, a);
        if ia.is_ion && self.scorer.ignore_ions {
            return None;
        }
        let pa0 = mv.fine.pos(c, 0)[slot];
        if (0..nf).any(|f| mv.fine.pos(c, f)[slot].dist_sq(pa0) > SHARE_TOL * SHARE_TOL) {
            return None;
        }
        let pa: Vec<Vec3> = (0..nf).map(|f| mv.fine.pos(c, f)[slot]).collect();
        let pr = p.probe.probe_radius;
        let nearby = self.max_vdw + w.info[a as usize].vdw_radius + 2.0 * pr;
        let excl: &[u32] = self.exclude.get(&a).map(|v| v.as_slice()).unwrap_or(&[]);
        let min_occ = p.min_occupancy;
        // the own Mover's other atoms must be no target in any pose
        for oslot in 0..mv.n_moved {
            let b = mv.atoms[oslot];
            if b == a || self.mover_atom_deleted(m, oslot, coarse) || excl.contains(&b) || w.occ[b as usize].abs() < min_occ {
                continue;
            }
            let ib = self.mover_atom_info(m, oslot, coarse, b);
            if (0..nf).any(|f| {
                let d = mv.fine.pos(c, f)[oslot].dist(pa[f]);
                (1e-5..=nearby).contains(&d) && d <= ia.vdw_radius + ib.vdw_radius + 2.0 * pr
            }) {
                return None;
            }
        }
        // excluded atoms: the own Mover's per pose, the others once
        let reach = self.dot_reach.get(&a).copied().unwrap_or(f64::INFINITY) + EXCLUSION_MARGIN;
        let covers = |pa: Vec3, pe: Vec3, re: f64| pe.dist_sq(pa) < (reach + re) * (reach + re);
        buf.excl.clear();
        let mut trims = vec![Vec::new(); nf];
        for &e in excl {
            match self.dyn_of.get(&e) {
                Some(&(em, eslot)) if em == m => {
                    if self.mover_atom_deleted(m, eslot as usize, coarse) {
                        continue;
                    }
                    let re = self.mover_atom_info(m, eslot as usize, coarse, e).vdw_radius;
                    for (f, tr) in trims.iter_mut().enumerate() {
                        let pe = mv.fine.pos(c, f)[eslot as usize];
                        if covers(pa[f], pe, re) {
                            tr.push((pe, re));
                        }
                    }
                }
                _ => {
                    let (pe, deleted, ie) = self.atom_state(st, e);
                    if !deleted && covers(pa[0], pe, ie.vdw_radius) {
                        buf.excl.push((pe, ie.vdw_radius));
                    }
                }
            }
        }
        // targets: static, then the other clique Movers' (as in score_atom_ext)
        buf.targets.clear();
        self.static_targets(a, pa[0], &ia, nearby, excl, &mut buf.targets);
        for &om in self.atom_movers.get(&a).map(|v| v.as_slice()).unwrap_or(&[]) {
            let li = self.local[om as usize];
            if om == m || li == NONE || st.movers.get(li as usize) != Some(&om) {
                continue;
            }
            let cfg = st.cfg[li as usize];
            let omv = &self.movers[om as usize];
            for oslot in 0..omv.n_moved {
                let b = omv.atoms[oslot];
                if b == a || self.mover_atom_deleted(om, oslot, cfg) {
                    continue;
                }
                let pb = self.mover_atom_pos(om, oslot, cfg);
                let d = pb.dist(pa[0]);
                if d < 1e-5 || d > nearby || excl.contains(&b) || w.occ[b as usize].abs() < min_occ {
                    continue;
                }
                let ib = self.mover_atom_info(om, oslot, cfg, b);
                if d <= ia.vdw_radius + ib.vdw_radius + 2.0 * pr {
                    buf.targets.push(Target { pos: pb, info: ib });
                }
            }
        }
        let s = &self.scorer;
        prepare_targets(s, &ia, pa[0], &buf.targets, pr, &mut buf.prepared);
        let mut adds = Vec::new();
        fold_dot_lanes(s, &ia, pa[0], &self.dots_soa[&a], &buf.excl, &buf.prepared, &mut buf.lanes, |k, r| {
            if let Some((kind, val)) = dot_addend(s, r) {
                adds.push(DotAdd { k: k as u32, kind, val });
            }
        });
        Some(FineShared { adds, pa, trims })
    }

    /// Score one Mover atom under the clique state (`OptimizerC::scoreAtom`).
    pub(super) fn score_atom(&self, st: &CliqueState, a: u32, buf: &mut ScoreBuf) -> ScoreDotsResult {
        self.score_atom_ext(st, a, buf, None)
    }

    /// `score_atom` with Reduce2's trimmed-dot cache keyed by (atom, location
    /// index) when `trim` is given (compat mode).
    pub(super) fn score_atom_ext(
        &self,
        st: &CliqueState,
        a: u32,
        buf: &mut ScoreBuf,
        trim: Option<(&std::cell::RefCell<FxHashMap<(u32, u32), Vec<DotPair>>>, u32)>,
    ) -> ScoreDotsResult {
        let w = self.w;
        if w.occ[a as usize] < self.p.min_occupancy {
            return ScoreDotsResult::default();
        }
        let (pa, _, ia) = self.atom_state(st, a);
        let pr = self.p.probe.probe_radius;
        let r_live = w.info[a as usize].vdw_radius;
        let nearby = self.max_vdw + r_live + 2.0 * pr;
        let excl: &[u32] = self.exclude.get(&a).map(|v| v.as_slice()).unwrap_or(&[]);
        let min_occ = self.p.min_occupancy;
        buf.targets.clear();
        buf.excl.clear();
        // excluded atoms' current geometry (for trimming); one farther than
        // the dots reach plus its radius cannot cover a dot
        let fixed = !self.p.compat;
        let reach = self.dot_reach.get(&a).copied().unwrap_or(f64::INFINITY) + EXCLUSION_MARGIN;
        for &e in excl {
            let (pe, deleted, ie) = self.atom_state(st, e);
            if fixed && deleted {
                continue;
            }
            if pe.dist_sq(pa) >= (reach + ie.vdw_radius) * (reach + ie.vdw_radius) {
                continue;
            }
            buf.excl.push((pe, ie.vdw_radius));
        }
        // static neighbors
        self.static_targets(a, pa, &ia, nearby, excl, &mut buf.targets);
        // dynamic neighbors from interacting Movers
        if let Some(ms) = self.atom_movers.get(&a) {
            for &m in ms {
                let li = self.local[m as usize];
                if li == NONE || st.movers.get(li as usize) != Some(&m) {
                    continue;
                }
                let cfg = st.cfg[li as usize];
                let mv = &self.movers[m as usize];
                for slot in 0..mv.n_moved {
                    let b = mv.atoms[slot];
                    if b == a || self.mover_atom_deleted(m, slot, cfg) {
                        continue;
                    }
                    let pb = self.mover_atom_pos(m, slot, cfg);
                    let d = pb.dist(pa);
                    if d < 1e-5 || d > nearby {
                        continue;
                    }
                    if excl.contains(&b) || w.occ[b as usize].abs() < min_occ {
                        continue;
                    }
                    let ib = self.mover_atom_info(m, slot, cfg, b);
                    if d <= ia.vdw_radius + ib.vdw_radius + 2.0 * pr {
                        buf.targets.push(Target { pos: pb, info: ib });
                    }
                }
            }
        }
        // trim dots inside excluded atoms
        if let Some((cache, loc)) = trim {
            if let Some(d) = cache.borrow().get(&(a, loc)) {
                return self.scorer_score(pa, &ia, d, &buf.targets, &mut buf.prepared);
            }
        }
        if trim.is_none() {
            // all dots in coordinate arrays: trimming and scoring in one pass
            let s = &self.scorer;
            let mut ret = ScoreDotsResult::default();
            if ia.is_ion && s.ignore_ions {
                return ret;
            }
            prepare_targets(s, &ia, pa, &buf.targets, pr, &mut buf.prepared);
            score_dot_lanes(s, &ia, pa, &self.dots_soa[&a], &buf.excl, &buf.prepared, &mut buf.lanes, &mut ret);
            let density = self.p.probe.density;
            ret.bump /= density;
            ret.hbond /= density;
            ret.attract /= density;
            return ret;
        }
        let dots = &self.dots[&a];
        buf.dots.clear();
        'dot: for dp in dots.iter() {
            let p = pa + dp.d;
            for &(pe, re) in &buf.excl {
                if p.dist_sq(pe) < re * re {
                    continue 'dot;
                }
            }
            buf.dots.push(*dp);
        }
        if let Some((cache, loc)) = trim {
            cache.borrow_mut().insert((a, loc), buf.dots.clone());
        }
        self.scorer_score(pa, &ia, &buf.dots, &buf.targets, &mut buf.prepared)
    }

    fn scorer_score(&self, pa: Vec3, ia: &AtomInfo, dots: &[DotPair], targets: &[Target], prepared: &mut Vec<PTarget>) -> ScoreDotsResult {
        // Same arithmetic as DotScorer::score_dots with precomputed probe offsets.
        let pr = self.p.probe.probe_radius;
        let density = self.p.probe.density;
        let s = &self.scorer;
        let mut ret = ScoreDotsResult::default();
        if ia.is_ion && s.ignore_ions {
            return ret;
        }
        prepare_targets(s, ia, pa, targets, pr, prepared);
        for dp in dots {
            let dot_abs = pa + dp.d;
            let probe_loc = pa + dp.probe;
            accumulate_dot(s, check_dot(s, ia, dot_abs, probe_loc, prepared), &mut ret);
        }
        ret.bump /= density;
        ret.hbond /= density;
        ret.attract /= density;
        ret
    }

    /// Sum of the scores of a Mover's scored atoms under the clique state.
    pub(super) fn score_mover_atoms(&self, st: &CliqueState, m: u32, buf: &mut ScoreBuf, skip_deleted: bool) -> (f64, bool) {
        self.score_mover_atoms_in(st, m, buf, skip_deleted, &[], 0)
    }

    /// `score_mover_atoms` in fine pose `f`, with each atom scored from its
    /// shared evaluation `shared[slot]` when present.
    pub(super) fn score_mover_atoms_in(&self, st: &CliqueState, m: u32, buf: &mut ScoreBuf, skip_deleted: bool, shared: &[Option<FineShared>], f: usize) -> (f64, bool) {
        let mv = &self.movers[m as usize];
        let cfg = self.cfg_of(st, m);
        let density = self.p.probe.density;
        let mut tot = 0.0;
        let mut bad = false;
        for slot in 0..mv.n_moved {
            if skip_deleted && self.mover_atom_deleted(m, slot, cfg) {
                continue;
            }
            let a = mv.atoms[slot];
            let r = match shared.get(slot).and_then(Option::as_ref) {
                Some(sh) => {
                    let mut r = sum_addends(&sh.adds, &self.dots_soa[&a], sh.pa[f], &sh.trims[f]);
                    r.bump /= density;
                    r.hbond /= density;
                    r.attract /= density;
                    r
                }
                None => self.score_atom_ext(st, a, buf, None),
            };
            tot += r.total();
            bad |= r.has_bad_bump;
        }
        (tot, bad)
    }

    pub(super) fn pref(&self, m: u32, cfg: Cfg) -> f64 {
        let mv = &self.movers[m as usize];
        self.p.preference_magnitude
            * match cfg {
                Cfg::Coarse(c) => mv.coarse_pref[c as usize],
                Cfg::Fine(c, f) => mv.fine.pref(c as usize, f as usize),
            }
    }
}

/// A target prepared for one source atom state: the ion and alternate
/// filters applied, and what a dot test needs precomputed (the probe reach and
/// the annular test's surface point and kissing distance), with the same
/// arithmetic as the per-dot test.
#[derive(Clone, Copy)]
pub(super) struct PTarget {
    pos: Vec3,
    info: AtomInfo,
    prb2: f64,
    surface: Vec3,
    kiss: f64,
}

pub(super) fn prepare_targets(s: &DotScorer, src: &AtomInfo, src_pos: Vec3, targets: &[Target], probe_radius: f64, out: &mut Vec<PTarget>) {
    out.clear();
    let src_r = src.vdw_radius;
    for b in targets {
        let bi = &b.info;
        if bi.is_ion && s.ignore_ions {
            continue;
        }
        if !compatible_alts(src.alt, bi.alt) {
            continue;
        }
        let prb = bi.vdw_radius + probe_radius;
        let tv = b.pos - src_pos;
        let surface = tv / tv.length() * src_r + src_pos;
        let tr = bi.vdw_radius;
        let kiss = 2.0 * src_r * (tr * probe_radius / ((src_r + tr) * (src_r + probe_radius))).sqrt();
        out.push(PTarget { pos: b.pos, info: *bi, prb2: prb * prb, surface, kiss });
    }
}

/// Running state of `check_dot` over a dot's targets: the target with the
/// smallest gap so far (first one wins ties) and its hydrogen-bond status.
#[derive(Clone, Copy)]
pub(super) struct DotFold {
    best_gap: f64,
    best: Option<PTarget>,
    is_hbond: bool,
    too_close: bool,
    hb_min: f64,
    cause_dummy: bool,
}

impl DotFold {
    const START: DotFold =
        DotFold { best_gap: 1e100, best: None, is_hbond: false, too_close: false, hb_min: 0.0, cause_dummy: false };

    #[inline]
    pub(super) fn step(&mut self, s: &DotScorer, src: &AtomInfo, dot_abs: Vec3, probe_loc: Vec3, b: &PTarget) {
        let bi = &b.info;
        if probe_loc.dist_sq(b.pos) > b.prb2 {
            return;
        }
        let gap = dot_abs.dist(b.pos) - bi.vdw_radius;
        if gap < self.best_gap {
            let cs = src.charge as i32;
            let cb = bi.charge as i32;
            let both = cs != 0 && cb != 0;
            let comp = both && cs * cb < 0;
            let could = (src.is_donor && bi.is_acceptor) || (src.is_acceptor && bi.is_donor);
            if could && (!both || comp) {
                self.is_hbond = true;
                self.hb_min = if both { s.max_charged_h_overlap } else { s.max_regular_h_overlap };
                self.too_close = gap < -self.hb_min;
            } else {
                if src.is_dummy_hydrogen || bi.is_dummy_hydrogen {
                    return;
                }
                self.is_hbond = false;
                self.too_close = false;
            }
            self.cause_dummy = bi.is_dummy_hydrogen;
            self.best_gap = gap;
            self.best = Some(*b);
        }
    }

    /// (overlap type, gap, overlap, annular) of the dot, or None without a target.
    #[inline]
    pub(super) fn finish(&self, s: &DotScorer, src: &AtomInfo, dot_abs: Vec3) -> Option<(OverlapType, f64, f64, bool)> {
        let b = self.best?;
        Some(finish_dot(s, src, self.best_gap, self.is_hbond, self.too_close, self.hb_min, self.cause_dummy, &b, dot_abs))
    }
}

/// The end of `check_dot` for a dot whose nearest target is `b`.
#[allow(clippy::too_many_arguments)]
#[inline(always)]
fn finish_dot(
    s: &DotScorer,
    src: &AtomInfo,
    best_gap: f64,
    is_hbond: bool,
    too_close: bool,
    hb_min: f64,
    cause_dummy: bool,
    b: &PTarget,
    dot_abs: Vec3,
) -> (OverlapType, f64, f64, bool) {
    let mut gap = best_gap;
    let overlap;
    let mut ot;
    if gap > 0.0 {
        overlap = 0.0;
        ot = if s.weak_hbonds && is_hbond { OverlapType::HydrogenBond } else { OverlapType::NoOverlap };
    } else if is_hbond {
        if too_close {
            gap += hb_min;
            overlap = -0.5 * gap;
            ot = OverlapType::Clash;
        } else {
            overlap = -0.5 * gap;
            ot = OverlapType::HydrogenBond;
        }
    } else {
        overlap = -0.5 * gap;
        ot = OverlapType::Clash;
    }
    // annularDots: dot2srcCenter > kissEdge2bullsEye
    let d2s = (b.surface - dot_abs).length();
    let annular = d2s > b.kiss;
    if (src.is_dummy_hydrogen || cause_dummy) && (too_close || ot != OverlapType::HydrogenBond) {
        ot = OverlapType::Ignore;
    }
    (ot, gap, overlap, annular)
}

/// `DotScorer::check_dot` with a precomputed probe location.
#[inline(always)]
fn check_dot(s: &DotScorer, src: &AtomInfo, dot_abs: Vec3, probe_loc: Vec3, targets: &[PTarget]) -> Option<(OverlapType, f64, f64, bool)> {
    let mut f = DotFold::START;
    for b in targets {
        f.step(s, src, dot_abs, probe_loc, b);
    }
    f.finish(s, src, dot_abs)
}

/// Dots per chunk in `score_dot_lanes`.
const LANES: usize = 4;

/// The dot offsets of one dot sphere as coordinate arrays (dot and probe),
/// padded to whole chunks.
pub(super) struct DotsSoa {
    n: usize,
    pub(super) d: [Vec<f64>; 3],
    p: [Vec<f64>; 3],
    /// Per chunk: center and radius of a sphere holding its dot offsets, and
    /// of one holding its probe offsets.
    bounds: Vec<(Vec3, f64, Vec3, f64)>,
    /// The probe spheres of `bounds` as coordinate arrays and radii.
    probe_bounds: ChunkSpheres,
    /// The dot spheres of `bounds` as coordinate arrays and radii.
    dot_bounds: ChunkSpheres,
}

/// One bounding sphere per chunk, as coordinate arrays and radii.
struct ChunkSpheres {
    c: [Vec<f64>; 3],
    r: Vec<f64>,
}

impl ChunkSpheres {
    fn new(spheres: impl Iterator<Item = (Vec3, f64)>) -> ChunkSpheres {
        let mut out = ChunkSpheres { c: [Vec::new(), Vec::new(), Vec::new()], r: Vec::new() };
        for (c, r) in spheres {
            out.c[0].push(c.x);
            out.c[1].push(c.y);
            out.c[2].push(c.z);
            out.r.push(r);
        }
        out
    }

    /// Set bit `bit` of word `word` (of `nw` per chunk) in `masks` for each
    /// chunk whose sphere comes within `r` of `q`, a point given relative to
    /// the atom center. The spheres carry `CHUNK_MARGIN`, which keeps this a
    /// superset of the chunks that an exact test of their dots can pass.
    #[inline]
    fn mark(&self, q: Vec3, r: f64, masks: &mut [u64], nw: usize, word: usize, bit: u32) {
        let [cx, cy, cz] = &self.c;
        let n = self.r.len();
        let (cx, cy, cz, cr) = (&cx[..n], &cy[..n], &cz[..n], &self.r[..n]);
        if nw == 1 {
            let masks = &mut masks[..n];
            for c in 0..n {
                let (x, y, z) = (cx[c] - q.x, cy[c] - q.y, cz[c] - q.z);
                let reach = r + cr[c];
                masks[c] |= ((x * x + y * y + z * z <= reach * reach) as u64) << bit;
            }
        } else {
            for c in 0..n {
                let (x, y, z) = (cx[c] - q.x, cy[c] - q.y, cz[c] - q.z);
                let reach = r + cr[c];
                masks[c * nw + word] |= ((x * x + y * y + z * z <= reach * reach) as u64) << bit;
            }
        }
    }
}

/// Bit masks of chunk spheres reached: `nw` words of 64 bits per chunk, bit
/// `i` for item `i`.
#[derive(Default)]
struct ChunkMasks {
    nw: usize,
    bits: Vec<u64>,
}

impl ChunkMasks {
    fn reset(&mut self, nc: usize, items: usize) {
        self.nw = items.div_ceil(64).max(1);
        self.bits.clear();
        self.bits.resize(nc * self.nw, 0);
    }

    #[inline]
    fn mark(&mut self, spheres: &ChunkSpheres, q: Vec3, r: f64, item: usize) {
        spheres.mark(q, r, &mut self.bits, self.nw, item / 64, (item % 64) as u32);
    }

    #[inline(always)]
    fn chunk(&self, c: usize) -> &[u64] {
        &self.bits[c * self.nw..(c + 1) * self.nw]
    }
}

/// Call `f` with each item whose bit is set in `words`, in order.
#[inline(always)]
fn for_bits(words: &[u64], mut f: impl FnMut(usize)) {
    for (wi, &w) in words.iter().enumerate() {
        let mut bits = w;
        while bits != 0 {
            f(wi * 64 + bits.trailing_zeros() as usize);
            bits &= bits - 1;
        }
    }
}

/// The debugging switches read from the environment: REDUCE3_ATOMWISE
/// (score atoms one state at a time) and REDUCE3_DEBUG_ATOM (print the
/// coarse state scores of the singleton Mover of one atom). They are read
/// once.
pub(super) struct DebugSwitches {
    pub(super) atomwise: bool,
    pub(super) debug_atom: Option<String>,
}

pub(super) fn debug_switches() -> &'static DebugSwitches {
    static SWITCHES: std::sync::OnceLock<DebugSwitches> = std::sync::OnceLock::new();
    SWITCHES.get_or_init(|| DebugSwitches {
        atomwise: std::env::var_os("REDUCE3_ATOMWISE").is_some(),
        debug_atom: std::env::var("REDUCE3_DEBUG_ATOM").ok(),
    })
}

/// Slack on the chunk bounds (far above rounding error).
const CHUNK_MARGIN: f64 = 1e-6;

impl DotsSoa {
    pub(super) fn new(dots: &[DotPair]) -> DotsSoa {
        let padded = dots.len().div_ceil(LANES) * LANES;
        let col = |f: &dyn Fn(&DotPair) -> f64| {
            let mut v: Vec<f64> = dots.iter().map(f).collect();
            v.resize(padded, 0.0);
            v
        };
        let sphere = |pts: &[Vec3]| {
            let c = pts.iter().fold(Vec3::ZERO, |a, &q| a + q) / pts.len() as f64;
            (c, pts.iter().map(|q| q.dist(c)).fold(0.0, f64::max))
        };
        let bounds: Vec<(Vec3, f64, Vec3, f64)> = dots
            .chunks(LANES)
            .map(|ch| {
                let (cd, rd) = sphere(&ch.iter().map(|x| x.d).collect::<Vec<_>>());
                let (cp, rp) = sphere(&ch.iter().map(|x| x.probe).collect::<Vec<_>>());
                (cd, rd + CHUNK_MARGIN, cp, rp + CHUNK_MARGIN)
            })
            .collect();
        let probe_bounds = ChunkSpheres::new(bounds.iter().map(|b| (b.2, b.3)));
        let dot_bounds = ChunkSpheres::new(bounds.iter().map(|b| (b.0, b.1)));
        DotsSoa {
            n: dots.len(),
            d: [col(&|x| x.d.x), col(&|x| x.d.y), col(&|x| x.d.z)],
            p: [col(&|x| x.probe.x), col(&|x| x.probe.y), col(&|x| x.probe.z)],
            bounds,
            probe_bounds,
            dot_bounds,
        }
    }

    pub(super) fn dot(&self, k: usize) -> DotPair {
        DotPair { d: v3(self.d[0][k], self.d[1][k], self.d[2][k]), probe: v3(self.p[0][k], self.p[1][k], self.p[2][k]) }
    }
}

/// What `DotFold::step` needs of a target for one source atom.
#[derive(Clone, Copy)]
struct LaneTarget {
    pos: Vec3,
    prb2: f64,
    /// Probe reach (for the chunk cull).
    prb: f64,
    vdw: f64,
    hb: bool,
    hbm: f64,
    /// Fold state bits of the target: hydrogen bond, dummy cause, index.
    base: u64,
}

const LANE_NO_TARGET: u64 = u64::MAX;

/// Scratch space of `score_dot_lanes`, kept between calls.
#[derive(Default)]
pub(super) struct LaneBufs {
    targets: Vec<LaneTarget>,
    /// Per chunk, the targets that can reach its probes.
    hits: ChunkMasks,
}

/// Fill `hits` with the chunks each of `lt` can reach from the atom at `pa`.
fn target_hits(dots: &DotsSoa, pa: Vec3, lt: &[LaneTarget], hits: &mut ChunkMasks) {
    hits.reset(dots.bounds.len(), lt.len());
    for (ti, t) in lt.iter().enumerate() {
        hits.mark(&dots.probe_bounds, t.pos - pa, t.prb, ti);
    }
}

/// The per-target constants of `DotFold::step` for `src`, in order; a target
/// the step leaves alone for every dot is dropped.
fn lane_targets(s: &DotScorer, src: &AtomInfo, targets: &[PTarget], lt: &mut Vec<LaneTarget>) {
    lt.clear();
    let cs = src.charge as i32;
    for (t, b) in targets.iter().enumerate() {
        let bi = &b.info;
        let cb = bi.charge as i32;
        let both = cs != 0 && cb != 0;
        let comp = both && cs * cb < 0;
        let could = (src.is_donor && bi.is_acceptor) || (src.is_acceptor && bi.is_donor);
        let hb = could && (!both || comp);
        if !hb && (src.is_dummy_hydrogen || bi.is_dummy_hydrogen) {
            continue;
        }
        let hbm = if both { s.max_charged_h_overlap } else { s.max_regular_h_overlap };
        let base = ((t as u64) << 8) | hb as u64 | ((bi.is_dummy_hydrogen as u64) << 2);
        lt.push(LaneTarget { pos: b.pos, prb2: b.prb2, prb: b.prb2.sqrt(), vdw: bi.vdw_radius, hb, hbm, base });
    }
}

/// The dots and probes of one chunk with the atom at `pa`.
struct Chunk {
    a: [[f64; LANES]; 3],
    q: [[f64; LANES]; 3],
}

impl Chunk {
    #[inline(always)]
    fn at(pa: Vec3, dots: &DotsSoa, start: usize) -> Chunk {
        let mut c = Chunk { a: [[0.0; LANES]; 3], q: [[0.0; LANES]; 3] };
        let pc = [pa.x, pa.y, pa.z];
        for ax in 0..3 {
            for l in 0..LANES {
                c.a[ax][l] = pc[ax] + dots.d[ax][start + l];
                c.q[ax][l] = pc[ax] + dots.p[ax][start + l];
            }
        }
        c
    }

    /// Lanes whose dot lies within (strictly) `r2` of `c`.
    #[inline(always)]
    fn dot_within(&self, c: Vec3, r2: f64) -> [bool; LANES] {
        let mut out = [false; LANES];
        for (l, o) in out.iter_mut().enumerate() {
            let (x, y, z) = (self.a[0][l] - c.x, self.a[1][l] - c.y, self.a[2][l] - c.z);
            *o = x * x + y * y + z * z < r2;
        }
        out
    }

    /// Lanes whose probe lies within `r2` of `c`.
    #[inline(always)]
    fn probe_within(&self, c: Vec3, r2: f64) -> [bool; LANES] {
        let mut out = [false; LANES];
        for (l, o) in out.iter_mut().enumerate() {
            let (x, y, z) = (self.q[0][l] - c.x, self.q[1][l] - c.y, self.q[2][l] - c.z);
            *o = x * x + y * y + z * z <= r2;
        }
        out
    }

    /// Live lanes of the first `m`: those outside every atom of `excl`. Only
    /// atoms that reach the chunk's dot sphere (center `dc`, radius `rd`) can
    /// remove one.
    #[inline(always)]
    fn live(&self, m: usize, excl: &[(Vec3, f64)], dc: Vec3, rd: f64) -> [bool; LANES] {
        let mut live = [false; LANES];
        for (l, v) in live.iter_mut().enumerate() {
            *v = l < m;
        }
        for &(pe, re) in excl {
            let reach = re + rd;
            if dc.dist_sq(pe) < reach * reach {
                let inside = self.dot_within(pe, re * re);
                for l in 0..LANES {
                    live[l] &= !inside[l];
                }
            }
        }
        live
    }

    #[inline(always)]
    fn dot(&self, l: usize) -> Vec3 {
        v3(self.a[0][l], self.a[1][l], self.a[2][l])
    }
}

/// `DotFold` of each lane of a chunk, as the target index and state bits of
/// `LaneTarget::base`.
#[derive(Clone, Copy)]
struct LaneFold {
    bg: [f64; LANES],
    st: [u64; LANES],
    hm: [f64; LANES],
}

impl LaneFold {
    const START: LaneFold = LaneFold { bg: [1e100; LANES], st: [LANE_NO_TARGET; LANES], hm: [0.0; LANES] };

    /// `DotFold::step` on the live lanes.
    #[inline(always)]
    fn step(&mut self, c: &Chunk, live: &[bool; LANES], t: &LaneTarget) {
        let mut within = c.probe_within(t.pos, t.prb2);
        for l in 0..LANES {
            within[l] &= live[l];
        }
        if !within.iter().any(|&x| x) {
            return;
        }
        let neg_hbm = -t.hbm;
        for l in 0..LANES {
            let (rx, ry, rz) = (c.a[0][l] - t.pos.x, c.a[1][l] - t.pos.y, c.a[2][l] - t.pos.z);
            let gap = (rx * rx + ry * ry + rz * rz).sqrt() - t.vdw;
            let upd = within[l] & (gap < self.bg[l]);
            let new_state = t.base | (((t.hb & (gap < neg_hbm)) as u64) << 1);
            self.bg[l] = if upd { gap } else { self.bg[l] };
            self.st[l] = if upd { new_state } else { self.st[l] };
            self.hm[l] = if upd & t.hb { t.hbm } else { self.hm[l] };
        }
    }

    /// `DotFold::finish` for lane `l`, which has a target (`targets` as
    /// indexed by `base`).
    #[inline(always)]
    fn finish(&self, l: usize, s: &DotScorer, src: &AtomInfo, targets: &[PTarget], dot_abs: Vec3) -> (OverlapType, f64, f64, bool) {
        let st = self.st[l];
        let b = &targets[(st >> 8) as usize];
        finish_dot(s, src, self.bg[l], st & 1 != 0, st & 2 != 0, self.hm[l], st & 4 != 0, b, dot_abs)
    }

    /// The `DotFold` of lane `l` (`targets` as indexed by `base`).
    #[inline(always)]
    fn fold(&self, l: usize, targets: &[PTarget]) -> DotFold {
        let st = self.st[l];
        if st == LANE_NO_TARGET {
            return DotFold::START;
        }
        DotFold {
            best_gap: self.bg[l],
            best: Some(targets[(st >> 8) as usize]),
            is_hbond: st & 1 != 0,
            too_close: st & 2 != 0,
            hb_min: self.hm[l],
            cause_dummy: st & 4 != 0,
        }
    }
}

/// Score the dots of one atom state against `targets` in order, as
/// `check_dot` and `accumulate_dot` do one dot at a time: dots inside an atom
/// of `excl` are skipped and the rest accumulate into `ret` in dot order.
/// Dots go through in chunks of `LANES`, each test over a chunk at once (the
/// compiler vectorizes the fixed-size loops) with the scalar arithmetic, so
/// the results are the same bits; each dot's fold keeps the target order.
#[allow(clippy::too_many_arguments)]
pub(super) fn score_dot_lanes(
    s: &DotScorer,
    src: &AtomInfo,
    pa: Vec3,
    dots: &DotsSoa,
    excl: &[(Vec3, f64)],
    targets: &[PTarget],
    lanes: &mut LaneBufs,
    ret: &mut ScoreDotsResult,
) {
    fold_dot_lanes(s, src, pa, dots, excl, targets, lanes, |_, r| accumulate_dot(s, r, ret));
}

/// The chunked fold of `score_dot_lanes`, handing each dot that finds a
/// target to `emit` (dot index and result), in dot order.
#[allow(clippy::too_many_arguments)]
#[inline]
pub(super) fn fold_dot_lanes(
    s: &DotScorer,
    src: &AtomInfo,
    pa: Vec3,
    dots: &DotsSoa,
    excl: &[(Vec3, f64)],
    targets: &[PTarget],
    lanes: &mut LaneBufs,
    mut emit: impl FnMut(usize, Option<(OverlapType, f64, f64, bool)>),
) {
    let LaneBufs { targets: lt, hits } = lanes;
    lane_targets(s, src, targets, lt);
    if lt.is_empty() {
        return;
    }
    // a chunk with no target in reach of its probes cannot score
    target_hits(dots, pa, lt, hits);
    for (ci, &(cd, rd, _, _)) in dots.bounds.iter().enumerate() {
        let near_t = hits.chunk(ci);
        if near_t.iter().all(|&w| w == 0) {
            continue;
        }
        let start = ci * LANES;
        let m = (dots.n - start).min(LANES);
        let c = Chunk::at(pa, dots, start);
        let live = c.live(m, excl, pa + cd, rd);
        if !live.iter().any(|&x| x) {
            continue;
        }
        let mut f = LaneFold::START;
        for_bits(near_t, |ti| f.step(&c, &live, &lt[ti]));
        // finish and accumulate in dot order
        for l in 0..m {
            if live[l] && f.st[l] != LANE_NO_TARGET {
                emit(start + l, Some(f.finish(l, s, src, targets, c.dot(l))));
            }
        }
    }
}

/// Add one dot's result to the running sums (the body of `score_dots`).
#[inline]
pub(super) fn accumulate_dot(s: &DotScorer, r: Option<(OverlapType, f64, f64, bool)>, ret: &mut ScoreDotsResult) {
    if let Some((kind, val)) = dot_addend(s, r) {
        add_to(ret, kind, val);
    }
}

/// Which sum of `ScoreDotsResult` a dot adds to.
const ADD_ATTRACT: u8 = 0;
const ADD_BUMP: u8 = 1;
/// A bump that also sets `has_bad_bump`.
const ADD_BAD_BUMP: u8 = 2;
const ADD_HBOND: u8 = 3;

/// The sum one dot's result adds to and the amount, as `accumulate_dot`
/// adds it, or None when it adds nothing.
#[inline]
pub(super) fn dot_addend(s: &DotScorer, r: Option<(OverlapType, f64, f64, bool)>) -> Option<(u8, f64)> {
    let (ot, gap, overlap, annular) = r?;
    match s.interaction_type(ot, gap, true) {
        InteractionType::Invalid => None,
        InteractionType::WideContact | InteractionType::CloseContact | InteractionType::WeakHydrogenBond => {
            if annular {
                return None;
            }
            let sg = gap / s.gap_scale;
            Some((ADD_ATTRACT, (-sg * sg).exp()))
        }
        InteractionType::SmallOverlap | InteractionType::Bump => Some((ADD_BUMP, -s.bump_weight * overlap)),
        InteractionType::BadBump => Some((ADD_BAD_BUMP, -s.bump_weight * overlap)),
        InteractionType::StandardHydrogenBond => Some((ADD_HBOND, s.hbond_weight * overlap)),
    }
}

#[inline]
pub(super) fn add_to(ret: &mut ScoreDotsResult, kind: u8, val: f64) {
    match kind {
        ADD_ATTRACT => ret.attract += val,
        ADD_HBOND => ret.hbond += val,
        _ => {
            ret.bump += val;
            ret.has_bad_bump |= kind == ADD_BAD_BUMP;
        }
    }
}

/// A Mover atom's dots evaluated once for all of its fine poses: the addends
/// of the dots that find a target, and per pose the atom's place and the own
/// Mover's excluded atoms that can trim them.
pub(super) struct FineShared {
    adds: Vec<DotAdd>,
    pa: Vec<Vec3>,
    trims: Vec<Vec<(Vec3, f64)>>,
}

/// One dot's addend to an atom score: dot index, sum and amount.
#[derive(Clone, Copy)]
pub(super) struct DotAdd {
    pub(super) k: u32,
    pub(super) kind: u8,
    pub(super) val: f64,
}

/// An atom's score at `pa` from its dot addends, leaving out the dots inside
/// an atom of `trim`, added in dot order as `score_dot_lanes` adds them.
pub(super) fn sum_addends(adds: &[DotAdd], dots: &DotsSoa, pa: Vec3, trim: &[(Vec3, f64)]) -> ScoreDotsResult {
    let mut ret = ScoreDotsResult::default();
    let [dx, dy, dz] = &dots.d;
    for ad in adds {
        let k = ad.k as usize;
        let q = v3(pa.x + dx[k], pa.y + dy[k], pa.z + dz[k]);
        if !trim.iter().any(|&(pe, re)| q.dist_sq(pe) < re * re) {
            add_to(&mut ret, ad.kind, ad.val);
        }
    }
    ret
}

/// Largest distance between an atom's positions in several poses for one
/// evaluation of its dots to serve them all: far above the rounding error of
/// a pose that should leave the atom in place, and far below any distance
/// that moves a score by more than rounding.
pub(super) const SHARE_TOL: f64 = 1e-9;

/// A sphere through which another Mover of the clique reaches a dot: one of
/// its excluded atoms in some state, tested on the dot (`probe` false, strict
/// test), or one of its targets in some state, tested on the dot's probe.
#[derive(Clone, Copy)]
pub(super) struct DepSphere {
    /// Bit of the Mover in the atom's list of other Movers.
    pub(super) k: u32,
    pub(super) probe: bool,
    pub(super) c: Vec3,
    pub(super) r: f64,
    pub(super) r2: f64,
}

/// Scratch space of `classify_dots`, kept between calls.
#[derive(Default)]
pub(super) struct ClassifyBufs {
    pub(super) lanes: LaneBufs,
    /// Per chunk, the spheres that can reach it.
    shits: ChunkMasks,
}

/// Sort the dots of one atom state for `atom_factors_dotwise`, in dot order,
/// with the chunked tests of `fold_dot_lanes`. `targets` holds the static
/// targets (the first `n_static`) and then the own Mover's. A dot inside an
/// atom of `excl` drops out. A dot that no sphere of `spheres` reaches goes
/// to `base` with its result over all of `targets`; any other goes to `group`
/// with the bits (`DepSphere::k`) of the Movers that reach it and its fold
/// over the static targets alone.
#[allow(clippy::too_many_arguments)]
pub(super) fn classify_dots(
    s: &DotScorer,
    src: &AtomInfo,
    pa: Vec3,
    dots: &DotsSoa,
    excl: &[(Vec3, f64)],
    targets: &[PTarget],
    n_static: usize,
    spheres: &[DepSphere],
    bufs: &mut ClassifyBufs,
    mut base: impl FnMut(usize, Option<(OverlapType, f64, f64, bool)>),
    mut group: impl FnMut(u64, usize, DotFold),
) {
    let ClassifyBufs { lanes: LaneBufs { targets: lt, hits }, shits } = bufs;
    lane_targets(s, src, targets, lt);
    // the chunks each target and each sphere can reach
    target_hits(dots, pa, lt, hits);
    shits.reset(dots.bounds.len(), spheres.len());
    for (si, sp) in spheres.iter().enumerate() {
        let b = if sp.probe { &dots.probe_bounds } else { &dots.dot_bounds };
        shits.mark(b, sp.c - pa, sp.r, si);
    }
    for (ci, &(cd, rd, _, _)) in dots.bounds.iter().enumerate() {
        // a chunk that meets neither a target nor another Mover adds nothing
        let (near_t, near_s) = (hits.chunk(ci), shits.chunk(ci));
        if near_t.iter().chain(near_s).all(|&w| w == 0) {
            continue;
        }
        let start = ci * LANES;
        let m = (dots.n - start).min(LANES);
        let c = Chunk::at(pa, dots, start);
        let live = c.live(m, excl, pa + cd, rd);
        if !live.iter().any(|&x| x) {
            continue;
        }
        let mut dep = [0u64; LANES];
        for_bits(near_s, |si| {
            let sp = &spheres[si];
            let hit = if sp.probe { c.probe_within(sp.c, sp.r2) } else { c.dot_within(sp.c, sp.r2) };
            for l in 0..LANES {
                dep[l] |= ((hit[l] & live[l]) as u64) << sp.k;
            }
        });
        // the static targets, then the own Mover's
        let mut f = LaneFold::START;
        let mut static_fold = f;
        let mut own_started = false;
        for_bits(near_t, |ti| {
            let t = &lt[ti];
            if !own_started && (t.base >> 8) as usize >= n_static {
                static_fold = f;
                own_started = true;
            }
            f.step(&c, &live, t);
        });
        if !own_started {
            static_fold = f;
        }
        for l in 0..m {
            if !live[l] {
                continue;
            }
            if dep[l] == 0 {
                if f.st[l] != LANE_NO_TARGET {
                    base(start + l, Some(f.finish(l, s, src, targets, c.dot(l))));
                }
            } else {
                group(dep[l], start + l, static_fold.fold(l, targets));
            }
        }
    }
}

#[derive(Default)]
pub(super) struct ScoreBuf {
    targets: Vec<Target>,
    prepared: Vec<PTarget>,
    lanes: LaneBufs,
    excl: Vec<(Vec3, f64)>,
    dots: Vec<DotPair>,
}
