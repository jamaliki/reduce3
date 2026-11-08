//! Compat mode: a faithful port of Reduce2's OptimizerC (sequential; vertex-cut
//! clique search; trimmed-dot and atom-score caches with the same keys and
//! evaluation order). Used to verify Reduce3 against Reduce2 exactly.

use super::*;

type TrimCache = std::cell::RefCell<FxHashMap<(u32, u32), Vec<DotPair>>>;

struct CompatOpt<'c, 'a> {
    ctx: &'c Ctx<'a>,
    st: CliqueState<'c>,
    coarse_loc: Vec<usize>,
    high: Vec<f64>,
    trim: TrimCache,
    score_cache: Option<FxHashMap<(u32, SmallVec<[u16; 6]>), f64>>,
    buf: ScoreBuf,
    calculated: usize,
    cached: usize,
    edges: &'c [(usize, usize)],
}

impl<'c, 'a> CompatOpt<'c, 'a> {
    fn mag(&self) -> f64 {
        self.ctx.p.preference_magnitude
    }

    fn score_atom(&mut self, a: u32, loc: u32) -> f64 {
        self.calculated += 1;
        self.ctx.score_atom_ext(&self.st, a, &mut self.buf, Some((&self.trim, loc))).total()
    }

    fn score_atom_cached(&mut self, a: u32, loc: u32) -> f64 {
        let key: SmallVec<[u16; 6]> = match self.ctx.atom_movers.get(&a) {
            Some(ms) => ms.iter().map(|&m| self.coarse_loc[m as usize] as u16).collect(),
            None => SmallVec::new(),
        };
        if let Some(&v) = self.score_cache.as_ref().unwrap().get(&(a, key.clone())) {
            self.cached += 1;
            return v;
        }
        let v = self.score_atom(a, loc);
        self.score_cache.as_mut().unwrap().insert((a, key), v);
        v
    }

    /// `scorePosition(states, index, offset)` for coarse states.
    fn score_position(&mut self, m: usize, index: usize, offset: u32) -> f64 {
        let mv = &self.ctx.movers[m];
        let mut ret = 0.0;
        for slot in 0..mv.n_moved {
            let d = in_state(&mv.coarse_del, index);
            if slot < d.len() && d[slot] {
                continue;
            }
            let a = mv.atoms[slot];
            ret += if self.score_cache.is_some() {
                self.score_atom_cached(a, index as u32 + offset)
            } else {
                self.score_atom(a, index as u32 + offset)
            };
        }
        ret
    }

    fn set(&mut self, m: usize, c: usize) {
        self.st.cfg[m] = Cfg::Coarse(c as u16);
    }

    fn all_states(num: &[usize]) -> Vec<Vec<usize>> {
        // generateAllStates: the first entry varies fastest
        let total: usize = num.iter().product();
        let mut out = Vec::with_capacity(total);
        let mut cur = vec![0usize; num.len()];
        for _ in 0..total {
            out.push(cur.clone());
            for k in 0..num.len() {
                cur[k] += 1;
                if cur[k] < num[k] {
                    break;
                }
                cur[k] = 0;
            }
        }
        out
    }

    fn brute_force(&mut self, movers: &[usize]) -> f64 {
        let num: Vec<usize> = movers.iter().map(|&m| self.ctx.movers[m].num_coarse()).collect();
        let mut best = -1e100;
        let mut best_state: Option<Vec<usize>> = None;
        let mag = self.mag();
        for cur in Self::all_states(&num) {
            for (k, &m) in movers.iter().enumerate() {
                if self.coarse_loc[m] != cur[k] {
                    self.set(m, cur[k]);
                    self.coarse_loc[m] = cur[k];
                }
            }
            let mut score = 0.0;
            for (k, &m) in movers.iter().enumerate() {
                score += mag * self.ctx.movers[m].coarse_pref[cur[k]];
                score += self.score_position(m, cur[k], 0);
            }
            if score > best || best_state.is_none() {
                best = score;
                best_state = Some(cur.clone());
            }
        }
        self.finish(movers, &best_state.unwrap())
    }

    fn finish(&mut self, movers: &[usize], best: &[usize]) -> f64 {
        for (k, &m) in movers.iter().enumerate() {
            self.set(m, best[k]);
            self.coarse_loc[m] = best[k];
        }
        let mag = self.mag();
        let mut ret = 0.0;
        for (k, &m) in movers.iter().enumerate() {
            let mut my = mag * self.ctx.movers[m].coarse_pref[best[k]];
            my += self.score_position(m, best[k], 0);
            self.high[m] = my;
            ret += my;
        }
        ret
    }

    fn components(verts: &[usize], edges: &[(usize, usize)]) -> Vec<Vec<usize>> {
        // connected components numbered by first vertex (boost DFS order)
        let mut comp_of: FxHashMap<usize, usize> = FxHashMap::default();
        let mut comps: Vec<Vec<usize>> = Vec::new();
        for &v in verts {
            if comp_of.contains_key(&v) {
                continue;
            }
            let id = comps.len();
            let mut stack = vec![v];
            comp_of.insert(v, id);
            let mut members = vec![];
            while let Some(x) = stack.pop() {
                members.push(x);
                for &(a, b) in edges {
                    let y = if a == x { b } else if b == x { a } else { continue };
                    if verts.contains(&y) && !comp_of.contains_key(&y) {
                        comp_of.insert(y, id);
                        stack.push(y);
                    }
                }
            }
            comps.push(members);
        }
        // members listed in vertex order
        for c in comps.iter_mut() {
            c.sort_by_key(|x| verts.iter().position(|y| y == x).unwrap());
        }
        comps
    }

    fn n_choose_m(n: usize, m: usize) -> Vec<Vec<usize>> {
        let mut result = Vec::new();
        let mut idx: Vec<usize> = (1..=m).collect();
        loop {
            result.push(idx.iter().map(|i| i - 1).collect());
            let mut i = m as isize - 1;
            while i >= 0 && idx[i as usize] == n - m + i as usize + 1 {
                i -= 1;
            }
            if i < 0 {
                break;
            }
            idx[i as usize] += 1;
            for j in (i as usize + 1)..m {
                idx[j] = idx[j - 1] + 1;
            }
        }
        result
    }

    fn vertex_cut(&mut self, movers: &[usize]) -> f64 {
        if movers.len() <= 2 {
            return self.brute_force(movers);
        }
        let edges: Vec<(usize, usize)> =
            self.edges.iter().copied().filter(|(a, b)| movers.contains(a) && movers.contains(b)).collect();
        // findVertexCut
        let n = movers.len();
        let mut cut: Vec<usize> = Vec::new();
        let mut keep: Vec<usize> = movers.to_vec();
        'search: for size in 0..n {
            for set in Self::n_choose_m(n, size) {
                let k: Vec<usize> = (0..n).filter(|i| !set.contains(i)).map(|i| movers[i]).collect();
                if Self::components(&k, &edges).len() > 1 {
                    cut = set.iter().map(|&i| movers[i]).collect();
                    keep = k;
                    break 'search;
                }
            }
        }
        if cut.is_empty() {
            return self.brute_force(movers);
        }
        let num: Vec<usize> = cut.iter().map(|&m| self.ctx.movers[m].num_coarse()).collect();
        let mut best = -1e100;
        let mut best_state: Option<Vec<usize>> = None;
        let mag = self.mag();
        let comps = Self::components(&keep, &edges);
        for cur in Self::all_states(&num) {
            for (k, &m) in cut.iter().enumerate() {
                if self.coarse_loc[m] != cur[k] {
                    self.set(m, cur[k]);
                    self.coarse_loc[m] = cur[k];
                }
            }
            let mut score = 0.0;
            for comp in &comps {
                score += self.vertex_cut(comp);
            }
            for (k, &m) in cut.iter().enumerate() {
                score += mag * self.ctx.movers[m].coarse_pref[cur[k]];
                score += self.score_position(m, cur[k], 0);
            }
            if score > best || best_state.is_none() {
                best = score;
                best_state = Some(movers.iter().map(|&m| self.coarse_loc[m]).collect());
            }
        }
        self.finish(movers, &best_state.unwrap())
    }
}

/// Run Reduce2's optimizer exactly (compat mode). Returns per-component results.
pub(super) fn compat_optimize(ctx: &Ctx, components: &[Vec<u32>], edges: &[(usize, usize)], v: i32) -> Vec<CliqueResult> {
    let nm = ctx.movers.len();
    let all: Vec<u32> = (0..nm as u32).collect();
    let mut co = CompatOpt {
        ctx,
        st: CliqueState { movers: &all, cfg: vec![Cfg::Coarse(0); nm] },
        coarse_loc: vec![0; nm],
        high: vec![0.0; nm],
        trim: std::cell::RefCell::new(FxHashMap::default()),
        score_cache: None,
        buf: ScoreBuf::default(),
        calculated: 0,
        cached: 0,
        edges,
    };
    let mag = ctx.p.preference_magnitude;
    // Initialize
    let mut initial = vec![0.0; nm];
    for m in 0..nm {
        let mut s = mag * ctx.movers[m].coarse_pref[0];
        s += co.score_position(m, 0, 0);
        co.high[m] = s;
        initial[m] = s;
    }
    let mut coarse_best = vec![0.0; components.len()];
    // singletons first, then groups (the original's order)
    for (ci, comp) in components.iter().enumerate() {
        if comp.len() != 1 {
            continue;
        }
        let m = comp[0] as usize;
        let n = ctx.movers[m].num_coarse();
        let mut scores: Vec<f64> = ctx.movers[m].coarse_pref.iter().map(|p| p * mag).collect();
        for i in 0..n {
            co.set(m, i);
            scores[i] += co.score_position(m, i, 0);
        }
        let mut max = scores[0];
        let mut mi = 0;
        for i in 1..n {
            if scores[i] > max {
                max = scores[i];
                mi = i;
            }
        }
        co.set(m, mi);
        co.coarse_loc[m] = mi;
        co.high[m] = max;
        coarse_best[ci] = max;
    }
    for (ci, comp) in components.iter().enumerate() {
        if comp.len() == 1 {
            continue;
        }
        co.score_cache = Some(FxHashMap::default());
        let movers: Vec<usize> = comp.iter().map(|&m| m as usize).collect();
        coarse_best[ci] = co.vertex_cut(&movers);
        co.score_cache = None;
    }
    // Fine optimization in global Mover order
    let mut fine_info: Vec<(u32, String)> = Vec::new();
    for m in 0..nm {
        let mv = &ctx.movers[m];
        let c = co.coarse_loc[m];
        let initial_score = co.score_position(m, c, 0);
        let nf = mv.fine.count();
        let mut msg = String::new();
        if nf > 0 {
            let mut scores: Vec<f64> = mv.fine.prefs(c).iter().map(|p| p * mag).collect();
            for i in 0..nf {
                co.st.cfg[m] = Cfg::Fine(c as u16, i as u16);
                let mut s = 0.0;
                for slot in 0..mv.n_moved {
                    s += co.score_atom(mv.atoms[slot], 100000 + i as u32);
                }
                scores[i] += s;
            }
            let mut max = scores[0];
            let mut mi = 0;
            for i in 1..nf {
                if scores[i] > max {
                    max = scores[i];
                    mi = i;
                }
            }
            if max > co.high[m] {
                if v >= 3 {
                    msg += &format!(
                        "   Setting Mover to fine orientation {}, max score = {:.2} (coarse score {:.2})\n",
                        mi, max, initial_score
                    );
                }
                co.st.cfg[m] = Cfg::Fine(c as u16, mi as u16);
                co.high[m] = max;
            } else {
                co.st.cfg[m] = Cfg::Coarse(c as u16);
                if v >= 3 {
                    msg += "   Leaving Mover at coarse orientation\n";
                }
            }
        }
        fine_info.push((m as u32, msg));
    }
    // Flip annotations, in report order does not matter for state (reset after)
    let p = ctx.p;
    let mut annot = vec![String::new(); nm];
    for m in 0..nm {
        let mv = &ctx.movers[m];
        if !mv.kind.is_flip() {
            annot[m] = " .".into();
            continue;
        }
        let n = mv.num_coarse();
        let fin = co.st.cfg[m].coarse();
        let other = (fin + n / 2) % n;
        let saved = co.st.cfg[m];
        co.st.cfg[m] = Cfg::Coarse(other as u16);
        let (os, ob) = ctx.score_mover_atoms(&co.st, m as u32, &mut co.buf, false);
        let os = os + ctx.pref(m as u32, co.st.cfg[m]);
        co.st.cfg[m] = Cfg::Coarse(fin as u16);
        let (fs, fb) = ctx.score_mover_atoms(&co.st, m as u32, &mut co.buf, false);
        let fs = fs + ctx.pref(m as u32, co.st.cfg[m]);
        co.st.cfg[m] = saved;
        let desc = mv.pose_description(fin, None, !p.skip_bond_fixup);
        annot[m] = if ob && fb {
            " BothClash".into()
        } else if desc.contains("Unflipped") && os > fs && os - fs <= p.non_flip_preference {
            " Uncertain".into()
        } else {
            " .".into()
        };
    }
    let mut out = Vec::with_capacity(components.len());
    for (ci, comp) in components.iter().enumerate() {
        out.push(CliqueResult {
            initial: comp.iter().map(|&m| initial[m as usize]).collect(),
            cfg: comp.iter().map(|&m| co.st.cfg[m as usize]).collect(),
            high: comp.iter().map(|&m| co.high[m as usize]).collect(),
            coarse_best: coarse_best[ci],
            coarse_info: String::new(),
            fine_info: comp.iter().map(|&m| fine_info[m as usize].clone()).collect(),
            annot: comp.iter().map(|&m| annot[m as usize].clone()).collect(),
            calculated: if ci == 0 { co.calculated } else { 0 },
            cached: if ci == 0 { co.cached } else { 0 },
        });
    }
    out
}
