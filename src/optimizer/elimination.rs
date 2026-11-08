//! Exact clique optimization by variable elimination

pub(super) struct Factor {
    pub(super) scope: Vec<usize>,
    pub(super) dims: Vec<usize>,
    pub(super) table: Vec<f64>,
}

impl Factor {
    pub(super) fn size(dims: &[usize]) -> usize {
        dims.iter().product()
    }
}

/// Index of an assignment (`vals` indexed by local variable) in a factor table.
#[inline]
pub(super) fn table_index(scope: &[usize], dims: &[usize], vals: &[usize]) -> usize {
    let mut idx = 0;
    for (k, &v) in scope.iter().enumerate() {
        idx = idx * dims[k] + vals[v];
    }
    idx
}

/// Maximize the sum of factors; returns the best assignment.
pub(super) fn variable_elimination(nvars: usize, doms: &[usize], mut factors: Vec<Factor>) -> Vec<usize> {
    let mut eliminated = vec![false; nvars];
    // (var, remaining scope, dims, argmax table)
    let mut trace: Vec<(usize, Vec<usize>, Vec<usize>, Vec<u16>)> = Vec::new();
    for _ in 0..nvars {
        // choose variable producing the smallest combined table
        let mut best_v = usize::MAX;
        let mut best_cost = usize::MAX;
        for v in 0..nvars {
            if eliminated[v] {
                continue;
            }
            let mut uni: Vec<usize> = vec![v];
            for f in &factors {
                if f.scope.contains(&v) {
                    for &s in &f.scope {
                        if !uni.contains(&s) {
                            uni.push(s);
                        }
                    }
                }
            }
            let cost: usize = uni.iter().map(|&s| doms[s]).fold(1usize, |a, b| a.saturating_mul(b));
            if cost < best_cost {
                best_cost = cost;
                best_v = v;
            }
        }
        let v = best_v;
        eliminated[v] = true;
        let (with_v, rest): (Vec<Factor>, Vec<Factor>) = factors.into_iter().partition(|f| f.scope.contains(&v));
        factors = rest;
        let mut uni: Vec<usize> = Vec::new();
        for f in &with_v {
            for &s in &f.scope {
                if s != v && !uni.contains(&s) {
                    uni.push(s);
                }
            }
        }
        uni.sort_unstable();
        let udims: Vec<usize> = uni.iter().map(|&s| doms[s]).collect();
        let n = Factor::size(&udims);
        let mut table = vec![f64::NEG_INFINITY; n];
        let mut arg = vec![0u16; n];
        let mut vals = vec![0usize; nvars];
        for idx in 0..n {
            // decode idx into vals for the union scope
            let mut r = idx;
            for k in (0..uni.len()).rev() {
                vals[uni[k]] = r % udims[k];
                r /= udims[k];
            }
            let mut best = f64::NEG_INFINITY;
            let mut bi = 0u16;
            for sv in 0..doms[v] {
                vals[v] = sv;
                let mut s = 0.0;
                for f in &with_v {
                    s += f.table[table_index(&f.scope, &f.dims, &vals)];
                }
                if s > best || sv == 0 {
                    if sv == 0 || s > best {
                        best = s;
                        bi = sv as u16;
                    }
                }
            }
            table[idx] = best;
            arg[idx] = bi;
        }
        trace.push((v, uni.clone(), udims.clone(), arg));
        factors.push(Factor { scope: uni, dims: udims, table });
    }
    let mut assign = vec![0usize; nvars];
    for (v, scope, dims, arg) in trace.into_iter().rev() {
        let idx = table_index(&scope, &dims, &assign);
        assign[v] = arg[idx] as usize;
    }
    assign
}
