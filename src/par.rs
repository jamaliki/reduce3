//! Parallelism within one model. Batch runs process many models at once and
//! run each one sequentially, which skips rayon's splitting overhead; the
//! results are the same either way.

use rayon::prelude::*;
use std::cell::Cell;

thread_local! {
    static SEQUENTIAL: Cell<bool> = const { Cell::new(false) };
}

/// Whether work on this thread runs sequentially.
pub fn sequential() -> bool {
    SEQUENTIAL.with(|s| s.get())
}

/// Run `f` with this thread's work on one model kept sequential.
pub fn run_sequential<R>(f: impl FnOnce() -> R) -> R {
    let prev = SEQUENTIAL.with(|s| s.replace(true));
    let r = f();
    SEQUENTIAL.with(|s| s.set(prev));
    r
}

/// `items.par_iter().map(f).collect()`, or the sequential equivalent.
pub fn map_collect<T, R, C, F>(items: &[T], f: F) -> C
where
    T: Sync,
    R: Send,
    F: Fn(&T) -> R + Sync + Send,
    C: FromIterator<R> + FromParallelIterator<R>,
{
    if sequential() {
        items.iter().map(f).collect()
    } else {
        items.par_iter().map(f).collect()
    }
}

/// `items.par_iter().flat_map_iter(f).collect()`, or the sequential equivalent.
pub fn flat_map_collect<T, R, I, F>(items: &[T], f: F) -> Vec<R>
where
    T: Sync,
    R: Send,
    I: Iterator<Item = R>,
    F: Fn(&T) -> I + Sync + Send,
{
    if sequential() {
        items.iter().flat_map(f).collect()
    } else {
        items.par_iter().flat_map_iter(f).collect()
    }
}

/// An unstable sort; callers' keys cover whole elements, so the order is unique.
pub fn sort_unstable_by_key<T: Send, K: Ord, F: Fn(&T) -> K + Sync>(v: &mut [T], key: F) {
    if sequential() {
        v.sort_unstable_by_key(key);
    } else {
        v.par_sort_unstable_by_key(key);
    }
}
