//! Placement and optimization of Movers (port of mmtbx/reduce/Optimizers.py,
//! Optimizers.cpp and InteractionGraph.*).
//!
//! The objective is exactly Reduce2's: for every Mover, the preference energy
//! of its state plus the Probe dot score of each of its (non-deleted) atoms
//! against everything else. Cliques of interacting Movers are optimized
//! exactly with variable elimination over cached per-atom score tables
//! instead of Reduce2's exhaustive/vertex-cut search, and independent cliques
//! are processed in parallel.

use crate::geom::*;
use crate::movers::{self, in_state, Mover, MoverKind, SingleHOptions};
use crate::probe::*;
use crate::world::World;
use rustc_hash::{FxHashMap, FxHashSet};
use smallvec::SmallVec;
use std::fmt::Write as _;
use std::sync::Arc;
use std::time::Instant;

#[cfg(test)]
mod lane_agreement_tests;

mod compat;
mod dots;
mod driver;
mod elimination;
mod neighbours;
mod placement;

use compat::*;
use dots::*;
use driver::*;
use elimination::*;
use neighbours::*;
use placement::*;
pub use driver::{optimize, optimize_debug, ConformerInput, DebugMovers};

/// Behavior switches: `compat` reproduces Reduce2 quirks that are bugs.
#[derive(Clone, Debug)]
pub struct OptParams {
    pub probe: ProbeParams,
    pub add_flip_movers: bool,
    pub alt_id: Option<String>,
    pub bonded_neighbor_depth: usize,
    pub use_neutron_distances: bool,
    pub min_occupancy: f64,
    pub preference_magnitude: f64,
    pub non_flip_preference: f64,
    /// Fixed mode: preference for a hydroxyl hydrogen on a planar atom (a
    /// phenol, an enol, a carboxylic acid) to lie in its plane.
    pub planar_hydroxyl_preference: f64,
    /// Fixed mode: further preference for an acid's hydroxyl hydrogen to be
    /// syn to the carbonyl oxygen rather than anti.
    pub acid_syn_preference: f64,
    pub skip_bond_fixup: bool,
    /// Fixed mode: place the hydrogen on the alpha carbon of an Asn, Gln or His
    /// again once its flip's fix-up has docked the side chain (fix 27 in the
    /// README). Off keeps Reduce2's position, for comparisons; compat mode
    /// never places it again.
    pub reride_flipped_alpha_h: bool,
    /// Fixed mode: the closest, in A, that a hydroxyl, thiol or ammonium
    /// hydrogen may be placed to a heavy atom that is not its parent, when
    /// any turn about its bond keeps that clear (fix 29 in the README).
    /// Reduce2 aims such a hydrogen straight at a nearby acceptor, which
    /// leaves it as close as 1.3 A. 0 keeps Reduce2's poses; compat mode
    /// never applies it.
    pub rotor_clearance: f64,
    pub flip_states: String,
    pub verbosity: i32,
    pub compat: bool,
    /// Build the report text. Without it the optimizer skips the scores and
    /// strings that only the report shows; the chosen states, deletions and
    /// fix-ups stay the same.
    pub report: bool,
}

impl Default for OptParams {
    fn default() -> Self {
        OptParams {
            probe: ProbeParams::reduce2_defaults(),
            // Reduce2 defaults to False; flips are on unless asked otherwise.
            add_flip_movers: true,
            alt_id: None,
            bonded_neighbor_depth: 4,
            use_neutron_distances: false,
            min_occupancy: 0.02,
            preference_magnitude: 1.0,
            non_flip_preference: 0.5,
            planar_hydroxyl_preference: 1.0,
            acid_syn_preference: 1.0,
            skip_bond_fixup: false,
            reride_flipped_alpha_h: true,
            rotor_clearance: 1.5,
            flip_states: String::new(),
            verbosity: 2,
            compat: false,
            report: true,
        }
    }
}

#[derive(Clone, Debug)]
pub struct FlipMoverState {
    pub mover_type: String,
    pub model_num: i64,
    pub alt_id: String,
    pub any_alt: bool,
    pub chain: String,
    pub res_name: String,
    pub res_id: i32,
    pub icode: String,
    pub flipped: bool,
    pub fixed_up: bool,
}

pub fn parse_flip_states(s: &str, compat: bool) -> Vec<FlipMoverState> {
    let mut ret = Vec::new();
    for state in s.split(',') {
        let words: Vec<&str> = state.split_whitespace().collect();
        if words.len() != 6 && words.len() != 7 {
            continue;
        }
        let model_num: i64 = words[0].parse().unwrap_or(1);
        let mut alt = words[1].to_ascii_lowercase();
        let any_alt = alt == ".";
        if any_alt {
            alt.clear();
        }
        let chain = if compat { words[2].to_ascii_uppercase() } else { words[2].to_string() };
        let res_name = words[3].to_ascii_uppercase();
        let t = if res_name == "HIS" { "HisFlip" } else { "AmideFlip" };
        let rid = words[4];
        let (res_id, icode) = match rid.parse::<i32>() {
            Ok(v) => (v, String::new()),
            Err(_) => (rid[..rid.len() - 1].parse().unwrap_or(0), rid[rid.len() - 1..].to_string()),
        };
        let flipped = words[5] == "Flipped";
        let fixed_up = words.len() == 7 && words[6] == "AnglesAdjusted";
        ret.push(FlipMoverState {
            mover_type: t.into(),
            // compat: the original adds one here and one more when comparing.
            model_num: if compat { model_num + 1 } else { model_num },
            alt_id: alt,
            any_alt,
            chain,
            res_name,
            res_id,
            icode,
            flipped,
            fixed_up,
        });
    }
    ret
}

fn find_flip_state<'a>(w: &World, a: u32, states: &'a [FlipMoverState], compat: bool, n_models: usize) -> Option<&'a FlipMoverState> {
    let l = &w.labels[a as usize];
    for fs in states {
        let model_ok = if compat {
            // The original compares a string model id against an integer, so
            // only files without MODEL records ever match.
            l.model_id.is_empty()
        } else {
            (l.model_index as i64 + 1 == fs.model_num) || (n_models == 1 && fs.model_num == 1)
        };
        let alt = l.altloc.trim();
        let alt_ok = if compat {
            alt.is_empty() || alt.to_ascii_lowercase() == fs.alt_id.to_ascii_lowercase()
        } else {
            fs.any_alt || alt.is_empty() || alt.to_ascii_lowercase() == fs.alt_id
        };
        if model_ok && *l.chain == *fs.chain && l.resseq == fs.res_id && alt_ok && l.icode.trim() == fs.icode.trim() {
            return Some(fs);
        }
    }
    None
}

/// Flip performed (by Mover or lock-down), for reporting/Flipkins.
#[derive(Clone, Debug)]
pub struct FlippedMoverInfo {
    pub alt: String,
    pub base_atom: u32,
}

#[derive(Default)]
pub struct OptOutput {
    pub info: String,
    pub warnings: String,
    pub hydrogens_to_delete: Vec<u32>,
    pub amide_flips: Vec<FlippedMoverInfo>,
    pub his_flips: Vec<FlippedMoverInfo>,
    /// Alpha carbons about which a flip's fix-up docked the side chain, by
    /// Mover or lock-down, in the order applied.
    pub docked_alpha_carbons: Vec<u32>,
    pub num_calculated: usize,
    pub num_cached: usize,
}

fn vcheck(verbosity: i32, level: i32, msg: &str) -> String {
    if verbosity >= level { format!("{}{}", " ".repeat(level as usize), msg) } else { String::new() }
}

struct Timer {
    t: Instant,
    verbosity: i32,
}
impl Timer {
    fn new(v: i32) -> Self {
        Timer { t: Instant::now(), verbosity: v }
    }
    fn report(&mut self, msg: &str) -> String {
        let d = self.t.elapsed().as_secs_f64();
        self.t = Instant::now();
        vcheck(self.verbosity, 2, &format!("Time to {}: {:.3}\n", msg, d))
    }
}

/// Python `str(float)`.
pub fn py_float(v: f64) -> String {
    let s = format!("{:?}", v);
    s
}

/// `getAtomsWithinNBonds` from probe Helpers.
pub fn atoms_within_n_bonds(w: &World, atom: u32, probe_rad: f64, n: usize, non_h_n: usize) -> Vec<u32> {
    let a_loc = w.pos[atom as usize];
    let a_rad = w.info[atom as usize].vdw_radius;
    let atom_is_h = w.is_h(atom);
    let n = if w.info[atom as usize].is_dummy_hydrogen { 1 } else { n };
    let mut set: Vec<u32> = vec![atom];
    for i in 0..n {
        // the atoms found before this round
        for k in 0..set.len() {
            let a = set[k];
            for &nb in &w.bonded[a as usize] {
                if i < non_h_n || atom_is_h || w.is_h(nb) {
                    if w.compatible(atom, nb) {
                        let d = w.pos[nb as usize].dist(a_loc);
                        if d <= a_rad + w.info[nb as usize].vdw_radius + 2.0 * probe_rad && !set.contains(&nb) {
                            set.push(nb);
                        }
                    }
                }
            }
        }
    }
    set.retain(|&x| x != atom);
    set
}

/// `fixupExplicitDonors`.
fn fixup_explicit_donors(w: &mut World, atoms: &[u32]) {
    for &a in atoms {
        if w.is_h(a) {
            for k in 0..w.bonded[a as usize].len() {
                let n = w.bonded[a as usize][k];
                let e = w.elem(n);
                if e == "N" || e == "O" || e == "S" {
                    w.info[a as usize].is_donor = true;
                    w.info[n as usize].is_donor = false;
                }
            }
        } else {
            let e = w.elem(a);
            if e == "N" || e == "O" || e == "S" {
                w.info[a as usize].is_donor = false;
            }
        }
    }
}

pub use crate::atomtypes::is_aromatic_acceptor;
