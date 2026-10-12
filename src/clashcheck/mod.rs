//! All-atom Probe2 contacts on an explicitly typed chemical graph.
//!
//! The orchestration adapts mmtbx/programs/probe2.py and probe/Helpers.py,
//! Copyright 2021-2023 Richardson Lab at Duke, under Apache-2.0. Existing
//! Probe scoring retains the CCTBX notices in the repository NOTICE file.
//! This module never infers missing chemical types or alters real atom sites.

mod evaluate;
mod neighbors;
mod prepare;

use crate::geom::Vec3;
use crate::probe::{InteractionType, ProbeParams};
use crate::world::{AtomLabels, World};

/// Probe2 selection and reporting options used by clashscore2.
#[derive(Clone, Debug)]
pub struct Options {
    pub probe: ProbeParams,
    pub nuclear: bool,
    pub occupancy_fraction: f64,
    pub minimum_occupancy: f64,
    pub bonded_depth: usize,
    pub include_water_water: bool,
    pub minimum_water_hydrogen_occupancy: f64,
    pub maximum_water_hydrogen_b: f64,
    /// Retain all contact classes; false retains clashes and hydrogen bonds.
    pub report_vdws: bool,
}

impl Default for Options {
    fn default() -> Self {
        let mut probe = ProbeParams::reduce2_defaults();
        probe.bump_weight = 10.0;
        probe.hydrogen_bond_weight = 4.0;
        Self {
            probe,
            nuclear: false,
            occupancy_fraction: 0.1,
            minimum_occupancy: 0.02,
            bonded_depth: 4,
            include_water_water: false,
            minimum_water_hydrogen_occupancy: 0.25,
            maximum_water_hydrogen_b: 80.0,
            report_vdws: false,
        }
    }
}

impl Options {
    /// Match clashscore2's actual per-model occupancy convention. It passes
    /// the raw maximum occupancy to a helper that divides that value by 100.
    pub fn for_clashscore(world: &World, nuclear: bool) -> Self {
        let maximum = world.occ.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        Self {
            nuclear,
            occupancy_fraction: (maximum / 100.0).min(0.1),
            ..Self::default()
        }
    }
}

/// Directed contacts of one class, aggregated over surface dots.
#[derive(Clone, Debug)]
pub struct Contact {
    pub source: u32,
    pub target: u32,
    pub interaction_type: InteractionType,
    pub dot_count: usize,
    pub atom_gap: f64,
    pub min_dot_gap: f64,
    /// Position of the dot with the lowest dot gap.
    pub location: Vec3,
    pub score_over_density: f64,
}

/// A unique physical atom pair that survives hydrogen-bond suppression.
#[derive(Clone, Debug)]
pub struct Clash {
    pub first: u32,
    pub second: u32,
    /// Rounded Probe2 dot gap, as consumed by clashscore2.
    pub overlap: f64,
    pub location: Vec3,
    pub max_b_factor: f64,
}

/// Complete contact results and labels for real and temporary phantom atoms.
pub struct Report {
    pub atoms: Vec<AtomLabels>,
    pub n_real: usize,
    pub n_selected_nonwater_atoms: usize,
    pub phantom_count: usize,
    pub contacts: Vec<Contact>,
    pub clashes: Vec<Clash>,
    pub clashscore: f64,
}

/// Validate chemical coverage supplied by the caller, prepare temporary atoms,
/// and evaluate contacts. Missing typing or symmetry must fail before this call.
pub fn run(mut world: World, options: &Options) -> Result<Report, String> {
    let n = world.len();
    if n == 0
        || n != world.n_real
        || world.info.len() != n
        || world.bonded.len() != n
        || world.labels.len() != n
        || world.occ.len() != n
        || world.b.len() != n
    {
        return Err("Probe inputs must cover every real atom without existing phantoms".into());
    }
    if world
        .labels
        .iter()
        .any(|a| a.model_index != world.labels[0].model_index)
    {
        return Err("Probe model-specific scores require one selected hierarchy model".into());
    }
    if options.occupancy_fraction < 0.0
        || !options.occupancy_fraction.is_finite()
        || options.minimum_occupancy < 0.0
        || !options.minimum_occupancy.is_finite()
        || options.probe.probe_radius < 0.0
        || !options.probe.probe_radius.is_finite()
        || options.probe.density <= 0.0
        || !options.probe.density.is_finite()
        || !options.minimum_water_hydrogen_occupancy.is_finite()
        || !options.maximum_water_hydrogen_b.is_finite()
        || [
            options.probe.worse_clash_cutoff,
            options.probe.clash_cutoff,
            options.probe.contact_cutoff,
            options.probe.uncharged_hydrogen_cutoff,
            options.probe.charged_hydrogen_cutoff,
            options.probe.bump_weight,
            options.probe.hydrogen_bond_weight,
            options.probe.gap_weight,
        ]
        .iter()
        .any(|value| !value.is_finite() || *value < 0.0)
    {
        return Err("Probe selection and sphere parameters must be finite and valid".into());
    }
    for i in 0..n {
        let p = world.pos[i];
        if !p.x.is_finite()
            || !p.y.is_finite()
            || !p.z.is_finite()
            || !world.occ[i].is_finite()
            || !world.b[i].is_finite()
            || !world.info[i].vdw_radius.is_finite()
            || world.info[i].vdw_radius <= 0.0
            || world.bonded[i].iter().any(|&j| j as usize >= n)
        {
            return Err("Probe atoms need finite sites, occupancies, radii and valid bonds".into());
        }
        let dots = 4.0
            * std::f64::consts::PI
            * options.probe.density
            * world.info[i].vdw_radius
            * world.info[i].vdw_radius;
        if !dots.is_finite() || dots > 1_000_000.0 {
            return Err("Probe sphere exceeds the supported numerical resource range".into());
        }
    }
    let ignored = prepare::prepare(&mut world, options)?;
    evaluate::evaluate(world, options, &ignored)
}
