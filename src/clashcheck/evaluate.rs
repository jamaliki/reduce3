use super::{Clash, Contact, Options, Report};
use crate::optimizer::atoms_within_n_bonds;
use crate::probe::{DotScorer, DotSphereCache, InteractionType, OverlapType, Target, trim_dots};
use crate::world::World;
use rustc_hash::FxHashSet;
use std::collections::BTreeMap;

fn interaction_index(kind: InteractionType) -> u8 {
    match kind {
        InteractionType::WideContact => 0,
        InteractionType::CloseContact => 1,
        InteractionType::WeakHydrogenBond => 2,
        InteractionType::SmallOverlap => 3,
        InteractionType::Bump => 4,
        InteractionType::BadBump => 5,
        InteractionType::StandardHydrogenBond => 6,
        InteractionType::Invalid => 7,
    }
}

// Probe2 writes three decimal places before clashscore2 tests -0.4.
fn reported_gap(value: f64) -> Result<f64, String> {
    // Formatting rounds the binary value directly; scaling before rounding can
    // change half-decimal ties compared with Probe2's Python formatter.
    format!("{value:.3}")
        .parse()
        .map_err(|_| "Probe gap serialization failed".into())
}

fn reported_location(position: crate::geom::Vec3) -> Result<crate::geom::Vec3, String> {
    Ok(crate::geom::v3(
        reported_gap(position.x)?,
        reported_gap(position.y)?,
        reported_gap(position.z)?,
    ))
}

#[cfg(test)]
mod tests {
    use super::reported_gap;

    #[test]
    fn decimal_output_matches_python_format_on_binary_neighbors() {
        // CPython format(value, ".3f") oracle, including both neighbors of
        // the -0.4005 clash boundary. Scaling first would change this tie.
        for (value, expected) in [
            (0.4005_f64, [0.400_f64, 0.401, 0.401]),
            (-0.4005, [-0.401, -0.401, -0.400]),
            (1.2345, [1.234, 1.234, 1.235]),
            (-1.2345, [-1.235, -1.234, -1.234]),
            (2.6755, [2.675, 2.675, 2.676]),
            (0.0005, [0.000, 0.001, 0.001]),
            (-0.0005, [-0.001, -0.001, -0.000]),
            (1234.5675, [1234.567, 1234.568, 1234.568]),
        ] {
            let bits = value.to_bits();
            let neighbors = if value.is_sign_negative() {
                [f64::from_bits(bits + 1), value, f64::from_bits(bits - 1)]
            } else {
                [f64::from_bits(bits - 1), value, f64::from_bits(bits + 1)]
            };
            for (input, output) in neighbors.into_iter().zip(expected) {
                assert_eq!(reported_gap(input).unwrap().to_bits(), output.to_bits());
            }
        }
    }
}

pub(super) fn evaluate(
    world: World,
    options: &Options,
    ignored: &[bool],
) -> Result<Report, String> {
    let query = super::neighbors::Query::new(&world.pos)?;
    let max_radius = world
        .info
        .iter()
        .fold(1.0_f64, |r, info| r.max(info.vdw_radius));
    let max_range = 2.0 * (max_radius + options.probe.probe_radius);
    let mut sources: Vec<_> = (0..world.n_real as u32)
        .filter(|&i| {
            !world.labels[i as usize].is_water && world.occ[i as usize] > options.occupancy_fraction
        })
        .collect();
    let n_selected = sources.len();
    sources.sort_by_key(|&i| world.atom_id_string(i));
    let scorer = DotScorer::new(&options.probe);
    let mut spheres = DotSphereCache::new(options.probe.density);
    let mut trimmed = Vec::new();
    let mut contacts = BTreeMap::<(u32, u32, u8), Contact>::new();
    for source in sources {
        let si = source as usize;
        if ignored[si] || world.occ[si] < options.minimum_occupancy {
            continue;
        }
        let mut nearby = Vec::new();
        query.for_each_within(world.pos[si], 0.00001, max_range, |j, _, d_sq| {
            let ti = j as usize;
            let cutoff = world.info[si].vdw_radius
                + world.info[ti].vdw_radius
                + 2.0 * options.probe.probe_radius;
            if d_sq <= cutoff * cutoff
                && !ignored[ti]
                && world.occ[ti] >= options.minimum_occupancy
                && world.labels[si].model_index == world.labels[ti].model_index
                && world.compatible(source, j)
                && !(world.labels[si].is_water
                    && world.labels[ti].is_water
                    && (!options.include_water_water || world.labels[si].ag == world.labels[ti].ag))
            {
                nearby.push(j);
            }
        });
        nearby.sort_unstable();
        let excluded = atoms_within_n_bonds(
            &world,
            source,
            options.probe.probe_radius,
            options.bonded_depth,
            3,
        );
        nearby.retain(|j| !excluded.contains(j));
        let targets: Vec<_> = nearby
            .iter()
            .map(|&j| Target {
                pos: world.pos[j as usize],
                info: world.info[j as usize],
            })
            .collect();
        let exclusions: Vec<_> = excluded
            .iter()
            .map(|&j| (world.pos[j as usize], world.info[j as usize].vdw_radius))
            .collect();
        let dots = spheres.get(world.info[si].vdw_radius);
        trim_dots(world.pos[si], &dots, &exclusions, &mut trimmed);
        for &offset in &trimmed {
            let Some(dot) = scorer.check_dot(
                world.pos[si],
                &world.info[si],
                offset,
                &targets,
                options.probe.probe_radius,
            ) else {
                continue;
            };
            if dot.overlap_type == OverlapType::Ignore || (dot.gap > 0.0 && dot.annular) {
                continue;
            }
            let target = nearby[dot.target_index];
            let ti = target as usize;
            if world.occ[ti] <= options.occupancy_fraction {
                continue;
            }
            if dot.overlap_type != OverlapType::HydrogenBond
                && (matches!(world.info[si].alt, 0 | b' ')
                    || matches!(world.info[ti].alt, 0 | b' '))
                && world.occ[si] + world.occ[ti] <= 1.0
            {
                continue;
            }
            let kind = scorer.interaction_type(dot.overlap_type, dot.gap, false);
            if !options.report_vdws && dot.overlap_type == OverlapType::NoOverlap {
                continue;
            }
            let location = world.pos[si] + offset;
            let atom_gap = world.pos[si].dist(world.pos[ti])
                - world.info[si].vdw_radius
                - world.info[ti].vdw_radius;
            let spike_length = dot.overlap.abs();
            let score = match kind {
                InteractionType::WideContact => {
                    (-(dot.gap / options.probe.gap_weight).powi(2)).exp()
                }
                InteractionType::SmallOverlap
                | InteractionType::Bump
                | InteractionType::BadBump
                | InteractionType::WeakHydrogenBond => -options.probe.bump_weight * spike_length,
                _ => options.probe.hydrogen_bond_weight * spike_length,
            } / options.probe.density;
            let item = contacts
                .entry((source, target, interaction_index(kind)))
                .or_insert(Contact {
                    source,
                    target,
                    interaction_type: kind,
                    dot_count: 0,
                    atom_gap,
                    min_dot_gap: dot.gap,
                    location,
                    score_over_density: 0.0,
                });
            item.dot_count += 1;
            item.score_over_density += score;
            if dot.gap < item.min_dot_gap {
                item.min_dot_gap = dot.gap;
                item.location = location;
            }
        }
    }
    let contacts: Vec<_> = contacts.into_values().collect();
    let mut hbonds = FxHashSet::default();
    let mut clashes = BTreeMap::<(u32, u32), Clash>::new();
    for item in &contacts {
        let pair = (item.source.min(item.target), item.source.max(item.target));
        if item.interaction_type == InteractionType::StandardHydrogenBond {
            hbonds.insert(pair);
        } else if item.interaction_type == InteractionType::Bump {
            let overlap = reported_gap(item.min_dot_gap)?;
            if overlap > -0.4 {
                continue;
            }
            let entry = clashes.entry(pair).or_insert(Clash {
                first: pair.0,
                second: pair.1,
                overlap,
                location: reported_location(item.location)?,
                max_b_factor: format!(
                    "{:.2}",
                    world.b[pair.0 as usize].max(world.b[pair.1 as usize])
                )
                .parse()
                .map_err(|_| "Probe B-factor serialization failed")?,
            });
            if overlap < entry.overlap {
                entry.overlap = overlap;
                entry.location = reported_location(item.location)?;
            }
        }
    }
    clashes.retain(|pair, _| !hbonds.contains(pair));
    let clashes: Vec<_> = clashes.into_values().collect();
    let clashscore = if n_selected > 0 {
        1000.0 * clashes.len() as f64 / n_selected as f64
    } else {
        0.0
    };
    Ok(Report {
        phantom_count: world.len() - world.n_real,
        atoms: world.labels,
        n_real: world.n_real,
        n_selected_nonwater_atoms: n_selected,
        contacts,
        clashes,
        clashscore,
    })
}
