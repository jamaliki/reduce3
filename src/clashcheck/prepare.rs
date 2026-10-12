use super::Options;
use crate::atomtypes::is_aromatic_acceptor;
use crate::probe::AtomInfo;
use crate::world::World;

pub(super) fn prepare(world: &mut World, options: &Options) -> Result<Vec<bool>, String> {
    let query = super::neighbors::Query::new(&world.pos)?;
    let mut ignored = vec![false; world.n_real];
    let mut waters: Vec<_> = (0..world.n_real as u32)
        .filter(|&i| {
            world.labels[i as usize].is_water && world.occ[i as usize] > options.occupancy_fraction
        })
        .collect();
    waters.sort_by_key(|&i| world.atom_id_string(i));
    for i in waters {
        let index = i as usize;
        if world.is_h(i) {
            ignored[index] = world.occ[index] < options.minimum_water_hydrogen_occupancy
                || world.b[index] > options.maximum_water_hydrogen_b;
        } else if world.elem(i) == "O" {
            world.info[index].is_donor = false;
            world.info[index].is_acceptor = true;
            if !world.bonded[index].is_empty() {
                continue;
            }
            let (radius, distance) = if options.nuclear {
                (1.0, 0.98)
            } else {
                (1.05, 0.84)
            };
            let mut nearby = Vec::new();
            query.for_each_within(world.pos[index], 0.001, 4.0, |j, _, _| nearby.push(j));
            nearby.sort_unstable();
            let mut candidates: Vec<(u32, f64)> = Vec::new();
            for j in nearby {
                let target = j as usize;
                if !world.compatible(i, j)
                    || world.labels[index].model_index != world.labels[target].model_index
                    || world.elem(j) == "H"
                    || !world.info[target].is_acceptor
                {
                    continue;
                }
                let overlap = world.pos[index].dist(world.pos[target])
                    - (radius + world.info[target].vdw_radius + distance);
                if overlap > -0.1 || world.occ[target] < 0.0 {
                    continue;
                }
                let mut found = false;
                if is_aromatic_acceptor(world.resname(j), world.name(j)) {
                    for item in &mut candidates {
                        if world.labels[item.0 as usize].ag == world.labels[target].ag
                            && is_aromatic_acceptor(world.resname(item.0), world.name(item.0))
                        {
                            if overlap < item.1 {
                                *item = (j, overlap);
                            }
                            found = true;
                            break;
                        }
                    }
                }
                if !found {
                    candidates.push((j, overlap));
                }
            }
            for (j, overlap) in candidates {
                let offset = world.pos[j as usize] - world.pos[index];
                let length = offset.length();
                if length == 0.0 {
                    continue;
                }
                let shift = distance + (-distance).max((overlap + 0.6).min(0.0));
                let info = AtomInfo {
                    vdw_radius: radius,
                    is_donor: true,
                    is_dummy_hydrogen: true,
                    ..AtomInfo::default()
                };
                world.add_phantom(i, world.pos[index] + offset * (shift / length), info);
                ignored.push(false);
            }
        }
    }
    for i in 0..world.n_real {
        if world.is_h(i as u32) {
            for &j in &world.bonded[i] {
                if matches!(world.elem(j), "N" | "O" | "S") {
                    world.info[i].is_donor = true;
                    world.info[j as usize].is_donor = false;
                }
            }
        } else if matches!(world.elem(i as u32), "N" | "O" | "S") {
            world.info[i].is_donor = false;
        }
    }
    for i in 0..world.n_real as u32 {
        if !world.is_positive_ion(i) {
            continue;
        }
        let radius = world.info[i as usize].vdw_radius;
        query.for_each_within(
            world.pos[i as usize],
            radius,
            3.25 + radius,
            |j, _, distance_sq| {
                let distance = distance_sq.sqrt();
                let expected = radius + world.info[j as usize].vdw_radius;
                if distance >= expected - 0.6 && distance <= expected + 0.2 {
                    world.bonded[i as usize].push(j);
                    world.bonded[j as usize].push(i);
                }
            },
        );
    }
    Ok(ignored)
}
