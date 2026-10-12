//! Probe2 orchestration tests on explicit, typed atoms.

use reduce3::clashcheck::{Options, run};
use reduce3::probe::AtomInfo;
use reduce3::world::World;

fn carbon_pair() -> World {
    let structure = reduce3::pdbio::read_pdb(concat!(
        "ATOM      1  C   GLY A   1       0.000   0.000   0.000  1.00  0.00           C\n",
        "ATOM      2  C   GLY A   2       2.000   0.000   0.000  1.00  0.00           C\n"
    ));
    let mut world = World::from_structure(&structure);
    assert_eq!(world.len(), 2);
    world.info = vec![
        AtomInfo {
            vdw_radius: 1.7,
            ..AtomInfo::default()
        };
        2
    ];
    world.bonded = vec![vec![]; 2];
    world
}

#[test]
fn directed_dots_reduce_to_one_physical_clash_pair() {
    let report = run(carbon_pair(), &Options::default()).unwrap();
    assert_eq!(report.clashes.len(), 1);
    assert_eq!(report.n_selected_nonwater_atoms, 2);
    assert_eq!(report.clashscore, 500.0);
    assert!(
        report
            .contacts
            .iter()
            .any(|contact| contact.source == 0 && contact.target == 1)
    );
    assert!(
        report
            .contacts
            .iter()
            .any(|contact| contact.source == 1 && contact.target == 0)
    );
    assert!(report.clashes[0].overlap < -1.0);
}

#[test]
fn bonded_atoms_and_incompatible_alternates_do_not_clash() {
    let mut world = carbon_pair();
    world.bonded = vec![vec![1], vec![0]];
    assert!(run(world, &Options::default()).unwrap().clashes.is_empty());
    let mut world = carbon_pair();
    world.labels[0].altloc = "A".into();
    world.labels[1].altloc = "B".into();
    world.info[0].alt = b'A';
    world.info[1].alt = b'B';
    assert!(run(world, &Options::default()).unwrap().clashes.is_empty());
}

#[test]
fn blank_alt_occupancy_rule_excludes_underfilled_pairs() {
    let mut world = carbon_pair();
    world.occ = vec![0.5, 0.5];
    assert!(run(world, &Options::default()).unwrap().clashes.is_empty());
}

#[test]
fn missing_typing_fails_closed() {
    let mut world = carbon_pair();
    world.info[1].vdw_radius = 0.0;
    assert!(run(world, &Options::default()).is_err());
}

#[test]
fn extreme_coordinates_use_sparse_cells_or_fail_closed() {
    let mut world = carbon_pair();
    world.pos[1] = reduce3::geom::v3(1e12, 1e12, 1e12);
    assert!(run(world, &Options::default()).unwrap().clashes.is_empty());
    let mut world = carbon_pair();
    world.pos[1].x = 1e200;
    assert!(run(world, &Options::default()).is_err());
}

#[test]
fn nonfinite_occupancy_bfactor_and_radius_fail_closed() {
    let mut world = carbon_pair();
    world.occ[0] = f64::NAN;
    assert!(run(world, &Options::default()).is_err());
    let mut world = carbon_pair();
    world.b[0] = f64::INFINITY;
    assert!(run(world, &Options::default()).is_err());
    let mut world = carbon_pair();
    world.info[0].vdw_radius = f64::NAN;
    assert!(run(world, &Options::default()).is_err());
}

#[test]
fn different_models_require_separate_scores() {
    let mut world = carbon_pair();
    world.labels[1].model_index = 1;
    assert!(run(world, &Options::default()).is_err());
}

#[test]
fn water_phantoms_are_temporary_and_do_not_hide_collisions() {
    let mut world = carbon_pair();
    world.labels[1].is_water = true;
    world.labels[1].element = "O".into();
    world.labels[1].resname = "HOH".into();
    world.info[1].vdw_radius = 1.52;
    world.info[0].is_acceptor = true;
    let report = run(world, &Options::default()).unwrap();
    assert_eq!(report.phantom_count, 1);
    assert_eq!(report.n_real, 2);
    assert_eq!(report.n_selected_nonwater_atoms, 1);
    assert_eq!(report.clashes.len(), 1);
    assert_eq!(report.clashscore, 1000.0);
}

#[test]
fn charged_hydrogen_bond_rules_preserve_excessive_overlap() {
    let mut world = carbon_pair();
    world.labels[0].element = "H".into();
    world.info[0].is_donor = true;
    world.info[0].charge = 1;
    world.labels[1].element = "C".into();
    world.info[1].is_acceptor = true;
    world.info[1].charge = -1;
    let report = run(world, &Options::default()).unwrap();
    assert!(
        report
            .contacts
            .iter()
            .any(|contact| contact.interaction_type
                == reduce3::probe::InteractionType::StandardHydrogenBond)
    );
    assert!(report.clashes.is_empty());
}

#[test]
fn fourth_bond_exclusion_requires_a_hydrogen_endpoint() {
    fn chain(hydrogen_end: bool) -> World {
        let mut world = carbon_pair();
        for i in 2..5 {
            let mut label = world.labels[1].clone();
            label.ag = i as u32;
            world.labels.push(label);
            world.pos.push(reduce3::geom::v3(i as f64 * 0.5, 0.0, 0.0));
            world.info.push(world.info[1]);
            world.occ.push(1.0);
            world.b.push(0.0);
            world.bonded.push(vec![]);
            world.riding.push(None);
        }
        world.pos[1].x = 0.5;
        world.n_real = 5;
        for i in 0..4 {
            world.bonded[i].push((i + 1) as u32);
            world.bonded[i + 1].push(i as u32);
        }
        if hydrogen_end {
            world.labels[4].element = "H".into();
        }
        world
    }
    let ordinary = chain(false);
    let excluded = reduce3::optimizer::atoms_within_n_bonds(&ordinary, 0, 0.25, 4, 3);
    assert!(!excluded.contains(&4));
    let hydrogen = chain(true);
    let excluded = reduce3::optimizer::atoms_within_n_bonds(&hydrogen, 0, 0.25, 4, 3);
    assert!(excluded.contains(&4));
}

#[test]
fn wrapper_occupancy_convention_uses_raw_maximum_divided_by_100() {
    let mut world = carbon_pair();
    world.occ[1] = 0.05;
    world.labels[1].altloc = "A".into();
    world.info[1].alt = b'A';
    let options = Options::for_clashscore(&world, false);
    assert_eq!(options.occupancy_fraction, 0.01);
    let report = run(world, &options).unwrap();
    assert_eq!(report.n_selected_nonwater_atoms, 2);
    assert_eq!(report.clashes.len(), 1);
}
