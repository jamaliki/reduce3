//! Fixed-mode clearance guard for rotor hydrogens (fix 29 in the README).
//!
//! A hydroxyl, thiol or ammonium hydrogen that Reduce2 aims straight at a
//! nearby acceptor can end up 1.3 A from it. The score does not stop it: the
//! hydrogen's own sphere trims the dots of its parent atom that overlap the
//! acceptor, which hides the clash between the two heavy atoms, and Probe
//! softens a hydrogen bond that is too short. So no pose that puts a rotor
//! hydrogen closer than `min_dist` to a heavy atom is offered, wherever a
//! pose that keeps clear exists.
//!
//! The guard is a property of the poses, not of the score, so the ordinary
//! scoring and everything that is not a rotor stay as they are.

use super::RotatorMeta;
use crate::geom::{Vec3, rotate_deg_axis_dir};
use crate::probe::SpatialGrid;
use crate::world::World;
use smallvec::SmallVec;

/// The heavy atoms a rotor hydrogen must keep clear of.
pub struct Clearance {
    grid: SpatialGrid,
    min_dist: f64,
}

impl Clearance {
    /// Index the heavy atoms among `atoms` at their current positions. This
    /// includes waters that the scoring ignores for their occupancy or B
    /// factor: a hydrogen has no more room there than anywhere else. `None`
    /// for a limit that is not positive.
    pub fn new(w: &World, atoms: &[u32], min_dist: f64) -> Option<Clearance> {
        if min_dist.is_nan() || min_dist <= 0.0 {
            return None;
        }
        let points: Vec<(u32, Vec3)> = atoms
            .iter()
            .copied()
            .filter(|&a| (a as usize) < w.n_real && !w.is_h(a))
            .map(|a| (a, w.pos[a as usize]))
            .collect();
        Some(Clearance { grid: SpatialGrid::new(&points, 3.0), min_dist })
    }

    /// Whether a heavy atom outside `skip` lies within the limit of `at`.
    fn blocks(&self, skip: &[u32], at: Vec3) -> bool {
        let limit = self.min_dist * self.min_dist;
        let mut blocked = false;
        self.grid.for_each_within(at, 0.0, self.min_dist, |id, _, d2| {
            blocked |= d2 < limit && !skip.contains(&id);
        });
        blocked
    }
}

/// The hydrogen's bonded neighbors and their neighbors: its parent and the
/// atom behind it. Their distance to the hydrogen is fixed by the geometry.
fn bonded_within_two(w: &World, h: u32) -> SmallVec<[u32; 8]> {
    let mut near = SmallVec::new();
    for &n in &w.bonded[h as usize] {
        near.push(n);
        near.extend(w.bonded[n as usize].iter().copied());
    }
    near
}

/// Remove the coarse angles of `rot` that leave a hydrogen of `hydrogens` (at
/// its start position in `w`) too close to a heavy atom, or that have a fine
/// angle that does, and add the nearest angle that is clear on either side of
/// each one removed, so a hydrogen can still turn up to the limit.
///
/// `period` is the angle after which the hydrogens fill each other's places
/// (360 for one hydrogen, 120 for an ammonium); `admissible(angle, fine)`
/// holds the limits that apply to a pose besides this one (the planar
/// hydroxyl limit).
///
/// Nothing changes when no coarse angle is blocked, nor when no turn at all
/// is clear (checked at each degree, first with the fine angles and then
/// without them, which are dropped when that is what leaves a turn): such a
/// rotor is scored as Reduce2 scores it.
pub(super) fn constrain(
    w: &World,
    clear: &Clearance,
    hydrogens: &[u32],
    rot: &mut RotatorMeta,
    period: f64,
    admissible: &dyn Fn(f64, &[f64]) -> bool,
) {
    let skip: Vec<SmallVec<[u32; 8]>> = hydrogens.iter().map(|&h| bonded_within_two(w, h)).collect();
    let (origin, dir) = (rot.axis_origin, rot.axis_dir);
    let clear_at = |angle: f64| {
        hydrogens.iter().zip(&skip).all(|(&h, skip)| !clear.blocks(skip, rotate_deg_axis_dir(origin, dir, w.pos[h as usize], angle)))
    };
    let legal = |angle: f64, fine: &[f64]| {
        admissible(angle, fine) && std::iter::once(0.0).chain(fine.iter().copied()).all(|offset| clear_at(angle + offset))
    };
    let wrap = |angle: f64| (angle + period / 2.0).rem_euclid(period) - period / 2.0;

    let coarse = rot.coarse_angles.clone();
    let fine = rot.fine_angles.clone();
    if coarse.iter().all(|&a| legal(a, &fine)) {
        return;
    }
    for fine_set in [&fine[..], &[]] {
        let mut angles: Vec<f64> = coarse.iter().copied().filter(|&a| legal(a, fine_set)).collect();
        for &blocked in coarse.iter().filter(|&&a| !legal(a, fine_set)) {
            for direction in [-1.0, 1.0] {
                let nearest = (1..=(period / 2.0) as i32)
                    .map(|step| wrap(blocked + direction * step as f64))
                    .find(|&candidate| legal(candidate, fine_set));
                if let Some(candidate) = nearest {
                    if !angles.iter().any(|&a| wrap(a - candidate).abs() < 1e-6) {
                        angles.push(candidate);
                    }
                }
            }
        }
            if !angles.is_empty() {
            if fine_set.is_empty() {
                rot.fine_angles.clear();
            }
            rot.coarse_angles = angles;
            return;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::v3;
    use crate::movers::{RidingRef, SingleHOptions, single_hydrogen_rotator};
    use crate::probe::AtomInfo;
    use crate::world::AtomLabels;

    /// A methanol-like hydroxyl: H (0), O (1), C (2) with three methyl
    /// neighbors (3-5), then one acceptor oxygen for each torsion in
    /// `acceptors`, 2.0 A from O at that torsion from the hydrogen.
    fn world(acceptors: &[f64]) -> World {
        let down = 70.5_f64.to_radians();
        let tetra = |phi: f64, len: f64| v3(len * down.cos(), len * down.sin() * phi.to_radians().cos(), len * down.sin() * phi.to_radians().sin());
        let h = tetra(0.0, 0.96);
        let mut pos = vec![v3(-h.x, h.y, h.z), v3(0.0, 0.0, 0.0), v3(1.43, 0.0, 0.0)];
        pos.extend([0.0, 120.0, 240.0].map(|phi| v3(1.43, 0.0, 0.0) + tetra(phi, 1.52)));
        let hydrogen = pos[0];
        pos.extend(acceptors.iter().map(|&t| rotate_deg_axis_dir(Vec3::ZERO, v3(1.0, 0.0, 0.0), hydrogen.normalize() * 2.0, t)));
        let n = pos.len();
        let label = |i: usize| {
            let (name, element) = match i {
                0 => ("HG", "H"),
                1 => ("OG", "O"),
                2..=5 => ("C", "C"),
                _ => ("OW", "O"),
            };
            AtomLabels { name: name.into(), element: element.into(), resname: "SER".into(), chain: "A".into(), ..Default::default() }
        };
        let mut bonded = vec![vec![]; n];
        for (a, b) in [(0, 1), (1, 2), (2, 3), (2, 4), (2, 5)] {
            bonded[a].push(b as u32);
            bonded[b].push(a as u32);
        }
        let mut riding = vec![None; n];
        riding[0] = Some(RidingRef { n: 0, a2: 3 });
        World {
            labels: (0..n).map(label).collect(),
            pos,
            occ: vec![1.0; n],
            b: vec![10.0; n],
            info: vec![AtomInfo::default(); n],
            bonded,
            riding,
            n_real: n,
        }
    }

    fn rotor(w: &mut World, clearance: Option<&Clearance>) -> crate::movers::Mover {
        let acceptors: Vec<u32> = (6..w.len() as u32).collect();
        let options = SingleHOptions {
            circular_angle_spacing: true,
            any_partner_valence: true,
            constrain_planar_hydroxyls: true,
            planar_preference: 0.0,
            syn_preference: 0.0,
            clearance,
        };
        single_hydrogen_rotator(w, 0, &acceptors, &[], &options).unwrap()
    }

    /// The hydrogen's distance to the nearest acceptor in each pose.
    fn gaps(w: &World, m: &crate::movers::Mover) -> Vec<f64> {
        m.all_poses().map(|pose| (6..w.len()).map(|a| pose[0].dist(w.pos[a])).fold(f64::INFINITY, f64::min)).collect()
    }

    fn clearance(w: &World) -> Clearance {
        Clearance::new(w, &(0..w.len() as u32).collect::<Vec<_>>(), 1.5).unwrap()
    }

    #[test]
    fn without_the_guard_a_pose_aims_at_the_acceptor() {
        let mut w = world(&[0.0]);
        let m = rotor(&mut w, None);
        assert!(gaps(&w, &m).into_iter().any(|d| d < 1.3));
    }

    #[test]
    fn no_pose_is_blocked_and_the_poses_run_up_to_the_limit() {
        let mut w = world(&[0.0]);
        let clear = clearance(&w);
        let m = rotor(&mut w, Some(&clear));
        let gaps = gaps(&w, &m);
        assert!(gaps.iter().all(|&d| d >= 1.5 - 1e-9), "{gaps:?}");
        // the poses nearest the acceptor are the clear ones next to its blocked arc
        assert!(gaps.iter().any(|&d| d < 1.6), "{gaps:?}");
        assert!(m.fine.count() > 0);
        // the hydroxyl's own atoms are where they were
        assert!(m.coarse_pos.iter().all(|pose| pose[1].dist(w.pos[1]) < 1e-9));
    }

    #[test]
    fn a_rotor_with_nothing_in_reach_is_not_touched() {
        let far = || {
            let mut w = world(&[0.0]);
            w.pos[6] = v3(-1.0, 3.0, 0.0);
            w
        };
        let clear = clearance(&far());
        let plain = rotor(&mut far(), None);
        let guarded = rotor(&mut far(), Some(&clear));
        assert_eq!(plain.coarse_pos, guarded.coarse_pos);
        assert_eq!(plain.fine.iter().collect::<Vec<_>>(), guarded.fine.iter().collect::<Vec<_>>());
    }

    #[test]
    fn with_no_clear_turn_the_rotor_is_left_as_it_was() {
        let all_round = [0.0, 60.0, 120.0, 180.0, 240.0, 300.0];
        let mut w = world(&all_round);
        let clear = clearance(&w);
        let plain = rotor(&mut world(&all_round), None);
        let guarded = rotor(&mut w, Some(&clear));
        assert!(gaps(&w, &guarded).into_iter().all(|d| d < 1.5));
        assert_eq!(plain.coarse_pos, guarded.coarse_pos);
        assert_eq!(plain.fine.iter().collect::<Vec<_>>(), guarded.fine.iter().collect::<Vec<_>>());
    }

    #[test]
    fn a_clear_arc_narrower_than_the_fine_angles_is_searched_without_them() {
        // four acceptors 100 degrees apart leave a gap of a few degrees
        // between the arcs they block, in three places
        let spaced = [0.0, 100.0, 200.0, 300.0];
        let mut w = world(&spaced);
        let clear = clearance(&w);
        let plain = rotor(&mut world(&spaced), None);
        assert!(plain.fine.count() > 0);
        let guarded = rotor(&mut w, Some(&clear));
        assert_eq!(guarded.fine.count(), 0);
        assert!(guarded.num_coarse() >= 2);
        assert!(gaps(&w, &guarded).into_iter().all(|d| d >= 1.5 - 1e-9));
    }

    #[test]
    fn the_parent_and_the_atom_behind_it_are_not_obstacles() {
        let w = world(&[]);
        let clear = clearance(&w);
        // 0.96 A from O and about 2.0 A from C, every pose is closer than the limit to them
        assert!(clear.blocks(&[], w.pos[0]));
        assert!(!clear.blocks(&bonded_within_two(&w, 0), w.pos[0]));
    }
}
