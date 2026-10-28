//! The lane kernel (`score_dot_lanes`, as `score_atom_ext` calls it) against
//! the scalar Probe scorer of probe.rs (`trim_dots`, then
//! `DotScorer::score_dots`): the same sums to the bit on random source atoms,
//! targets and excluded atoms, with more than 64 targets in some (several
//! mask words per chunk) and under each scorer switch.

use super::*;

/// splitmix64, for a fixed sequence of cases.
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^ (z >> 31)
    }
    fn unit(&mut self) -> f64 {
        (self.next() >> 11) as f64 / (1u64 << 53) as f64
    }
    fn range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.unit()
    }
    fn below(&mut self, n: usize) -> usize {
        (self.next() % n as u64) as usize
    }
    fn chance(&mut self, p: f64) -> bool {
        self.unit() < p
    }
    fn direction(&mut self) -> Vec3 {
        loop {
            let v = v3(self.range(-1.0, 1.0), self.range(-1.0, 1.0), self.range(-1.0, 1.0));
            let l = v.length();
            if l > 0.1 && l <= 1.0 {
                return v / l;
            }
        }
    }
}

/// Probe radii of the atom types Reduce2 scores (polar H, H, N, O, C, S, ...).
const RADII: [f64; 8] = [1.0, 1.05, 1.17, 1.4, 1.5, 1.55, 1.65, 1.8];

fn atom_info(rng: &mut Rng) -> AtomInfo {
    AtomInfo {
        vdw_radius: RADII[rng.below(RADII.len())],
        is_acceptor: rng.chance(0.3),
        is_donor: rng.chance(0.3),
        is_dummy_hydrogen: rng.chance(0.05),
        is_ion: rng.chance(0.05),
        charge: [0, 0, 0, 1, -1][rng.below(5)],
        alt: [0, 0, b'A', b'B'][rng.below(4)],
    }
}

#[test]
fn lane_kernel_matches_scalar_probe_scorer() {
    let mut rng = Rng(20261028);
    // cases with a bump, an H-bond, a contact and a bad bump
    let (mut cases, mut seen) = (0, [0; 4]);
    for (weak_hbonds, ignore_ions) in [(false, false), (true, false), (false, true), (true, true)] {
        let mut p = ProbeParams::reduce2_defaults();
        p.allow_weak_hydrogen_bonds = weak_hbonds;
        p.ignore_ion_interactions = ignore_ions;
        let s = DotScorer::new(&p);
        let (pr, density) = (p.probe_radius, p.density);
        let (mut lanes, mut prepared, mut trimmed) = (LaneBufs::default(), Vec::new(), Vec::new());
        for _ in 0..400 {
            let src = atom_info(&mut rng);
            let pa = v3(rng.range(-60.0, 60.0), rng.range(-60.0, 60.0), rng.range(-60.0, 60.0));
            // the dots as `optimize_structure` builds them
            let sphere = dot_sphere(src.vdw_radius, density);
            let pairs: Vec<DotPair> = sphere
                .iter()
                .map(|&d| {
                    let l = d.length();
                    DotPair { d, probe: if l > 0.0 { d / l * (l + pr) } else { Vec3::ZERO } }
                })
                .collect();
            let soa = DotsSoa::new(&pairs);
            let n_targets = [0, 1, 4, 17, 63, 64, 65, 150][rng.below(8)];
            let targets: Vec<Target> = (0..n_targets)
                .map(|_| {
                    let info = atom_info(&mut rng);
                    let reach = src.vdw_radius + info.vdw_radius + 2.0 * pr;
                    Target { pos: pa + rng.direction() * rng.range(0.05, reach + 0.3), info }
                })
                .collect();
            let excl: Vec<(Vec3, f64)> =
                (0..rng.below(5)).map(|_| (pa + rng.direction() * rng.range(0.9, 1.7), RADII[rng.below(RADII.len())])).collect();

            trim_dots(pa, &sphere, &excl, &mut trimmed);
            let want = s.score_dots(pa, &src, &trimmed, &targets, pr, density, false);

            let mut got = ScoreDotsResult::default();
            if !(src.is_ion && s.ignore_ions) {
                prepare_targets(&s, &src, pa, &targets, pr, &mut prepared);
                score_dot_lanes(&s, &src, pa, &soa, &excl, &prepared, &mut lanes, &mut got);
                got.bump /= density;
                got.hbond /= density;
                got.attract /= density;
            }
            let bits = |r: &ScoreDotsResult| (r.bump.to_bits(), r.hbond.to_bits(), r.attract.to_bits(), r.has_bad_bump);
            assert_eq!(
                bits(&got),
                bits(&want),
                "case {cases} (weak H-bonds {weak_hbonds}, ignore ions {ignore_ions}, {n_targets} targets): lanes {:?} {:?} {:?}, scalar {:?} {:?} {:?}",
                got.bump, got.hbond, got.attract, want.bump, want.hbond, want.attract
            );
            for (n, hit) in seen.iter_mut().zip([want.bump != 0.0, want.hbond != 0.0, want.attract != 0.0, want.has_bad_bump]) {
                *n += hit as usize;
            }
            cases += 1;
        }
    }
    // the cases exercise every kind of score, not just empty sums
    assert!(seen.iter().all(|&n| n * 10 > cases), "of {cases} cases, too few score each kind: {seen:?}");
}
