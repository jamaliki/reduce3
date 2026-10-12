# Native Probe2 contact validation

`clashcheck::run` checks one selected hierarchy model. It consumes a temporary
`World` with the real atom sites, labels, occupancies, B factors, chemical radii,
donor/acceptor flags, charges and full covalent adjacency. The caller must first
resolve chemical typing and reject unsupported symmetry. The scorer never uses
an element-radius guess to claim coverage.

`Options::for_clashscore` preserves CCTBX 2026.9's actual occupancy convention:
the wrapper passes the raw maximum occupancy to a helper that divides it by
100 and caps the result at 0.1. Source atoms exclude water. Target atoms include
water; unselected atoms can still block surface dots.

The preparation step adds temporary water phantom hydrogens, fixes explicit
hydrogen donor flags and adds ionic neighbor bonds. Contact evaluation reuses
the existing Probe dot spheres and scoring rules, with bond-depth exclusions,
alternate compatibility, occupancy filters and occlusion. Each directed
contact keeps its dot count, minimum dot gap and score sum. Clash reports join
the two directions, suppress pairs with hydrogen bonds and preserve Probe2's
rounded gap, location and B-factor values. The clashscore denominator includes
all selected non-water atoms, even when they produce no dots.

The spatial index uses the existing dense grid for compact bounds and sparse
buckets for large disconnected bounds. Nonfinite inputs, missing radii,
multiple hierarchy models and unsupported numerical ranges return errors.
These errors require the reference path or an incomplete report. They do not
show that a structure passed the checks.

## Reference tests

The fixture tool uses CCTBX's captured initial chemical typing and covalent
bonds, then independently runs the native preparation and dot calculation:

```sh
cargo build --no-default-features --features refcheck --example probe_fixture
python tools/check_probe_reference.py \
  --reference-dir /absolute/reference/capture/directory \
  --native-bin /absolute/reduce3/target/debug/examples/probe_fixture \
  --pattern '*-probe.json' --pattern '*-cif-no-cdl.json' \
  --pattern '*-nuclear-no-cdl.json' --pattern '*post-probe.json'
```

The comparison checks atom-pair IDs, gaps, locations, B factors, clashscores,
directed contact counts and minimum dot gaps against standalone Probe2 output.
This isolates the scorer. A separate end-to-end test must check the native
chemical graph and initial atom information before claiming full replacement.
Use one model per audit view; do not merge model-specific denominators.

`tools/check_probe_charge_oracle.py` compares direct CCTBX dot scoring on
synthetic neutral, opposite-charge and same-charge pairs. It checks the regular
and charged overlap boundaries and positive-gap contacts. The positive donor
hydrogen is deliberate test input: CCTBX's full model interpretation clears
that formal charge before Probe, so this oracle does not prove full-model
chemical typing.

Do not infer a speedup from fixture runtime. Build first, then run paired
benchmarks serially, with warmup, alternating order, repeated measurements and
laptop load records. The current archive scan must remain unchanged until a
separate promotion decision.

The orchestration adapts CCTBX's Apache-2.0 Probe2 implementation, Copyright
2021-2023 Richardson Lab at Duke. The repository's `LICENSE`,
`LICENSE-cctbx.txt` and `NOTICE` retain the applicable license notices.
