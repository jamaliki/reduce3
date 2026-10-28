# Golden outputs

`tests/golden.rs` runs fixed mode (the default) on the 301 structures here and
compares each result with the golden file for the platform it runs on.

## What is recorded

`fixed-<arch>-<os>.tsv` has one row per structure:

| Column | Contents |
|---|---|
| `status` | `ok`, or `error` if the run failed (then `hydrogen_sha256` hashes the message) |
| `hydrogens` | number of hydrogens in the output |
| `hydrogen_sha256` | sha256 of every hydrogen's label and coordinates as written (`{:.3}`, 1e-3 Å), in output order |
| `movers` | number of Movers in the report |
| `decision_sha256` | sha256 of the report's decisions: each Mover (type and place) and its final pose, without the scores |
| `heavy_sha256` | sha256 of every other atom's label and coordinates as written (flips move heavy atoms) |

Any change in a value fails the test. `determinism_across_thread_counts`
checks that one rayon thread, many threads and batch mode (several models at
once, each sequential) give the same output to the bit: every coordinate,
occupancy and B, and the whole description apart from its timings.

## Platforms

Fixed mode is not bit-identical across platforms. `riding::cx`,
`ccdrestraints` and `cell` use fused multiply-adds on aarch64 only, and
trigonometry comes from the platform's libm. So each platform has its own golden file:

- `fixed-aarch64-macos.tsv`: Apple silicon.
- `fixed-x86_64-linux.tsv`: Sandpit Tokyo compute nodes (Ubuntu 24.04, glibc 2.39), where the corpora are built.

A platform without a golden file fails the test. Bless one there to add it.

## chem_data

The outputs depend on the monomer libraries. `chem_data.tsv` pins chem_data
2026.9 (the conda package `chem_data-2026.9-pyhf50f6f0_0`, sha256 in the
file). It lists the 559 files that this test and the rest of the suite read,
each with its sha256. The list was found by logging every file `MonLib`
reads. The golden test checks those files before it runs. To fetch them
without the rest of chem_data (1.3 GB download, 5 MB kept):

```bash
python3 tools/golden/fetch_chem_data.py tests/golden/chem_data.tsv chem_data
REDUCE3_CHEM_DATA=chem_data REDUCE3_REQUIRE_CHEM_DATA=1 cargo test
```

Without chem_data, the tests that need it skip. Set
`REDUCE3_REQUIRE_CHEM_DATA=1`, as CI does, to make them fail instead. Some
tests already fail without it (`expect`).

## Running and blessing

```bash
cargo test --release --test golden                       # both tests
REDUCE3_GOLDEN_DUMP=/tmp/dump cargo test --release --test golden fixed_mode
REDUCE3_GOLDEN_BLESS=1 cargo test --release --test golden fixed_mode
```

- The golden test runs in debug builds too (about 15–30 s). The determinism test is ignored in debug builds.
- `REDUCE3_GOLDEN_DUMP` writes each structure's hydrogens, heavy atoms (as written, then exact) and report, so two dumps can be diffed.
- `REDUCE3_GOLDEN_BLESS` rewrites this platform's file. Bless only after a deliberate change, and review the diff. A change that moves output on one platform usually moves it on the others, so bless every platform's file in the same commit.

## The structures

`structures.tsv` lists each structure's class, where it came from and the
sha256 of its original file and of its input text. They are a seeded random
sample (`tools/golden/select_structures.py`, seed 20261028) of the local
mirrors of the PDB and AlphaFold DB v6:

| Class | Count |
|---|---|
| AlphaFold DB models | 60 |
| proteins with ligands | 55 |
| alternate locations | 35 |
| metals | 35 |
| proteins only | 30 |
| nucleic acids | 30 |
| covalent links (`_struct_conn` covale) | 20 |
| modified residues | 15 |
| metal clusters (Fe-S, haem, ...) | 10 |
| NMR ensembles | 10 |
| a large assembly (9epc, 65,302 heavy atoms) | 1 |

Each entry has at most 4,000 heavy atoms and 8,000 atom rows, except the large
assembly. An entry can fit several classes but is counted once.

`tools/golden/prepare_inputs.py` trims each file to the categories and items
Reduce3 reads (`_atom_site`, `_atom_site_anisotrop`, `_cell`, `_symmetry`,
`_space_group`, `_struct_conn`), keeping every value as written. This needs
checking. Before the inputs were committed, golden runs on the untrimmed
originals gave identical dumps (coordinates, reports and scores) for all 301
structures. A first trim that dropped `_atom_site_anisotrop` changed one
rotor. Anisotropic ADPs set their atoms' B, and a water's B decides whether
it gets phantom hydrogens.

PDB entries are CC0. AlphaFold DB models are CC BY 4.0 (AlphaFold DB,
EMBL-EBI and Google DeepMind).
