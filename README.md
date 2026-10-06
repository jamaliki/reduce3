# Reduce3

Reduce3 is a Rust reimplementation of cctbx **Reduce2** (`mmtbx.reduce2`). It adds hydrogens to a
macromolecular model and optimizes the rotatable and flippable groups (OH/SH/NH3+/methyl rotations,
Asn/Gln/His flips) by Probe dot scoring. It covers the whole pipeline:

* restraint interpretation (monomer library, links, modifications, disulfides, zinc coordination,
  automatic linking such as N-glycosylation and glycosidic bonds),
* riding-hydrogen placement,
* Probe contact scoring, the Movers and the clique optimizer,
* PDB and mmCIF input and output, with the same formatting as iotbx.

It has two modes:

* **fixed** (default) corrects the Reduce2 bugs listed below and optimizes every clique exactly, in
  parallel.
* **compat** (`--compat`) reproduces Reduce2 exactly, bugs included. It exists to validate the port.

## Build

```bash
cargo build --release
```

The binary is `target/release/reduce3`. It needs the cctbx `chem_data` directory (the `geostd` and
`mon_lib` monomer libraries plus `chemical_components`). Reduce3 looks for it in `--chem-data DIR`,
then `$REDUCE3_CHEM_DATA`, `$CHEM_DATA`, and finally the active conda environment.

## Usage

```bash
reduce3 model.pdb                              # writes modelH.pdb and modelH.txt
reduce3 model.cif add_flip_movers=True         # writes modelFH.cif and modelFH.txt
reduce3 --compat model.pdb -o out.pdb          # behave exactly like Reduce2
reduce3 modelH.pdb approach=optimize           # optimize existing hydrogens
reduce3 modelH.pdb approach=remove             # strip hydrogens
```

Reduce2's `name=value` parameters work the same way, with the same defaults: `approach`,
`add_flip_movers`, `n_terminal_charge`, `keep_existing_H`, `exclude_water`,
`use_neutron_distances`, `preference_magnitude`, `non_flip_preference`, `skip_bond_fix_up`,
`set_flip_states`, `model_id`, `alt_id`, `bonded_neighbor_depth`, `verbosity`,
`stop_on_any_missing_hydrogen`, `ignore_missing_restraints`, `output.filename`,
`output.description_file_name`, `output.write_files`, and all `probe.*` scoring parameters.
Extra options: `--threads N`, `-q`. `reduce3 --help` lists everything.

## Library use

Reduce3 is also a Rust library. It can work on a model that another program has already
parsed, without writing or re-reading text:

* implement `cifsource::CifSource` (categories as tables of values) for the program's parsed
  mmCIF data block, and call `reduce3::run_cif(&block, &monlib, &params)`;
* stream the result into the program's own document type by implementing `cifsource::CifSink`
  and calling `mmcif::write_cif(&output.structure, &mut sink)`.

The mmCIF reader and writer use the same two traits, so this path builds exactly the model a
file would and emits exactly the items and loops `reduce3` writes. A source that already holds
parsed numbers can hand them over through `CifTable::number` instead of text.

## Validation against Reduce2 (cctbx 2026.9)

The reference data come from the original Reduce2 run on 13 structures: 1a28, 1crn, 1d3z (NMR,
10 models), 1dfu, 1ehz (tRNA), 1ubq, 1xso, 2oob, 3gfh, 3vyk (glycans), 4fen (RNA with cobalt
hexammine), 6oge (glycoprotein, 22k atoms with H) and 7c31 (anisotropic B). Each was run with and
without flips.

| Check | Result in compat mode |
|---|---|
| Output model (PDB), 26 runs | **byte-identical** in all 26 |
| Output model (mmCIF): 4fen (mmCIF input), 7c31, 1d3z, 1ehz | **byte-identical** in all 4 |
| Description file (scores, Mover report, deletions) | identical apart from timing lines, the header, and the order of the deletion list (Reduce2 prints it in Python set order) |
| `approach=optimize` and `approach=remove` (1crn) | byte-identical |
| Hydrogen placement vs. dumps, 13 structures | same atoms, names, riding types; coordinates equal to 0.0000 Å, about 98% of them bit-identical |
| Optimizer on Reduce2's own intermediate state (26 dumps) | same report, coordinates and deletions in all 26 |

Bit-exact agreement needed two things beyond porting the code:

* cctbx moves every site through fractional coordinates and back (twice) while processing the
  model, and anisotropic ADPs go through U\*. Reduce3 repeats that round-off using the
  symmetry-averaged unit cell cctbx derives.
* The arm64 cctbx build fuses multiply-adds in its C++ vector code. Compat mode uses the same fused
  forms on arm64 so that exact score ties break the same way.

`tools/compare_outputs.sh` reruns the end-to-end comparison. The developer subcommands
`refcheck`, `hcheck` and `wcheck` (built with `--features refcheck`) compare against JSON dumps
written by `reference/harness/dump_ref.py`. [reference/README.md](reference/README.md) explains
how to recreate the reference data, and [docs/INTERPRETATION_SPEC.md](docs/INTERPRETATION_SPEC.md)
documents the cctbx interpretation behavior that Reduce3 reproduces.

## Performance

These are wall-clock times on an Apple-silicon Mac, including file I/O. Reduce3 is the best of 3
runs; Reduce2 is a single run, and its time includes about 1.5 s of Python/cctbx start-up.

| Structure | Atoms in | Reduce2 | Reduce3 fixed | Reduce3 compat | Speed-up (fixed) |
|---|---:|---:|---:|---:|---:|
| 1crn | 327 | 1.73 s | 0.013 s | 0.014 s | 133× |
| 1ubq | 660 | 1.93 s | 0.017 s | 0.020 s | 112× |
| 7c31 | 1,532 | 2.53 s | 0.022 s | 0.026 s | 117× |
| 1ehz | 1,821 | 2.92 s | 0.028 s | 0.035 s | 105× |
| 4fen | 2,084 | 17.27 s | 0.036 s | 1.89 s | 484× |
| 1xso | 2,541 | 5.86 s | 0.041 s | 0.050 s | 145× |
| 3gfh | 3,291 | 4.76 s | 0.042 s | 0.051 s | 114× |
| 1a28 | 4,262 | 6.22 s | 0.054 s | 0.068 s | 116× |
| 6oge | 11,494 | 12.73 s | 0.120 s | 0.172 s | 106× |
| 1d3z (10 models) | 12,310 | 12.48 s | 0.078 s | 0.067 s | 161× |
| 3j3q (HIV capsid, mmCIF) | 2,440,800 | not run | 24 s (about 10 GB peak memory) | | |

Where the speed comes from:

* Fixed mode solves each clique of interacting Movers exactly by variable elimination, and the
  cliques run in parallel. The score tables are decomposed per Probe dot: a dot that no other
  Mover can reach is scored once per state, and the remaining dots are grouped by the Movers that
  reach them. For 4fen, whose cobalt hexammine forms a 7-Mover clique, this took the optimizer from
  1.6 s to 0.03 s.
* Compat mode runs a direct port of Reduce2's `OptimizerC` (vertex cuts, brute force, caches), so
  its cost follows Reduce2's search. That is why 4fen takes 1.9 s there.
* Interpretation, riding placement and scoring are flat-array code with spatial grids, so every
  stage is linear in the number of atoms.

## Reduce2 bugs fixed in the default mode

`--compat` keeps all of these, for validation.

**Program flow**

1. Only the last model of a multi-model file was optimized: the Optimizer's `modelIndex=0` default
   becomes index −1. All models are now optimized (in parallel).
2. `model_id` selected the wrong model: the 1-based id was used as a 0-based index, and an
   out-of-range id silently kept every model and used the first. It is now 1-based, and an
   out-of-range id is an error.
3. `use_neutron_distances` never reached the Optimizer, so phantom water hydrogens and the report
   always used X-ray values. The setting now applies there too.
4. Hydrogens deleted for one alternate conformation were forgotten when the next alternate ran.
   Each alternate now keeps its own deletions.

**Hydrogen placement**

5. Amide and guanidinium NH2 names (HD21/HD22, HE21/HE22, HH11/HH12) came from the periodic image
   of the dictionary torsion nearest a temporary hydrogen position. In the reference outputs this
   swapped the names for 15 of 197 groups in 6oge, 2/32 in 1xso, 1/40 in 3gfh and 1/16 in 7c31.
   Fixed mode uses the dictionary torsion, and every group now follows the CCD/dictionary
   convention.
6. `flip_symmetric_amino_acids` flipped every alternate conformation of a residue as soon as one
   of them needed it. Each atom group is now tested on its own geometry.
7. Modified amino acids had their hydrogens renamed to standard amino-acid names. That made them
   "unexpected", and they were silently dropped. Only standard and D amino acids are renamed now.
8. In automatic linking, the per-residue-pair link bookkeeping was keyed by a slice of
   `atom.id_str()`. For multi-model files that slice contains part of the model id and the atom
   name, so the limits applied per atom pair. It is now keyed by model and atom groups.

**Optimizer and Movers**

9. The spatial query kept atoms filed in the grid cell of their position when last inserted.
   Atoms moved by Movers or by methyl staggering could then be missed by neighbor searches. Queries
   now use current positions.
10. Waters later rejected for low occupancy or high B were still used to aim phantom hydrogens.
    They are now removed first.
11. Phantom hydrogens were added to waters with bonded neighbors, although the code's own comment
    says to leave those alone. Such waters are now skipped.
12. Phantom hydrogens were not counted as potential contacts when Movers were placed. They now are.
13. Atoms marked for deletion during placement still acted as scoring targets and still trimmed
    dots. They no longer do either.
14. OptimizerC's clique decomposition (vertex cuts with greedy sub-searches) is not guaranteed to
    find the best combination. Variable elimination is exact.
15. The fine-optimization step compared against a stale coarse score. The report re-parsed the
    printed (rounded) initial score, and it did not give each Mover's score in the final
    configuration. All three are fixed.
16. Flip annotations scored deleted atoms, and were computed for His flips whose flipping was
    disabled. Both are fixed.
17. The Asn/Gln metal check tested only the unflipped oxygen position, and used the nitrogen's
    radius. It now tests the oxygen in both positions with its own radius, and locks the flip so the
    oxygen faces the ion.
18. Single-hydrogen rotators measured angular spacing without wrapping around 360°. It is now
    measured on the circle.
19. `_rotateHingeDock` raised a math domain error when rounding pushed a cosine just outside
    [−1, 1], and the Mover was silently dropped. The value is now clamped.
20. His flips recorded the wrong atom for the flip report because a loop variable was reused. The
    NE2 is recorded, as for every other flip.
21. `set_flip_states` never matched in multi-model files and was off by two for model numbers: one
    was added at parsing and another when comparing, and a string model id was compared with an
    integer. Chain ids were also upper-cased, and the `.` "any altloc" wildcard was ignored. All
    of these are fixed.

## Known limitations

* Crystallographic symmetry contacts are not considered. Reduce2 makes disulfides across symmetry
  operators and snaps atoms on special positions; Reduce3 makes neither. Neither case occurs in the
  test set.
* User-supplied restraint CIF files are supported by the library code but not yet exposed on the
  command line.
* Reduce2 options that only produce side files are accepted and ignored: `comparison_file`,
  `output.flipkin_directory`, `output.clique_outline_file_name`, `output.print_atom_info` and
  `profile`.
* Compat mode's bit-for-bit agreement assumes the reference cctbx was built for arm64 (fused
  multiply-add). On other platforms compat mode uses plain arithmetic, which is what x86-64 builds
  normally do.

## Source layout

| File | Contents |
|---|---|
| `src/pipeline.rs` | program flow (`Program.run`) |
| `src/hplace.rs`, `src/riding.rs` | `place_hydrogens`, riding connectivity and parameterization |
| `src/interp.rs`, `src/autolink.rs`, `src/monlib.rs`, `src/names.rs` | restraint interpretation, automatic links, monomer library, atom-name mapping |
| `src/atominfo.rs`, `src/probe.rs` | `getExtraAtomInfo`, Probe dot scoring |
| `src/movers.rs`, `src/optimizer.rs` | Movers, optimizer (exact VE and the compat OptimizerC port) |
| `src/pdbio.rs`, `src/mmcif.rs`, `src/cif.rs`, `src/model.rs` | I/O and the iotbx-style hierarchy |
| `src/cifsource.rs`, `src/lib.rs` | library interface for already-parsed CIF data |
| `src/cell.rs` | cctbx unit-cell arithmetic used by compat mode |
| `tools/` | table generators and the end-to-end comparison script |
| `reference/` | scripts that produce the Reduce2 reference data |

## License

Reduce3 is licensed under the [Apache License 2.0](LICENSE). It is a port of cctbx code, and
the parts derived from cctbx also carry the [cctbx license](LICENSE-cctbx.txt); see
[NOTICE](NOTICE).
