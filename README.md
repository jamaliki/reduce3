# Reduce3

Reduce3 is a Rust reimplementation of cctbx **Reduce2** (`mmtbx.reduce2`). It adds hydrogens to a
macromolecular model and optimizes the rotatable and flippable groups (OH/SH/NH3+/methyl rotations,
Asn/Gln/His flips) by Probe dot scoring. It covers the whole pipeline:

* restraint interpretation (monomer library, links, modifications, disulfides including those to
  symmetry copies, zinc and iron-sulfur cluster coordination, automatic linking such as
  N-glycosylation and glycosidic bonds),
* Reduce2's fallback for residues that only the wwPDB chemical component dictionary describes
  (restraints built from the CCD entry, as Reduce2 does through RDKit, but without needing RDKit),
* riding-hydrogen placement,
* Probe contact scoring, the Movers and the clique optimizer,
* PDB and mmCIF input and output, with the same formatting as iotbx.

It has two modes:

* **fixed** (default) corrects the Reduce2 bugs listed below, goes on where Reduce2 gives up (see
  "Where Reduce2 gives up"), and optimizes the cliques in parallel, exactly unless a clique is too
  dense to search (see "Dense cliques" below).
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
Extra options: `--threads N`, `-q`. `reduce3 --help` lists everything. In fixed mode a residue
without restraints does not stop the run (so `ignore_missing_restraints` has no effect there);
`stop_on_any_missing_hydrogen=True` makes it stop.

## Library use

Reduce3 is also a Rust library. It can work on a model that another program has already
parsed, without writing or re-reading text:

* implement `cifsource::CifSource` (categories as tables of values) for the program's parsed
  mmCIF data block, and call `reduce3::run_cif(&block, &monlib, &params)`;
* stream the result into the program's own document type by implementing `cifsource::CifSink`
  and calling `mmcif::write_cif(&output.structure, &mut sink)` (Reduce2's layout), or
  `mmcif::write_cif_preserving(&output.structure, &block, code, &mut sink)` to write the model back
  into its source block with every other category kept.

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
| Residues only the CCD describes: 1fdo (6MO, Fe4S4 cluster), 2atz (DGT, disulfide to a symmetry copy), 3fx8 (FE2), with and without flips | **byte-identical** in all 6 (Reduce2 run with RDKit) |
| Restraints Reduce2 builds from the CCD, per entry (`reference/harness/dump_ccd_restraints.py`) | identical acceptance, values and order for all 3,093 CCD entries that neither library describes (2,215 accepted, 878 rejected) |

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

**Restraints and output**

22. Residues with CCD-built restraints get the bond and angle values of the GeoStd low-pH or
    neutron variant of their dictionary, but cctbx applied them only to the first residue of each
    name. Every residue gets them now, and in neutron mode the bonds it changes take neutron
    distances, like the bonds it adds.
23. The output dropped everything but the coordinates: Reduce2's mmCIF has none of the input's
    other categories (`_struct_conn`, `_entity`, the sequence schemes, ...), and its PDB has no
    SSBOND, LINK or CONECT records. Fixed mode writes mmCIF input back into its own data block,
    with every category kept, `_atom_site` rebuilt with the input's items and label identifiers
    (new hydrogens take their residue's), atom ids renumbered and `_atom_site_anisotrop` and
    `_atom_type` updated. PDB output keeps SSBOND and LINK and renumbers CONECT.

## Where Reduce2 gives up

Fixed mode carries on in these cases; compat mode stops or loses data exactly as Reduce2 does.

* **CCD entries RDKit rejects.** Reduce2 builds restraints for a residue that only the CCD describes
  through RDKit, and gets none when RDKit cannot read the entry: an unknown formal charge (the Fe
  of A1IW2 in 9hpx, the Ru of RU7 in 5v4h), an aromatic or delocalized bond order, a valence RDKit
  rejects, or missing coordinates. Fixed mode builds them from the CCD geometry in all these cases,
  taking model coordinates where ideal ones are missing and leaving out only restraints it cannot
  measure.
* **Hydrogen bond lengths for those residues.** Reduce2 makes every bond 0.9 times its CCD length,
  which assumes the CCD has neutron-length X-H bonds; entries with X-ray-length coordinates come out
  near 0.87 A. Fixed mode gives bonds to hydrogen GeoStd's X-ray and neutron lengths for the parent
  element and its bond count (`src/h_distances.rs`; C-H 0.97/0.93, N-H 0.86, O-H 0.85 A for X-ray),
  so neutron runs get neutron lengths too.
* **Atom types for those residues.** Reduce2 leaves their atoms untyped, so Probe treats them all
  alike: no hydrogen-bond donors or acceptors, element radii, no polar-hydrogen radius. Fixed mode
  gives them CCP4/GeoStd energy types from the CCD chemistry (element, charge, bonded hydrogens,
  bond orders, aromatic and conjugated rings), so donors, acceptors and radii follow GeoStd. Over
  the 47,921 residues both GeoStd and the CCD describe, the types give the same Probe properties
  as GeoStd's for 97.6% of 2.2 million atoms (`typecheck`; much of the rest is GeoStd typing some
  entries inconsistently). Like GeoStd, fixed mode also builds them at physiological protonation:
  carboxylic, phosphoric and sulfonic acids lose their acidic hydrogen, which the CCD's neutral form
  carries. The hydrogen sets then agree with GeoStd's for 99.5% of a million hydrogens.
* **Atoms of unknown element** (element X, such as UNX in 1h0h and 4iio). Reduce2 deletes them.
  Fixed mode keeps them unchanged; they get no hydrogens and take no part in scoring.
* **Models with no site for a hydrogen** (C-alpha or phosphate traces). Reduce2 stops ("It was not
  possible to place any H atoms"). Fixed mode writes the model unchanged and says so in the report.
* **Residues no dictionary describes** (UNL, or a code missing from the CCD). Reduce2 stops
  ("Restraints were not found"), and with `ignore_missing_restraints=True` deletes their input
  hydrogens. Fixed mode reports them, keeps their input hydrogens (not scored), and places the
  hydrogens of everything else.
* **Hydrogens oriented by an atom on the bond axis.** When the atom that orients a one-neighbor
  group lies on the parent bond axis (cobalt hexammine on a crystallographic two-fold in 9ciy),
  Reduce2 divides by zero. Fixed mode orients the group from another neighbor and optimizes it.
* **Rotatable groups Reduce2 leaves alone.** Reduce2 counts the neighbors of every alternate
  conformation together, so an atom shared by two alternates seems to have too many bonds:
  alternate Ser/Thr hydroxyls get no rotator, alternate Thr/Val methyls are not staggered and a
  histidine split at CG is not flipped. It also rejects an OH whose partner atom has other than two
  or three further bonds (S-OH, P-OH, metal hydroxides), an NH3 whose partner has fewer than three
  (metal ammines), and groups whose dictionary names no reference atom (methanol). Fixed mode gives
  each alternate its own neighbors and builds these rotators.
* **Dense cliques.** Reduce2 searches each clique exhaustively. Rings of Thr/Ser hydroxyls from
  four chains meeting on a channel axis (6een) make cliques of 10 to 18 interacting Movers that
  neither Reduce2 nor exact search finishes. Fixed mode searches a clique exactly when its tables
  stay within a work budget, and otherwise by block coordinate ascent (each Mover, then each
  touching pair, given the others). On 3,509 cliques that exact search can solve, the ascent
  finds the same optimum for 3,497.

## Known limitations

* Symmetry is used for disulfides and, in fixed mode, zinc coordination (a Cys or His binding the
  zinc of a neighboring copy). Reduce2 also bonds atoms to symmetry copies of other metal sites
  through cctbx's asymmetric-unit mappings, snaps atoms on special positions, and skips iron-sulfur
  cluster coordination entirely when a symmetry copy comes within 3.5 A of a cluster; Reduce3 does
  none of these, so compat mode differs from Reduce2 where they occur.
* With `keep_existing_H=True`, Reduce2 applies the low-pH/neutron dictionary values to a CCD-built
  residue only if one of its hydrogens was added (input hydrogens of PDB files have a padded
  element field it does not recognize). Reduce3 does not track that padding and applies them
  whenever the residue has hydrogens.
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
| `src/ccdrestraints.rs`, `src/rdkit_valence.rs`, `src/h_distances.rs` | restraints built from the CCD (Reduce2's RDKit fallback and the fixed-mode builder), RDKit's valence verdicts, GeoStd X-H lengths |
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
