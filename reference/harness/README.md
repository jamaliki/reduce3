# Reduce2 reference dump harness

`dump_ref.py` runs the original cctbx `mmtbx.reduce2` (v3.1.0, `approach=add`, default
parameters) step by step and writes its intermediate state as JSON. Use it to compare the
Rust port (Reduce3) against the original.

## Files

| Path | What |
|---|---|
| `harness/dump_ref.py` | The harness. `reference/env/bin/python harness/dump_ref.py <model> [--flips] [--neutron] [--out-dir D] [--id ID] [--no-pdb] [--indent N]` |
| `harness/run_all.sh` | Runs every `reference/pdbs/*.pdb` with and without `--flips` (sequentially; set `JOBS=n` to run in parallel), then runs `compare_out.py`. Per-run logs go to `dumps/logs/`. |
| `harness/compare_out.py` | Compares `final_atoms` in each dump with `reference/out/<id>H.pdb` / `<id>FH.pdb` (atoms matched by model/chain/resseq/icode/altloc/resname/name, default tolerance 0.001 Å). It also byte-compares the harness PDB with the reduce2 PDB. |
| `dumps/<id>.json`, `dumps/<id>_flips.json` | Dumps (`_flips` = `add_flip_movers=True`; `--neutron` adds `_neutron`). |
| `dumps/<id>[_flips].pdb` | The model text reduce2 would write (`model_as_pdb(output_cs=...)`). |

### How it reproduces reduce2

- The `reduce2.Program` object is built with the same `CCTBXParser` that `mmtbx.reduce2`
  uses. That gives the same PHIL defaults, including the probe overrides
  `bump_weight=100` and `hydrogen_bond_weight=40`, and the same DataManager options.
- The harness then runs the steps of `Program.run()` itself:
  1. `data_manager.get_model()`
  2. select `~element X`
  3. `get_output_crystal_symmetry`
  4. `add_crystal_symmetry_if_necessary(box_cushion=5)`
  5. reduce2's own `_AddHydrogens()`, which calls `reduce_hydrogen.place_hydrogens`
  6. `Optimizers.Optimizer(...)`, with exactly the arguments reduce2 passes
  7. delete `getHydrogensToDelete()`
  8. `_ReinterpretModel(False)`
- Some data only exists inside `Optimizer.__init__`. The harness reads it through
  monkeypatches that only observe and are removed afterwards. They never change
  arguments or return values. The patched calls are:
  - `Helpers.getBondedNeighborLists`
  - `Helpers.getExtraAtomInfo`
  - `Optimizer._PlaceMovers`
  - `Optimizers.OptimizerC`
  - `Movers.MoverTetrahedralMethylRotator`
- **Verification.** All 26 runs (13 structures, with and without flips) write a PDB that is
  byte-identical to the one in `reference/out`. All final coordinates agree within 0.001 Å; the
  largest difference, 0.0009 Å, is PDB rounding. `--neutron --flips` on 1crn is also
  byte-identical to a fresh `mmtbx.reduce2 use_neutron_distances=True add_flip_movers=True`
  run.

## JSON schema (`schema_version` 1)

### Conventions

- **i_seq.** Unless a field says otherwise, `i_seq` means the atom index right after H
  placement: the model the Optimizer sees. In `after_h_placement`, `i_seq == list index`.
- **xyz.** Coordinates are `[x, y, z]` floats at full double precision.
- **Non-finite floats.** These are written as `null`.
- **Atom names.** Names, resseq, icode and altloc are the raw hierarchy strings, so they keep
  their padding (for example `' N  '` and `'   1'`).

### Atom record

Used by `after_h_placement` and `final_atoms`:

```
{i_seq, model_id, chain_id, resseq, resseq_int, icode, resname, altloc, name, element,
 charge, xyz, occ, b, hetero}
```

`final_atoms` entries add `i_seq_after_h_placement`. Their own `i_seq` is the index after
deletion and reinterpretation.

### Top-level keys

| Key | Content |
|---|---|
| `id`, `input_file`, `options`, `reduce2_version` | Run identification. `options = {add_flip_movers, use_neutron_distances}`. |
| `params` | The reduce2 PHIL values used, including `probe.*`. `optimizer_defaults_not_passed_by_reduce2` lists the Optimizer keyword defaults that reduce2 does not pass: `modelIndex=0`, `useNeutronDistances=False`, `minOccupancy=0.02`. |
| `symmetry` | `unit_cell` and `space_group` used for the calculation (possibly a box), and `output_cs_unit_cell` (written to CRYST1, or null). |
| `h_placement` | Counters and labels from the `place_hydrogens` object: `n_H_initial`, `n_H_final`, `no_H_placed_mlq`, `site_labels_*`, `sl_removed` (for example SS-bond HG). |
| `after_h_placement` | (1) Every atom right after `_AddHydrogens()`, before any Optimizer change. |
| `restraints` | (2) From `model.get_restraints_manager().geometry` at that point. See below. |
| `riding_existed_before_optimizer`, `riding` | (3) `riding_h_manager.h_parameterization`: one entry per H, `{index, htype, ih, a0, a1, a2, a3, a, b, h, n, disth}`. This is the manager the Optimizer used. It is always the one built by `place_hydrogens` with `use_ideal_dihedral=True`. |
| `atom_info` | (4) One entry per atom. See below. |
| `initial_extra_atom_info_warnings` | The warnings string from `getExtraAtomInfo`, for example carbonyl radius overrides and aromatic acceptors. |
| `optimizer_runs` | One entry per `_PlaceMovers` call, that is, per optimized model index × alternate, in execution order. See below. |
| `movers` | (5) Every Mover from every run. See below. |
| `staggered_methyls` | Methyls that `MoverTetrahedralMethylRotator` re-staggered during placement. These are not Movers, but they do move coordinates. `{carbon, carbon_label, atoms, xyz_before, xyz_after, axis_origin, axis_dir, offset}` |
| `after_optimization_before_deletion` | `{i_seq, xyz}` for every atom after the Optimizer, before H deletion. |
| `info_text`, `warnings_text` | (6) The full `opt.getInfo()` (verbosity 2) and `opt.getWarnings()` strings. |
| `hydrogens_to_delete` | (6) `{i_seq, label, id_str}`, sorted by i_seq. `label` is reduce2's own format, for example `chain A HIS 68 HE2`. |
| `final_atoms` | (7) Atoms reduce2 writes, after deletion and `_ReinterpretModel(False)`. |
| `timings` | (8) Seconds. See below. |
| `program_log` | Parser, validate and H-placement log text. |

### `restraints`

`counts` holds totals. Each list:

| List | Entries |
|---|---|
| `bonds` | All simple bond proxies: `{i, j, distance_ideal, weight, slack, origin_id}` |
| `bonds_asu` | All ASU bond proxies: `{i, j, j_sym, distance_ideal, weight, slack, origin_id}` |
| `angles_with_h` | Angle proxies with any H/D: `{i, j, k, angle_ideal, weight, origin_id}` |
| `dihedrals_with_h_end` | Dihedral proxies with H/D at i or l: `{i, j, k, l, angle_ideal, periodicity, weight, origin_id, alt_angle_ideals?}` |
| `planarities_with_h` | Planarity proxies containing any H/D: `{i_seqs, weights, origin_id}` |

### `atom_info[]`

Each entry has:

- `i_seq`
- `energy_type`: from `model._type_energies`
- `model_vdw_radius`: `model.get_specific_vdw_radius(i, False)`
- `model_h_bond_type`: `model.get_specific_h_bond_type(i)`
- `model_ion_radius`: ions only
- `bonded`: the neighbor i_seqs from the Optimizer's `getBondedNeighborLists` call, in list
  order. Order matters, because Movers use `[0]`.
- `initial`: the ExtraAtomInfo right after `Helpers.getExtraAtomInfo`
- `final`: the ExtraAtomInfo from `Optimizer._extraAtomInfo` after the Optimizer finished

`initial` and `final` both have the form
`{vdwRadius, isAcceptor, isDonor, isDummyHydrogen, isIon, charge, altLoc}`.

`final` differs from `initial` for these reasons:

- `fixupExplicitDonors` clears `isDonor` on heavy atoms that carry explicit H.
- Water O is forced to acceptor and non-donor.
- HisFlip state or fix-up and ionic lock-down change some atoms.

### `optimizer_runs[]`

Each entry has:

| Field | Content |
|---|---|
| `run`, `alt`, `model_id` | Run identification. |
| `report_model_index` | The model index printed in the report. It is `-1`; see the quirks. |
| `n_atoms_considered` | Includes phantoms. |
| `conformer_atom_i_seqs` | Atoms considered in this run. |
| `phantom_hydrogens[]` | `{i_seq (> max model i_seq), xyz, occ, b, parent_i_seq, extra_info}` for water phantoms. |
| `extra_info_before_placement[]`, `xyz_before_placement[]` | State right before Mover placement, after phantoms and `fixupExplicitDonors`. Per conformer atom. |
| `delete_atoms_from_placement` | H deleted unconditionally during placement: HIS ionic lock-down. |
| `mover_indices` | Indices into `movers`. |
| `report` | The parsed `BEGIN REPORT ... END REPORT` block. See below. |
| `num_calculated_atoms`, `num_cached_atoms` | Counters from OptimizerC. |

`report` has this shape:

```
{model_index, alt, raw_lines[], groups[{kind: set|singletons, size?, initial_total?, final_total?, raw_line}],
 entries[{info, initial_score, final_score, pose, group_index, angle_deg?, flipped?, flag, raw_line}]}
```

### `movers[]`

Each entry has:

| Field | Content |
|---|---|
| `run`, `alt`, `index_in_run` | Insertion order, matching the "Added Mover… N" lines. |
| `class` | For example `MoverSingleHydrogenRotator`, `MoverNH3Rotator`, `MoverAromaticMethylRotator`, `MoverAmideFlip` or `MoverHisFlip`. |
| `info_at_placement`, `info` | The `_moverInfo` string, without and with ` Initial score: x`. |
| `in_final_optimizer_movers` | Whether the Mover is in `Optimizer._movers` at the end. |
| `atoms`, `atom_labels` | i_seqs of `CoarsePositions().atoms`. |
| `coarse_positions[c][j]` | Position of `atoms[j]` in coarse state `c`. It can be shorter than `atoms`: flips move only the first 5 or 9 atoms. |
| `preference_energies[c]` | Per coarse state. |
| `coarse_extra_infos[c][j]`, `coarse_delete_mes[c][j]` | Empty lists for rotators and amide flips; filled for HisFlip. |
| `num_fine_positions[c]` | Fine positions per coarse state. |
| `fixups[c]` | `FixUp(c)` as `{atoms, positions, extra_infos, delete_mes}`. `atoms` is always listed, even when `positions` is empty. |
| `rotator` | Rotators only: `{axis_origin, axis_dir, offset, coarse_range, coarse_step_degrees, fine_step_degrees, do_fine_rotations, has_preference_function, preferred_orientation_scale, coarse_angles, fine_angles}` |
| `flip` | Flips only: `{non_flip_preference, enabled_flip_states?, enable_fixup?}` |
| `final` | `{coarse_index, fine_index (-1 = none), score}`, read directly from that run's `OptimizerC`: `GetCoarseLocation`, `GetFineLocation`, `GetHighScore`. |
| `report` | The matching parsed report entry, which also keeps `raw_line`. |

### `timings`

| Key | Meaning |
|---|---|
| `h_placement_s` | `_AddHydrogens` wall time. |
| `optimization_s` | Optimizer construction wall time, including hook overhead. |
| `optimization_hook_overhead_s`, `optimization_minus_hook_overhead_s` | The hook overhead, and the optimization time with it subtracted. |
| Other keys | Load, dump, delete and reinterpret times. |

## Pipeline quirks worth knowing when porting

- **Only the last model is optimized.** reduce2 does not pass `modelIndex`, so the default
  `0` gives `range(-1, 0)`.
  - For multi-model files such as 1d3z (10 models), only the last model gets Movers,
    methyl staggering and phantom H.
  - The other models keep their riding-placed H.
  - The report prints "Model -1".
- **Mover lists are reset per alternate.** `Optimizer._movers` and `_moverInfo` are
  reassigned on every alternate. Alternates run in reverse order, for example B then A. At
  the end, only the last alternate's Movers remain on the Optimizer. The harness captures
  every run.
- **`useNeutronDistances` is never passed to the Optimizer.** With
  `use_neutron_distances=True`, H placement and interpretation use neutron distances.
  Phantom water H and `getExtraAtomInfo` still use the X-ray settings: 0.84 Å bond and 1.05
  Å radius.
- **Mover constructors move atoms.** Rotators rotate their H to a canonical start
  orientation. `MoverTetrahedralMethylRotator` is built only to stagger methyls. HIS ionic
  lock-down sets coordinates, sets acceptor flags and deletes H even without `--flips`
  (1xso: 14 H).
