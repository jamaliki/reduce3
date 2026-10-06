# cctbx PDB interpretation as used by Reduce2 hydrogen placement: implementation spec

Status: written from the installed reference environment
(`ref/env/lib/python3.12/site-packages`, `cctbx_base-2026.9`). It was checked by
running that Python on `ref/pdbs/*` (scripts and raw dumps are in
`ref/spec_work/`, see §11).

Path abbreviations used for citations (`file:line`). All paths are relative to
`site-packages/` unless noted otherwise.

| abbrev | file |
|---|---|
| PI | `mmtbx/monomer_library/pdb_interpretation.py` |
| SRV | `mmtbx/monomer_library/server.py` |
| CT | `mmtbx/monomer_library/cif_types.py` |
| LM | `mmtbx/monomer_library/linking_mixins.py` |
| LU | `mmtbx/monomer_library/linking_utils.py` |
| LS | `mmtbx/monomer_library/linking_setup.py` |
| RH | `mmtbx/hydrogens/reduce_hydrogen.py` (the **env** copy; `ref/cctbx_project` has a slightly newer one, see §10) |
| CONN / PARAM / HH | `mmtbx/hydrogens/connectivity.py`, `parameterization.py`, `hydrogens.h` |
| RIDING | `mmtbx/hydrogens/riding.py` |
| MODEL | `mmtbx/model/model.py` |
| IP | `iotbx/pdb/__init__.py` |
| ANI | `iotbx/pdb/atom_name_interpretation.py` |
| CRN | `iotbx/pdb/common_residue_names.h` (+ `modified_aa_names.h`, `modified_rna_dna_names.h`) |
| HIER | `iotbx/pdb/hierarchy.py`; `HIER.cpp` = `ref/cctbx_project/iotbx/pdb/hierarchy.cpp` (the env ships no `.cpp`) |
| GR | `cctbx/geometry_restraints/__init__.py`; `GR/utils.h`, `GR/dihedral.h` |
| RSU / RSB | `mmtbx/ligands/ready_set_utils.py`, `mmtbx/ligands/ready_set_basics.py` |
| MCL | `mmtbx/conformation_dependent_library/mcl.py`, `metal_coordination_library.py`, `mcl_sf4_coordination.py` |

---------------------------------------------------------------------------

## 0. Pipeline context and effective parameters

### 0.1 What `place_hydrogens.run()` (RH:938) feeds into interpretation

These steps happen in order:

1. If the model has no crystal symmetry, call `shift_and_box_model(shift_model=False)`.
   This makes a P1 box and leaves coordinates unchanged.
2. Delete atoms with element `X`. Then delete all existing H/D unless `keep_existing_H`.
3. Run `add_missing_H_atoms_at_bogus_position` (RH:1404). For every atom group:
   * Skip blank-altloc atom groups when the residue group has **more than 2** atom groups.
     This is the case blank + A + B. The H then go into A and B only, and each copy includes the backbone H.
   * Skip water when `exclude_water` is set. This is the Reduce2 default (`exclude_water=True`).
   * `mon_lib_query` (RH:424): call `mon_lib_srv.get_comp_comp_id_and_atom_name_interpretation(resname, names)` (see §1).
     If that fails, fall back to the CCD + RDKit (`get_h_restraints`, RH:284). The fallback writes a throwaway
     `auto_<RES>` restraint object with no `type_energy`.
   * Collect the expected H: the dictionary `atom_dict()` entries with `type_symbol=="H"`.
     First, if `test_for_peptide` (CT:456) is true, remove the free-form-only H found by `_terminal_h` (RH:398).
     These are all H on OXT, plus every N-bound H except one. The one kept is `H` or `D` if the dictionary has
     it, otherwise the first by atom order.
   * For classes common/modified/d amino acid **only**, apply the v2 to v3 rename hack (RH:1499-1506). If a pair
     `(X1,X2)` from `HA*, HB*, HG*, HD*, HE*, HG1*` is present and both have `type_energy=='HCH2'`, then X3 is
     added and X1 is removed. The model therefore gets v3 names `HB2/HB3`. The dictionary uses v2 `HB1/HB2`.
   * missing = expected − names already present. The order is Python `set` order, but atoms are re-sorted later.
   * Each missing H is created with element `"H"`, name `(' '+n).ljust(4)` when `len(n)<4`, and the `hetero`
     and `segid` of the first atom in the group. The **bogus position** of every new H in an atom group is
     `mean(xyz of the atoms in that atom group before addition) + (0.5,0.5,0.5)`. All new H in one group sit on
     that single point. This point matters later (§8.4).
4. `place_n_terminal_propeller` (§4.1). Default `n_terminal_charge='residue_one'`.
5. `pdb_hierarchy.sort_atoms_in_place()` (§6.9), then `atoms().reset_serial()`.
6. Build a new `mmtbx.model.manager(pdb_hierarchy=…, restraint_objects=ro)` and call
   `model.process(make_restraints=True, pdb_interpretation_params=get_reduce_pdb_interpretation_params(neutron))`.
   That call is everything in §1-§7. Inside it:
   `mmtbx.utils.process_pdb_file_srv` (mmtbx/utils/__init__.py:236) builds a **fresh** `server.server()`
   and a fresh `server.ener_lib(use_neutron_distances=neutron)` (lines 261-266).
   Restraint objects (user CIFs and the `auto_*` ones) are loaded into both with `process_cif_object`.
   Then `pdb_interpretation.process(... strict_conflict_handling=False, force_symmetry=True, substitute_non_crystallographic_unit_cell_if_necessary=True)`
   runs, and finally `geometry_restraints_manager(assume_hydrogens_all_missing = not has_hd)`.
7. RH-level post-processing, outside this spec's core but listed for completeness:
   `add_link_h_restraints`, riding setup and idealization (§8), removal of unparameterized H,
   `exclude_H_on_links` (uses origin ids, §5.6), `exclude_H_on_esterified_O`, and `name_prochiral_h`.

### 0.2 Effective `pdb_interpretation` parameters

These are from `get_reduce_pdb_interpretation_params` (RH:485) on top of the defaults in PI:224-575. The values
were printed by the reference Python.

| parameter | value | note |
|---|---|---|
| restraints_library.cdl | **False** | reduce turns CDL off |
| restraints_library.mcl | True | Zn and Fe-S coordination, §5.3 |
| omega_cdl / rdl / hpdl / cdl_nucleotides.enable | False | |
| sort_atoms | True | |
| flip_symmetric_amino_acids | **True** | changes coordinates in place, §6.10 |
| use_ncs_to_build_restraints | True | speed only, §6.11 |
| use_neutron_distances | as requested (Reduce2 `use_neutron_distances`, default False) | |
| automatic_linking.link_metals | **Auto** | behaves as "off" for automatic linking; MCL still runs, §5 |
| link_residues / link_carbohydrates / link_ligands | True / True / True | |
| link_amino_acid_rna_dna / link_small_molecules | False / False | |
| metal_coordination_cutoff | 3.0 | |
| amino_acid_bond_cutoff | 1.9 | |
| inter_residue_bond_cutoff | 2.2 | |
| buffer_for_second_row_elements | 0.5 | |
| carbohydrate_bond_cutoff / ligand_bond_cutoff / small_molecule_bond_cutoff | 1.99 / 1.99 / 1.98 | |
| exclude_hydrogens_from_bonding_decisions | True | effectively a no-op, §5.4 |
| link_distance_cutoff | 3.0 | polymer link and chain-break cutoff |
| disulfide_distance_cutoff / exclusion_distance_cutoff | 3.0 / 3.0 | |
| add_angle_and_dihedral_restraints_for_disulfides | True | |
| dihedral_function_type | determined_by_sign_of_periodicity | |
| chir_volume_esd | 0.2 | |
| peptide_link | cis_threshold 45, apply_all_trans False, discard_omega False, discard_psi_phi **True**, apply_peptide_plane False, omega_esd_override None | |
| rna_sugar_pucker_analysis.enable | True | bond_min 1.2, bond_max 1.8 |
| c_beta_restraints | True | C-beta improper dihedrals, no H |
| secondary_structure / ramachandran / reference_coordinate | off | |
| clash_guard.nonbonded_distance_threshold | None | |
| disable_uc_volume_vs_n_atoms_check, proceed_with_excessive_length_bonds, allow_polymer_cross_special_position | True | |

---------------------------------------------------------------------------

## 1. Monomer-library lookup and atom-name mapping

### 1.1 Library files and environment variables (SRV:16-107)

`find_mon_lib_file(relative_path_components)` (SRV:27) searches in this order:

1. `$MMTBX_CCP4_MONOMER_LIB/<components>` if that environment variable is set. If it is set but the file is
   not there, a request for `geostd_list.cif` returns `None` without searching further.
2. The relative roots, in order: `chem_data/geostd`, `chem_data/mon_lib`, `mon_lib`, `geostd`,
   `ext_ref_files/mon_lib`. Each is resolved with `libtbx.env.find_in_repositories`. The order is **reversed**
   when the request is `mon_lib_list.cif`.
3. `$CLIBD_MON/<components>`.

Resolved in this environment:

| file | path used |
|---|---|
| link/mod/synonym list | `chem_data/mon_lib/list/mon_lib_list.cif` (the copy in `geostd/list/` is **not** used) |
| geostd list (merged on top) | `chem_data/geostd/list/geostd_list.cif` |
| energy library | `chem_data/geostd/ener_lib.cif` (**not** `mon_lib/ener_lib.cif`; `geostd_ener_lib.cif` is not merged because the code after the early `return` at SRV:96 is dead) |
| RNA/DNA chain links and mods | `chem_data/geostd/rna_dna/{chain_link_rna2p,chain_link_rna3p,mod_rna2p,mod_rna3p,mod_rna2p_pur,mod_rna3p_pur,mod_rna2p_pyr,mod_rna3p_pyr}.cif` (SRV:494) |
| CCD (Reduce fallback, saccharide typing, N-terminal test) | `chem_data/chemical_components/<lower first char>/data_<CODE>.cif` (`mmtbx/chemical_components/__init__.py:58-105`) |

The two energy libraries agree on every non-metal type. They differ only in the metal and ion rows (for example MG
`vdw` is 1.39 in geostd and 1.73 in mon_lib). Use the geostd values.

### 1.2 Server construction (SRV:405-506)

* Read `mon_lib_list.cif`, then `merge_and_overwrite_cifs(geostd_list, mon_lib_list)` (SRV:234):
  * A geostd data block whose name does not end in `_list` or `energy` **replaces** the mon_lib block with the
    same name, or is added if absent. Blocks in geostd_list.cif: `link_NAG-ASN`, `link_ALPHA*`/`BETA*`,
    `link_ACE_C-N`, `link_SS`, `link_SSRAD`, `link_PEPTIDE-PLANE`, `link_PRE/POST-BETA-TRANS`,
    `link_ASP_CG-ANY_N`, `mod_5*END`, `mod_3*END`, `mod_COOH`, `mod_ACID-ASP/GLU`, `mod_NH3`, `mod_NH1`,
    `mod_NH1NOTPRO`, `mod_NH2`, `mod_NH2NOTPRO`, `mod_NH2N`, `mod_CF-CBH`, `mod_CF-COH`. `mod_COO` comes from
    mon_lib.
  * For `*_list` blocks, rows are **appended** to the mon_lib loop. Later rows win when the list is turned into
    a dict (SRV:127-170), so geostd `chem_link` and `chem_mod` rows override.
* `convert_all(..., skip_comp_list=True)`. This fills:
  * `comp_synonym_list_dict` (alternative id → standard id). All 25 entries: BOX→BEZ, CEG→CEG-b-D, DAL→ALA-D,
    DLE→LEU-D, DOD→HOH, DPN→PHE-D, DPR→PRO, DTH→THR, DTR→TRP-D, DTY→TYR, DVA→VAL-D, FCA→FUC-a-D, FCB→FUC-b-D,
    FUC→FUC-a-L, GAL→GAL-b-D, GCU→GCU-b-D, GLC→GLC-b-D, H2O→HOH, MAN→MAN-b-D, NAG→NAG-b-D, OH2→HOH, SUL→SO4,
    WAT→HOH, ZN1→ZN, ZN2→ZN.
  * `comp_synonym_atom_list_dict[comp][alt_atom] = atom`. Examples: ILE `CD→CD1`; HOH `OW,OH2,OW0→O`,
    `D1→H1`, `D2→H2`; ALA `HN1→H`.
  * `link_link_id_dict`, `link_link_id_list` (list order matters, §3.3), and `mod_mod_id_dict`.
  * The comp cache `comp_comp_id_dict` starts **empty**. Dictionaries are loaded lazily.
* `process_geostd_rna_dna()` loads the 8 rna_dna files. This defines links `rna2p`/`rna3p` and the mods
  `rna2p, rna3p, rna2p_pur, …`.

### 1.3 `get_comp_comp_id_direct(comp_id)` (SRV:544-678)

1. `comp_id = comp_id.strip().upper()`. Return `None` if it is empty.
2. If it is cached in `comp_comp_id_dict`, return the cached entry. This includes dictionaries pre-loaded from
   restraint objects, which therefore **win over the files**.
3. `std = comp_synonym_list_dict.get(comp_id, "")`.
4. `find_file` (SRV:592). The user-supplied directory is not used by Reduce. For `i_pass in (0,1)` and for
   `trial in (std, comp_id)` (empty entries skipped):
   * geostd: `chem_data/geostd/<trial[0].lower()>/data_<trial>.cif`. Pass 1 is a case-insensitive scan, and it
     **asserts** (crashes) if anything matches there.
   * mon_lib: `chem_data/mon_lib/<trial[0].lower()>/<trial>.cif`. Windows device names (CON, PRN, AUX, NUL,
     COM1…) become `<X>_<X>.cif`. Pass 1 is a case-insensitive scan.

   So the precedence is geostd(std) > mon_lib(std) > geostd(id) > mon_lib(id), and then the case-insensitive pass.
5. `process_cif(file)` caches **every** `comp_*` block in the file under its `chem_comp.id.upper()`. The lookup
   order is then `cache[std]` (and alias `cache[comp_id]=…`), then `cache[comp_id]`. If neither is found, raise
   `Sorry`.

Observed sources: the 20 amino acids, MSE, HOH, A/C/G/U and AD/CD/GD/TD all come from geostd. Single-atom ions
(MG, ZN, NA, CL, CA, K, MN, CU, HG) come from `mon_lib/<x>/<X>.cif`.

`normalize_atom_ids_in_place` (CT:217) is effectively a no-op. The renaming code is commented out, so dictionary
atom ids are used exactly as written. geostd standard residues use: protein **v2 H names** (`HB1/HB2`,
`HA1/HA2` for GLY, `HG11/HG12` for ILE); nucleic acids **v3** names (`OP1`, `C1'`, `H5''`). DNA thymine uses
`C5M/H5M1-3`. Some modified nucleotides (for example PSU) still use v2 names (`O1P`, `O5*`).

### 1.4 Residue name to dictionary (`residue_name_plus_atom_names_interpreter`, IP:497-560)

Called through `SRV.get_comp_comp_id_and_atom_name_interpretation` (SRV:680-710) with
`return_mon_lib_dna_name=True`.

1. `work = resname.strip().upper()`. An empty name gives `(None, None)`.
2. **D amino acids.** `three_letter_l_given_three_letter_d`: DAL DAR DAS DCY DGL DGN DHI DIL DLE DLY DPN DPR DSG
   DSN DTH DTR DTY DVA MED map to the L name. The L name is used for interpretation and `d_aa_residue_name` is
   recorded.
3. If a protein interpreter exists for `work` (20 standard amino acids + MSE, ANI:494-516), then
   `ani = interpreter.match_atom_names(names)` (§1.6).
4. Otherwise, RNA/DNA (IP:113-141): `rna_dna_reference_residue_names` maps
   `A,C,G → ?A,?C,?G`; `U→U`; `T→DT`; `+X` the same; `DA,DC,DG,DT` as is; CNS `ADE,CYT,GUA,URI,THY`;
   `AD,CD,GD,TD → DA…DT`.
   * `?X` becomes **RNA if the residue has `O2'` (or `HO2'`, `2HO*` and other aliases), else DNA**.
     Verified: an `A` without `O2'` becomes DNA (`AD`).
   * If the interpretation has unexpected names and the residue is a single atom named like
     `AD,A,CD,C,GD,G,TD,U` → `(None, None)`.
   * For CNS names with unexpected atoms the interpretation is dropped. `translate_cns_dna_rna_residue_names=None`
     is the default.
   * `work` = interpretation residue name, then mapped to the mon_lib DNA name:
     `A,C,G,U` unchanged; `DA→AD, DC→CD, DG→GD, DT→TD`.
5. `SRV:698`: if `d_aa_rn` is set **and** the D dictionary is **already cached**, return
   `(D dict, None)` with no interpretation. Otherwise return `(get_comp_comp_id_direct(work), ani)`.

   In the fresh per-`process()` server the cache is empty, so the first D residue gets the **L dictionary plus
   mod `PEPT-D`** (CA chirality becomes positive), plus `DIL_chir_02_both` or `DTH_chir_02_both` for DIL/DTH
   (PI:1283-1292). Whether later D residues of the same type also go this way depends on whether something loaded
   the D dictionary into the cache. In practice nothing does during `process()`, so they all do.

### 1.5 Residue classes and per-residue processing in `build_chain_proxies`

See §2 for classes. A residue with `monomer is None` is an "unknown residue" (§1.9).

### 1.6 Atom-name mapping (`monomer_mapping._get_mappings`, PI:1368-1459)

Input:
* `atom_names_given`: names with **all spaces removed**, in hierarchy order. Atoms with element `Q` are ignored
  (PI:1311).
* `mon_lib_names` (from `ani.mon_lib_names()`, or `None`).
* `atom_dict` = the current monomer's `atom_dict()`, after any mods.

**Protein interpreters** (ANI). Each residue has patterns where `h` stands for H or D, so D names are accepted.
Synonyms include `OT1→O`, `OT2/OC→OXT`, `hC→hXT`, `hN→h`, `1hN/1hT/h0A→1h` and similar. Every `1hX` pattern
also accepts `hX1`, and vice versa (`alternative_hydrogen_pattern`).

A `mutually_exclusive_pairs` triple `(1hB,2hB,3hB)` converts v3 to v2: **if the third member is present, then
2nd→1st and 3rd→2nd; otherwise identity**. Residue-specific synonyms: ILE `CD→CD1`, `1hD→1hD1`…; MSE `SED→SE`;
LEU/CYS `1hG→hG`; PHE `1hZ→hZ`.

`mon_lib_name = pattern.upper()`. A leading digit is rotated to the end (`1HB → HB1`).
Unmatched names give `None`.

Full mapping as applied to the names Reduce writes (from `spec_work/data/aa_h_name_mapping.txt`):

```
GLY  HA2->HA1 HA3->HA2           ILE HG12->HG11 HG13->HG12
LEU/PHE/TRP/SER/TYR/CYS/ASN/HIS/ASP  HB2->HB1 HB3->HB2
MET/MSE/GLN/GLU  HB2->HB1 HB3->HB2 HG2->HG1 HG3->HG2
PRO  HB2->HB1 HB3->HB2 HG2->HG1 HG3->HG2 HD2->HD1 HD3->HD2
LYS  + HD2->HD1 HD3->HD2 HE2->HE1 HE3->HE2      ARG + HD2->HD1 HD3->HD2
ALA/VAL/THR  identity
```

Other examples (verified): `HN→H`, `1HB→HB1`, `D/DA/DB2/DB3/DG → H/HA/HB1/HB2/HG`, `1HG2→HG21`,
N-terminal `H1/H2/H3` (or `HT1/HT2`, `1H/2H/3H`) → `H1/H2/H3`, `OT1→O`, `OT2→OXT`, `HXT` kept.

**RNA/DNA interpreter** (`rna_dna_atom_names_interpretation`, IP:457-495). Model names map to reference v3 names,
and then to **mon_lib v2 names** via `rna_dna_atom_names_reference_to_mon_lib_translation_dict` (IP:162-218):

* `C1'→C1*`, `OP1→O1P`, `OP2→O2P`, `OP3→O3T`, `H5'→H5*1`, `H5''→H5*2`, `C7→C5M`, `H71-73→H5M1-3`,
  `HO2'→HO2*`, `HO3'→HO3*`, `HO5'→HO5*`.
* `H2'` → `H2*` for RNA, `H2*1` for DNA. `H2''→H2*2`. `HOP3 → None`.

The alias table (IP:271-447) accepts v2 forms such as `O1P`, `C1*`, `1H5*`, `2HO*`, `H5T` and D variants.

**Per-atom algorithm** (in hierarchy order):

1. `replace_primes = False` if `ani` is set. Otherwise it is `True` if the residue looks like RNA/DNA
   (`residue_analysis` OK or dictionary classified RNA/DNA). Otherwise it is
   `(#names containing "'" > 0 and #names containing "*" == 0)`.
2. `handle_case_insensitive`: if every dictionary id is all upper or all lower case, upper-case the given name.
3. `atom_name = mon_lib_names[i]` if not None, else the given name.
4. If `atom_name` is not in `atom_dict`, build auto-synonyms:
   * Digit rotation: a leading digit moves to the end, or a trailing digit moves to the front.
   * If `replace_primes`, replace `'` with `*` and add that name (and its digit rotation).
   * Take the first candidate in `atom_dict`. If none:
     * prepend the given name;
     * try `monomer.hydrogen_deuterium_aliases()` (`"D"+id[1:] → id` for every dictionary H whose id starts
       with "H", CT:291);
     * then `comp_synonym_atom_list_dict[chem_comp.id]`;
     * else keep the given name.
5. If it is still not found and the residue is RNA/DNA: use `rna_dna_atom_names_backbone_aliases`
   (IP:449-455; only aliases whose reference name is a **backbone** atom: `C1'…O5'`, `OP1-3`, `P`, the sugar H,
   `HO2'`, `HO3'`, `HO5'`, `HOP3`). Map `name → reference → the cif atom id with the same reference`.

   This is how v2 mon_lib names (`O1P`, `C5*`, `H5*1`, `H2*2`, `O3T`) are mapped **back** to the v3 geostd ids
   (`OP1`, `C5'`, `H5'`, `H2''`, `OP3`). Base atoms are not aliased. They must match directly, for example
   `C5M/H5M1` in TD.
6. Record `atom_names_mappings[model_name.strip()] = atom_name` when they differ.

   The first atom with a given mapped name is **expected** if the name is in `atom_dict`, else **unexpected**.
   Later atoms with the same mapped name are **duplicates**.
7. For peptides with `ani is None` only: `_rename_ot1_ot2` and `_auto_alias_h_h1`. If the dictionary has `H1` and
   neither `H1` nor `D1` is expected, then an unexpected `H` (or `D`, but not both) becomes `H1`.
8. `_set_missing_atoms`: the dictionary atoms that were not matched, split into H and non-H.

The bond/angle/dihedral/chirality/plane builders retry ids that are not found (§6.1).

### 1.7 Applying modifications (`CT.comp_comp_id.apply_mod`, CT:426; `monomer_mapping.apply_mod`, PI:1672)

* Each mod is applied at most once per residue (`chem_mod_ids`).
* `atom_list`:
  * `add` appends `chem_comp_atom(new_atom_id, new_type_symbol, new_type_energy, …)`;
  * `delete` removes the atom **and every bond, angle, torsion, chirality or plane row that mentions it**.
    An unknown id is recorded in `skipped_deletes`;
  * `change` overwrites non-empty `new_*` fields, with a rename propagated to every restraint row.
* `bond/angle/tor/chir/plane` lists: `add` appends, `delete` removes matching rows, `change` overwrites the
  non-empty `new_*` values. Torsion matching requires the exact ordered 4-tuple. Angle matching accepts
  1↔3 swapped.
* If the mod `name` contains "terminus", the residue is flagged `is_terminus`.
* After every mod, `_get_mappings()` is **rerun** on the new dictionary.

### 1.8 Mods applied by `resolve_unexpected` (PI:1539-1660)

This runs only when `incomplete_info is None`. Peptide residues that contain only CA, only N/CA/C, only backbone,
or only backbone+CB, and RNA residues with only P, skip it. See §4 for the terminal rules. It also covers:
* `HXT` → `COOH`; else `OXT` → `COO`;
* `HC`/`DC` without `OC` (common/modified amino acid only) → `CF-COH`; `HBC` → `CF-CBH`;
* GLU `HE2` → `ACID-GLU`; ASP `HD2` → `ACID-ASP`.

### 1.9 Atoms and residues that are not in the dictionary

* **Unexpected or duplicate atoms** get **no** restraints of any kind, so they are isolated in the bond graph.
  Their energy and H-bond type is the Python value `False`, stored as the **string `'False'`** by
  `conformer_i_seq.convert` (PI:2739).
  The nonbonded energy registry leaves them `''` (unknown). Their scattering type comes from the element column.
  For H atoms this also triggers the pH fix-up at PI:6193-6230. That looks for
  `geostd/<x>/data_<RES>_neutron.cif`, then `data_<RES>_pH_low.cif` and may add bond and angle proxies and energy types for
  those H in the nonbonded registry only (`type_energies` stays `'False'`).
* **Dictionary atoms with `type_energy None`**, for example Reduce's CCD/RDKit `auto_*` dictionaries: type
  `'None'`, H-bond type `'None'`.
* **Unknown residue** (`monomer is None`, PI:2911-2946):
  * If it is a single-atom residue and `ad_hoc_single_atom_residue` (PI:66) recognises the element (list at
    PI:39-43: `ZN CA MG CL NA MN K FE CU CD HG NI CO BR XE SR CS PT BA TL PB SM AU RB YB LI KR MO LU CR OS GD TB LA F AR AG HO GA CE W SE RU RE PR IR EU AL V TE SB PD U I S`,
    or NH3/CH4 special cases) **and** `ener_lib` has that element type, then energy type = element and H-bond
    type from ener_lib.
  * Otherwise every atom gets energy type `''` and H-bond type `'N'`.
  * The residue always counts as a chain break.
* With `stop_for_unknowns=False` (Reduce2's `stop_on_any_missing_hydrogen` default), unknown nonbonded types are
  tolerated (mmtbx/utils/__init__.py:358-359).

### 1.10 Energy and scattering types (PI:2886-2909 and `type_symbol_registry_base`, PI:985-1060)

For each residue with a monomer, `type_energies[i_seq] = atom_dict[mapped_name].type_energy` and
`type_h_bonds[i_seq] = ener_lib.lib_atom[type_energy].hb_type`, or `None` if that type is missing. An atom that is
neither unexpected nor found in `atom_dict` after mapping (for example an element-`Q` atom, which is ignored by the
mapper) raises `Sorry`.

`type_energies` is a dict updated per conformer, so the last conformer wins for shared atoms. It is converted to a
`flex.std_string` in i_seq order: `False→'False'`, `None→'None'`.

The nonbonded registry (`assign_from_monomer_mapping`) with `strict_conflict_handling=False`:
* first assignment wins;
* a later different symbol: if the expected-atom counts are equal, raise `Sorry`; otherwise the residue with
  **more** expected atoms wins;
* charge comes from the atom's charge field. A single-atom residue whose name equals its element gets charge 99.

The scattering type registry starts from the element column. `D` is kept when the dictionary says `H`. A
conflicting element is overridden silently only for the trusted standard amino acids (PI:905).

---------------------------------------------------------------------------

## 2. Residue classes

### 2.1 `iotbx.pdb.common_residue_names_get_class(name, consider_ccp4_mon_lib_rna_dna=False)` (CRN, C++)

* `padded = name`, **left-padded with spaces to 3 characters** if shorter. Longer names are not padded. The
  comparison is exact and **case-sensitive** (`'hoh'` → `other`). The first match in this order wins:

| class | members (exact padded strings) |
|---|---|
| `common_amino_acid` | GLY ALA VAL LEU ILE MET **MSE** PRO PHE TRP SER THR ASN GLN TYR CYS LYS ARG HIS ASP GLU |
| `d_amino_acid` | DAL DAR DAS DCY DGL DGN DHI DIL DLE DLY DPN DPR DSG DSN DTH DTR DTY DVA MED |
| `modified_amino_acid` | 2311 names in `modified_aa_names.h` (`spec_work/data/modified_aa_names.txt`) |
| `common_rna_dna` | `"A  ","C  ","G  ","T  ","U  ","  A".."  U"," A ".." U ","+A ".."+U "," +A".." +U","DA ","DC ","DG ","DT "," DA"," DC"," DG"," DT","ADE","CYT","GUA","THY","URI"` |
| `modified_rna_dna` | 834 names in `modified_rna_dna_names.h` (`data/modified_rna_dna_names.txt`) |
| `ccp4_mon_lib_rna_dna` | only if the flag is set: `AD,CD,GD,TD` and `Ad,Cd,Gd,Td`, each as `"XX "` or `" XX"` |
| `common_water` | HOH H2O OH2 DOD OD2 WAT |
| `common_small_molecule` | GOL PO4 SO4 |
| `common_saccharide` | NAG NDG MAN BMA FUC FUL BGC GLC |
| `common_element` | ZN CA MG CL NA MN K FE CU CD HG NI CO BR XE SR CS PT BA TL PB SM AU RB YB LI KR MO LU CR OS GD TB LA F AR AG HO GA CE W SE RU RE PR IR EU AL V TE SB PD, in padded forms `"XX "`, `" XX"` (one-letter: `"K  "," K ","  K"`, and the same for F, W, V) |
| `other` | everything else |

Quirk: the generated lists contain entries shorter than 3 characters (1 in the amino-acid list, `MA`; 40 in the
RNA/DNA list, for example `0A`, `DI`, `DU`, `PU`). The query is padded but the set entries are not, so these
entries **can never match**. Verified: `get_class('0A') == 'other'`.

### 2.2 `comp_comp_id.set_classification` (CT:456-499) and the monomer classification used for polymer linking

* `test_for_peptide`: the dictionary has N, CA, C, O **and** bonds `CA-N`, `C-CA` and a `C-O` bond whose
  `type != "single"` (`coval` counts). Result `"peptide"`.
* `test_for_rna_dna` (`iotbx/pdb/rna_dna_detection.classification`): needs P, OP1/OP2 (or O1P/O2P), O5' (or O5*)
  and all of C5' C4' O4' C3' O3' C2' C1', plus all of the required bonds `OP1-P OP2-P O5'-P C1'-C2' C2'-C3'
  C3'-C4' C3'-O3' C4'-C5' C4'-O4' C1'-O4' C5'-O5'`. A bond `C2'-O2'` makes it RNA, otherwise DNA. The suffix is
  `v2` for star names and `_mixed` for mixed names.
* Exactly one positive gives that class. Otherwise `water` if the atoms are exactly {H1,H2,O}, else
  `"undetermined"`.
* `mm.is_rna_dna` also becomes True when `iotbx.pdb.rna_dna_detection.residue_analysis(atoms)` reports no
  problems (PI:1325-1335). This covers modified nucleotides with any dictionary.

### 2.3 `linking_utils.get_classes(atom)` (LU:225-320), used only by automatic linking

1. `gc = get_class(resname, consider_ccp4_mon_lib_rna_dna = (atom group has >1 atom))`. `UNK` is treated as
   `common_amino_acid`. Then `modified_amino_acid → "other"` and `modified_rna_dna → "other"`.
2. The flag `common_saccharide` is set if `gc == 'common_saccharide'`, or otherwise if the CCD `_chem_comp.type`
   is one of `SACCHARIDE`, `D-/L-SACCHARIDE[, ALPHA/BETA LINKING]` (LU:21-28). Residues in
   `one_letter_given_three_letter` count as L-peptide and HOH as non-polymer without a CCD lookup. This flag is
   checked first.
3. `important_only` = the first class found in the order `common_saccharide, common_water, common_element,
   common_small_molecule, common_amino_acid, common_rna_dna, ccp4_mon_lib_rna_dna, other, …`, after
   `_filter_for_metal`:
   * a `common_element` whose element is in the metal list becomes `"metal"`;
   * an `other` single-atom residue with a metal element becomes `"metal"`.

   The metal list (LS:21-26): `ZN CA MG NA MN K FE CU CD HG NI CO SR CS PT BA TL PB SM AU RB YB LI MO LU CR OS GD TB LA AG HO GA CE W RU RE PR IR EU AL V PD U SB SE TE`.
4. `other` residues with ≥4 atoms named `C, CA, N, O, OXT` also set the `uncommon_amino_acid` flag.

### 2.4 `test_for_peptide` in Reduce

RH:1482 uses the **unmodified** dictionary to decide whether to drop the free-form terminal H (§0.1).

---------------------------------------------------------------------------

## 3. Polymer linking between consecutive residues (`build_chain_proxies`, PI:2748-3290)

### 3.1 Iteration

For each model, each chain and each conformer (`chain.conformers()`: one per altloc, or `""`), the code walks
`conformer.residues()` in order. Each residue holds the blank atoms plus that altloc's atoms. Settings:
`i_conformer` = index of the altloc in sorted order, with blank = 0 (PI:3721); `is_first_conformer_in_chain =
(j_conformer==0)`.

The registry tables are re-initialised per **model** only (PI:3878).

For residue *i* (`mm`) with previous residue `prev_mm`:

* `mm.monomer is None` → unknown residue (§1.9), `n_chain_breaks += 1`, no link.
* `prev_mm is not None and not residue.link_to_previous` → chain break. `link_to_previous` is False only for the first
  residue after a PDB `BREAK` record (`construct_hierarchy.cpp:148`); models built from a hierarchy keep the flag.
* otherwise, if `prev_mm.monomer` exists: `prev_mm.lib_link = get_lib_link(prev_mm, mm)` (§3.2). `None` gives
  "Not linked" and no restraints.
* If there is a link (PI:2961-3063):
  1. `mod_id_1` / `mod_id_2` of the chem_link are applied to prev/this residue, followed by
     `resolve_unexpected()`, but only if `within_linking_cutoff()` is true. Quirk: that function checks **only the
     first** link bond, and only that its distance is ≤ `link_distance_cutoff`. A missing atom returns 999,
     which is truthy. **None of the peptide or RNA chain links have mods.**
  2. Link bonds go through `add_bond_proxies(..., sites_cart, distance_cutoff=link_distance_cutoff=3.0)`. A bond
     longer than 3.0 Å is **not added**, and its sorted i_seq pair goes into `broken_bond_i_seq_pairs`, which
     counts as a chain break. Any link angle, dihedral, chirality or plane containing both atoms of a broken pair
     is skipped (`involves_broken_bonds`, PI:847). Link bonds always use `value_dist`, never the neutron value.
     Origin id 0.
  3. Link angles, dihedrals (with `chem_link_id`, `sites_cart`, `peptide_link_params`; §3.4), chiralities and
     planes are added. The order is bonds, angles, dihedrals, chir, planes, **before** the residue's own
     restraints.
* Then the residue's own restraints are added in the order: scattering and energy types, bonds, angles,
  dihedrals, chiralities, planes. A CYS SG (chem_comp.id == "CYS") is recorded for disulfides (PI:3149).

### 3.2 `get_lib_link(m_i, m_j)` (PI:1903-2000)

1. Either residue is water → `None`.
2. **Both residues peptide** (`monomer.is_peptide()`) → `get_lib_link_peptide` (PI:1893):
   `NMTRANS` if `m_j` has an expected atom `CN`; else `PTRANS` if `m_j.monomer.chem_comp.id == "PRO"`
   (D-PRO is interpreted as PRO so it qualifies; HYP and other analogues do not); else `TRANS`.
   **CIS, PCIS and NMCIS are never chosen here.** Only the label and the omega ideal change later (§3.4).
3. **Both residues RNA/DNA** (`mm.is_rna_dna or monomer.is_rna_dna()`) → `rna2p` if `m_i.is_rna2p`, else
   `rna3p`. DNA, and RNA whose pucker could not be analysed, always get `rna3p` (`is_rna2p` stays None). The
   mon_lib `p` link is never used for chains.
4. Otherwise, a generic search over `link_link_id_list` in file order:
   * skip links named in `non_chain_links=("SS-bridge",)`;
   * skip links whose comp, mod and group are all empty;
   * `link_match_one` (PI:1798) normalises the comp group (`L-peptide, D-peptide → "peptide"`;
     `DNA, RNA → "DNA/RNA"`). A link side matches if its `comp_id` is empty or equal (case-insensitive), **and**
     its `group_comp` is empty or equal, or the residue group is empty while the comp_id matched with length > 0.
   * For each match, count `n_unresolved_bonds` and `n_unresolved_angles` against the two `monomer_atom_dict`s.
   * `matches.sort()` uses the non-transitive `link_match.__lt__` (PI:1867). It returns True as soon as **any**
     one of these is better: fewer unresolved bonds, fewer unresolved angles, longer comp_id match 1, 2, longer
     group match 1, 2. CPython sorts fewer than 64 items with binary insertion sort, so reproduce that exactly.
   * Take `matches[0]`. Require `is_proper_match` (some comp_id or group length > 0), else `None`.

   Examples: `ACE_C-N` (ACE→peptide), `NH2_CTERM` (peptide→NH2), `NME_N-C`, `FOR_*`. These are
   **origin_id 0** chain links.

### 3.3 Peptide and nucleic chain link contents

These are the merged definitions as loaded. Full dump: `spec_work/data/links_and_mods_merged.txt`.

```
TRANS  (grp peptide -> grp peptide)   [CIS identical except omega 0.0]
  bond  1C-2N 1.329 (0.014)
  angle 1O-1C-2N 123.0(1.6)  1CA-1C-2N 116.2(2.0)  1C-2N-2H 124.3(3.0)  1C-2N-2CA 121.7(1.8)
  tor   psi 1N-1CA-1C-2N 160 esd30 per2 [discarded]   omega 1CA-1C-2N-2CA 180 esd5 per0
        phi 1C-2N-2CA-2C 60 esd20 per3 [discarded]
  plane plane1: 1CA 1C 1O 2N (0.02)     plane2: 1C 2N 2CA 2H (0.02)
PTRANS (peptide -> PRO)
  bond 1C-2N 1.341(0.016); angles O-C-N 123.0(1.6) CA-C-N 116.9(1.5) C-N-CD 125.0(4.1) C-N-CA 122.6(5.0)
  plane2: 1C 2N 2CA 2CD (0.05); no H anywhere
NMTRANS (peptide -> residue with CN): as TRANS with 2CN instead of 2H; plane2 esd 0.05
rna3p / rna2p (geostd/rna_dna):
  bond 1O3'-2P 1.607(0.012); angles O3'-P-O5' 104.0(1.9), O3'-P-OP1 108(3), O3'-P-OP2 108(3), C3'-O3'-P 119.7(1.2)
  tor epsilon C4'-C3'-O3'-P (-140 | 2p:-110) esd35 per1; zeta C3'-O3'-P-O5' (172 | 2p:145) esd30 per3;
      alpha O3'-P-O5'-C5' (300 | 2p:165) esd20 per3
```

**Restraints involving H that come from links:**
* TRANS and CIS: angle `C(i-1)-N-H 124.3±3.0`, and `plane2` (C, N, CA, H). If H is absent the plane has three
  atoms and is dropped silently (planes need ≥4 resolved atoms).
* NMTRANS: the same with CN.
* rna2p/rna3p: none.

**Links never delete atoms**, because no chain link carries a mod. The amide H of every residue comes from the
residue dictionary. The free-amine extra H are never added by Reduce (§0.1 `_terminal_h`). Dictionary torsion
`Var_01 C-CA-N-H 170 esd180 per72` is present in most geostd amino acids (§6.3).

### 3.4 Omega, cis/trans, psi/phi (PI:2242-2417)

* `chem_link_id ∈ {TRANS,PTRANS,NMTRANS,CIS,PCIS,NMCIS}` and `tor.id ∈ {psi,phi}` → skipped
  (`discard_psi_phi=True`).
* `tor.id=="omega"` and `chem_link_id ∈ {TRANS,PTRANS,NMTRANS}`: compute the model dihedral. If
  `|angle_delta_deg(model, 180, per0→1)| > 180 − 45` (that is, |ω| < 45°), set the label to CIS/PCIS/NMCIS and
  the proxy `angle_ideal = 0`. Every other restraint stays the TRANS/PTRANS one. The PCIS-specific angle values in
  mon_lib are therefore never used.
* `periodicity = tor.period` (`determined_by_sign_of_periodicity`). `period ≤ 0` means harmonic in cctbx.

### 3.5 RNA sugar pucker → mods (PI:1325-1366; `rna_sugar_pucker_analysis.py`)

This applies only when the residue is not a peptide dictionary and `residue_analysis` reports no problems.
DNA is detected but returns before the pucker step (`if not ra1.is_rna: return`, where `is_rna = O2' present`).

Analysis:
* `delta = dihedral(C5',C4',C3',O3')`, normalised to [0,360).
* `P(i+1)` is used only if it is within 1.8 Å of O3'.
* `c1p_outbound` = the closest N or C (no prime) to C1' within 1.463+0.5 Å.
* All the sugar bonds checked must be within [1.2, 1.8] Å, otherwise no decision.
* `is_2p` = (distance from P to the C1'-outbound line < 2.9 Å). If P is absent, use (distance from O3' to that
  line < 2.4 Å).

Primary mod by `modernize_rna_resname`:
* `A,G → rna2p_pur / rna3p_pur`;
* `C,U → rna2p_pyr / rna3p_pyr`;
* anything else (including modified nucleotides that pass `residue_analysis`) → `rna2p / rna3p`.

These mods change only heavy-atom bonds, angles and torsions. They do change the heavy angles that riding uses as
`angle_a1a0a2` for sugar H (§8). `is_rna2p` also chooses the chain link.

---------------------------------------------------------------------------

## 4. Terminal modifications

### 4.1 Reduce side: N-terminal H (RH:1367-1402; RSU:58-197; RSB:9-60)

`place_n_terminal_propeller`, for every model and chain:
1. `rgs = chain.residue_groups()[0]`. With `residue_one`, skip unless `rgs.resseq_as_int()==1`.
   `first_in_chain` takes any first residue. `no_charge` skips the step entirely.
2. Find `N`: the first atom group with an atom named `N`. If there is none, skip.
3. `bonds_in_restraints(N, exclude_hydrogens=True)` reads the **CCD** file of the residue and lists N's non-H
   bond partners. If there are ≥2 (≥3 for PRO), skip: N is already substituted.
4. For common, modified and D amino acids, delete the atom named `H` from every atom group of the residue group.
5. `add_n_terminal_hydrogens_to_residue_group(rgs)` with `bonds=None` and `retain_original_hydrogens=True`. For
   each altloc combination that has N, CA and C (`generate_atom_group_atom_names`; blank atoms are merged into
   each altloc):
   * `proton = 'D'` if the atom group contains only D among its H, otherwise `'H'`.
   * dihedral = 120°, or `dihedral(H,N,CA,C)` if an `H` still exists.
   * `rh3 = construct_xyz(N, 1.0 Å, CA, 109.5°, C, dihedral)`. This gives three points
     `rn + 1.0*(sinα(cos(φ+k·120°)e1 + sin(φ+k·120°)e2) − cosα e0)`, k = 0, 1, 2, where
     `e0=unit(N−CA)`, `e1=unit((C−CA) − ((C−CA)·e0)e0)`, `e2=e0×e1`.
   * If ≥3 of `H,H1,H2,H3,HT1,HT2` (or the D names) already exist, do nothing.
   * Names are `H1`, `H2`, `H3` (or `D1-3`). **PRO skips `H1`**, so PRO gets `H2` at `rh3[0]` and `H3` at
     `rh3[1]`. Others get `H1→rh3[0]`, `H2→rh3[1]`, `H3→rh3[2]`. An existing name is skipped. Occupancy and B
     are copied from N, hetero from the first atom of the group. The new atoms are appended to **the atom group
     of N**.
   These coordinates are later replaced by riding (§8). They only affect `check_propeller_order`.

### 4.2 Interpretation side (`resolve_unexpected`, PI:1539-1660)

**Peptide**, when `ani` is set (standard residues). `u` = unexpected atoms renamed to their mon_lib names:
* 3 of `H1/H2/H3` unexpected and the dictionary has `H` → **mod NH3** ("NH3-terminus"): N type `NT3`; add `H1,H2,H3`
  (type `HNT3`); delete `H` (and through it `Var_01`, `CA-N-H` and so on); bonds `N-H1/2/3 0.89` (neutron 1.04),
  esd 0.02; `N-CA` becomes 1.491; angles `H-N-H` and `H-N-CA` 109.47±3.0. No torsions.
* 2 of them → PRO: **NH2** ("NH2-terminus_for_proline": `HN1,HN2` type `HNH2`, bonds 0.96 (n 1.04), angles
  `HN-N-CA` and `HN-N-CD` 109.47, `HN1-N-HN2` 109.47). Non-PRO: **NH2NOTPRO** (the same without the CD angles).
  `mon_lib_names` H1/H2/H3 are then renamed in order to `HN1, HN2`.
* 1 of them → PRO: **NH1** (`HN`; `HN-N-CA` and `HN-N-CD` 120). Non-PRO: **NH1NOTPRO**. Renamed to `HN`.

When `ani is None`, the same mods are triggered by name: `H1|1H & H2|2H & H3|3H` (needs `H` in the dictionary) →
NH3; `HN1|1HN & HN2|2HN` → NH2; `HN` → NH1; any one of H1/H2/H3 → NH3; any HN1/HN2 → NH2.

C terminus:
* `HXT` → **COOH** (`O`, `OXT`, `HXT`; C=O 1.185 double, C-OXT 1.33, OXT-HXT 0.948 (n 1.02); angle
  `C-OXT-HXT 108`; torsion `O-C-OXT-HXT 0 per0`; plane `oxt` with HXT).
* else `OXT` (unexpected; geostd standard amino acids have no OXT) → **COO**: O and OXT type `OC`, both C-O
  1.231 `deloc`, angles CA-C-O and CA-C-OXT 121, O-C-OXT 118, torsion `N-CA-C-OXT 160 per2`, plane
  `CA C O OXT`. Reduce never adds HXT.

GLU `HE2` → ACID-GLU, ASP `HD2` → ACID-ASP. These come from mon_lib names, so a v3 `HD2` maps to `HD2`.

**RNA/DNA**, when `ani` is set:
* `ani.have_op3_or_hop3` → **p5*END**: adds `O3T`(OP) and `HOP3`, bond `O3T-P 1.48`, plus angles written with v2
  names that fail on v3 dictionaries. The model `OP3` maps to `O3T` because `OP3→O3T`.
* elif not `ani.have_phosphate` → **5*END**: deletes P, OP1, OP2; O5' type OH1; adds `HO5'` with bond 0.84
  (n 1.02) and angle `C5'-O5'-HO5' 120`.
* `ani.have_ho3prime` → **3*END**: O3' type OH1, adds `HO3'`, torsion `hh C4'-C3'-O3'-HO3' 0 per3`.
* Without `ani`: `O3T` unexpected → p5*END; no P, OP1, OP2 expected → 5*END; `HO3*` unexpected → 3*END.

**What Reduce actually gets** (verified on 1ehz and 4fen). geostd A/C/G/U/AD/… dictionaries have **no**
OP3, HOP3, HO5' or HO3', and Reduce adds only dictionary H. So 5' residues that carry OP3 get p5*END (no HOP3 placed),
5' residues with P but no OP3 get no mod, 5'-OH residues (no P) get 5*END (**no HO5'** placed), and 3' ends get no mod (**no HO3'**). The `HO2'` hydroxyl on RNA
is placed (dictionary `hh2 C1'-C2'-O2'-HO2' 0 per2`).

**D residues**: PEPT-D (`chir CA N CB C positiv`), plus `DIL_chir_02_both` / `DTH_chir_02_both`.

---------------------------------------------------------------------------

## 5. Automatic linking and other non-polymer bonds

### 5.1 Overview

All of this runs in `construct_geometry_restraints_manager` (PI:5508), after the chain proxies, in this order:

disulfides (5.2) → geometry edits (none) → `process_nonbonded_for_links` (5.4) → grm build → enol-peptide fix
(PI:5944; sets C-N to 1.27 for enol peptides only; irrelevant here) → in `process.geometry_restraints_manager`:
MCL (5.3) → C-beta dihedrals.

**PDB `LINK`/`SSBOND` records and mmCIF `struct_conn` are never read.** Only `apply_cif_link` params do that, and
Reduce sets none. In any case the Reduce model is built from a bare hierarchy (`model_input=None`).

### 5.2 Disulfides (PI:3940-3990, 4701-4772, 5549-5730)

* Candidates: the SG atom of every residue whose dictionary `chem_comp.id == "CYS"`. DCY is interpreted as CYS
  and counts.
* **Exclusions**: drop any SG that is within `exclusion_distance_cutoff=3.0 Å` of an atom whose element (from
  `determine_chemical_element_simple`) is **not** one of `H D T S O P N C SE`, typically a metal.
* Pairs are SG pairs within `disulfide_distance_cutoff=3.0 Å` from a nonbonded search that respects symmetry,
  model and conformer exclusions (different altloc conformers never pair).
* The bond comes from `link_SS`: **2.031 Å**, esd 0.02, `origin_id = SS BOND (1)`. It goes into
  `bond_params_table` and `bond_asu_table`, so symmetry partners are allowed.
* Only if the symmetry operator is `x,y,z`, for each altloc combination of the CA/CB atoms (with blank fallback):
  * angles `CB-SG-SG'` and `SG-SG'-CB'` **104.2±2.1**;
  * dihedrals `ss CB-SG-SG'-CB' 93 alt[-86] esd10 per1`, `chi2_1 CA-CB-SG-SG' 79 alt[183,-73] esd20 per1`,
    `chi2_2 SG-SG'-CB'-CA' -79 alt[-183,73]`;
  * all with origin SS BOND, added with `add_if_not_duplicated`.
* **Nothing is deleted.** CYS `HG` keeps its dictionary restraints: `SG-HG 1.20` (n 1.20), `CB-SG-HG 109±5`,
  `chi2 CA-CB-SG-HG 180 esd15 per3`. HG is removed later by Reduce's `exclude_H_on_links` because SG takes part in
  an origin≠0 bond.

### 5.3 Metal Coordination Library (PI:6274-6282; MCL)

This runs because `mcl=True` and `link_metals in [Auto, True]`. It works on the nonbonded proxies of the built grm.

* **Zn tetrahedral** (`metal_coordination_library.get_metal_coordination_proxies`): for nonbonded pairs within
  3.0 Å where one atom is **named** `ZN`, the partner's name must appear in the Zn database (`SG`, `ND1`, `NE2`).
  At most one partner per residue.
  * Fewer than 4 partners: bonds 2.30 Å, σ 0.03.
  * Exactly 4 partners: ideals by (number of CYS, number of HIS), for example (4,0) `ZN-SG 2.330(0.029)`,
    `SG-ZN-SG 109.45(5.46)`; (3,1), (2,2), (1,3) see source. Angles are added too.
  * Origin **metal coordination (3)**.
* **Fe-S clusters** (`SF4`, `F3S`, `FES`): an Fe within 3.5 Å of an S or N of another residue, first contact per
  residue. Bonds from tables (for example SF4-CYS `FE-S 2.268(0.034)`) and angles. Origin metal coordination.
  Any symmetry contact aborts the whole step.
* Verified on 1xso: `ZN-ND1 HIS ×3` at 2.30. The CU ions are not linked.

### 5.4 `process_nonbonded_for_links` (LM:301-1171)

1. `max_bonded_cutoff = max(3.0, 1.9, 1.99, 1.99, 1.98, 2.2+0.5) = 3.0`. Nonbonded pair search over all atoms
   (simple + asu), iterated in `sorted_value_proxies_generator(by_value="delta")` order.
2. `skip_if_longer` (LS:118-168), called with `amino_acid=1.9, rna_dna=3.4, intra_residue=inter_residue_bond_cutoff=2.2, saccharide=1.99, metal=3.0, sulfur=2.5, other=2.0` (squared):
   * (aa,aa) 1.9; (aa,other) 2.2; (rna,rna) 3.4; (rna,metal) 3.0; (rna,other) 2.0; (ccp4rna,other) 2.0;
   * (aa,sacch) 1.99; (sacch,sacch) 1.99; (element,water) 3.0; (other,other) 2.0; (sulfur,sulfur) 2.5;
   * (aa,rna) 1.9; (rna,small) 2.0; (aa,small) 1.9·2.0; (small,other) 2.0·2.0; (metal,metal) 3.0;
   * (metal, other|water|aa) 3.0.
3. Per pair `(i, j, distance, sym_op)`:
   * Skip if i and j are already bonded (`bond_asu_table.contains(i, j, j_sym)` for j_sym 0 or 1).
   * Skip if **either atom is H or D**. H are never linked, so `exclude_hydrogens_from_bonding_decisions` has no
     further effect.
   * Skip exclude selections (none here).
   * Skip if either residue is SF4, F3S or FES; or if both are in {ZN, CYS}; or both in {ZN, HIS}. MCL handles
     those.
   * Skip if both atoms are in the same residue group.
   * Skip different non-blank altlocs.
   * Class tests: `link_ligands` is True, so `other` is allowed. `link_small_molecules` is False, so
     `common_small_molecule` (GOL, PO4, SO4) never links.
     * aa/d_aa/uncommon_aa pairs: `possible_cyclic_peptide` (C–N names, one atom in the first residue group of the
       chain and the other in the last, same chain id) → `use_only_bond_cutoff`.
     * With a symmetry op, skip amino acid–saccharide pairs.
   * Skip if `i` was already custom-bonded (Misc. or metal) by this routine to some atom `t` that is covalently
     bonded to `j` (shell 0), or the same with `i` and `j` swapped. This is the "bonded atoms can't link to the
     same atom" rule (LM:601-611).
   * **`is_atom_pair_linked`** (LU:396-544):
     * false if either element is in `H D F CL BR I AT HE NE AR KR XE`, or the pair is O–O;
     * `class = adjust_class(important_only)`: `common_element` → `metal` or `ion`, `other` + metal element →
       `metal`; both sulfur-class (S in aa or other) → `sulfur`;
     * lookups in `skip_if_both` return false: (water,water), (aa,water), (sacch,water), (water,other);
     * `limit = skip_if_longer[sorted classes]`. `+0.5²` is added when either element is not in
       `["LI","BE","B","C","N","O","F"]`. (Correction, verified by tracing: by the time `model.process()` links,
       iotbx has stripped the elements to `"C"`, so the buffer applies only to second-row and heavier
       elements; an earlier version of this note said it always applied.)
     * then the rules in order: aa+rna, rna+small, aa+small, other+small → true; sulfur+sulfur → true;
       saccharide in lookup → `d² ≤ 1.99²` (metal: 3.0²); `common_element` → true if a metal (else `Sorry`);
       `metal` → false for HIS CE1, CD2 and CB, otherwise returns `link_metals` (Auto is truthy); aa+aa needs an
       N–C element pair; `d_amino_acid` O → false; any `other` → true; else false.
   * If both partners are not `common_element`, `check_valence` (LU:823): an **O** with exactly 2 atoms within
     1.8 Å in its own atom group (H at bogus positions count; same model, chain, resseq and resname) is
     rejected.
   * `class_key = sorted(important_only pair)`:
     * **`if link_metals != True and "metal" in class_key: continue`**. With Auto, links with a metal-class
       residue (single-atom metal ions, common_element metals) are **not made**. A metal inside a multi-atom
       `other` residue (for example HEM FE) has `important_only=='other'`, so it **is** linked as a `Misc. bond`.
       This follows from the code; no HEM case was tested.
     * **`link_residues` is True → `if len(key)>2: continue`**: any candidate with a symmetry operator is
       skipped. No symmetry auto-links are made at all.
     * Carbohydrate rules: `link_carbohydrates` True. aa–rna is excluded because `link_amino_acid_rna_dna` is
       False.
   * Per-atom and per-residue-pair limits:
     * `maximum_per_atom_links` = 1 for classes saccharide, rna_dna, amino_acid and other. A second link for the
       same atom is allowed only by `_may_link_again`: distance ≤ 2.0, an element with `neutral_valence`
       (N:3), and fewer existing heavy neighbours + links than the valence.
     * `maximum_inter_residue_links[class_key]` defaults to 1: (element,other) 8, (metal,other) 6, (aa,metal) 2,
       (sacch,metal) 3, (rna,metal) 2, (rna,rna) 5.
   * Skip if any bond already exists between the two residues' atom groups.
   * **Link resolution**:
     1. `is_atom_group_pair_linked`: a link id `"RES1-RES2"` or `"RES2-RES1"` exists (for example `NAG-ASN`) →
        `_apply_link_using_proxies` with `origin_id = linking_class['link_<KEY>']`. If the key is not in the
        origin table (user links), **skip**.
     2. else `is_atom_pair_linked_tuple`: two peptide-class residues with N–C names → **`TRANS`** link (origin
        `link_TRANS` 79); or `RES_ATOM-ANY_ATOM` keys (for example `ASP_CG-ANY_N`).
     3. else `process_atom_groups_for_linking_single_link`: the key is `RES:ATOM-RES:ATOM`, or a glyco key
        (`NAG-ASN`, `ALPHA1-4`, `BETA1-4` … from the anomeric-carbon chirality).
        * Glyco links → `glyco_utils.apply_glyco_link_using_proxies_and_atoms`: bond **1.439 Å σ0.02**, angles,
          chirality. Origin `link_<ALPHA|BETA…>` or `glycosidic custom (5)`.
        * A known link id → `_apply_link_using_proxies`.
        * Otherwise `check_for_peptide_links` (aa ↔ other residue with C, N, O) → TRANS or SS.
        * Otherwise a **custom bond**: origin **`Misc. bond (10)`**, or **`metal coordination (3)`** if an
          `important_only` class is `metal`. Ideal from `bondlength_defaults.run(atom1, atom2)` (element-pair or
          metal tables; fallback 2.3), weight `1/σ²` (σ default 0.02), slack from the table.
   * `_apply_link_using_proxies` (LM:23-228):
     * Atoms are found **by exact atom name** in the two atom groups, so link rows that mention H work only if
       that H exists by that name. If the two groups are more than 3 Å apart for the first bond they are swapped.
     * Bonds use **`value_dist` only (never neutron)**.
     * Angles, dihedrals (periodicity from the link), chiralities (`volume_ideal ±2.4`, weight 25) and planes
       (≥4 atoms) are added with `add_if_not_duplicated` and **the same origin_id**.
     * **Link mods (for example `DEL-HD22`, `DEL-O1`, `DEL-HG`) are not applied.** No atom is deleted. The H stays
       with its dictionary restraints until Reduce's `exclude_H_on_links`.

### 5.5 Verified examples

* **3vyk**:
  * glycan `O3 MMA – C1 MAN` `link_ALPHA1-3 (19)`, `O4 MMA – C1 NAG` `link_BETA1-4 (27)`, `link_ALPHA1-6 (21)`,
    `link_BETA1-2 (25)`, all 1.439 Å. The glyco angle list includes `O(i)-C1-H1`.
  * CA ion LINK records are ignored, and no Ca–O link is made.
  * The Reduce step then removes `HO3/HO4/HO6` of MMA and `HO2` of MAN.
* **6oge**: `ND2 ASN – C1 NAG` `link_NAG-ASN (59)` at 1.439. Reduce removes `HD22`.

### 5.6 origin_id table and which bonds have origin_id ≠ 0

The table is `cctbx.geometry_restraints.linking_class` (`auto_linking_types.py`); the full table is in
`spec_work/data/origin_ids.tsv`:

```
0 covalent geometry   1 SS BOND   2 hydrogen bonds   3 metal coordination   4 edits
5 glycosidic custom   6 basepair stacking   7 basepair parallelity   8 side-chain parallelity
9 basepair planarity  10 Misc. bond   11 User supplied cif_link   12 solvent network
13 reference hydrogen bonds  14 C-beta  15 chi angles  16.. link_<KEY> in sorted key order:
16 link_ACE_C-N 17 link_AHT-ALA 18 link_ALPHA1-2 19 link_ALPHA1-3 20 link_ALPHA1-4 21 link_ALPHA1-6
22 link_ALPHA2-3 23 link_ALPHA2-6 24 link_ASP_CG-ANY_N 25 link_BETA1-2 26 link_BETA1-3 27 link_BETA1-4
28 link_BETA1-6 29 link_BETA2-3 … 35 link_CIS … 43 link_FE-CYS … 53 link_MAN-ASN … 59 link_NAG-ASN
… 62 link_NH2_CTERM 63 link_NMCIS 65 link_NMTRANS 66 link_PCIS 70 link_PTRANS 72 link_SS 79 link_TRANS
83 link_ZN-CYS 84 link_gap 85 link_p 86 link_symmetry
```

| source of bond | origin_id |
|---|---|
| residue dictionary bonds (including mods), polymer chain links (TRANS, PTRANS, NMTRANS, rna2p, rna3p, ACE_C-N, NH2_CTERM…) built in `build_chain_proxies` | **0** |
| disulfide | 1 (SS BOND) |
| MCL Zn and Fe-S | 3 |
| automatic named link | `link_<KEY>` (16-86), for example TRANS for cyclic or "other"-residue peptide bonds = **79** |
| glycosidic | `link_ALPHA*/BETA*` or 5 |
| automatic custom bond | 10 (Misc. bond) or 3 (metal) |
| apply_cif_link | 11 (not used) |

Note: dihedral proxies from `c_beta` carry origin 14 but involve no H. `exclude_H_on_links` only looks at bond
proxies (simple + asu) with origin ≠ 0.

---------------------------------------------------------------------------

## 6. How proxies are generated

### 6.1 General rules (PI:2088-2630)

For each restraint row of a residue (`m_j = None`) or a link (atoms carry comp ids 1 or 2):
1. Look up the atom ids in `monomer_atom_dict`. On failure retry:
   * bonds, dihedrals, chiralities, planes: all ids with `'`→`*`, then `*`→`'`;
   * angles: `convert_v3_to_v2` (`'→*, OP1→O1P, OP2→O2P`), then `convert_v2_to_v3`.

   The row is **mutated in place**. If still not found → `corrupt_monomer_library_definitions`, skip.
2. Map ids to expected atoms. If any is missing → skip. The counter is `unresolved_hydrogen` if any of the
   restraint's atoms is a dictionary H (`type_symbol=="H"`), else `unresolved_non_hydrogen`.

   **For planes**, a missing atom just drops that atom (counter only). The plane is kept if at least 4 atoms are
   resolved, otherwise discarded.
3. No ideal value, or esd None/0 → `undefined`, skip. Exception: dihedrals with esd None/0 become `const`
   proxies (§6.3).
4. Discard angles, dihedrals and chiralities that involve atoms on special positions (`special_position_dict`).
   For planes, only the special atom is dropped.
5. Skip anything spanning a broken link bond pair (§3.1).
6. Register through the proxy registry (§6.6).

### 6.2 Bonds (PI:2088-2166)

* `distance_ideal = value_dist`, or `value_dist_neutron` if `use_neutron_distances` and that field is non-empty.
  This applies to residue bonds only; link calls do not pass the flag (PI:2986).
* `weight = 1/value_dist_esd²`. The **X-ray esd is always used**.
* In conformers after the first, bonds whose atoms are all blank-altloc are skipped
  (`already_assigned_to_first_conformer`).
* geostd `_chem_comp_bond.value_dist_neutron` column: SER `N-H 0.86/1.02`, `CA-HA 0.97/1.09`, `OG-HG 0.84/0.98`.
  Mod bonds carry `new_value_dist_neutron` (NH3 `N-H 0.89/1.04`).

### 6.3 Dihedrals (PI:2242-2417)

* **Every** `_chem_comp_tor` row is used. There is **no filtering by id** (`chi*`, `hh*`, `Var_*`, `CONST_*`,
  `nu*` are all treated the same), and no "esd 180 means ignore" rule. Only `psi`, `phi` and `omega` of the
  peptide chain links get special handling (§3.4).
* `value_angle None` → skip.
* `value_angle_esd None or 0` → **const dihedral proxy** (weight 0, periodicity 0) in the separate
  `const_dihedral_proxies` list. Examples: HIS ring `CB-CG-ND1-CE1 180` and the nucleic-base `CONST_xx` rows.
* Otherwise `dihedral_proxy(angle_ideal=value_angle, weight=1/esd², periodicity=period, alt_angle_ideals=…)`.
  `alt_value_angle` is a comma-separated string; unparseable → None.
* Geostd `Var_01 C-CA-N-H 170 esd 180 period 72` (SER, THR, CYS, ASN, HIS…) becomes a normal proxy with weight
  3.1e-5 and periodicity 72. It is **resolvable only when the residue keeps its amide `H`**: not with NH3, not for
  PRO.
* Typical H torsions: SER `hh1 CA-CB-OG-HG 180 esd30 per3`; THR `hh1 CA-CB-OG1-HG1 180/30/3`, `hh2 CA-CB-CG2-HG23 60/30/3`; LYS `CD-CE-NZ-HZ3 60/30/3`; TYR `CE1-CZ-OH-HH 180/30/2`; ASN `CB-CG-ND2-HD22 180/30/2`; ARG `NE-CZ-NH1-HH12`, `NE-CZ-NH2-HH22 180/20/2`; CYS `chi2 CA-CB-SG-HG 180/15/3`; HIS ring H `… 180 esd5 per0`; PRO `CA-N-CD-HD2/HD3`, `CD-N-CA-HA` per3; RNA `hh2 C1'-C2'-O2'-HO2' 0/40/2`, G `hh1 N3-C2-N2-H22 0/40/1`.

### 6.4 Angles, chiralities, planes, parallelity

* Angle: `weight = 1/esd²`.
* **Chirality** (PI:2419-2515):
  * `volume_sign[:4].lower() ∈ {posi, nega, both}`, else skip.
  * `volume_ideal` from the residue: `uctbx.unit_cell([b1,b2,b3, ang(2,c,3), ang(1,c,3), ang(1,c,2)]).volume()`,
    negated for "neg". It uses the dictionary bonds centre–atom (x-ray `value_dist`) and angles. If any of those
    are missing, the result is None and the row is skipped.
  * For links: `lib_link.get_chir_volume_ideal(m_i, m_j)`.
  * `both_signs = (sign=="both")`, `weight = 1/0.2²`.
* Plane: each atom's weight is `1/dist_esd²`. A plane is kept only with ≥4 atoms.
* No parallelity proxies are generated (`add_parallelity_proxies` is never called by default).

### 6.5 C-beta dihedrals

These come from `mmtbx.geometry_restraints.c_beta`, added after grm construction. Each has two proxies
(`N-C-CA-CB`, `C-N-CA-CB`), harmonic, origin `C-beta`. No H is involved, but they do appear in
`dihedral_proxies`, which riding loops over.

### 6.6 Registries and canonical ordering (GR:120-400, `GR/dihedral.h:117-150`)

* **Bond**: `sort_i_seqs` (i<j). The table is keyed `[i][j]`.
* **Angle**: `sort_i_seqs` swaps the ends so that `i0 < i2`. Keyed by `(i1, (i0, i2))`.
* **Dihedral and const dihedral**: `sort_i_seqs` swaps the ends if `i0 > i3` **and independently** swaps the
  middle if `i1 > i2`. **Each swap negates `angle_ideal` and every `alt_angle_ideal`.** This is geometrically
  exact, because a torsion depends only on the projections of the end atoms. Keyed by `(i0, (i1, i2, i3))`.
  Example: `Var_01 (C,CA,N,H)` with `CA > N` is stored as `(C, N, CA, H)` with ideal `-170`.
* **Chirality**: the centre stays first. The other three i_seqs are sorted ascending by pairwise swaps, and each
  swap negates `volume_ideal` unless `both_signs`. Keyed by `(i0, (i1, i2, i3))`.
* **Planarity**: `sort_i_seqs` sorts the atoms by i_seq, with the weights carried along. Keyed by
  `(i0, tuple(i1…))`. So the "first H of a plane" used by riding (§8.2 step 4) is the H with the lowest i_seq.
* `process`:
  * new → append;
  * existing with equal ideal (within 1e-6, using `angle_delta_deg` with periodicity for angles and dihedrals)
    and equal weight → `is_new = False`. If that happens and any atom is not blank-altloc, raise `Sorry`
    ("Duplicate … restraints").
  * different values → `_handle_conflict`: if the two sources have the same number of expected atoms → `Sorry`
    (one exception: names exactly `[' H2 ',' N  ']`, print only); else keep the proxy from the source with
    **more** expected atoms.
* Proxy order in the final grm:
  * angles: for each residue, the link angles to the previous residue, then the residue's own angles;
    then SS angles, then auto-link angles, then MCL angles;
  * dihedrals: the same pattern, then auto-link, then C-beta.

  Riding depends on this order (§8.2).

### 6.7 Alternate conformations

Restraints are built separately for each conformer (blank + altloc X). Effects:
* Atoms that only exist in altloc X bond only to blank atoms and X atoms. Cross-altloc proxies never arise.
* Restraints on blank atoms only are built once for bonds (first conformer) and deduplicated for the other types.
* Different mods or links per conformer that change a shared proxy go through the conflict rule above.
* Nonbonded and auto-linking searches use `conformer_indices`: blank = 0; altlocs are numbered 1, 2, … in sorted
  order whether or not a blank exists (PI:3721-3732). Different non-zero indices never interact.
* Reduce puts all H of a partially-altloc residue into the altloc groups (§0.1). Each conformer then has its own H
  set, bonded to the blank or altloc heavy atoms.

### 6.8 Missing atoms

Covered in §6.1. A missing heavy atom silently drops every restraint that needs it. H atoms still get their bond,
angle and so on as long as their own partners exist. Riding later refuses H whose parent lost a dictionary heavy
neighbour (§8.2 step 7).

### 6.9 Atom order (`hierarchy.sort_atoms_in_place`; HIER.cpp:1335-1492)

* Atoms are sorted **within each atom group** with `std::sort`, which is not stable, using
  `sorting_general(a1, a2, table)`.
* Names are compared after `strip()` and `'*'→"'"`. "H" means the stripped name **starts with `H`**; this is not
  the element.
* `pos` = index of the upper-cased name in the table, or `table.size` if absent.
  * If `pos1 == pos2`: non-H before H, then raw 4-character names compared lexicographically.
  * Otherwise: non-H before H, then by `pos`.
* Table choice by class of the residue name:
  * `common_rna_dna` or `modified_rna_dna`: `big_na` if the group has an `N9` atom, else `small_na`;
  * everything else, **including ligands**: the amino-acid table (`N CA C O CB NB OB SB CB1 … OXT H H1 1H H2 2H H3 3H HA HA1 … HH33`).
* Consequences: ligand heavy atoms not in the table are sorted alphabetically by padded name, after any names
  that are in the table, and H go last. 1crn THR 1 ends up as `N CA C O CB OG1 CG2 H1 H2 H3 HA HB HG1 HG21 HG22 HG23`.
* This runs before `process` (RH) and again inside `model.process` (sort_atoms=True), then `reset_i_seq`. The
  resulting i_seq order drives proxy canonicalisation and riding.

### 6.10 `flip_symmetric_amino_acids` (HIER:2470-2601) changes the model

This runs at the start of `build_all_chain_proxies` (PI:3637) on the model's own hierarchy. For each residue group,
`flip_it` is decided from the first atom group that has the needed atoms:
* ARG `|dihedral(CD,NE,CZ,NH1)| > 90`;
* ASP `|dih(CA,CB,CG,OD1)| > 90`;
* GLU `|dih(CB,CG,CD,OE1)| > 90`;
* PHE and TYR `|dih(CA,CB,CG,CD1)| > 90`;
* VAL chirality `(CB,CA,CG1,CG2)`, LEU `(CG,CB,CD1,CD2)`: `|−2.5 − V| > 2`.

If true, **swap `xyz` and `b`** of the pairs (ARG NH1/NH2, HH11/HH21, HH12/HH22; ASP OD1/OD2; GLU OE1/OE2;
PHE/TYR CD1/CD2, CE1/CE2, HD1/HD2, HE1/HE2; VAL CG1/CG2 + HG1x/HG2x; LEU CD1/CD2 + HDxx; plus D variants) in
**all** atom groups of the residue group. A pair with only one member present means no flip. All H sit at the
same bogus point at this time, so only the heavy-atom coordinates change in practice.

### 6.11 Other

* The NCS shortcut (PI:3816-3870) is only used if the copies match to within 0.01 Å RMSD. It copies the master's
  proxies, which gives the same result. It can be ignored.
* `merge_atoms_at_end_to_residues` and `format_correction_for_H`: the latter re-pads H names shorter than 4
  characters that are left-justified, to `' %-3s'`, unless they start with a digit.
* Excessive bond lengths are allowed (`proceed_with_excessive_length_bonds`).

---------------------------------------------------------------------------

## 7. Energy types → ExtraAtomInfo

### 7.1 Model accessors (MODEL:2695-2735)

```
get_specific_vdw_radius(i_seq, vdw_radius_without_H=False):
    t = model._type_energies[i_seq]                 # from acp.type_energies (§1.10)
    e = model.get_ener_lib()                         # the ener_lib built by process():
                                                     # server.ener_lib(use_neutron_distances=flag)
    return e.lib_atom[t].vdwh_radius if vdw_radius_without_H else e.lib_atom[t].vdw_radius
get_specific_ion_radius(i_seq) -> e.lib_atom[t].ion_radius
get_specific_h_bond_type(i_seq) -> model._type_h_bonds[i_seq]   # 'A','B','D','N','H', or 'False'/'None'
```

* Reduce2 calls `get_specific_vdw_radius(i, False)`, which reads the **`vdw_radius`** column, not `vdwh_radius`.
* With `use_neutron_distances=True`, `ener_lib.convert_lib_atom` (SRV:777) **replaces `vdw_radius` by
  `vdw_radius_neutron` when present**. This applies only to H types: 1.22→1.17 for the C-H types and 1.05→1.00 for
  polar H. Ion radii are unchanged.
* After `model.select`, `_type_energies` and `_type_h_bonds` are selected too (MODEL:4137).
* **No energy type**: `t` is `''` (unknown residue), `'False'` (unexpected atom) or `'None'` (dictionary without
  `type_energy`, for example Reduce `auto_*`). Then `e.lib_atom[t]` raises **KeyError**.
  `mmtbx.probe.Helpers.getExtraAtomInfo` (Helpers.py:359-507) catches it. It sets the radius to
  `ener_lib.lib_atom[element].vdw_radius` if an entry named after the element exists (`C`, `N`, `O`, `S`, `H`, …),
  else leaves it at 0. The H-bond type then stays at its default (not donor, not acceptor). An H-bond type
  `'False'`, `'None'` or `''` gives a warning and no flags.
* Probe-side overrides applied after the lookup are **not** part of interpretation, listed here only so they are
  not confused with it:
  * `element_is_ion()` (element-based: metals and F/CL/BR/I) → `ion_radius`;
  * carbonyl `C`, plus ASP/ASN `CG` and GLU/GLN `CD` → 1.65;
  * polar H → 1.05;
  * aromatic acceptor table;
  * HET N with no H → acceptor.

### 7.2 Table: types used by the 20 amino acids, MSE, water, nucleic acids, and the auto mods

Source: `chem_data/geostd/ener_lib.cif`. Generated by `spec_work/energy_table.py` → `dumps/energy_table.txt`.
`hb`: A acceptor, B both, D donor, H hydrogen, N none.

| type | el | hb | vdw_radius | vdw (neutron) | vdwh_radius | ion_radius | used by |
|---|---|---|---|---|---|---|---|
| C | C | N | 1.70 | 1.70 | 1.75 | – | backbone C, ASN CG, GLN/GLU CD, ASP CG, ARG CZ |
| CH1 | C | N | 1.70 | 1.70 | 1.95 | – | CA, CB (VAL ILE THR), LEU CG, sugar C1'–C4' |
| CH2 | C | N | 1.70 | 1.70 | 1.92 | – | CH2 groups, GLY CA, C5', DNA C2' |
| CH3 | C | N | 1.70 | 1.70 | 1.94 | – | methyls, TD C5M |
| CR15 | C | N | 1.75 | 1.75 | 1.74 | – | HIS CD2/CE1, TRP CD1, purine C8 |
| CR16 | C | N | 1.75 | 1.75 | 1.82 | – | aromatic CH (PHE, TYR, TRP), pyrimidine C5/C6, A C2 |
| CR5 | C | N | 1.75 | 1.75 | 1.74 | – | HIS CG, TRP CG |
| CR56 | C | N | 1.75 | 1.75 | 1.74 | – | TRP CD2/CE2, purine C4/C5 |
| CR6 | C | N | 1.75 | 1.75 | 1.74 | – | PHE/TYR CG, TYR CZ, base C2/C4/C6 |
| N | N | N | 1.55 | 1.55 | 1.60 | 1.32 | PRO N |
| NH1 | N | D | 1.55 | 1.55 | 1.60 | 1.32 | backbone N, NH1 mods |
| NH2 | N | D | 1.55 | 1.55 | 1.60 | 1.32 | ASN ND2, GLN NE2, base NH2, NH2 mods |
| NT3 | N | D | 1.55 | 1.55 | 1.60 | 1.32 | LYS NZ, NH3 mod (N-terminal N) |
| NC1 / NC2 | N | D | 1.55 | 1.55 | 1.60 | 1.32 | ARG NE / NH1, NH2 |
| NR15 | N | D | 1.55 | 1.55 | 1.60 | 1.32 | HIS ND1, NE2; TRP NE1 |
| NR16 | N | D | 1.55 | 1.55 | 1.60 | 1.32 | G N1, U/T N3 |
| NR5 | N | N | 1.55 | 1.55 | 1.60 | 1.32 | purine N9 |
| NR6 | N | N | 1.55 | 1.55 | 1.60 | 1.32 | pyrimidine N1 |
| NRD5 / NRD6 | N | A | 1.55 | 1.55 | 1.60 | 1.32 | purine N7 / N1, N3 |
| O | O | A | 1.40 | 1.40 | 1.52 | 1.28 | backbone O, ASN OD1, GLN OE1, base O |
| OC | O | A | 1.40 | 1.40 | 1.52 | 1.28 | ASP OD*, GLU OE*, COO/COOH O, OXT |
| OH1 | O | B | 1.40 | 1.40 | 1.52 | 1.28 | SER OG, THR OG1, TYR OH, O2', 3*END/5*END O |
| OH2 | O | B | 1.40 | 1.40 | 1.52 | 1.28 | water O |
| O2 | O | A | 1.40 | 1.40 | 1.52 | 1.28 | O4' |
| OC2 | O | A | 1.40 | 1.40 | 1.52 | 1.28 | O5', O3' |
| OP | O | A | 1.40 | 1.40 | 1.52 | 1.28 | OP1, OP2, O3T |
| P | P | N | 1.80 | 1.80 | 1.88 | 0.59 | P |
| S | S | A | 1.80 | 1.80 | 1.88 | 0.40 | CYS SG, MET SD |
| SE | SE | N | 1.90 | 1.90 | – | 0.42 | MSE SE |
| H | H | N | 1.22 | 1.17 | – | – | **CYS HG** (not polar in ener_lib), CF-COH/CBH |
| HCH1/HCH2/HCH3 | H | N | 1.22 | 1.17 | – | – | all aliphatic H |
| HCR5 / HCR6 | H | N | 1.05 | 1.00 | – | – | aromatic CH H |
| HNH1 | H | H | 1.05 | 1.00 | – | – | amide H, NH1 mods |
| HNH2 | H | H | 1.05 | 1.00 | – | – | NH2 H |
| HNT3 | H | H | 1.05 | 1.00 | – | – | LYS HZ*, N-terminal H1-3 |
| HNC1 / HNC2 | H | H | 1.05 | 1.00 | – | – | ARG HE / HH** |
| HNR5 / HNR6 | H | H | 1.05 | 1.00 | – | – | HIS HD1/HE2, TRP HE1 / G H1, U H3 |
| HOH1 | H | H | 1.05 | 1.00 | – | – | SER HG, THR HG1, TYR HH, HO2', COOH HXT, ASP HD2 / GLU HE2 |
| HOH2 | H | H | 1.05 | 1.00 | – | – | water H |

Per-residue atom → type maps for all of these residues are in `dumps/energy_table.txt`. Common ions (geostd
values; `get_specific_ion_radius` in parentheses): MG 1.39 (0.65), ZN 1.45 (0.71), NA 1.69 (0.95), K 2.07
(1.33), CA 1.73 (0.99), MN 1.54 (0.80), FE 1.48 (0.74), CU 1.46 (0.72), CD 1.65 (0.91), HG 1.74 (1.00),
NI 1.40 (0.66), CO 1.44 (0.70), CL 1.75 (1.67, hb A), BR 1.85 (0.73), I 1.98 (0.56), F 1.47 (1.19, hb B). All
others are in `data/ener_lib_lib_atom.txt`. Ion types are hb N except CL (A) and F (B).

---------------------------------------------------------------------------

## 8. Riding-H manager

### 8.1 Setup

`model.setup_riding_h_manager(use_ideal_dihedral=True)` (MODEL:3011) creates
`riding.manager(pdb_hierarchy, grm, use_ideal_dihedral=True, use_ideal_bonds_angles=True, ignore_h_with_dof=False, mon_lib_srv)`
(RIDING:11-50). That builds `connectivity.determine_connectivity` (CONN:40) from **the grm proxies** and the
**current coordinates**, with H still at their bogus or construct_xyz positions. Then
`parameterization.manager` (PARAM:9). Since `idealize=True`, it immediately calls `idealize_h_riding()`
(MODEL:3025), which in turn calls `apply_new_H_positions` (HH:126).

H atoms are recognised with `atom.element_is_hydrogen()` (H or D).

### 8.2 Connectivity (CONN:47-106)

1. **First neighbours** (CONN:108): loop over `grm.get_all_bond_proxies(sites)[0]` (simple proxies only;
   symmetry/asu bonds are ignored). For each bond with exactly one H, `a0 = {iseq: parent, dist_ideal:
   proxy.distance_ideal}`, using the first bond seen. `fsc0` = full simple connectivity of shell 0. An H with more
   than one neighbour goes into `double_H`.
2. `count_H`: any H without an entry is "slipped" and later gets `number_non_h_neighbors=0`, so it is not
   parameterized.
3. **Second neighbours** (CONN:138): for each angle proxy in grm order (`get_all_angle_proxies`) with the H at
   one end and the parent `i_seqs[1]` in the parent set, record `(ih, parent, other)` → `angle_ideal`. An
   inconsistent parent raises `Sorry` unless the H is in `double_H`.
4. `process_plane_proxies`: for each planarity proxy containing H, `plane_h[first_H] = other_H`. **Only the first
   H of each plane is "in plane".**
5. `process_second_neighbors`:
   * Blank-altloc second neighbours are kept. Altloc ones are reduced by `process_alternate_neighbors`: keep
     singletons; for same-named atoms keep the highest occupancy, and only one altloc.
   * H second neighbours become `h1`, `h2` (in order); the others become `a1`, `a2`, `a3` in **Python `set`
     iteration order of the i_seqs**. The final positions do not depend on that order (§8.5).
   * `number_h_neighbors` and `number_non_h_neighbors` are set.
   * For H with exactly one non-H second neighbour, record `a0a1_dict[parent] += [a1]`.
6. **Third neighbours from dihedral proxies** (CONN:373-441):
   * First `const_dihedral_proxies`: if H is at `i1` → `b1 = (i4, angle_ideal)`; if at `i4` →
     `b1 = (i1, −angle_ideal)`.
   * Then `dihedral_proxies` (residue, link, SS, auto-link, C-beta; last write wins). If H is at `i1` or `i4`:
     `model = dihedral(i1..i4)` with the **current (bogus) coordinates**, and
     `b1.dihedral_ideal = model + angle_delta_deg(model, angle_ideal, periodicity)`. If the dihedral is undefined
     (collinear) the loop returns early for this proxy.
   * The stored proxy order is used **as stored** (§6.6), and riding later interprets the value as
     `dihedral(H, a0, a1, b1)`. That is exact only when the stored order is `(b1, a1, a0, H)` or
     `(H, a0, a1, b1)`. For `Var_01` stored as `(C, N, CA, H)` it is not, which gives the mirror quirk in §8.4.
   * `assign_b1_for_H_atom_groups`: the sibling H (h1, h2) get `b1 = {iseq: i_third}` with **no** ideal, unless
     their own proxy is seen later.

   `angle_delta_deg(a1, a2, p)` (`GR/utils.h:42`):
   ```
   half = 180/max(1,|p|); d = fmod(a2-a1, 2*half); if d < -half: d += 2*half elif d > half: d -= 2*half
   ```
7. **Angles at the parent and fallback third neighbours** (CONN:210-307):
   * `parent_angles[a0][(x,z)]` = all non-H X-a0-Z angles.
   * For H with exactly one non-H neighbour `a1`, the angle proxies centred on `a1` with one end the parent give
     the raw third neighbours, in **angle-proxy order**.
   * `assign_a0_angles`: if ≥2 non-H neighbours, `angle_a1a0a2` must exist; with 3, also `angle_a2a0a3` and
     `angle_a3a0a1`. Otherwise the H is reset to unparameterized.
   * For 1-non-H-neighbour H **without** a `b1` iseq:
     * no non-H third neighbour → unparameterized;
     * the parent has fewer heavy neighbours in the model than in its dictionary (`_parent_lost_heavy_neighbor`,
       deduplicated by name and altloc-aware) → unparameterized;
     * else `b1 = first third neighbour` and:
       * if `number_h_neighbors == 2` (NH3- or CH3-like with no torsion; for example N-terminal H1-H3):
         `check_for_plane_proxy`, which does nothing here because h2 exists. Then
         `dihedral_ideal = +60° if ih < h1.iseq else −60°`;
       * else (solo H, or one H sibling): `dihedral_ideal = 180°`.
8. `add_slipped`. Set `is_in_plane = ih in plane_h`.

### 8.3 Parameterization (PARAM:27-460) and position formulas (HH:52-121)

`use_ideal_bonds_angles=True`, so `disth = a0.dist_ideal` (the bond proxy ideal, neutron or x-ray) and all angles
are the proxy ideals. Loop over H in **i_seq order**, skipping H already parameterized.

Notation used below. Unit vectors are taken from the **current** heavy-atom coordinates:
`u10 = unit(r1−r0)`, `u20 = unit(r2−r0)`, `u30 = unit(r3−r0)`.

| condition | htype | coefficients (PARAM) | position (HH `compute_h_position`) |
|---|---|---|---|
| 2 non-H neighbours (`process_2_neighbors`) | see below | `α0 = ∠a1a0a2`, `α1 = ∠H-a0-a1`, `α2 = ∠H-a0-a2` (ideal); `c_i = cos α_i`; `a = (c1−c0c2)/(1−c0²)`, `b = (c2−c0c1)/(1−c0²)`; `sumang = α0+α1+α2` | |
| … `abs(sumang−2π) < 0.05` | `flat_2neigbs` | h = 0 | `r0 + d·unit(a·u10 + b·u20)` |
| … not flat, 1 H sibling | `2tetra` (both H) | `h = ½·∠(H-a0-H')` in radians, negated if `(u10×u20)·uhH < 0`; the **sibling gets −h** | `d0 = unit(a·u10 + b·u20)`, `v0 = unit(u10×u20)`; `r0 + d(cos h·d0 + sin h·v0)` |
| … not flat, no H sibling | `2neigbs` (or `flat_2neigbs` if `is_in_plane`) | `h = sqrt(1−c1²−c2²−c0²+2c0c1c2)/sin α0`, sign as above | `r0 + d·unit(a·u10 + b·u20 + h·v0)` |
| … `sumang > 2π+0.05` and the root < 0 | none | not parameterized (`unk_ideal_list`) | |
| 3 non-H, 0 H (`process_3_neighbors`) | `3neigbs` | Cramer: `D = [[1,w12,w13],[w12,1,w23],[w13,w23,1]]`, `w_ij = cos(ideal ∠a_i-a0-a_j)`, `c_i = cos(ideal ∠H-a0-a_i)`; `a = det(Dx)/det D` and so on (x, y, z columns replaced by c) | `r0 + d·unit(a·u10 + b·u20 + h·u30)` |
| 1 non-H, 0 or 2 H (`process_1_neighbor`) | `alg1b` (0 H) / `prop` (2 H) | `α = ideal ∠H-a0-a1`; `φ = b1.dihedral_ideal` (radians); falls back to the model dihedral if no ideal. For `prop`, the parameterization is built from whichever of {ih, h1, h2} has `dihedral_ideal` (ih first); that H gets `n=0` and the two siblings `n=1,2` in `check_propeller_order` order | `u1 = unit(r0−r1)`, `u2 = unit((rb1−r1) − ((rb1−r1)·u1)u1)`, `u3 = u1×u2`, `φ' = φ + n·2π/3`; `r0 + d(sin α(cos φ'·u2 + sin φ'·u3) − cos α·u1)`. This gives `dihedral(H, a0, a1, b1) = φ'` exactly (verified numerically) |
| 1 non-H, 1 H (`process_1_neighbor_type_arg`) | `alg1a` (both H) | needs `dihedral_ideal` on ih or h1; if `107 < ∠(H-a0-H') < 111` → not parameterized; the H with the ideal gets `φ`, the other `φ+π` | same formula as alg1b, n=0 |
| anything else | none | | |

* `check_propeller_order(a0, a1, ih, h1, h2)`: keep `(h1, h2)` if `((rH−r0)×(rh2−r0))·(r1−r0) ≥ 0`, else swap.
  When all three H sit on the same bogus point the cross product is 0, so the order is kept: `h1, h2` are in
  second-neighbour (angle-proxy) order.
* `check_if_atoms_superposed`: any pair H/a0 or a_i/a0 closer than 0.001 Å raises `Sorry`. Every H must therefore
  be ≥0.001 Å from its parent at setup time; the bogus +0.5 offset exists for this reason.
* `apply_new_H_positions` (HH:126): runs `compute_h_position` for each parameterized H in i_seq order and writes
  the result immediately. Every formula uses only a0, a1, a2/b1 and a3. These are non-H, except that `b1` from a
  dihedral proxy could in principle be an H, in which case the update order matters.
* After riding, `idealize_h_riding` re-idealizes waters with 2 H: O–H 0.85 (0.98 for neutron/electron scattering
  tables), H-O-H 103.91°, keeping the plane and bisector. Not relevant with `exclude_water=True`.

### 8.4 Where the bogus H position leaks into the result

Reproduce these exactly:
1. **Periodic dihedral choice.** `b1.dihedral_ideal` is the periodic image of the dictionary value closest to the
   bogus-position dihedral.
   * 1crn THR 1 `HG1` (dict 180, per 3): model 74.25° → **60°**.
   * `HG23` (dict 60, per 3) → **−60°**.
   * SER 6 `HG` → −60°; TYR `HH` (per 2) → 0° or 180°; ASN `HD22` (per 2) → 0°.
2. **2tetra side.** The first H of a CH2 pair in i_seq order goes to the side of the a1-a0-a2 plane where the
   bogus point lies; its sibling goes to the other side. RH later renames the CH2 H stereochemically
   (`name_prochiral_h`, outside this spec).
3. **Var_01 quirk at chain starts without a propeller** (first residue not numbered 1, or after a break), where
   the amide `H` has only CA as heavy neighbour:
   * The proxy is stored as `(C, N, CA, H)` with −170, per 72.
   * `dihedral_ideal = model(C,N,CA,Hbogus)` rounded to the 5° grid of −170 (half period 2.5°).
   * It is then used as `dihedral(H, N, CA, C)`. The result is the **mirror** of the bogus H dihedral.
   * Verified on 1crn without residue 1: bogus `dih(H,N,CA,C) = +3.18°` → placed at `−5.0°`.
4. **N-terminal NH3 with no torsion**: b1 = the first third neighbour from angle-proxy order, which is **CB** for
   THR 1 (the THR dictionary lists `N-CA-CB` before `N-CA-C`; for GLY the first one is `N-CA-C`, so C). `H1` gets φ = +60°; `H2` and `H3`
   follow at +180° and +300° (verified: `dih(H1,N,CA,CB) = 60.0°`, H2 180.0°, H3 −60.0°).

### 8.5 Order independence

`3neigbs`, `2tetra`, `2neigbs` and `flat_2neigbs` give the same position for any permutation of a1/a2/a3, because
the h-sign test and `v0` flip together. So the Python `set` ordering of a1–a3 does not matter. What does matter:
* angle-proxy order: third-neighbour fallback, and h1/h2, which decides which name gets n=1 or n=2;
* dihedral-proxy order: last write wins;
* the bogus coordinates.

---------------------------------------------------------------------------

## 9. Worked examples (verified; raw dumps in `spec_work/dumps/`)

### 9.1 1crn THR 1, N-terminal (`dumps/1crn_res1.txt`)

* Mods: NH3 (`Modifications used: {'NH3': 1}`). N type `NT3`/hb D; H1-3 `HNT3`/hb H, radius 1.05.
* Bonds: `N-H1/H2/H3 0.89` (neutron 1.04) σ0.02; `N-CA 1.491`.
* Angles: `CA-N-H1/H2/H3 109.47 ±3`; `H1-N-H2`, `H1-N-H3`, `H2-N-H3 109.47 ±3`. No torsion and no plane on H1-3.
  `Var_01` was deleted along with `H`.
* Riding: H1/H2/H3 are `prop`, `a0=N a1=CA b1=CB`, α=1.9106 (109.47°), φ=+60° (fallback), n=0,1,2, d=0.89.
  Final `dih(H1,N,CA,CB) = 60°`, H2 180°, H3 −60°.
* The other THR 1 H:
  * HA: `3neigbs (N,C,CB)`, d 0.97.
  * HB: `3neigbs (CA,OG1,CG2)`.
  * HG1: `alg1b a0=OG1 a1=CB b1=CA`, α 110°, φ 60°, d 0.84.
  * HG21-23: `prop a0=CG2 a1=CB b1=CA`, φ −60° (HG23 n=0, HG22 n=1, HG21 n=2).

### 9.2 Peptide H (1crn CYS 3 `H`)

* Bond `N-H 0.86` (n 1.02).
* Angles: residue `CA-N-H 114.0±3`; link (TRANS) `C(2)-N-H 124.3±3`; link `C(2)-N-CA 121.7`.
* Link plane2 `C(2) N CA H` σ0.02. Dihedral `Var_01` stored `(C,N,CA,H) −170 per72`, harmless here.
* Riding: `flat_2neigbs a0=N a1=C(2) a2=CA`. The sum is 360°, so `a = −1.0737`, `b = −0.9710`, d = 0.86.

### 9.3 1crn SER 6 `HG`

* `OG-HG 0.84` (n 0.98); `CB-OG-HG 110±3`; `hh1 CA-CB-OG-HG 180 esd30 per3`.
* OG type `OH1`/hb B; HG `HOH1`/hb H.
* Riding: `alg1b a0=OG a1=CB b1=CA`, φ = −60° (from the bogus point), d 0.84.

### 9.4 1crn CYS 3 SG in disulfide 3–40

* `SG-SG(40) 2.031 σ0.02 origin SS BOND`.
* Angles `CB-SG-SG' 104.2±2.1` (both sides).
* Dihedrals `CB-SG-SG'-CB' 93 alt[−86] per1`, `CA-CB-SG-SG' 79 alt[183,−73]`, and the mirror on 40.
* HG still has `SG-HG 1.20`, `CB-SG-HG 109±5`, `CA-CB-SG-HG 180 esd15 per3`, and is placed by riding
  (`alg1b`, φ 60°, d 1.2). Reduce then deletes it because SG has an origin ≠ 0 bond. All six CYS HG in 1crn are
  removed.

### 9.5 1ubq HIS 68 ring (`dumps/1ubq_his68.txt`)

* The dictionary has both HD1 and HE2 (HIS+). Types: ND1/NE2 `NR15`/hb D; HD1/HE2 `HNR5`/hb H 1.05;
  CD2/CE1 `CR15`; HD2/HE1 `HCR5`/hb N 1.05.
* Bonds `ND1-HD1 0.86`, `NE2-HE2 0.86`, `CD2-HD2 0.93`, `CE1-HE1 0.93`.
* Angles: each ring H has two angles with its ring neighbours, for example `CG-ND1-HD1` and `CE1-ND1-HD1` both
  125.35.
* Dihedrals `CD2-CG-ND1-HD1 180 esd5 per0` (one per ring H).
* Const dihedrals on ring heavy atoms.
* One 10-atom plane `CB CG ND1 CD2 CE1 NE2 HD1 HD2 HE1 HE2`.
* Riding: all ring H `flat_2neigbs`.
* With Zn (1xso), MCL adds `ND1-ZN 2.30 origin 3`; Reduce keeps HD1 (HIS exception for metal coordination).

### 9.6 Others

* ASN 46 C-terminus (COO, OXT `OC`, plane `CA C O OXT`): `HD21/HD22` are `alg1a` (φ = 0 for HD22, π for HD21).
* ARG 10: `HH11-22` `alg1a`, `HE` `flat_2neigbs`.
* PRO 5: all CH2 `2tetra`, HA `3neigbs`; `PTRANS` link, no H.
* RNA G 1 (1ehz):
  * p5*END maps `OP3`→`O3T`; `rna3p_pur`;
  * H5'/H5'' `2tetra`, H1'-H4' `3neigbs`, HO2' `alg1b` (φ 0), H8/H1 `flat_2neigbs`, H21/H22 `alg1a`.
* Synthetic DNA (`dumps/dna3*`): `DT` → TD dictionary; model `C7` maps to `C5M`. **Reduce names the methyl H
  `H5M1/H5M2/H5M3`** (dictionary names), not H71-73. They are `prop`.

---------------------------------------------------------------------------

## 10. Caveats and things that are easy to get wrong

* The env `RH` differs from `ref/cctbx_project/mmtbx/hydrogens/reduce_hydrogen.py` in two places:
  1. The env applies the HB1/HB2 → HB2/HB3 rename hack to **modified** amino acids as well. Their dictionaries
     usually have no protein interpreter, so the renamed H becomes **unexpected** (no restraints), is not
     parameterized, and is dropped.
  2. The water-selection line moved.

  This spec follows the env.
* `link_metals=Auto` is **not** `True`. Automatic metal links are off, but MCL Zn and Fe-S runs, and metals inside
  multi-atom "other" residues can still get `Misc. bond`.
* Symmetry-related automatic links are never made under these parameters. Disulfides across symmetry **are**
  made (bond only, no angles or dihedrals).
* The second-row buffer (+0.25 Å² on the squared cutoffs) applies only to second-row and heavier elements: elements are already stripped when linking runs.
* `check_valence` counts bogus-placed H within 1.8 Å of an O. This could rarely veto a link.
* The process-time server is fresh, so D amino acids get L dictionary + PEPT-D (§1.4).
* `get_class` is case-sensitive and the short modified-name entries can never match (§2.1).
* `?X` nucleotides without O2' are treated as DNA (§1.4).
* `Var_01`'s weight is about 0 but it still drives riding at chain starts (§8.4).
* `type_energies` strings `'False'`, `'None'` and `''` mean "no type" (KeyError in `get_specific_vdw_radius`).

## 11. Files in `ref/spec_work/`

* `dump_proxies.py <pdb> "<selection>" [neutron]`: raw proxies in grm order with **raw** ideals (`.geo` prints
  model+delta for periodic dihedrals instead), origin ids, types, riding connectivity and parameterization, and
  final coordinates. It hooks `place_hydrogens.add_link_h_restraints` and `setup_riding_h_manager`.
* `dump_restraints.py`: the full `.geo` plus per-atom type and radius TSV and the riding TSV.
* `run_with_log.py`: the interpretation log for Reduce's `process()` call, all bonds and angles with origin ≠ 0,
  and the H removed on links.
* `dump_links_mods.py`: merged link and mod definitions. `energy_table.py`: §7.2.
* `data/`: `modified_aa_names.txt`, `modified_rna_dna_names.txt`, `origin_ids.tsv`,
  `links_and_mods_merged.txt`, `ener_lib_lib_atom.txt`, `aa_h_name_mapping.txt`.
* `dumps/`: outputs for 1crn (residues 1, 3/4/40, 5, 6, 10/14, 46; neutron), 1crn without residue 1, 1ubq
  (HIS 68, LYS 6, TYR 59), 1ehz (residues 1-2), 4fen, 1xso, 3vyk, 6oge, and synthetic DNA.
