# Native validation chemistry

`interp::validation::build` creates weighted restraints on the supplied atom
order and coordinates. It does not prepare hydrogens or change the source
structure. The caller must use `MonLib::load_validation` and check every
`Coverage` record before treating a check as complete.

The validation loader starts with `mon_lib/list/mon_lib_list.cif`, then applies
the GeoStd library. The older preparation loader keeps its prior merged-list
choice. This separation matters because the effective libraries can assign
different targets to the same linked N-H bond.

## Restraints and interpretation

The graph contains bond and angle targets with weights `1 / sigma^2`, signed
chiral-volume targets, per-atom plane weights, hydrogen flags, atom paths,
typed formal charges, dictionary provenance and the covalent graph used by
Probe. Atom paths, coordinates, names, elements, altlocs and hydrogen flags
must match the input structure.

Protein interpretation supports terminal modifications, D-amino acid
modifications, cis/trans peptide links, disulfides and the protein backbone
CDL v1.2 tables. The builder emits each residue's restraints before applying
its outgoing peptide-link modification, as CCTBX does. Atom energy types also
come from the earlier residue interpretation. Applying all modifications
before emitting restraints changes targets and therefore does not match the
reference.

Nuclear mode uses the nuclear bond target when the library supplies one, but
keeps the ordinary bond ESD. CCTBX's `add_bond_proxies` follows that rule. The
library reader retains nuclear ESD metadata for provenance; validation does
not use it as a weight.

Chiral targets derive from the dictionary's three bond lengths and three
angles. The default volume sigma is 0.2. The graph records both-signs and
pseudochiral flags for the evaluator. Planes require at least four supplied
sites. SS targets use the effective `SS` link, including its 104.2-degree
CB-SG-SG angle with sigma 2.1.

The CDL exporter preserves the reference table's double precision values and
records the source and binary SHA-256 digests. The lookup uses nearest
ten-degree bins, Python's ties-to-even rounding and angle wrapping. It also
matches the reference adjustment that makes C-N-CA, C-N-H and CA-N-H targets
sum to 360 degrees. Table and accepted-ion-label exports use CCTBX 2026.9;
`CCTBX_LICENSE.txt` contains its license.

## Supplied and optional hydrogens

The graph validates supplied water hydrogen atoms. Absent water hydrogens do
not make an oxygen-only water incomplete; the Probe stage creates its own
temporary water phantoms. For neutral HIS, the supplied HD1 or HE2 proton
allows the other ring proton to be absent. Missing both remains incomplete.
Disulfide CYS may omit HG. Other missing expected hydrogen or restraint sites
keep the relevant check incomplete. These rules never infer that a supplied
hydrogen has a valid parent when its topology is unknown.

`typed_charges` matches the formal charge left by CCTBX scattering-label
normalization. For example, the reference neutralizes H1+ but retains O1-.
Probe must use these typed charges; reading the raw source charge instead can
change its hydrogen-bond overlap threshold. This normalization does not
change input atoms or legacy hydrogen preparation.

## Coverage limits

The graph reports unsupported chemistry instead of substituting guessed
targets, weights or radii. Current unsupported cases include ASU symmetry
bonds, multiple models in one view, nucleotide pucker/CDL, enol peptide
targets, unknown or CCD-derived components, cross-residue chiral links,
unresolved explicit links, and automatic metal or glycosidic links without a
validated weighted dictionary graph. Missing energy types, conflicting
proxies and missing nonoptional sites also prevent completion.

CCTBX can swap coordinates and B factors under symmetric residue labels
before interpreting them. The native builder preserves sites and detects
inputs needing those ARG/ASP/GLU/PHE/TYR/VAL/LEU swaps, then reports them
unsupported. This is separate from weighted graph parity on unchanged sites.
The caller must also reject input metadata that its parser drops, including
custom restraints and nonidentity `struct_conn` symmetry, or route it to the
reference implementation.

These checks establish the tested native domain. They do not establish
coverage of every CCTBX interpretation option or the entire PDB. The local
archive scanner and its immutable criteria remain unchanged.

## Reference tests

Synthetic tests check effective-library selection, sequential outgoing and
incoming modifications, water-only inputs and both neutral HIS states.
Unit tests check CDL binning, signed chiral targets, formal-charge
normalization and symmetric-label swap detection.

Independent full 1CRN nuclear-mode captures match all 652 bond, 1183 angle,
56 chirality and 101 plane targets and weights. The five-residue crambin
capture matches all 66 bonds, 119 angles, seven chiralities and eight planes.
These counts include the disulfide and CDL cases in those fixtures. They are
reference-parity results, not performance measurements.

The 2160-atom protein case also matches graph targets, weights, flags,
initial Probe typing, clash pairs and scores. CCTBX rounds some input
hydrogen coordinates from five to three decimals, changing sites by at most
0.0005 A, while the native path preserves them. Full parity therefore uses
an explicit reference-processed coordinate fixture. Its manifest records
the changed sites and original input hash; the raw-input discrepancy remains
recorded rather than being treated as floating-point agreement.
