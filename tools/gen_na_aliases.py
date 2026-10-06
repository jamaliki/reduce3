"""Print the Rust table for iotbx.pdb.rna_dna_atom_names_backbone_aliases."""
import iotbx.pdb as p
for a, r in sorted(p.rna_dna_atom_names_backbone_aliases.items()):
  print('    ("%s", "%s"),' % (a, r.strip()))
