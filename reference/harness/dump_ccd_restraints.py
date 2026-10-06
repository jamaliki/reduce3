#!/usr/bin/env python
"""
dump_ccd_restraints.py -- write the restraints Reduce2 builds for residues that only the
CCD describes (mmtbx.hydrogens.reduce_hydrogen.get_h_restraints(resname, strict=False)),
as JSON lines, so Reduce3's port can be compared with them.

Usage (with the reference cctbx Python, which needs RDKit):
  reference/env/bin/python reference/harness/dump_ccd_restraints.py 6MO FE2 DGT > ccd.jsonl
  reference/env/bin/python reference/harness/dump_ccd_restraints.py --sample 500 > ccd.jsonl

--sample N takes N entries, in sorted order with a fixed stride, from the CCD entries
that neither mon_lib nor geostd describe.

Each line is {"id": ..., "ok": bool, "atoms": [[id, type_symbol], ...],
"bonds": [[a1, a2, type, value_dist], ...], "angles": [[a1, a2, a3, value], ...],
"tors": [[id, a1, a2, a3, a4, value], ...]}, with the values get_h_restraints sets
(non-finite ones, such as the torsion of three collinear atoms, as strings).
"ok" is false when get_h_restraints returns None.
"""
from __future__ import absolute_import, division, print_function

import io
import json
import math
import os
import sys
from contextlib import redirect_stdout

from mmtbx.hydrogens.reduce_hydrogen import get_h_restraints


def chem_data():
  import libtbx.env_config  # noqa: F401
  from mmtbx.chemical_components import get_cif_filename
  return os.path.dirname(os.path.dirname(os.path.dirname(get_cif_filename('ALA'))))


def sample(n):
  root = chem_data()
  ccd = os.path.join(root, 'chemical_components')
  ids = []
  for letter in sorted(os.listdir(ccd)):
    d = os.path.join(ccd, letter)
    if not os.path.isdir(d):
      continue
    for f in sorted(os.listdir(d)):
      if not (f.startswith('data_') and f.endswith('.cif')):
        continue
      code = f[5:-4]
      first = code[0].lower()
      if os.path.exists(os.path.join(root, 'geostd', first, 'data_%s.cif' % code)):
        continue
      if os.path.exists(os.path.join(root, 'mon_lib', first, '%s.cif' % code)):
        continue
      ids.append(code)
  stride = max(1, len(ids) // n)
  return ids[::stride][:n]


def number(v):
  return v if not isinstance(v, float) or math.isfinite(v) else str(v)


def dump(code):
  quiet = io.StringIO()
  try:
    with redirect_stdout(quiet):
      md = get_h_restraints(code, strict=False)
  except Exception as e:  # Reduce2 would stop here
    return {'id': code, 'ok': False, 'error': repr(e)}
  if md is None:
    return {'id': code, 'ok': False}
  return {
    'id': code,
    'ok': True,
    'group': md.chem_comp.group,
    'atoms': [[a.atom_id, a.type_symbol] for a in md.atom_list],
    'bonds': [[b.atom_id_1, b.atom_id_2, b.type, number(b.value_dist)] for b in md.bond_list],
    'angles': [[a.atom_id_1, a.atom_id_2, a.atom_id_3, number(a.value_angle)]
               for a in md.angle_list],
    'tors': [[t.id, t.atom_id_1, t.atom_id_2, t.atom_id_3, t.atom_id_4, number(t.value_angle)]
             for t in md.tor_list],
  }


def main(argv):
  if argv[:1] == ['--sample']:
    codes = sample(int(argv[1]))
  else:
    codes = argv
  for code in codes:
    print(json.dumps(dump(code)))
    sys.stdout.flush()


if __name__ == '__main__':
  main(sys.argv[1:])
