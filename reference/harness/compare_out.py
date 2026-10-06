#!/usr/bin/env python
"""
compare_out.py -- check harness dumps against the files written by the real mmtbx.reduce2.

For each reference/dumps/<id>[_flips].json it compares `final_atoms` with
reference/out/<id>H.pdb (or <id>FH.pdb for _flips), matching atoms by
(model, chain, resseq, icode, altloc, resname, name), and reports the max coordinate
deviation plus missing/extra atoms.  It also does a byte comparison of the
harness-written reference/dumps/<id>[_flips].pdb with the reduce2 output file.

Usage: reference/env/bin/python reference/harness/compare_out.py [id ...] [--tol 0.001]
"""
from __future__ import print_function
import argparse
import glob
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REF = os.path.normpath(os.path.join(HERE, '..'))
DUMPS = os.path.join(REF, 'dumps')
OUT = os.path.join(REF, 'out')


def parse_pdb(path):
  atoms = {}
  dups = 0
  model = ''
  for line in open(path):
    rec = line[:6]
    if rec == 'MODEL ':
      model = line[10:14].strip()
    elif rec in ('ATOM  ', 'HETATM'):
      key = (model, line[20:22].strip(), line[22:26].strip(), line[26:27].strip(),
             line[16:17].strip(), line[17:20].strip(), line[12:16].strip())
      xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
      if key in atoms:
        dups += 1
      atoms[key] = xyz
  return atoms, dups


def dump_atoms(final_atoms):
  atoms = {}
  for a in final_atoms:
    key = (a['model_id'].strip(), a['chain_id'].strip(), a['resseq'].strip(), a['icode'].strip(),
           a['altloc'].strip(), a['resname'].strip(), a['name'].strip())
    atoms[key] = a['xyz']
  return atoms


def compare(dump_path, tol):
  base = os.path.basename(dump_path)[:-5]
  flips = base.endswith('_flips')
  sid = base[:-6] if flips else base
  if '_neutron' in sid:
    return None
  ref_pdb = os.path.join(OUT, sid + ('FH' if flips else 'H') + '.pdb')
  if not os.path.exists(ref_pdb):
    return {'name': base, 'status': 'no reference output'}
  d = json.load(open(dump_path))
  mine = dump_atoms(d['final_atoms'])
  ref, dups = parse_pdb(ref_pdb)
  missing = [k for k in ref if k not in mine]
  extra = [k for k in mine if k not in ref]
  maxdev = 0.0
  worst = None
  nbad = 0
  for k, r in ref.items():
    if k in mine:
      m = mine[k]
      dev = math.sqrt(sum((m[i] - r[i]) ** 2 for i in range(3)))
      if dev > tol:
        nbad += 1
      if dev > maxdev:
        maxdev, worst = dev, k
  my_pdb = os.path.join(DUMPS, base + '.pdb')
  identical = None
  if os.path.exists(my_pdb):
    identical = open(my_pdb, 'rb').read() == open(ref_pdb, 'rb').read()
  ok = (not missing) and (not extra) and nbad == 0
  return {'name': base, 'status': 'OK' if ok else 'MISMATCH', 'n_ref': len(ref), 'n_dump': len(mine),
          'missing': len(missing), 'extra': len(extra), 'n_over_tol': nbad, 'max_dev': maxdev,
          'worst': worst, 'pdb_identical': identical, 'ref_duplicate_keys': dups,
          'missing_examples': missing[:5], 'extra_examples': extra[:5]}


def main():
  ap = argparse.ArgumentParser()
  ap.add_argument('ids', nargs='*')
  ap.add_argument('--tol', type=float, default=0.001)
  a = ap.parse_args()
  paths = sorted(glob.glob(os.path.join(DUMPS, '*.json')))
  if a.ids:
    paths = [p for p in paths if os.path.basename(p).split('_')[0].split('.')[0] in a.ids]
  allok = True
  for p in paths:
    r = compare(p, a.tol)
    if r is None:
      continue
    if r['status'] != 'OK':
      allok = False
    if 'n_ref' in r:
      print('{name:14s} {status:8s} ref={n_ref:6d} dump={n_dump:6d} missing={missing} extra={extra} '
            'over_tol={n_over_tol} max_dev={max_dev:.4f} pdb_identical={pdb_identical}'.format(**r))
      if r['status'] != 'OK':
        print('    worst', r['worst'], 'missing', r['missing_examples'], 'extra', r['extra_examples'])
    else:
      print('{name:14s} {status}'.format(**r))
  sys.exit(0 if allok else 1)


if __name__ == '__main__':
  main()
