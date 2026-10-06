#!/usr/bin/env python
"""
dump_ref.py -- run the ORIGINAL cctbx mmtbx.reduce2 pipeline step by step and dump
intermediate data as JSON, so that a port (Reduce3, Rust) can be compared against it.

Usage (with the reference cctbx Python):
  reference/env/bin/python reference/harness/dump_ref.py reference/pdbs/1crn.pdb            # -> reference/dumps/1crn.json
  reference/env/bin/python reference/harness/dump_ref.py reference/pdbs/1crn.pdb --flips    # -> reference/dumps/1crn_flips.json
  reference/env/bin/python reference/harness/dump_ref.py reference/pdbs/1crn.pdb --neutron  # -> reference/dumps/1crn_neutron.json

The flow below mirrors mmtbx/programs/reduce2.py Program.run() for approach=add with
default parameters (model_id=None, alt_id=None, ...).  The Program object itself is
built through the same CCTBXParser machinery that `mmtbx.reduce2` uses, so the PHIL
parameters (including the probe overrides bump_weight=100, hydrogen_bond_weight=40) and
DataManager options are identical, and reduce2's own _AddHydrogens() / _ReinterpretModel()
methods are called.  Observation-only hooks are installed around the Optimizer
construction to capture data that only exists inside Optimizer.__init__ (initial
ExtraAtomInfo, bonded-neighbor lists, Movers per alternate, OptimizerC final states,
phantom water hydrogens, methyl staggering).  The hooks never change arguments or
return values.

See README.md next to this file for the JSON schema.
"""
from __future__ import absolute_import, division, print_function

import argparse
import io
import json
import math
import os
import re
import sys
import time
import traceback

from libtbx.utils import multi_out
from iotbx.cli_parser import CCTBXParser
from scitbx.array_family import flex

import mmtbx_probe_ext as probeExt
from mmtbx.programs import reduce2
from mmtbx.reduce import Optimizers, Movers
from mmtbx.probe import Helpers
from mmtbx.hydrogens import reduce_hydrogen

SCHEMA_VERSION = 1

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_OUT_DIR = os.path.normpath(os.path.join(HERE, '..', 'dumps'))


# ------------------------------------------------------------------------------
# Small serialisation helpers

def fl(x):
  """float -> JSON-safe float (non-finite -> None)."""
  if x is None:
    return None
  x = float(x)
  return x if math.isfinite(x) else None


def vec(p):
  """Anything vec3-like (tuple, list, scitbx.matrix.col, flex element) -> [x, y, z]."""
  if hasattr(p, 'elems'):
    p = p.elems
  return [fl(p[0]), fl(p[1]), fl(p[2])]


def eai_dict(e):
  """probe ExtraAtomInfo -> dict."""
  return {
    'vdwRadius': fl(e.vdwRadius),
    'isAcceptor': bool(e.isAcceptor),
    'isDonor': bool(e.isDonor),
    'isDummyHydrogen': bool(e.isDummyHydrogen),
    'isIon': bool(e.isIon),
    'charge': int(e.charge),
    'altLoc': str(e.altLoc),
  }


def atom_record(a):
  """Hierarchy atom -> dict with the fields reduce2 writes."""
  ag = a.parent()
  rg = ag.parent()
  ch = rg.parent()
  md = ch.parent()
  return {
    'i_seq': int(a.i_seq),
    'model_id': md.id,
    'chain_id': ch.id,
    'resseq': rg.resseq,
    'resseq_int': rg.resseq_as_int(),
    'icode': rg.icode,
    'resname': ag.resname,
    'altloc': ag.altloc,
    'name': a.name,
    'element': a.element,
    'charge': a.charge,
    'xyz': vec(a.xyz),
    'occ': fl(a.occ),
    'b': fl(a.b),
    'hetero': bool(a.hetero),
  }


def atom_label(a):
  """Short human-readable label, e.g. 'A  12 ASN ND2' (with altloc / model when present)."""
  try:
    ag = a.parent()
    rg = ag.parent()
    ch = rg.parent()
    md = ch.parent()
    s = '{} {}{} {}{} {}'.format(ch.id, rg.resseq.strip(), rg.icode.strip(),
                                 ag.altloc.strip(), ag.resname.strip(), a.name.strip())
    if md.id.strip():
      s = 'model {} '.format(md.id.strip()) + s
    return s
  except Exception:
    return 'phantom i_seq {}'.format(a.i_seq)


def reduce2_res_label(a):
  """Same string reduce2 uses when it reports deleted hydrogens."""
  aName = a.name.strip().upper()
  chainID = a.parent().parent().parent().id
  resName = a.parent().resname.strip().upper()
  resID = str(a.parent().parent().resseq_as_int())
  altLoc = a.parent().altloc
  insertionCode = a.parent().parent().icode.strip()
  return "chain " + str(chainID) + " " + altLoc + resName + " " + resID + insertionCode + " " + aName


# ------------------------------------------------------------------------------
# Restraints / riding / model-level atom typing

def dump_restraints(model):
  grm = model.get_restraints_manager().geometry
  carts = flex.vec3_double()
  for a in model.get_atoms():
    carts.append(a.xyz)
  bps, asu = grm.get_all_bond_proxies(sites_cart=carts)
  hd = model.get_hd_selection()

  bonds = []
  for p in bps:
    i, j = p.i_seqs
    bonds.append({'i': int(i), 'j': int(j), 'distance_ideal': fl(p.distance_ideal),
                  'weight': fl(p.weight), 'slack': fl(p.slack), 'origin_id': int(p.origin_id)})
  bonds_asu = []
  for p in asu:
    bonds_asu.append({'i': int(p.i_seq), 'j': int(p.j_seq), 'j_sym': int(p.j_sym),
                      'distance_ideal': fl(p.distance_ideal), 'weight': fl(p.weight),
                      'slack': fl(p.slack), 'origin_id': int(p.origin_id)})

  angles = []
  for p in grm.angle_proxies:
    i, j, k = p.i_seqs
    if hd[i] or hd[j] or hd[k]:
      angles.append({'i': int(i), 'j': int(j), 'k': int(k), 'angle_ideal': fl(p.angle_ideal),
                     'weight': fl(p.weight), 'origin_id': int(p.origin_id)})

  dihedrals = []
  for p in grm.dihedral_proxies:
    i, j, k, l = p.i_seqs
    if hd[i] or hd[l]:
      d = {'i': int(i), 'j': int(j), 'k': int(k), 'l': int(l), 'angle_ideal': fl(p.angle_ideal),
           'periodicity': int(p.periodicity), 'weight': fl(p.weight),
           'origin_id': int(p.origin_id)}
      alts = list(p.alt_angle_ideals) if p.alt_angle_ideals is not None else []
      if alts:
        d['alt_angle_ideals'] = [fl(x) for x in alts]
      dihedrals.append(d)

  planarities = []
  for p in grm.planarity_proxies:
    seqs = list(p.i_seqs)
    if any(hd[s] for s in seqs):
      planarities.append({'i_seqs': [int(s) for s in seqs],
                          'weights': [fl(w) for w in p.weights],
                          'origin_id': int(p.origin_id)})

  return {
    'counts': {
      'bond_simple_total': len(bps), 'bond_asu_total': len(asu),
      'angle_total': grm.angle_proxies.size() if grm.angle_proxies is not None else 0,
      'angle_with_h': len(angles),
      'dihedral_total': grm.dihedral_proxies.size() if grm.dihedral_proxies is not None else 0,
      'dihedral_with_h_end': len(dihedrals),
      'planarity_total': grm.planarity_proxies.size() if grm.planarity_proxies is not None else 0,
      'planarity_with_h': len(planarities),
    },
    'bonds': bonds,
    'bonds_asu': bonds_asu,
    'angles_with_h': angles,
    'dihedrals_with_h_end': dihedrals,
    'planarities_with_h': planarities,
  }


def dump_riding(model):
  rm = model.get_riding_h_manager()
  if rm is None:
    return None
  out = []
  for idx, p in enumerate(rm.h_parameterization):
    if p is None:
      continue
    out.append({'index': idx, 'htype': p.htype, 'ih': int(p.ih), 'a0': int(p.a0),
                'a1': int(p.a1), 'a2': int(p.a2), 'a3': int(p.a3), 'a': fl(p.a),
                'b': fl(p.b), 'h': fl(p.h), 'n': int(p.n), 'disth': fl(p.disth)})
  return out


def model_atom_typing(model):
  """Per-atom values reported by the model itself (what getExtraAtomInfo reads)."""
  out = []
  te = getattr(model, '_type_energies', None)
  for a in model.get_atoms():
    i = a.i_seq
    d = {'i_seq': int(i)}
    try:
      d['energy_type'] = te[i] if te is not None else None
    except Exception as e:
      d['energy_type'] = None
    try:
      d['model_vdw_radius'] = fl(model.get_specific_vdw_radius(i, False))
    except Exception as e:
      d['model_vdw_radius'] = None
      d['model_vdw_radius_error'] = str(e)
    try:
      d['model_h_bond_type'] = model.get_specific_h_bond_type(i)
    except Exception as e:
      d['model_h_bond_type'] = None
      d['model_h_bond_type_error'] = str(e)
    if a.element_is_ion():
      try:
        d['model_ion_radius'] = fl(model.get_specific_ion_radius(i))
      except Exception:
        d['model_ion_radius'] = None
    out.append(d)
  return out


# ------------------------------------------------------------------------------
# Mover serialisation

def _fixup_dict(fu):
  return {
    # NOTE: atoms are listed even when positions are empty; Optimizer adds FixUp(0).atoms to
    # the set of Mover atoms for which it computes excluded-atom lists.
    'atoms': [int(a.i_seq) for a in fu.atoms],
    'positions': [vec(p) for p in fu.positions],
    'extra_infos': [eai_dict(e) for e in fu.extraInfos],
    'delete_mes': [bool(d) for d in fu.deleteMes],
  }


def mover_snapshot(m, info_at_placement):
  cp = m.CoarsePositions()
  n_coarse = len(cp.positions)
  d = {
    'class': type(m).__name__,
    'info_at_placement': info_at_placement,
    'atoms': [int(a.i_seq) for a in cp.atoms],
    'atom_labels': [atom_label(a) for a in cp.atoms],
    'coarse_positions': [[vec(p) for p in pos] for pos in cp.positions],
    'preference_energies': [fl(e) for e in cp.preferenceEnergies],
    'coarse_extra_infos': [[eai_dict(e) for e in row] for row in cp.extraInfos],
    'coarse_delete_mes': [[bool(x) for x in row] for row in cp.deleteMes],
    'num_fine_positions': [len(m.FinePositions(i).positions) for i in range(n_coarse)],
    'fixups': [_fixup_dict(m.FixUp(i)) for i in range(n_coarse)],
  }
  # Rotator internals (useful to compare angle bookkeeping directly).
  if isinstance(m, Movers._MoverRotator):
    d['rotator'] = {
      'axis_origin': vec(m._axis[0]),
      'axis_dir': vec(m._axis[1]),
      'offset': fl(m._offset),
      'coarse_range': fl(m._coarseRange),
      'coarse_step_degrees': fl(m._coarseStepDegrees),
      'fine_step_degrees': fl(m._fineStepDegrees),
      'do_fine_rotations': bool(m._doFineRotations),
      'has_preference_function': m._preferenceFunction is not None,
      'preferred_orientation_scale': fl(m._preferredOrientationScale),
      'coarse_angles': [fl(x) for x in m._coarseAngles],
      'fine_angles': [fl(x) for x in m._fineAngles],
    }
  flip = {}
  for attr, key in (('_nonFlipPreference', 'non_flip_preference'),
                    ('_enabledFlipStates', 'enabled_flip_states'),
                    ('_enableFixup', 'enable_fixup')):
    if hasattr(m, attr):
      v = getattr(m, attr)
      flip[key] = v if isinstance(v, (bool, int)) else fl(v)
  if flip:
    d['flip'] = flip
  return d


# ------------------------------------------------------------------------------
# Report parsing

_RE_BEGIN = re.compile(r"^\s*BEGIN REPORT: Model (-?\d+) Alt '(.*)':\s*$")
_RE_SET = re.compile(r"^\s*Set of (\d+) Movers:\s+Totals: initial score (\S+), final score (\S+)\s*$")
_RE_MOVER = re.compile(r"^\s*(?P<info>.*?) Initial score: (?P<initial>\S+) final score: (?P<final>\S+) pose (?P<pose>.*?)\s*$")
_RE_ANGLE = re.compile(r"^Angle (\S+) deg")


def parse_reports(info_text):
  """Parse every BEGIN REPORT ... END REPORT block of Optimizer.getInfo()."""
  reports = []
  cur = None
  group = None
  for line in info_text.split('\n'):
    m = _RE_BEGIN.match(line)
    if m:
      cur = {'model_index': int(m.group(1)), 'alt': m.group(2), 'raw_lines': [line],
             'groups': [], 'entries': []}
      group = None
      continue
    if cur is None:
      continue
    cur['raw_lines'].append(line)
    if line.strip() == 'END REPORT':
      reports.append(cur)
      cur = None
      continue
    m = _RE_SET.match(line)
    if m:
      group = {'kind': 'set', 'size': int(m.group(1)), 'initial_total': float(m.group(2)),
               'final_total': float(m.group(3)), 'raw_line': line}
      cur['groups'].append(group)
      continue
    if line.strip() == 'Singleton Movers:':
      group = {'kind': 'singletons', 'raw_line': line}
      cur['groups'].append(group)
      continue
    m = _RE_MOVER.match(line)
    if m:
      pose = m.group('pose')
      e = {'info': m.group('info'), 'initial_score': float(m.group('initial')),
           'final_score': float(m.group('final')), 'pose': pose,
           'group_index': len(cur['groups']) - 1, 'raw_line': line}
      words = pose.split()
      a = _RE_ANGLE.match(pose)
      if a:
        e['angle_deg'] = float(a.group(1))
      if words:
        e['flag'] = words[-1]          # '.', 'BothClash' or 'Uncertain'
      if words and words[0] in ('Flipped', 'Unflipped'):
        e['flipped'] = words[0] == 'Flipped'
      cur['entries'].append(e)
  return reports


# ------------------------------------------------------------------------------
# Observation hooks around Optimizer construction

class OptimizerHooks(object):
  """Installs monkeypatches that only observe; restores everything on exit."""

  def __init__(self, model):
    self.model = model
    self.max_i_seq = Helpers.getMaxISeq(model)
    self.overhead = 0.0
    self.bonded = None            # i_seq -> [i_seq] (in list order) from getBondedNeighborLists
    self.initial_eai = None       # i_seq -> dict, right after Helpers.getExtraAtomInfo
    self.initial_eai_warnings = None
    self.runs = []                # one per _PlaceMovers() call (= per model index / alternate)
    self.optcs = []               # OptimizerC instances, one per run
    self.staggered_methyls = []   # MoverTetrahedralMethylRotator constructions (not Movers)
    self._saved = []

  def _patch(self, obj, name, new):
    self._saved.append((obj, name, getattr(obj, name)))
    setattr(obj, name, new)

  def __enter__(self):
    hooks = self
    orig_bnl = Helpers.getBondedNeighborLists
    orig_eai = Helpers.getExtraAtomInfo
    orig_place = Optimizers.Optimizer._PlaceMovers
    orig_optc = Optimizers.OptimizerC
    orig_tet = Movers.MoverTetrahedralMethylRotator

    def bnl_hook(atoms, *a, **k):
      ret = orig_bnl(atoms, *a, **k)
      t0 = time.time()
      if hooks.bonded is None:
        hooks.bonded = {int(at.i_seq): [int(n.i_seq) for n in ret[at]] for at in atoms}
      hooks.overhead += time.time() - t0
      return ret

    def eai_hook(*a, **k):
      ret = orig_eai(*a, **k)
      t0 = time.time()
      if hooks.initial_eai is None:
        m = ret.extraAtomInfo
        hooks.initial_eai = {int(at.i_seq): eai_dict(m.getMappingFor(at))
                             for at in hooks.model.get_atoms()}
        hooks.initial_eai_warnings = ret.warnings
      hooks.overhead += time.time() - t0
      return ret

    def place_hook(opt, atoms, rotatableHydrogenIDs, bondedNeighborLists, hParameters,
                   addFlipMovers, alt):
      t0 = time.time()
      phantoms = []
      for at in atoms:
        if at.i_seq > hooks.max_i_seq:
          parents = bondedNeighborLists.get(at, [])
          phantoms.append({
            'i_seq': int(at.i_seq),
            'xyz': vec(at.xyz),
            'occ': fl(at.occ),
            'b': fl(at.b),
            'parent_i_seq': int(parents[0].i_seq) if parents else None,
            'extra_info': eai_dict(opt._extraAtomInfo.getMappingFor(at)),
          })
      # Donor/acceptor state right before Mover placement (after phantom placement and
      # Helpers.fixupExplicitDonors), for every non-phantom atom of this conformer.
      eai_before_place = {int(at.i_seq): eai_dict(opt._extraAtomInfo.getMappingFor(at))
                          for at in atoms if at.i_seq <= hooks.max_i_seq}
      xyz_before_place = {int(at.i_seq): vec(at.xyz) for at in atoms if at.i_seq <= hooks.max_i_seq}
      hooks.overhead += time.time() - t0

      ret = orig_place(opt, atoms, rotatableHydrogenIDs, bondedNeighborLists, hParameters,
                       addFlipMovers, alt)

      t0 = time.time()
      model_id = None
      for at in atoms:
        if at.i_seq <= hooks.max_i_seq:
          model_id = at.parent().parent().parent().parent().id
          break
      run = {
        'alt': alt,
        'model_id': model_id,
        'n_atoms_considered': len(atoms),
        'conformer_atom_i_seqs': [int(at.i_seq) for at in atoms if at.i_seq <= hooks.max_i_seq],
        'phantoms': phantoms,
        'extra_info_before_placement': eai_before_place,
        'xyz_before_placement': xyz_before_place,
        'delete_atoms_from_placement': sorted(int(at.i_seq) for at in ret),
        '_movers': list(opt._movers),            # object refs, resolved later
        '_moverInfo': opt._moverInfo,            # dict ref; values get ' Initial score: ...' later
        'mover_snapshots': [mover_snapshot(m, opt._moverInfo.get(m)) for m in opt._movers],
      }
      hooks.runs.append(run)
      hooks.overhead += time.time() - t0
      return ret

    def optc_factory(*a, **k):
      o = orig_optc(*a, **k)
      hooks.optcs.append(o)
      return o

    class TetCapture(orig_tet):
      def __init__(self, atom, bondedNeighborLists, hParameters, *a, **k):
        t0 = time.time()
        before = {int(atom.i_seq): vec(atom.xyz)}
        try:
          for n in bondedNeighborLists[atom]:
            before[int(n.i_seq)] = vec(n.xyz)
        except Exception:
          pass
        hooks.overhead += time.time() - t0
        orig_tet.__init__(self, atom, bondedNeighborLists, hParameters, *a, **k)
        t0 = time.time()
        hooks.staggered_methyls.append({
          'carbon': int(atom.i_seq),
          'carbon_label': atom_label(atom),
          'atoms': [int(x.i_seq) for x in self._atoms],
          'xyz_before': [before.get(int(x.i_seq)) for x in self._atoms],
          'xyz_after': [vec(x.xyz) for x in self._atoms],
          'axis_origin': vec(self._axis[0]),
          'axis_dir': vec(self._axis[1]),
          'offset': fl(self._offset),
        })
        hooks.overhead += time.time() - t0

    self._patch(Helpers, 'getBondedNeighborLists', bnl_hook)
    self._patch(Helpers, 'getExtraAtomInfo', eai_hook)
    self._patch(Optimizers.Optimizer, '_PlaceMovers', place_hook)
    self._patch(Optimizers, 'OptimizerC', optc_factory)
    self._patch(Movers, 'MoverTetrahedralMethylRotator', TetCapture)
    return self

  def __exit__(self, *exc):
    for obj, name, old in reversed(self._saved):
      setattr(obj, name, old)
    self._saved = []
    return False


class PlaceHydrogensCapture(object):
  """Wraps reduce_hydrogen.place_hydrogens so the placement object can be inspected."""

  def __init__(self):
    self.objs = []
    self._orig = None

  def __enter__(self):
    orig = reduce_hydrogen.place_hydrogens
    cap = self

    class Captured(orig):
      def __init__(self, *a, **k):
        orig.__init__(self, *a, **k)
        cap.objs.append(self)

    self._orig = orig
    reduce_hydrogen.place_hydrogens = Captured
    return self

  def __exit__(self, *exc):
    reduce_hydrogen.place_hydrogens = self._orig
    return False


# ------------------------------------------------------------------------------
# Main driver

def run(model_file, flips=False, neutron=False, out_dir=DEFAULT_OUT_DIR, model_tag=None,
        write_pdb=True, indent=None):
  model_file = os.path.abspath(model_file)
  base = model_tag or os.path.splitext(os.path.basename(model_file))[0]
  suffix = ('_flips' if flips else '') + ('_neutron' if neutron else '')
  out_json = os.path.join(out_dir, base + suffix + '.json')
  out_pdb = os.path.join(out_dir, base + suffix + '.pdb')

  timings = {}
  t_total = time.time()

  # --- Build the reduce2 Program exactly like `mmtbx.reduce2` (iotbx.cli_parser.run_program).
  log_buf = io.StringIO()
  logger = multi_out()
  logger.register('buffer', log_buf)
  # output.write_files=False: nothing is written by the Program object (run() is not called
  # anyway; the harness performs the run() steps itself below).
  phil_args = [model_file, 'output.write_files=False']
  if flips:
    phil_args.append('add_flip_movers=True')
  if neutron:
    phil_args.append('use_neutron_distances=True')
  parser = CCTBXParser(program_class=reduce2.Program, logger=logger)
  parser.parse_args(phil_args)
  prog = reduce2.Program(parser.data_manager, parser.working_phil.extract(),
                         master_phil=parser.master_phil, logger=logger)
  prog.validate()
  params = prog.params
  assert params.approach == 'add', 'harness only replicates approach=add'
  assert params.model_id is None, 'harness only replicates model_id=None'

  # --- Program.run(), approach=add ------------------------------------------------
  prog._bondedNeighborDepth = params.bonded_neighbor_depth

  t0 = time.time()
  prog.model = prog.data_manager.get_model()
  prog.model = prog.model.select(~prog.model.selection('element X'))
  prog._output_cs = reduce_hydrogen.get_output_crystal_symmetry(prog.data_manager.get_model())
  prog.model.add_crystal_symmetry_if_necessary(box_cushion=5)
  if prog.data_manager.has_restraints():
    prog.model.set_stop_for_unknowns(params.stop_on_any_missing_hydrogen)
    prog.model.process(make_restraints=False)
  timings['load_model_s'] = time.time() - t0
  input_cs = prog.model.crystal_symmetry()

  # Add hydrogens (reduce2's own method -> reduce_hydrogen.place_hydrogens)
  with PlaceHydrogensCapture() as phc:
    startAdd = time.time()
    prog._AddHydrogens()
    doneAdd = time.time()
  timings['h_placement_s'] = doneAdd - startAdd
  model = prog.model

  # ---- (1) atoms right after H placement
  t0 = time.time()
  after_h = [atom_record(a) for a in model.get_atoms()]
  riding_existed_before_optimizer = model.get_riding_h_manager() is not None
  # ---- (2) restraints at this point
  restraints = dump_restraints(model)
  # ---- (4, model part) typing reported by the model
  typing = model_atom_typing(model)
  ph = phc.objs[-1] if phc.objs else None
  h_placement_info = None
  if ph is not None:
    h_placement_info = {
      'n_H_initial': getattr(ph, 'n_H_initial', None),
      'n_H_final': getattr(ph, 'n_H_final', None),
      'no_H_placed_mlq': list(getattr(ph, 'no_H_placed_mlq', [])),
      'site_labels_disulfides': list(getattr(ph, 'site_labels_disulfides', [])),
      'site_labels_no_para': list(getattr(ph, 'site_labels_no_para', [])),
      'site_labels_tertiary_amide': list(getattr(ph, 'site_labels_tertiary_amide', [])),
      'site_labels_missing_neighbor': list(getattr(ph, 'site_labels_missing_neighbor', [])),
      'sl_removed': [str(x) for x in getattr(ph, 'sl_removed', [])],
      'auto_restraint_names': list(getattr(ph, 'auto_restraint_names', [])),
    }
  cs = model.crystal_symmetry()
  symmetry_info = {
    'unit_cell': list(cs.unit_cell().parameters()) if cs and cs.unit_cell() else None,
    'space_group': str(cs.space_group_info()) if cs and cs.space_group_info() else None,
    'output_cs_unit_cell': (list(prog._output_cs.unit_cell().parameters())
                            if prog._output_cs is not None else None),
  }
  timings['dump_after_h_s'] = time.time() - t0

  # --- Optimizer (same arguments reduce2 passes) -----------------------------------
  with OptimizerHooks(model) as hooks:
    startOpt = time.time()
    opt = Optimizers.Optimizer(params.probe, params.add_flip_movers,
      model, altID=params.alt_id,
      preferenceMagnitude=params.preference_magnitude,
      bondedNeighborDepth=prog._bondedNeighborDepth,
      nonFlipPreference=params.non_flip_preference,
      skipBondFixup=params.skip_bond_fix_up,
      flipStates=params.set_flip_states,
      verbosity=params.verbosity,
      cliqueOutlineFileName=params.output.clique_outline_file_name,
      fillAtomDump=params.output.print_atom_info)
    doneOpt = time.time()
  timings['optimization_s'] = doneOpt - startOpt
  timings['optimization_hook_overhead_s'] = hooks.overhead
  timings['optimization_minus_hook_overhead_s'] = (doneOpt - startOpt) - hooks.overhead

  t0 = time.time()
  warnings_text = opt.getWarnings()
  info_text = opt.getInfo()
  hToDelete = opt.getHydrogensToDelete()

  # ---- (3) riding parameterization actually used by the optimizer
  riding = dump_riding(model)

  # ---- (4) atom info
  final_eai = {}
  for a in model.get_atoms():
    final_eai[int(a.i_seq)] = eai_dict(opt._extraAtomInfo.getMappingFor(a))
  atom_info = []
  for d in typing:
    i = d['i_seq']
    rec = dict(d)
    rec['bonded'] = hooks.bonded.get(i, []) if hooks.bonded is not None else None
    rec['initial'] = hooks.initial_eai.get(i) if hooks.initial_eai is not None else None
    rec['final'] = final_eai.get(i)
    atom_info.append(rec)

  # ---- (5) movers, per optimizer run (model index x alternate)
  reports = parse_reports(info_text)
  final_movers = set(id(m) for m in getattr(opt, '_movers', []))
  movers_out = []
  runs_out = []
  for run_idx, run in enumerate(hooks.runs):
    optc = hooks.optcs[run_idx] if run_idx < len(hooks.optcs) else None
    report = reports[run_idx] if run_idx < len(reports) else None
    entries = list(report['entries']) if report else []
    used = [False] * len(entries)
    mover_indices = []
    for k, (m, snap) in enumerate(zip(run['_movers'], run['mover_snapshots'])):
      d = {'run': run_idx, 'alt': run['alt'], 'index_in_run': k}
      d.update(snap)
      info = run['_moverInfo'].get(m)
      d['info'] = info
      d['in_final_optimizer_movers'] = id(m) in final_movers
      if optc is not None:
        d['final'] = {'coarse_index': int(optc.GetCoarseLocation(m)),
                      'fine_index': int(optc.GetFineLocation(m)),
                      'score': fl(optc.GetHighScore(m))}
      # Match the report line '   <moverInfo> final score: X pose ...', where moverInfo is
      # the placement-time info string followed by ' Initial score: Y'.
      rep = None
      base_info = snap['info_at_placement']
      if base_info is not None:
        for ei, e in enumerate(entries):
          if not used[ei] and e['info'] == base_info:
            used[ei] = True
            rep = e
            break
      d['report'] = rep
      mover_indices.append(len(movers_out))
      movers_out.append(d)
    runs_out.append({
      'run': run_idx,
      'alt': run['alt'],
      'model_id': run['model_id'],
      'report_model_index': report['model_index'] if report else None,
      'n_atoms_considered': run['n_atoms_considered'],
      'conformer_atom_i_seqs': run['conformer_atom_i_seqs'],
      'phantom_hydrogens': run['phantoms'],
      'extra_info_before_placement': [dict(i_seq=i, **v) for i, v in
                                      sorted(run['extra_info_before_placement'].items())],
      'xyz_before_placement': [{'i_seq': i, 'xyz': v} for i, v in
                               sorted(run['xyz_before_placement'].items())],
      'delete_atoms_from_placement': run['delete_atoms_from_placement'],
      'mover_indices': mover_indices,
      'report': report,
      'num_calculated_atoms': int(optc.GetNumCalculatedAtoms()) if optc is not None else None,
      'num_cached_atoms': int(optc.GetNumCachedAtoms()) if optc is not None else None,
    })

  # ---- (6) hydrogens to delete
  h_del = [{'i_seq': int(a.i_seq), 'label': reduce2_res_label(a), 'id_str': a.id_str()}
           for a in sorted(hToDelete, key=lambda x: x.i_seq)]

  # Positions after optimization, before deletion (same i_seqs as after_h_placement)
  after_opt = [{'i_seq': int(a.i_seq), 'xyz': vec(a.xyz)} for a in model.get_atoms()]
  pre_delete_ids = {}
  for a in model.get_atoms():
    pre_delete_ids[a.memory_id()] = int(a.i_seq)
  timings['dump_after_opt_s'] = time.time() - t0

  # --- Deletion + reinterpretation, as reduce2 does ----------------------------------
  t0 = time.time()
  for a in hToDelete:
    a.parent().remove_atom(a)
  prog._ReinterpretModel(False)
  timings['delete_and_reinterpret_s'] = time.time() - t0

  # ---- (7) final atoms
  final_atoms = []
  for a in prog.model.get_atoms():
    r = atom_record(a)
    r['i_seq_after_h_placement'] = pre_delete_ids.get(a.memory_id())
    final_atoms.append(r)

  # Write the model text reduce2 would write (for byte-level comparison).
  pdb_text = None
  if write_pdb:
    output_cs = prog._output_cs is not None
    if output_cs:
      prog.model.set_unit_cell_crystal_symmetry(prog._output_cs)
    pdb_text = prog.model.model_as_pdb(output_cs=output_cs)

  timings['total_s'] = time.time() - t_total

  probe = params.probe
  out = {
    'schema_version': SCHEMA_VERSION,
    'id': base,
    'input_file': model_file,
    'options': {'add_flip_movers': bool(flips), 'use_neutron_distances': bool(neutron)},
    'reduce2_version': reduce2.version,
    'params': {
      'approach': params.approach,
      'n_terminal_charge': params.n_terminal_charge,
      'exclude_water': params.exclude_water,
      'keep_existing_H': params.keep_existing_H,
      'use_neutron_distances': params.use_neutron_distances,
      'preference_magnitude': params.preference_magnitude,
      'alt_id': params.alt_id,
      'model_id': params.model_id,
      'add_flip_movers': params.add_flip_movers,
      'non_flip_preference': params.non_flip_preference,
      'skip_bond_fix_up': params.skip_bond_fix_up,
      'set_flip_states': params.set_flip_states,
      'verbosity': params.verbosity,
      'bonded_neighbor_depth': params.bonded_neighbor_depth,
      'probe': {k: getattr(probe, k) for k in sorted(probe.__dict__) if not k.startswith('_')},
      # Optimizer.__init__ keyword arguments that reduce2 does NOT pass (defaults used):
      'optimizer_defaults_not_passed_by_reduce2': {
        'modelIndex': 0, 'useNeutronDistances': False, 'minOccupancy': 0.02},
    },
    'symmetry': symmetry_info,
    'h_placement': h_placement_info,
    'after_h_placement': after_h,
    'restraints': restraints,
    'riding_existed_before_optimizer': riding_existed_before_optimizer,
    'riding': riding,
    'atom_info': atom_info,
    'initial_extra_atom_info_warnings': hooks.initial_eai_warnings,
    'optimizer_runs': runs_out,
    'movers': movers_out,
    'staggered_methyls': hooks.staggered_methyls,
    'after_optimization_before_deletion': after_opt,
    'info_text': info_text,
    'warnings_text': warnings_text,
    'hydrogens_to_delete': h_del,
    'final_atoms': final_atoms,
    'timings': timings,
    'program_log': log_buf.getvalue(),
  }

  os.makedirs(out_dir, exist_ok=True)
  tmp = out_json + '.tmp'
  with open(tmp, 'w') as f:
    if indent:
      json.dump(out, f, indent=indent)
    else:
      json.dump(out, f, separators=(',', ':'))
  os.replace(tmp, out_json)
  if pdb_text is not None:
    with open(out_pdb, 'w') as f:
      f.write(pdb_text)
  return out_json, out


def main(argv=None):
  ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
  ap.add_argument('model_file', help='PDB or mmCIF model file')
  ap.add_argument('--flips', action='store_true', help='add_flip_movers=True')
  ap.add_argument('--neutron', action='store_true', help='use_neutron_distances=True')
  ap.add_argument('--out-dir', default=DEFAULT_OUT_DIR, help='output directory (default: reference/dumps)')
  ap.add_argument('--id', default=None, help='structure id used for file names (default: input basename)')
  ap.add_argument('--no-pdb', action='store_true', help='do not write the reduce2-equivalent .pdb next to the JSON')
  ap.add_argument('--indent', type=int, default=None, help='pretty-print JSON with this indent')
  a = ap.parse_args(argv)
  t0 = time.time()
  path, out = run(a.model_file, flips=a.flips, neutron=a.neutron, out_dir=a.out_dir,
                  model_tag=a.id, write_pdb=not a.no_pdb, indent=a.indent)
  print('wrote {} ({:.1f} MB) in {:.1f}s: {} atoms after H placement, {} final atoms, {} movers, '
        'H placement {:.2f}s, optimization {:.2f}s'.format(
          path, os.path.getsize(path) / 1e6, time.time() - t0, len(out['after_h_placement']),
          len(out['final_atoms']), len(out['movers']), out['timings']['h_placement_s'],
          out['timings']['optimization_s']))


if __name__ == '__main__':
  main()
