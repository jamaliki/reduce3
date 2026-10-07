#!/usr/bin/env python3
"""
neutron_benchmark.py -- compare Reduce3's hydrogens with the hydrogens and deuteriums that
neutron crystallography located, on the PDB's neutron structures.

Usage:
  tools/neutron_benchmark.py --pdb-dir DIR                       # fixed and compat mode
  tools/neutron_benchmark.py --pdb-dir DIR --config nopref="planar_hydroxyl_preference=0 acid_syn_preference=0"
  tools/neutron_benchmark.py --pdb-dir DIR --find                # rescan DIR, rewrite the entry list

DIR holds mmCIF files, either flat (1abc.cif or 1abc.cif.gz) or in the PDB's divided layout
(ab/1abc.cif.gz). The entries are listed in tools/neutron_entries.txt; entries missing from DIR
are skipped. Each --config NAME="ARGS" runs the reduce3 binary with those extra arguments
(default: fixed="" and compat="--compat"), and the report has one column per configuration.
Needs Python 3 with gemmi (pip install gemmi) and a release build of reduce3.

What is compared, residue by residue (residues without alternate conformations, deposited
hydrogens at occupancy 0.5 or more, deuterium named as hydrogen):
  * His: whether the ring is oriented as deposited (ND1 within 0.6 A of the deposited ND1),
    and where it is, whether the protonation (HID, HIE or HIP) matches. Histidines left with
    no ring hydrogen are also counted, over every His not bound to a metal.
  * Asn/Gln: whether the amide nitrogen is where the deposited one is (deposited amide
    hydrogens required, so the neutron data fixed the orientation).
  * Ser/Thr/Tyr hydroxyls: the difference between the deposited and the placed hydrogen's
    dihedral, split by what the deposited hydrogen bonds to (an N or O within 2.6 A).
  * Cys: whether the thiol hydrogen's presence matches, for cysteines on a metal and others.
  * Asp/Glu: deposited and placed carboxylic-acid hydrogens.
"""
import argparse
import collections
import gzip
import math
import os
import subprocess
import sys
from multiprocessing import Pool

import gemmi

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
ENTRIES = os.path.join(HERE, 'neutron_entries.txt')
HYDROXYL = {'SER': ('CA', 'CB', 'OG', 'HG'), 'THR': ('CA', 'CB', 'OG1', 'HG1'), 'TYR': ('CE1', 'CZ', 'OH', 'HH')}
ACID = {'ASP': (('OD1', 'OD2'), ('HD1', 'HD2')), 'GLU': (('OE1', 'OE2'), ('HE1', 'HE2'))}
WATERS = ('HOH', 'DOD', 'WAT')


# --- finding the entries ----------------------------------------------------------------

def is_neutron(path):
    """Whether the file's _exptl.method mentions neutron diffraction (read from its head)."""
    opener = gzip.open if path.endswith('.gz') else open
    try:
        with opener(path, 'rb') as f:
            return b'NEUTRON DIFFRACTION' in f.read(400000).upper()
    except OSError:
        return False


def mmcif_files(pdb_dir):
    for dirpath, _, names in os.walk(pdb_dir):
        for n in names:
            if n.endswith(('.cif', '.cif.gz')):
                yield os.path.join(dirpath, n)


def entry_path(pdb_dir, code):
    for p in (os.path.join(pdb_dir, code[1:3], code + '.cif.gz'), os.path.join(pdb_dir, code + '.cif.gz'),
              os.path.join(pdb_dir, code + '.cif')):
        if os.path.exists(p):
            return p
    return None


# --- comparison ---------------------------------------------------------------------------

def residues(st):
    """(chain, seqnum, icode) -> (resname, {atom name: atom}) for single-conformer residues,
    deuteriums named as hydrogens, keeping the fuller of an H/D pair on one site."""
    out = {}
    for ch in st[0]:
        for res in ch:
            if any(a.altloc != '\0' for a in res):
                continue
            atoms = {}
            for a in res:
                n = a.name
                if a.element.name in ('H', 'D') and n.startswith('D'):
                    n = 'H' + n[1:]
                if n not in atoms or a.occ > atoms[n].occ:
                    atoms[n] = a
            out[(ch.name, res.seqid.num, res.seqid.icode)] = (res.name, atoms)
    return out


def dihedral(a, b, c, d):
    return math.degrees(gemmi.calculate_dihedral(a.pos, b.pos, c.pos, d.pos))


def present(atoms, name):
    return name in atoms and atoms[name].occ >= 0.5


def compare(job):
    deposited_path, placed_path = job
    c, errors = collections.Counter(), collections.defaultdict(list)
    try:
        dep_st, out_st = gemmi.read_structure(deposited_path), gemmi.read_structure(placed_path)
    except (RuntimeError, ValueError):
        return c, errors
    dep, out = residues(dep_st), residues(out_st)
    dep_model = dep_st[0]
    ns = gemmi.NeighborSearch(dep_model, gemmi.UnitCell(), 5).populate()
    near = lambda pos, r: (m.to_cra(dep_model) for m in ns.find_atoms(pos, '\0', radius=r))

    for key, (rn, d) in dep.items():
        if key not in out or out[key][0] != rn:
            continue
        o = out[key][1]
        if rn == 'HIS' and {'ND1', 'NE2'} <= set(d) and {'ND1', 'NE2'} <= set(o):
            deposited = (present(d, 'HD1'), present(d, 'HE2'))
            if deposited != (False, False):
                c['his'] += 1
                if o['ND1'].pos.dist(d['ND1'].pos) <= 0.6:
                    c['his ring'] += 1
                    c['his state'] += ('HD1' in o, 'HE2' in o) == deposited
        if rn in ('ASN', 'GLN'):
            n, h = ('ND2', 'HD21') if rn == 'ASN' else ('NE2', 'HE21')
            if n in d and h in d and n in o:
                c[rn.lower()] += 1
                c[rn.lower() + ' agrees'] += o[n].pos.dist(d[n].pos) < 0.6
        if rn in HYDROXYL:
            names = HYDROXYL[rn]
            if all(x in d and x in o for x in names) and present(d, names[3]):
                e = abs((dihedral(*(d[x] for x in names)) - dihedral(*(o[x] for x in names)) + 180) % 360 - 180)
                errors[rn].append(e)
                c[rn.lower()] += 1
                c[rn.lower() + ' within 30'] += e <= 30
                partner, best = 'none', 2.6
                for cra in near(d[names[3]].pos, 2.6):
                    same = cra.chain.name == key[0] and cra.residue.seqid.num == key[1]
                    if cra.atom.element.name in ('N', 'O') and not same:
                        dist = cra.atom.pos.dist(d[names[3]].pos)
                        if dist < best:
                            best, partner = dist, 'water' if cra.residue.name in WATERS else 'protein'
                c['hydroxyl ' + partner] += 1
                c['hydroxyl ' + partner + ' within 30'] += e <= 30
        if rn == 'CYS' and 'SG' in d and 'SG' in o:
            group = 'cys metal' if any(cra.atom.element.is_metal for cra in near(d['SG'].pos, 2.9)) else 'cys other'
            c[group] += 1
            c[group + ' agrees'] += present(d, 'HG') == ('HG' in o)
        if rn in ACID and all(x in d and x in o for x in ACID[rn][0]):
            dep_h = any(present(d, h) for h in ACID[rn][1])
            out_h = any(h in o for h in ACID[rn][1])
            c['acid'] += 1
            c['acid deposited H'] += dep_h
            c['acid placed H'] += out_h
            c['acid both H'] += dep_h and out_h

    # histidines left with no ring hydrogen, over the whole output
    out_model = out_st[0]
    out_ns = gemmi.NeighborSearch(out_model, gemmi.UnitCell(), 5).populate()
    for ch in out_model:
        for res in ch:
            nd1, ne2 = res.find_atom('ND1', '*'), res.find_atom('NE2', '*')
            if res.name != 'HIS' or not (nd1 and ne2):
                continue
            if any(m.to_cra(out_model).atom.element.is_metal for n in (nd1, ne2)
                   for m in out_ns.find_atoms(n.pos, '\0', radius=2.6)):
                continue
            c['his no metal'] += 1
            c['his bare'] += not (res.find_atom('HD1', '*') or res.find_atom('HE2', '*'))
    return c, errors


def report(names, totals, errors):
    def pct(config, num, den):
        n, d = totals[config][num], totals[config][den]
        return f'{100 * n / d:.1f}% of {d}' if d else 'n/a'

    def median(config, rn):
        e = sorted(errors[config][rn])
        return f'{e[len(e) // 2]:.0f}°' if e else 'n/a'

    rows = [
        ('His ring oriented as deposited', lambda k: pct(k, 'his ring', 'his')),
        ('His protonation, where the ring is', lambda k: pct(k, 'his state', 'his ring')),
        ('His not on a metal with no ring H', lambda k: pct(k, 'his bare', 'his no metal')),
        ('Asn amide orientation', lambda k: pct(k, 'asn agrees', 'asn')),
        ('Gln amide orientation', lambda k: pct(k, 'gln agrees', 'gln')),
    ]
    for rn in ('SER', 'THR', 'TYR'):
        rows.append((f'{rn} hydroxyl H within 30° (median error)',
                     lambda k, rn=rn: f"{pct(k, rn.lower() + ' within 30', rn.lower())} ({median(k, rn)})"))
    for partner in ('protein', 'water', 'none'):
        label = {'protein': 'protein or ligand', 'water': 'water', 'none': 'nothing'}[partner]
        rows.append((f'Hydroxyl H within 30°, deposited H bonds to {label}',
                     lambda k, p=partner: pct(k, f'hydroxyl {p} within 30', f'hydroxyl {p}')))
    rows += [
        ('Cys on a metal: thiol H matches', lambda k: pct(k, 'cys metal agrees', 'cys metal')),
        ('Other Cys: thiol H matches', lambda k: pct(k, 'cys other agrees', 'cys other')),
        ('Asp/Glu protonated: deposited / placed / both',
         lambda k: f"{totals[k]['acid deposited H']} / {totals[k]['acid placed H']} / {totals[k]['acid both H']} of {totals[k]['acid']}"),
    ]
    print('| Measure | ' + ' | '.join(names) + ' |')
    print('|---|' + '---:|' * len(names))
    for label, cell in rows:
        print(f'| {label} | ' + ' | '.join(cell(k) for k in names) + ' |')


def placed_file(out_dir, code):
    for suffix in ('FH.cif', 'H.cif'):
        p = os.path.join(out_dir, code + suffix)
        if os.path.exists(p):
            return p
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--pdb-dir', required=True, help='directory of mmCIF files (flat or divided)')
    ap.add_argument('--entries', default=ENTRIES, help='entry list (default: tools/neutron_entries.txt)')
    ap.add_argument('--find', action='store_true', help='rescan --pdb-dir for neutron entries and rewrite --entries')
    ap.add_argument('--config', action='append', default=[], metavar='NAME=ARGS',
                    help='a configuration: reduce3 arguments to add (default: fixed="" and compat="--compat")')
    ap.add_argument('--reduce3', default=os.path.join(ROOT, 'target', 'release', 'reduce3'))
    ap.add_argument('--chem-data', help='passed to reduce3 as --chem-data')
    ap.add_argument('--work', default=os.path.join(os.environ.get('TMPDIR', '/tmp'), 'reduce3_neutron'))
    ap.add_argument('--jobs', type=int, default=os.cpu_count())
    args = ap.parse_args()

    if args.find:
        with Pool(args.jobs) as pool:
            files = list(mmcif_files(args.pdb_dir))
            found = sorted({os.path.basename(f).split('.')[0].lower()
                            for f, hit in zip(files, pool.map(is_neutron, files, chunksize=64)) if hit})
        with open(args.entries, 'w') as f:
            f.write('# PDB entries whose _exptl.method includes NEUTRON DIFFRACTION (joint X-ray/neutron included);\n'
                    '# regenerate with tools/neutron_benchmark.py --find.\n')
            f.write(''.join(code + '\n' for code in found))
        print(f'{len(found)} neutron entries among {len(files)} files, written to {args.entries}')
        return

    codes = [l.strip().lower() for l in open(args.entries) if l.strip() and not l.startswith('#')]
    inputs = {code: entry_path(args.pdb_dir, code) for code in codes}
    missing = [c for c, p in inputs.items() if p is None]
    inputs = {c: p for c, p in inputs.items() if p is not None}
    print(f'{len(inputs)} of {len(codes)} entries found in {args.pdb_dir}'
          + (f' (missing: {" ".join(missing[:10])}{" ..." if len(missing) > 10 else ""})' if missing else ''))
    os.makedirs(args.work, exist_ok=True)
    batch = os.path.join(args.work, 'inputs.txt')
    with open(batch, 'w') as f:
        f.write(''.join(p + '\n' for p in inputs.values()))

    configs = [tuple(c.split('=', 1)) if '=' in c else (c, '') for c in args.config] or [('fixed', ''), ('compat', '--compat')]
    totals, errors = {}, {}
    for name, extra in configs:
        out_dir = os.path.join(args.work, name)
        cmd = [args.reduce3, '--out-dir', out_dir, '--batch', batch, '--no-description', '--jobs', str(args.jobs)]
        if args.chem_data:
            cmd += ['--chem-data', args.chem_data]
        run = subprocess.run(cmd + extra.split(), capture_output=True, text=True)
        summary = [l for l in run.stderr.splitlines() if ' models, ' in l]
        print(f'{name}: {summary[-1] if summary else run.stderr.strip()[-300:]}')
        jobs = [(p, placed_file(out_dir, code)) for code, p in inputs.items()]
        jobs = [j for j in jobs if j[1]]
        totals[name], errors[name] = collections.Counter(), collections.defaultdict(list)
        with Pool(args.jobs) as pool:
            for c, e in pool.imap_unordered(compare, jobs, chunksize=4):
                totals[name].update(c)
                for rn, v in e.items():
                    errors[name][rn] += v
    print()
    report([n for n, _ in configs], totals, errors)


if __name__ == '__main__':
    sys.exit(main())
