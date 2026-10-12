"""Test charged contact rules against direct CCTBX Probe, without interpretation.

Synthetic +1 donor hydrogens exercise the charged cutoff. CCTBX's model
interpretation clears these charges, so this control does not test full
structure chemistry or the default validation wrapper.
"""

import argparse
import json
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path


def oracle(distance, donor_charge, acceptor_charge):
    import iotbx.pdb
    import mmtbx_probe_ext as probe

    root = iotbx.pdb.hierarchy.root()
    model = iotbx.pdb.hierarchy.model()
    root.append_model(model)
    atoms = []
    rows = []
    infos = []
    for index, (name, element, radius, charge, donor, acceptor, x) in enumerate([
        ("H", "H", 1.05, donor_charge, True, False, 0.0),
        ("O", "O", 1.52, acceptor_charge, False, True, distance),
    ]):
        chain = iotbx.pdb.hierarchy.chain(id=str(index))
        model.append_chain(chain)
        residue = iotbx.pdb.hierarchy.residue_group(resseq="1", icode="")
        chain.append_residue_group(residue)
        group = iotbx.pdb.hierarchy.atom_group(resname="LIG", altloc="")
        residue.append_atom_group(group)
        atom = iotbx.pdb.hierarchy.atom()
        atom.name = name
        atom.element = element
        atom.xyz = (x, 0.0, 0.0)
        atom.occ = 1.0
        atom.b = 10.0
        group.append_atom(atom)
        atoms.append(atom)
        infos.append(probe.ExtraAtomInfo(radius, acceptor, donor, False, False, charge))
        rows.append({"name": name, "element": element, "xyz": [x, 0.0, 0.0],
                     "occupancy": 1.0, "b_iso": 10.0, "hetero": True,
                     "model": "", "chain": str(index), "resseq": "1",
                     "component": "LIG", "icode": "", "altloc": ""})
    root.atoms().reset_i_seq()
    scorer = probe.DotScorer(probe.ExtraAtomInfoMap(atoms, infos))
    cache = probe.DotSphereCache(16.0)
    expected = {}
    for index, atom in enumerate(atoms):
        for dot in scorer.check_dots(atom, cache.get_sphere(infos[index].vdwRadius).dots(),
                                     .25, [atoms[1-index]], [], .5):
            kind = scorer.interaction_type(dot.overlapType, dot.gap, False)
            key = (str(index), str(1-index), scorer.interaction_type_short_name(kind))
            item = expected.setdefault(key, {"count": 0, "gap": dot.gap})
            item["count"] += 1
            item["gap"] = min(item["gap"], dot.gap)
    metadata = [{"vdw_radius": info.vdwRadius, "is_acceptor": info.isAcceptor,
                 "is_donor": info.isDonor, "is_dummy_hydrogen": False,
                 "is_ion": False, "charge": info.charge, "alt": ""} for info in infos]
    captured = {"nuclear": False, "report_vdws": True, "probe2": {"initial_atom_info_calls": [
        {"atoms": rows, "extra_atom_info": metadata, "bonds": []}]}}
    return captured, expected


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--native-bin", type=Path, required=True)
    args = parser.parse_args()
    cases = []
    for donor_charge, acceptor_charge in [(0, 0), (1, -1), (1, 1)]:
        for overlap in [-.3, -.1, .599999, .600001, .7, .799999, .800001, 1.2]:
            capture, expected = oracle(2.57 - overlap, donor_charge, acceptor_charge)
            with tempfile.TemporaryDirectory() as temporary:
                file = Path(temporary) / "synthetic-probe-input.json"
                file.write_text(json.dumps(capture))
                output = subprocess.run([str(args.native_bin), str(file)],
                                        check=True, capture_output=True, text=True)
                native = json.loads(output.stdout)
            actual = defaultdict(lambda: {"count": 0, "gap": float("inf")})
            for row in native["contacts"]:
                key = (row["source"][1], row["target"][1], row["type"])
                actual[key]["count"] += row["dot_count"]
                actual[key]["gap"] = min(actual[key]["gap"], row["min_dot_gap"])
            passed = expected.keys() == actual.keys() and all(
                row["count"] == actual[key]["count"] and
                abs(row["gap"] - actual[key]["gap"]) < 1e-12
                for key, row in expected.items())
            cases.append({"donor_charge": donor_charge, "acceptor_charge": acceptor_charge,
                          "overlap": overlap, "passed": passed})
    print(json.dumps({"complete": all(row["passed"] for row in cases), "cases": cases}, indent=2))
    return 0 if all(row["passed"] for row in cases) else 1


if __name__ == "__main__":
    raise SystemExit(main())
