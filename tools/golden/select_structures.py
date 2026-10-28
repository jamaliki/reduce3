#!/usr/bin/env python3
"""Choose the golden structures from local PDB and AlphaFold DB mirrors.

The choice is a seeded random sample, so the same mirrors give the same set.
Each candidate's `_atom_site` loop is read (first model only for the counts)
and the entry is filed under every class it belongs to; each class then takes
its quota in sample order, an entry counting once.

    python3 tools/golden/select_structures.py \
        --pdb ~/data/pdb --afdb ~/data/kaveh/afdb-v6/structures \
        --out selection.tsv

`prepare_inputs.py` then trims the chosen files into tests/golden/inputs.
Standard library only.
"""

import argparse
import gzip
import hashlib
import os
import random
import re
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

SEED = 20261028

# How many entries each class takes. An entry is filed under every class it
# fits; the classes are filled in this order.
QUOTAS = [
    ("large", 1),
    ("nmr", 10),
    ("nucleic", 30),
    ("metal_cluster", 10),
    ("metal", 35),
    ("covalent", 20),
    ("modified", 15),
    ("altloc", 35),
    ("ligand", 55),
    ("protein", 30),
    ("afdb", 60),
]

# Size caps for everything but "large": heavy atoms of the first model, and
# atom rows of all models (deposited hydrogens included).
MAX_ATOMS = 4000
MAX_ROWS = 8000
LARGE_ATOMS = (60000, 160000)

WATER = {"HOH", "DOD", "WAT"}
AMINO = set(
    "ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split()
)
NUCLEIC = set("A C G U I DA DC DG DT DI DU".split())
MODIFIED = set("MSE SEP TPO PTR HYP CSO CME KCX LLP MLY OCS CSD PCA SAH".split())
METALS = set(
    "LI NA K RB CS BE MG CA SR BA AL GA IN TL SN PB BI SC TI V CR MN FE CO NI CU ZN Y ZR NB MO "
    "TC RU RH PD AG CD HF TA W RE OS IR PT AU HG LA CE PR ND SM EU GD TB DY HO ER TM YB LU U".split()
)
CLUSTERS = {"SF4", "FES", "F3S", "HEM", "HEC", "CLF", "FS4", "CFM", "ICS", "MOS", "MGD", "B12", "COB", "CUA", "PCD"}
IONS_ONLY = METALS | {"F", "CL", "BR", "IOD", "SO4", "PO4", "NO3", "ACT", "GOL", "EDO", "PEG", "FMT", "UNX", "UNL"}

TOKEN = re.compile(r"'(?:[^']|'(?=\S))*'|\"[^\"]*\"|\S+")


def atom_site(text):
    """Header names and token rows of the `_atom_site` loop."""
    start = text.find("\n_atom_site.")
    if start < 0:
        return [], []
    lines = text[start + 1 :].split("\n")
    names, i = [], 0
    while i < len(lines) and lines[i].startswith("_atom_site."):
        names.append(lines[i].split()[0][len("_atom_site.") :])
        i += 1
    rows = []
    for line in lines[i:]:
        if not line or line.startswith("#") or line.startswith("_") or line.startswith("loop_"):
            break
        rows.append(TOKEN.findall(line))
    return names, rows


def struct_conn_types(text):
    """The `conn_type_id` values of `_struct_conn` (a loop or one item)."""
    if "_struct_conn.conn_type_id" not in text:
        return set()
    return {t for t in ("covale", "metalc", "disulf") if re.search(r"\b%s\b" % t, text)}


def classify(path):
    """Counts and classes of one entry, or None when it cannot be read."""
    try:
        data = Path(path).read_bytes()
        text = gzip.decompress(data).decode("utf-8", "replace")
    except (OSError, EOFError, gzip.BadGzipFile):
        return None
    names, rows = atom_site(text)
    if not rows:
        return None
    col = {n: k for k, n in enumerate(names)}
    need = ["group_PDB", "type_symbol", "label_alt_id", "auth_comp_id", "pdbx_PDB_model_num"]
    if any(n not in col for n in need):
        return None
    first_model = rows[0][col["pdbx_PDB_model_num"]]
    models, comps, metals, altloc, heavy, hydrogens = set(), set(), set(), False, 0, 0
    for r in rows:
        if len(r) != len(names):
            return None
        models.add(r[col["pdbx_PDB_model_num"]])
        if r[col["pdbx_PDB_model_num"]] != first_model:
            continue
        element = r[col["type_symbol"]].upper()
        if element in ("H", "D"):
            hydrogens += 1
            continue
        heavy += 1
        comp = r[col["auth_comp_id"]].strip("\"'").upper()
        comps.add(comp)
        if element in METALS:
            metals.add(element)
        if r[col["label_alt_id"]] not in (".", "?"):
            altloc = True
    method = re.search(r"_exptl\.method\s+['\"]?([A-Z ]+)", text)
    method = method.group(1).strip() if method else ""
    conn = struct_conn_types(text)
    classes = []
    nucleic = bool(comps & NUCLEIC)
    ligands = comps - WATER - AMINO - NUCLEIC - MODIFIED - IONS_ONLY
    if LARGE_ATOMS[0] <= heavy <= LARGE_ATOMS[1] and len(models) == 1:
        classes.append("large")
    if heavy <= MAX_ATOMS and len(rows) <= MAX_ROWS:
        if len(models) > 1 and "NMR" in method:
            classes.append("nmr")
        if nucleic:
            classes.append("nucleic")
        if comps & CLUSTERS:
            classes.append("metal_cluster")
        if metals:
            classes.append("metal")
        if "covale" in conn:
            classes.append("covalent")
        if comps & MODIFIED:
            classes.append("modified")
        if altloc:
            classes.append("altloc")
        if ligands:
            classes.append("ligand")
        if comps <= AMINO | WATER and len(models) == 1:
            classes.append("protein")
    return dict(
        heavy=heavy, hydrogens=hydrogens, rows=len(rows), models=len(models), method=method,
        metals=",".join(sorted(metals)), ligands=",".join(sorted(ligands)[:8]),
        classes=classes, sha256=hashlib.sha256(data).hexdigest(),
    )


def pdb_candidates(root, rng, n):
    entries = []
    for sub in sorted(os.listdir(root)):
        d = Path(root) / sub
        if d.is_dir():
            entries.extend(sorted(str(d / f) for f in os.listdir(d) if f.endswith(".cif.gz")))
    rng.shuffle(entries)
    return entries[:n]


def afdb_candidates(root, rng, n):
    subs = sorted(s for s in os.listdir(root) if (Path(root) / s).is_dir())
    out = []
    for sub in rng.sample(subs, min(len(subs), n)):
        files = sorted(f for f in os.listdir(Path(root) / sub) if f.endswith(".cif.gz"))
        out.append(str(Path(root) / sub / rng.choice(files)))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pdb", required=True, help="PDB mmCIF mirror, divided layout (xy/1xyz.cif.gz)")
    ap.add_argument("--afdb", required=True, help="AlphaFold DB mirror, divided layout (xy/AF-*.cif.gz)")
    ap.add_argument("--pdb-sample", type=int, default=6000, help="PDB entries to examine")
    ap.add_argument("--out", required=True)
    ap.add_argument("--workers", type=int, default=os.cpu_count())
    args = ap.parse_args()
    rng = random.Random(SEED)
    pdb = pdb_candidates(args.pdb, rng, args.pdb_sample)
    afdb = afdb_candidates(args.afdb, rng, 200)
    with ProcessPoolExecutor(args.workers) as pool:
        pdb_info = list(pool.map(classify, pdb, chunksize=16))
        afdb_info = list(pool.map(classify, afdb, chunksize=16))
    chosen, taken = [], set()
    for cls, quota in QUOTAS:
        pool_ = zip(afdb, afdb_info) if cls == "afdb" else zip(pdb, pdb_info)
        got = 0
        for path, info in pool_:
            if got == quota:
                break
            if info is None or path in taken:
                continue
            if cls == "afdb":
                if info["heavy"] > MAX_ATOMS or info["rows"] > MAX_ROWS:
                    continue
                info["classes"] = ["afdb"]
            elif cls not in info["classes"]:
                continue
            taken.add(path)
            chosen.append((cls, path, info))
            got += 1
        print(f"{cls}: {got} of {quota}", file=sys.stderr)
    with open(args.out, "w") as out:
        out.write("id\tclass\tclasses\theavy_atoms\tmodels\tmetals\tligands\tsource\tsource_sha256\n")
        for cls, path, info in chosen:
            ident = Path(path).name.removesuffix(".cif.gz")
            out.write("\t".join([
                ident, cls, ",".join(info["classes"]), str(info["heavy"]), str(info["models"]),
                info["metals"] or ".", info["ligands"] or ".", path, info["sha256"],
            ]) + "\n")


if __name__ == "__main__":
    main()
