#!/usr/bin/env python3
"""Trim the chosen structures into tests/golden/inputs and write structures.tsv.

Each mmCIF file keeps only the categories Reduce3 reads (`_atom_site`,
`_atom_site_anisotrop`, `_cell`, `_symmetry`, `_space_group`,
`_struct_conn`), and of the atom loops only the items it reads; every value
is kept as authored. (The anisotropic ADPs matter: they set the B of their
atoms, and a water's B decides whether it gets phantom hydrogens.) Check that
the trimmed inputs give the golden outputs of the originals (README.md)
before committing them. The golden test checks every input against the
sha256 of its text recorded in structures.tsv.

    python3 tools/golden/prepare_inputs.py selection.tsv tests/golden

Standard library only.
"""

import argparse
import gzip
import hashlib
import io
import re
from pathlib import Path

KEEP = {"_atom_site", "_atom_site_anisotrop", "_cell", "_symmetry", "_space_group", "_struct_conn"}

# The items of the atom loops that `mmcif::structure_from_cif` reads.
ITEMS = {
    "_atom_site": {
        "group_PDB", "id", "type_symbol", "label_atom_id", "label_alt_id", "label_asym_id",
        "auth_asym_id", "auth_comp_id", "auth_seq_id", "pdbx_PDB_ins_code", "Cartn_x", "Cartn_y",
        "Cartn_z", "occupancy", "B_iso_or_equiv", "pdbx_formal_charge", "pdbx_PDB_model_num",
        "auth_segid", "auth_break",
    },
    "_atom_site_anisotrop": {
        "id", "U[1][1]", "U[2][2]", "U[3][3]", "U[1][2]", "U[1][3]", "U[2][3]",
        "B[1][1]", "B[2][2]", "B[3][3]", "B[1][2]", "B[1][3]", "B[2][3]",
    },
}
# Fallbacks the reader uses when an item above is absent.
FALLBACK = {
    "_atom_site": {"label_atom_id": "auth_atom_id", "auth_asym_id": "label_asym_id",
                   "auth_comp_id": "label_comp_id", "auth_seq_id": "label_seq_id"},
}

TOKEN = re.compile(r"'(?:[^']|'(?=\S))*'|\"(?:[^\"]|\"(?=\S))*\"|\S+")


def select_items(category, body):
    """A loop of `category` with only the items the reader uses, one row per line."""
    lines = body.splitlines()
    tags = [l.split()[0] for l in lines if l.startswith(category + ".")]
    names = [t[len(category) + 1:] for t in tags]
    want = set(ITEMS[category])
    fallback = FALLBACK.get(category, {})
    want |= {fallback[n] for n in ITEMS[category] if n in fallback and n not in names}
    keep = [k for k, n in enumerate(names) if n in want]
    rows = [l for l in lines if l.strip() and l.strip() != "loop_" and not l.startswith("_")]
    out = ["loop_"] + [tags[k] for k in keep]
    for row in rows:
        values = TOKEN.findall(row)
        if len(values) != len(names):
            return body  # rows over several lines: keep the loop as it is
        out.append(" ".join(values[k] for k in keep))
    return "\n".join(out) + "\n"


def sections(text):
    """Split one mmCIF data block into (category, text) sections, in order.

    A section is a loop or a run of items of one category. Comment and blank
    lines go with the section before them.
    """
    lines = text.splitlines(keepends=True)
    header, out = [], []
    current = None  # [category, lines, is_loop, rows_started]
    in_text = False
    for line in lines:
        if in_text:
            current[1].append(line)
            if line.startswith(";"):
                in_text = False
            continue
        stripped = line.strip()
        if line.startswith(";"):
            current[1].append(line)
            in_text = True
            continue
        if not stripped or stripped.startswith("#"):
            (current[1] if current else header).append(line)
            continue
        if stripped.startswith("data_"):
            if current or out:
                raise ValueError("expected one data block")
            header.append(line)
            continue
        if stripped == "loop_":
            current = [None, [line], True, False]
            out.append(current)
            continue
        if stripped.startswith("_"):
            category = stripped.split()[0].split(".")[0]
            if current and current[2] and current[0] is None:
                current[0] = category
                current[1].append(line)
                continue
            if current and current[2] and not current[3] and current[0] == category:
                current[1].append(line)
                continue
            if current and not current[2] and current[0] == category:
                current[1].append(line)
                continue
            current = [category, [line], False, False]
            out.append(current)
            continue
        current[1].append(line)
        if current[2]:
            current[3] = True
    return "".join(header), [(c[0], "".join(c[1])) for c in out]


def without_trailing_comments(body):
    lines = body.splitlines(keepends=True)
    while lines and (not lines[-1].strip() or lines[-1].lstrip().startswith("#")):
        lines.pop()
    return "".join(lines)


def trim(text):
    header, parts = sections(text)
    name = header.strip().splitlines()[0]
    kept = []
    for category, body in parts:
        if category not in KEEP:
            continue
        body = without_trailing_comments(body)
        if category in ITEMS and body.lstrip().startswith("loop_"):
            body = select_items(category, body)
        kept.append(body)
    return name + "\n#\n" + "#\n".join(kept) + "#\n"


def write_gz(path, text):
    buf = io.BytesIO()
    with gzip.GzipFile(filename="", mode="wb", fileobj=buf, compresslevel=9, mtime=0) as gz:
        gz.write(text.encode())
    path.write_bytes(buf.getvalue())


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("selection", help="select_structures.py output")
    ap.add_argument("golden", type=Path, help="tests/golden")
    args = ap.parse_args()
    rows = [line.rstrip("\n").split("\t") for line in open(args.selection)]
    head, rows = rows[0], rows[1:]
    col = {n: k for k, n in enumerate(head)}
    inputs = args.golden / "inputs"
    inputs.mkdir(parents=True, exist_ok=True)
    with open(args.golden / "structures.tsv", "w") as out:
        out.write("id\tclass\tclasses\theavy_atoms\tmodels\tsource\tsource_sha256\tinput_sha256\n")
        for r in rows:
            ident, source = r[col["id"]], Path(r[col["source"]])
            data = source.read_bytes()
            if hashlib.sha256(data).hexdigest() != r[col["source_sha256"]]:
                raise SystemExit(f"{source}: changed since selection")
            text = trim(gzip.decompress(data).decode())
            write_gz(inputs / f"{ident}.cif.gz", text)
            origin = f"AFDB {ident}" if r[col["class"]] == "afdb" else f"PDB {ident}"
            out.write("\t".join([
                ident, r[col["class"]], r[col["classes"]], r[col["heavy_atoms"]], r[col["models"]],
                origin, r[col["source_sha256"]], hashlib.sha256(text.encode()).hexdigest(),
            ]) + "\n")


if __name__ == "__main__":
    main()
