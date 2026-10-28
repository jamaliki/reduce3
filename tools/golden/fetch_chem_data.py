#!/usr/bin/env python3
"""Fetch the pinned chem_data files that the tests read.

tests/golden/chem_data.tsv names a conda package (URL and sha256) and the
files of it the tests read, each with its sha256. This downloads the package
(about 1.3 GB), checks it, and writes those files under OUT, checking each:

    python3 tools/golden/fetch_chem_data.py tests/golden/chem_data.tsv OUT
    REDUCE3_CHEM_DATA=OUT REDUCE3_REQUIRE_CHEM_DATA=1 cargo test

`--package FILE` uses a package already downloaded. A `.conda` package is a
zip holding `pkg-*.tar.zst`; decompressing it needs Python 3.14's
`compression.zstd` or the `zstd` program. Otherwise standard library only.
"""

import argparse
import hashlib
import shutil
import subprocess
import sys
import tarfile
import tempfile
import urllib.request
import zipfile
from pathlib import Path

# Where the package installs chem_data (a noarch Python package).
PREFIX = "site-packages/chem_data/"


def read_pins(file):
    meta, files = {}, {}
    for line in open(file):
        line = line.rstrip("\n")
        if line.startswith("# ") and "\t" in line:
            key, value = line[2:].split("\t", 1)
            meta[key] = value
        elif line and not line.startswith("#") and line != "path\tsha256":
            path, sha = line.split("\t")
            files[path] = sha
    return meta["package"], meta["package_sha256"], files


def sha256_file(file):
    h = hashlib.sha256()
    with open(file, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def download(url, dest):
    print(f"downloading {url}", file=sys.stderr)
    with urllib.request.urlopen(url) as r, open(dest, "wb") as f:
        shutil.copyfileobj(r, f, 1 << 20)


def zstd_stream(raw):
    """A readable stream of the zstd-decompressed `raw` file object."""
    try:
        from compression import zstd  # Python 3.14

        return zstd.ZstdFile(raw), None
    except ImportError:
        proc = subprocess.Popen(["zstd", "-dc"], stdin=raw, stdout=subprocess.PIPE)
        return proc.stdout, proc


def extract(package, files, out):
    """Write `files` (paths under chem_data) from the package to `out`."""
    wanted = {PREFIX + p: p for p in files}
    found = set()
    with zipfile.ZipFile(package) as z:
        inner = [n for n in z.namelist() if n.startswith("pkg-") and n.endswith(".tar.zst")]
        if len(inner) != 1:
            raise SystemExit(f"{package}: expected one pkg-*.tar.zst, found {inner}")
        with tempfile.TemporaryDirectory() as tmp:
            # the zstd program reads a real file
            tar_zst = Path(tmp) / inner[0]
            with z.open(inner[0]) as src, open(tar_zst, "wb") as dst:
                shutil.copyfileobj(src, dst, 1 << 20)
            with open(tar_zst, "rb") as raw:
                stream, proc = zstd_stream(raw)
                with tarfile.open(fileobj=stream, mode="r|") as tar:
                    for member in tar:
                        path = wanted.get(member.name)
                        if path is None or not member.isfile():
                            continue
                        dest = out / path
                        dest.parent.mkdir(parents=True, exist_ok=True)
                        with tar.extractfile(member) as src, open(dest, "wb") as dst:
                            shutil.copyfileobj(src, dst)
                        found.add(path)
                if proc is not None and proc.wait() != 0:
                    raise SystemExit("zstd failed")
    missing = sorted(set(files) - found)
    if missing:
        raise SystemExit(f"{len(missing)} files are not in the package, e.g. {missing[:5]}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("pins", type=Path, help="tests/golden/chem_data.tsv")
    ap.add_argument("out", type=Path, help="chem_data directory to write")
    ap.add_argument("--package", type=Path, help="the package, already downloaded")
    args = ap.parse_args()
    url, package_sha, files = read_pins(args.pins)
    with tempfile.TemporaryDirectory() as tmp:
        package = args.package
        if package is None:
            package = Path(tmp) / url.rsplit("/", 1)[1]
            download(url, package)
        got = sha256_file(package)
        if got != package_sha:
            raise SystemExit(f"{package}: sha256 {got}, expected {package_sha}")
        extract(package, files, args.out)
    bad = [p for p, sha in files.items() if sha256_file(args.out / p) != sha]
    if bad:
        raise SystemExit(f"{len(bad)} files differ from their pins, e.g. {bad[:5]}")
    # MonLib::locate looks for these two directories
    for d in ("geostd", "mon_lib"):
        (args.out / d).mkdir(parents=True, exist_ok=True)
    print(f"wrote {len(files)} files to {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
