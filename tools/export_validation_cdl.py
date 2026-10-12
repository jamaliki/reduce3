"""Export pinned CCTBX CDL values without changing their precision."""

import hashlib
import importlib.metadata
import importlib
import json
import struct
from pathlib import Path

source = importlib.import_module("mmtbx.conformation_dependent_library.cdl_database")

groups = sorted(source.cdl_database)
root = Path(__file__).resolve().parents[1] / "src/interp/validation"
output = root / "cdl_targets.bin"
rows = bytearray()
for group in groups:
    for phi in range(-180, 180, 10):
        for psi in range(-180, 180, 10):
            record = source.cdl_database[group][(phi, psi)]
            assert len(record) == 26
            rows.extend(struct.pack("<24d", *record[2:]))
output.write_bytes(rows)
metadata = dict(
    source="CCTBX mmtbx.conformation_dependent_library.cdl_database",
    cctbx_version=importlib.metadata.version("cctbx-base"),
    database_version=source.version, groups=groups,
    source_sha256=hashlib.sha256(Path(source.__file__).read_bytes()).hexdigest(),
    binary_sha256=hashlib.sha256(rows).hexdigest(),
    encoding="8 groups x 36 phi bins x 36 psi bins x 24 little-endian f64 values",
    source_license="CCTBX BSD license; see CCTBX project LICENSE.txt",
)
(root / "cdl_targets.provenance.json").write_text(json.dumps(metadata, indent=2) + "\n")
print(json.dumps(metadata, indent=2))
