"""Export CCTBX's accepted ionic scattering labels for typed charge parity."""

from pathlib import Path
from cctbx.eltbx.xray_scattering import standard_labels_list

labels = sorted(label.upper() for label in standard_labels_list() if label.endswith(("+", "-")))
root = Path(__file__).resolve().parents[1]
body = "//! Ionic scattering labels from CCTBX 2026.9; see docs/CCTBX_LICENSE.txt.\n\n"
body += "pub(super) const IONIC_LABELS: &[&str] = &[\n"
body += "".join(f'    "{label}",\n' for label in labels)
body += "];\n"
(root / "src/interp/validation/charge_labels.rs").write_text(body)
