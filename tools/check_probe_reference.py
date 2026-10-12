"""Compare native atom-pair reports with captured standalone CCTBX Probe2 runs."""

import argparse
import json
import subprocess
from pathlib import Path


def reference_key(row, target=False):
    prefix = "target_" if target else ""
    return (row[prefix + "model_id"], row[prefix + "chain_id"],
            int(row[prefix + "resseq"]), row[prefix + "icode"].strip(),
            row[prefix + "resname"], row[prefix + "name"].strip(),
            row[prefix + "altloc"].strip())


def pair_key(first, second):
    return tuple(sorted((tuple(first), tuple(second))))


def contact_key(atom):
    return ("", atom["chainID"], int(atom["resID"]), atom["iCode"].strip(),
            atom["resName"], atom["atomName"].strip(), atom["alt"].strip())


def contact_summary(rows, native=False):
    summary = {}
    for row in rows:
        if native:
            key = (tuple(row["source"]), tuple(row["target"]), row["type"])
            count, gap = row["dot_count"], float(f'{row["min_dot_gap"]:.3f}')
        else:
            key = (contact_key(row["src"]), contact_key(row["target"]), row["type"])
            count, gap = row.get("dotCount", 1), row.get("dotGap", row["gap"])
        prior = summary.setdefault(key, {"count": 0, "gap": gap})
        prior["count"] += count
        prior["gap"] = min(prior["gap"], gap)
    return summary


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference-dir", type=Path, required=True)
    parser.add_argument("--native-bin", type=Path, required=True)
    parser.add_argument("--pattern", action="append")
    args = parser.parse_args()
    rows = []
    files = {file for pattern in (args.pattern or ["*-probe.json"])
             for file in args.reference_dir.glob(pattern)}
    for file in sorted(files):
        captured = json.loads(file.read_text())
        if not captured.get("complete"):
            raise ValueError(f"Incomplete reference capture: {file}")
        result = subprocess.run([str(args.native_bin), str(file)],
                                check=True, capture_output=True, text=True)
        native = json.loads(result.stdout)
        reference = captured["probe2"]
        expected = {pair_key(reference_key(row), reference_key(row, True)): row
                    for row in reference["pairs"]}
        actual = {pair_key(row["first"], row["second"]): row for row in native["pairs"]}
        missing = sorted(expected.keys() - actual.keys())
        extra = sorted(actual.keys() - expected.keys())
        overlap_errors = [{"pair": key, "reference": expected[key]["overlap"],
                           "native": actual[key]["overlap"]}
                          for key in expected.keys() & actual.keys()
                          if abs(expected[key]["overlap"] - actual[key]["overlap"]) > 1e-10]
        detail_errors = [{"pair": key, "reference": expected[key], "native": actual[key]}
                         for key in expected.keys() & actual.keys()
                         if expected[key]["xyz"] != actual[key]["location"] or
                         expected[key]["max_b_factor"] != actual[key]["max_b_factor"]]
        raw = json.loads(reference["program_runs"][0]["raw_output"])["flat_results"]
        expected_contacts = contact_summary(raw)
        actual_contacts = contact_summary(native["contacts"], native=True)
        contact_errors = [{"key": key, "reference": expected_contacts.get(key),
                           "native": actual_contacts.get(key)}
                          for key in expected_contacts.keys() | actual_contacts.keys()
                          if expected_contacts.get(key) != actual_contacts.get(key)]
        rows.append({"case": file.stem,
                     "passed": not missing and not extra and not overlap_errors and
                     not detail_errors and not contact_errors and
                     abs(native["clashscore"] - reference["clashscore"]) < 1e-10,
                     "reference_pair_count": len(expected), "native_pair_count": len(actual),
                     "reference_clashscore": reference["clashscore"],
                     "native_clashscore": native["clashscore"],
                     "phantom_count": native["phantom_count"],
                     "missing": missing, "extra": extra, "overlap_errors": overlap_errors,
                     "reference_contacts": len(expected_contacts),
                     "native_contacts": len(actual_contacts), "contact_errors": contact_errors,
                     "detail_errors": detail_errors})
    if not rows:
        raise ValueError("No captured Probe2 cases found")
    print(json.dumps({"complete": all(row["passed"] for row in rows), "cases": rows}, indent=2))
    return 0 if all(row["passed"] for row in rows) else 1


if __name__ == "__main__":
    raise SystemExit(main())
