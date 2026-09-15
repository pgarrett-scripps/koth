"""CLI/export regression: replicated weak candidates cannot invent cell signal.

Usage: python scripts/lfq_confidence_smoke.py target/release/koth_align
All inputs and outputs live in a temporary directory.
"""
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile


def table(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main():
    binary = str(Path(sys.argv[1]).resolve())
    with tempfile.TemporaryDirectory(prefix="lfq-confidence-smoke-") as tmp:
        root = Path(tmp)
        for run in ["A", "B"]:
            folder = root / "batch" / run
            folder.mkdir(parents=True)
            rows = []
            for i in range(200):
                mass = 1000.0 + 2 * i
                rows.append(dict(massCalib=mass, mz=mass / 2 + 1.007276466621,
                    charge=2, rtApex=5.0 + 0.4 * i, rtStart=0.0, rtEnd=100.0, im="",
                    intensityApex=1000.0, intensitySum=5000.0,
                    nIsotopes=3, nScans=5, cosine_score=0.5, isotope_score=0.6,
                    combined_score=0.3, ppm_error=0.0, neutron_offset=0,
                    theoretical_pattern="[0.6,0.3,0.1]", isotope_profile="[]", elution_profile="[]"))
            table(folder / "features.tsv", rows)
            # No measured hill signal, despite convincing cross-run coordinates.
            (folder / "hills.tsv").write_text("mz\n")
        config = root / "config.toml"
        config.write_text('[output]\nexport_long = true\n[lfq.consensus]\nmax_group_qvalue = 0.05\n')
        out = root / "out"
        subprocess.run([binary, str(root / "batch"), "--config", str(config), "--output", str(out)], check=True, capture_output=True)
        manifest = json.loads((out / "matrix_manifest.json").read_text())
        assert manifest["schema_version"] == 3
        assert manifest["n_consensus"] == 200
        for name in ["cells", "observations"]:
            entry = manifest[name]
            assert hashlib.sha256((out / entry["path"]).read_bytes()).hexdigest() == entry["sha256"]
        with (out / "lfq_matrix.long.tsv").open() as handle:
            cells = list(csv.DictReader(handle, delimiter="\t"))
        assert cells
        for cell in cells:
            assert None not in cell, "column/header mismatch"
            assert float(cell["group_qvalue"]) <= 0.05
            assert float(cell["lfq_intensity"]) == 0.0
            assert cell["lfq_status"] == "no_signal"
            assert int(cell["lfq_owned_samples"]) == 0
            assert int(cell["lfq_excluded_samples"]) == 0
            assert cell["lfq_competing_consensus_id"] == ""
            if cell["extraction_kind"] == "decoy":
                assert cell["lfq_q_value"] == ""
            else:
                assert float(cell["lfq_q_value"]) == 1.0
        print(json.dumps({"retained_weak_groups": 200, "cells": len(cells),
                          "invented_signal_cells": 0, "schema_version": 3}))


if __name__ == "__main__":
    main()
