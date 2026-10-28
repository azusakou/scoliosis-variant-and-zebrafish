#!/usr/bin/env python3
"""Reproduce core downstream results; network downloads are optional."""
import argparse
import datetime
import json
import shutil
import subprocess
from pathlib import Path

from io_utils import DATA_SHA256, PACKAGE_ROOT, load_payload, table_rows, write_csv
from sequence_annotations import check_candidate_annotations
from structure_analysis import (
    fetch_public_sources,
    reconstruct_structures,
    summarize_assemblies,
)
from variant_analysis import summarize_variants


def main():
    """Run independent analytical modules in a fresh output directory."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--fetch-public", choices=["assemblies", "all"])
    args = parser.parse_args()
    output = args.output.resolve()
    if output == PACKAGE_ROOT or PACKAGE_ROOT in output.parents:
        raise SystemExit("Choose a fresh output directory outside this package.")
    if output.exists() and any(output.iterdir()):
        raise SystemExit("Refusing to mix results into a nonempty output directory.")
    if not shutil.which("Rscript"):
        raise SystemExit(
            "Rscript is required (base R only); no R packages need installation."
        )

    data = load_payload("data.json.gz", DATA_SHA256)
    output.mkdir(parents=True, exist_ok=True)
    inputs = output / "inputs"
    inputs.mkdir()
    for name, table in data["tables"].items():
        write_csv(inputs / (name + ".csv"), table_rows(data, name), table["columns"])

    summary = {
        "scope": data["scope"],
        "started_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "molecular": summarize_variants(data, output),
        "assemblies": summarize_assemblies(data, output),
        "candidates": check_candidate_annotations(data, output),
        "coordinates": reconstruct_structures(output),
    }

    mechanics = output / "mechanics"
    mechanics.mkdir()
    with (output / "R_analysis.log").open("w") as log:
        subprocess.run(
            ["Rscript", str(PACKAGE_ROOT / "mechanics.R"), str(inputs), str(mechanics)],
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True,
        )
    if args.fetch_public:
        summary["public_files_retrieved_and_hash_verified"] = fetch_public_sources(
            data, output, args.fetch_public
        )
    summary["finished_utc"] = datetime.datetime.now(datetime.timezone.utc).isoformat()
    summary["public_upload_performed"] = False
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    print(
        "Completed core downstream reproduction. Frozen molecular predictions were reaggregated, not newly inferred."
    )


if __name__ == "__main__":
    main()
