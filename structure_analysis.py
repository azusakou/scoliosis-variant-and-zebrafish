"""Assembly summaries, coordinate reconstruction and optional checked downloads.

Frozen interface records remain distinct from newly retrievable public files;
this module does not authenticate historical atom coordinates or recalculate SASA.
"""

import collections
import hashlib
import statistics
import urllib.request
from io_utils import (
    COORDINATES_SHA256,
    load_payload,
    parse_boolean,
    table_rows,
    write_csv,
)


def summarize_assemblies(d, out):
    """Aggregate frozen interfaces and audit exact variant coverage, without atom-level recalculation."""
    interfaces = table_rows(d, "interfaces")
    primary = [r for r in interfaces if parse_boolean(r["primary_interface"])]
    panel_pairs = [r for r in primary if r["gene_a"] and r["gene_b"]]
    keep = [
        r
        for r in panel_pairs
        if all((len(set(r[k].split(";"))) == 1 for k in ["accession_a", "accession_b"]))
    ]
    by = collections.defaultdict(list)
    for r in keep:
        by[tuple(sorted([r["gene_a"], r["gene_b"]]))].append(r)
    edges = []
    for (a, b), rs in sorted(by.items()):
        v = [float(r["buried_surface_A2"]) for r in rs]
        edges.append(
            {
                "from": a,
                "to": b,
                "n_assemblies": len({r["pdb"] for r in rs}),
                "mean_BSA_A2": statistics.mean(v),
                "min_BSA_A2": min(v),
                "max_BSA_A2": max(v),
                "pdbs": ";".join(sorted({r["pdb"] for r in rs})),
            }
        )
    write_csv(out / "assembly_edges.csv", edges)
    threshold = []
    for area in [50, 100, 200, 500]:
        for distance in [4, 5, 6]:
            q = [
                r
                for r in keep
                if float(r["buried_surface_A2"]) >= area
                and int(r[f"contacts_{distance}A"]) >= 3
            ]
            pairs = {tuple(sorted([r["gene_a"], r["gene_b"]])) for r in q}
            threshold.append(
                {
                    "area_A2": area,
                    "distance_A": distance,
                    "n_chain_pairs": len(q),
                    "n_gene_edges": len(pairs),
                    "n_heterotypic_edges": sum((a != b for a, b in pairs)),
                }
            )
    write_csv(out / "assembly_threshold_sensitivity.csv", threshold)
    cov = {
        (r["gene"], int(r["position"]), r["aa"], r["uniprot_aa"]): r["pdbs"]
        for r in table_rows(d, "coverage")
    }
    coverage = []
    for r in table_rows(d, "variants_121"):
        p = r["canonical_position"]
        aa = r["canonical_variant"][:1]
        pdbs = cov.get((r["gene"], int(float(p)), aa, aa), "") if p else ""
        coverage.append(
            {
                "record_id": r["source_record_id"],
                "gene": r["gene"],
                "variant": r["protein_variant"],
                "canonical_variant": r["canonical_variant"],
                "covering_pdbs": pdbs,
                "n_covering_assemblies": len(pdbs.split(";")) if pdbs else 0,
            }
        )
    write_csv(out / "variant_assembly_coverage.csv", coverage)
    return {
        "panel": len(table_rows(d, "panel_35")),
        "assemblies": len(table_rows(d, "assemblies_70")),
        "contact_entries": len({r["pdb"] for r in interfaces}),
        "primary_chain_pairs": len(primary),
        "retained_chain_pairs": len(keep),
        "excluded_chimeric_pairs": len(panel_pairs) - len(keep),
        "protein_pairs": len(by),
        "heterotypic_pairs": sum((a != b for a, b in by)),
        "exactly_covered_variants": [
            r["gene"] + " " + r["variant"]
            for r in coverage
            if r["n_covering_assemblies"]
        ],
        "surface_max_relative_change": max(
            (
                abs(float(r["relative_change"]))
                for r in table_rows(d, "surface_resolution")
            )
        ),
    }


def reconstruct_structures(out):
    """Losslessly rebuild eight constructed candidate PDBs and verify their original SHA256 values."""
    p = load_payload("coordinates.json.gz", COORDINATES_SHA256)
    dest = out / "constructed_structures"
    dest.mkdir()
    verified = []
    for name, r in p["models"].items():
        a = p["bases"][r["gene"]].copy()
        for i, j, replacement in reversed(r["edits"]):
            a[i:j] = replacement
        b = (
            r["prefix"]
            + "".join((s[:6] + f"{i:5d}" + s[6:] for i, s in enumerate(a, 1)))
        ).encode()
        if hashlib.sha256(b).hexdigest() != r["sha256"]:
            raise ValueError("Coordinate reconstruction hash fails: " + name)
        (dest / name).write_bytes(b)
        verified.append(name)
    return {"losslessly_reconstructed_PDBs": verified, "count": len(verified)}


def fetch_public_sources(d, out, which):
    """Optionally retrieve public resources, stopping if any expected checksum differs."""
    dest = out / "public_sources"
    dest.mkdir()
    selected = [
        r
        for r in d["public_inputs"]
        if which == "all" or r["name"].startswith("assemblies/")
    ]
    for r in selected:
        p = dest / r["name"]
        p.parent.mkdir(parents=True, exist_ok=True)
        with urllib.request.urlopen(r["url"], timeout=60) as response:
            b = response.read()
        if hashlib.sha256(b).hexdigest() != r["sha256"]:
            raise ValueError(
                "Public source changed; exact snapshot not retrieved: " + r["url"]
            )
        p.write_bytes(b)
    return len(selected)
