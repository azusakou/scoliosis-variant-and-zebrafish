"""Candidate sequence correspondence and saved prediction-response checks.

Live ClinVar, gnomAD and AlphaMissense requests are not made by this replay.
DynaMut2 results are identity-checked, not resubmitted or treated as experiments.
"""

import hashlib
import json
import math
from io_utils import write_csv


def check_candidate_annotations(d, out):
    """Recompute sequence correspondences and verify the four completed saved DynaMut2 responses."""
    from Bio import Align
    from Bio.Align import substitution_matrices

    sequences = d["sequences"]
    alignments = []
    for s in sequences:
        if hashlib.md5(s["sequence"].encode()).hexdigest().upper() != s["md5"].upper():
            raise ValueError("Sequence hash mismatch: " + s["accession"])
    for gene, fishgenes, mutations in [
        ("FLNB", ["flnba", "flnbb"], ["L582P", "Y1294C"]),
        ("FNDC1", ["fndc1"], ["F108S", "Y1880C"]),
    ]:
        h = next((s for s in sequences if s["gene"] == gene))["sequence"]
        for fish in [s for s in sequences if s["gene"] in fishgenes]:
            f = fish["sequence"]
            for op, ex in [(-10.0, -0.5), (-12.0, -1.0)]:
                al = Align.PairwiseAligner()
                al.mode = "global"
                al.substitution_matrix = substitution_matrices.load("BLOSUM62")
                al.open_gap_score = op
                al.extend_gap_score = ex
                z = al.align(h, f)[0]
                mapping = {}
                for (a, b), (c, e) in zip(*z.aligned):
                    mapping.update(
                        {x + 1: y + 1 for x, y in zip(range(a, b), range(c, e))}
                    )
                for v in mutations:
                    p = int(v[1:-1])
                    q = mapping.get(p)
                    local = [
                        k
                        for k in range(max(1, p - 15), min(len(h), p + 15) + 1)
                        if k in mapping
                    ]
                    alignments.append(
                        {
                            "human_gene": gene,
                            "human_variant": v,
                            "fish_gene": fish["gene"],
                            "fish_accession": fish["accession"],
                            "mapped_fish_position": q,
                            "mapped_fish_residue": f[q - 1] if q else "-",
                            "site_conserved": q is not None and h[p - 1] == f[q - 1],
                            "local_31aa_identity_aligned": (
                                sum((h[k - 1] == f[mapping[k] - 1] for k in local))
                                / len(local)
                                if local
                                else None
                            ),
                            "open_gap_score": op,
                            "extend_gap_score": ex,
                        }
                    )
    write_csv(out / "fish_residue_alignments.csv", alignments)
    aa = {"LEU": "L", "PRO": "P", "TYR": "Y", "CYS": "C", "PHE": "F", "SER": "S"}
    results = []
    for r in d["dynamut2"]:
        q = json.loads(r["raw_completed_response"])
        v = r["mutation"]
        raw = float(q["prediction"])
        if (
            q["status"] != "DONE"
            or q["chain"] != "A"
            or int(q["position"]) != int(v[1:-1])
            or (aa[q["wild-type"]] != v[0])
            or (aa[q["mutant"]] != v[-1])
            or (not math.isfinite(raw))
        ):
            raise ValueError("Invalid completed DynaMut2 result")
        results.append(
            {
                "gene": r["gene"],
                "variant": v,
                "raw_ddG_kcal_mol": raw,
                "positive_destabilization_ddG_kcal_mol": -raw,
                "job_id": r["job_id"],
                "input_sha256": r["pdb_sha256"],
                "retrieved_utc": r["retrieved_utc"],
            }
        )
    write_csv(out / "dynamut2_results.csv", results)
    return {
        "DynaMut2_completed": len(results),
        "DynaMut2_above_1": [
            r["gene"] + " " + r["variant"]
            for r in results
            if r["positive_destabilization_ddG_kcal_mol"] > 1
        ],
        "sequence_correspondences": len(alignments),
        "annotation_scope": "Dated frozen gnomAD/ClinVar/AlphaMissense/source-CADD annotations exported as input tables; current database queries are not rerun.",
    }
