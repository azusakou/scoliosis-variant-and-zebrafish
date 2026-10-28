"""Frozen stability-score aggregation and exact threshold/null sensitivities.

These calculations do not rerun ThermoMPNN or EvoEF inference. NMR conformers
are averaged within variants before variant-level summaries are calculated.
"""

import collections
import math
import statistics
from io_utils import parse_boolean, table_rows, write_csv


def average_ranks(x):
    """Return average ranks, including tied observations."""
    order = sorted(range(len(x)), key=x.__getitem__)
    r = [None] * len(x)
    i = 0
    while i < len(x):
        j = i + 1
        while j < len(x) and x[order[j]] == x[order[i]]:
            j += 1
        for k in order[i:j]:
            r[k] = (i + j + 1) / 2
        i = j
    return r


def hypergeometric_probabilities(N, a, b):
    """Conditional overlap probabilities with both marginal list sizes fixed."""
    return {
        k: math.comb(a, k) * math.comb(N - a, b - k) / math.comb(N, b)
        for k in range(max(0, a + b - N), min(a, b) + 1)
    }


def summarize_variants(d, out):
    """Aggregate frozen model scores and recompute eligibility and overlap tests."""
    records = table_rows(d, "variants_121")
    scores = table_rows(d, "model_scores")
    by = collections.defaultdict(list)
    if len(records) != 121 or len({r["source_record_id"] for r in records}) != 121:
        raise ValueError("Variant inventory is incomplete or duplicated")
    for r in scores:
        if int(float(r["position"])) != int(r["canonical_variant"][1:-1]):
            raise ValueError("Canonical position mismatch")
        eligible = r["mapping_status"] in ("verified", "isoform_consensus") and (
            r["structural_confidence"].startswith("experimental NMR")
            or float(r["plddt"]) >= 70
        )
        if parse_boolean(r["primary_eligible"]) != eligible:
            raise ValueError("Eligibility does not match recorded rule")
        e, ep, shift = [
            float(r[k])
            for k in [
                "EvoEF_ddG_repaired_WT",
                "EvoEF_ddG_paired_local_WT",
                "local_WT_minus_repaired_WT",
            ]
        ]
        if not math.isclose(e, ep + shift, abs_tol=1e-09):
            raise ValueError("WT reference energy identity fails")
        by[r["gene"], r["variant"], r["canonical_variant"]].append(r)
    variants = []
    for key, rs in sorted(by.items()):
        estimates = [
            statistics.mean((float(r[k]) for r in rs))
            for k in [
                "ThermoMPNN_ddG_kcal_mol",
                "EvoEF_ddG_kcal_mol",
                "EvoEF_ddG_paired_local_WT",
            ]
        ]
        if not all((math.isfinite(x) for x in estimates)):
            raise ValueError("Nonfinite stability estimate")
        variants.append(
            dict(zip(["gene", "variant", "canonical_variant"], key))
            | dict(zip(["ThermoMPNN", "EvoEF", "EvoEF_paired_local_WT"], estimates))
            | {
                "n_models": len(rs),
                "primary_eligible": all(
                    (parse_boolean(r["primary_eligible"]) for r in rs)
                ),
            }
        )
    primary = [r for r in variants if r["primary_eligible"]]
    if len(variants) != 102 or len(primary) != 56:
        raise ValueError("Scored/primary denominator mismatch")
    sensitivity = []
    for cutoff in [0.5, 1.0, 1.5]:
        for reference, evo in [
            ("globally_repaired_WT", "EvoEF"),
            ("paired_local_WT", "EvoEF_paired_local_WT"),
        ]:
            N = len(primary)
            a = sum((r["ThermoMPNN"] > cutoff for r in primary))
            b = sum((r[evo] > cutoff for r in primary))
            sel = [r for r in primary if r["ThermoMPNN"] > cutoff and r[evo] > cutoff]
            null = hypergeometric_probabilities(N, a, b)
            strat = {0: 1.0}
            expected = 0
            for gene in sorted({r["gene"] for r in primary}):
                g = [r for r in primary if r["gene"] == gene]
                na = sum((r["ThermoMPNN"] > cutoff for r in g))
                nb = sum((r[evo] > cutoff for r in g))
                expected += na * nb / len(g)
                dist = hypergeometric_probabilities(len(g), na, nb)
                joint = collections.defaultdict(float)
                for x, px in strat.items():
                    for y, py in dist.items():
                        joint[x + y] += px * py
                strat = joint
            sensitivity.append(
                {
                    "cutoff_kcal_mol": cutoff,
                    "reference": reference,
                    "N": N,
                    "ThermoMPNN_positive": a,
                    "EvoEF_positive": b,
                    "observed_overlap": len(sel),
                    "expected_overlap_independent_label_null": a * b / N,
                    "overlap_tail_probability": sum(
                        (p for k, p in null.items() if k >= len(sel))
                    ),
                    "gene_stratified_expected_overlap": expected,
                    "gene_stratified_overlap_tail_probability": sum(
                        (p for k, p in strat.items() if k >= len(sel))
                    ),
                    "selected_variants": ";".join(
                        (r["gene"] + " " + r["variant"] for r in sel)
                    ),
                }
            )
    write_csv(out / "stability_summary_102.csv", variants)
    write_csv(out / "stability_threshold_sensitivity.csv", sensitivity)
    write_csv(
        out / "variant_selection_flow.csv",
        [
            {"selection_status": k, "n": v}
            for k, v in collections.Counter(
                (r["review_selection_status"] for r in records)
            ).items()
        ],
    )
    return {
        "records": 121,
        "scored": len(variants),
        "primary": len(primary),
        "model_rows": len(scores),
        "Spearman_default": statistics.correlation(
            average_ranks([r["ThermoMPNN"] for r in primary]),
            average_ranks([r["EvoEF"] for r in primary]),
        ),
    }
