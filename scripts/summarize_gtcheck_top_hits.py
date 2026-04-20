#!/usr/bin/env python3
"""
Summarize one or more bcftools gtcheck TSV outputs.

The script is intentionally species- and filename-agnostic. It reads gtcheck
files produced for any SNP-chip or marker-panel VCF and ranks the best panel
matches for each query sample. It uses only the Python standard library so it
can run on most HPC login nodes without extra package installation.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import math
import sys
from collections import defaultdict
from pathlib import Path
from typing import Iterable, Optional


TOP_HIT_COLUMNS = [
    "query_sample",
    "query_sample_basename",
    "panel_sample",
    "rank",
    "discordance",
    "discordance_per_site",
    "avg_neg_log_hwe",
    "sites_compared",
    "matching_genotypes",
    "match_fraction",
    "confidence_score",
    "source_file",
]


SAMPLE_SUMMARY_COLUMNS = [
    "query_sample",
    "query_sample_basename",
    "top_panel_sample",
    "top_discordance",
    "top_sites_compared",
    "top_match_fraction",
    "top_confidence_score",
    "second_panel_sample",
    "second_match_fraction",
    "match_fraction_gap_rank1_rank2",
    "max_sites_compared_any_hit",
    "top_sites_compared_fraction_of_panel",
    "max_sites_compared_fraction_of_panel",
    "query_vcf",
    "vcf_total_sites",
    "vcf_called_sites",
    "vcf_missing_sites",
    "vcf_partial_or_non_diploid_sites",
    "vcf_call_rate",
    "vcf_missing_rate",
    "vcf_het_sites",
    "vcf_hom_ref_sites",
    "vcf_hom_alt_sites",
    "vcf_het_rate_among_called",
    "query_vcf_qc_error",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Rank top matches from bcftools gtcheck output TSV files.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--input",
        nargs="+",
        default=["results/*.gtcheck.tsv"],
        help=(
            "One or more gtcheck TSV files. Shell globs are accepted when quoted, "
            "for example 'results/*.gtcheck.tsv'."
        ),
    )
    parser.add_argument(
        "-o",
        "--output",
        default="results/gtcheck_top_hits.tsv",
        help="Ranked top-hit TSV to write.",
    )
    parser.add_argument("-n", "--top-n", type=int, default=10, help="Hits per query sample.")
    parser.add_argument(
        "--min-sites",
        type=int,
        default=0,
        help="Minimum compared sites required for a hit to be ranked.",
    )
    parser.add_argument(
        "--keep-zero-sites",
        action="store_true",
        help="Keep gtcheck rows with zero compared sites. These rank last.",
    )
    parser.add_argument(
        "--sample-summary",
        default=None,
        help=(
            "Optional one-row-per-query output. Defaults to '<output stem>.sample_summary.tsv'. "
            "Use 'none' to disable."
        ),
    )
    parser.add_argument(
        "--panel-size",
        type=int,
        default=None,
        help="Total marker count in the panel, used only for coverage fractions.",
    )
    parser.add_argument(
        "--query-vcf-dir",
        default=None,
        help="Optional directory containing filtered single-sample query VCFs for QC metrics.",
    )
    parser.add_argument(
        "--query-vcf-suffix",
        default=".panel.filtered.diploid.vcf.gz",
        help="Suffix used to match a query sample basename to its query VCF.",
    )
    return parser.parse_args()


def expand_inputs(patterns: Iterable[str]) -> list[Path]:
    paths: list[Path] = []
    for pattern in patterns:
        if any(ch in pattern for ch in "*?["):
            matches = sorted(Path().glob(pattern))
        else:
            matches = [Path(pattern)]
        paths.extend(matches)
    return list(dict.fromkeys(paths))


def safe_float(value: str) -> Optional[float]:
    try:
        return float(value)
    except ValueError:
        return None


def safe_int(value: str) -> Optional[int]:
    try:
        return int(float(value))
    except ValueError:
        return None


def fmt(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        if math.isnan(value):
            return ""
        return f"{value:.8g}"
    return str(value)


def read_gtcheck(path: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    malformed = 0
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#") or line.startswith("INFO"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 7 or fields[0] != "DCv2":
                continue

            discordance = safe_int(fields[3])
            avg_neg_log_hwe = safe_float(fields[4])
            sites_compared = safe_int(fields[5])
            matching_genotypes = safe_int(fields[6])
            if None in {discordance, avg_neg_log_hwe, sites_compared, matching_genotypes}:
                malformed += 1
                continue

            assert discordance is not None
            assert avg_neg_log_hwe is not None
            assert sites_compared is not None
            assert matching_genotypes is not None

            if sites_compared > 0:
                discordance_per_site: Optional[float] = discordance / sites_compared
                match_fraction: Optional[float] = matching_genotypes / sites_compared
                confidence_score: Optional[float] = match_fraction * math.log10(sites_compared)
            else:
                discordance_per_site = None
                match_fraction = None
                confidence_score = None

            rows.append(
                {
                    "query_sample": fields[1],
                    "query_sample_basename": Path(fields[1]).name,
                    "panel_sample": fields[2],
                    "discordance": discordance,
                    "discordance_per_site": discordance_per_site,
                    "avg_neg_log_hwe": avg_neg_log_hwe,
                    "sites_compared": sites_compared,
                    "matching_genotypes": matching_genotypes,
                    "match_fraction": match_fraction,
                    "confidence_score": confidence_score,
                    "source_file": str(path),
                }
            )

    if malformed:
        print(f"WARNING: dropped {malformed} malformed row(s) from {path}", file=sys.stderr)
    if not rows:
        raise ValueError(f"No DCv2 rows found in {path}")
    return rows


def rank_hits(
    rows: list[dict[str, object]],
    top_n: int,
    min_sites: int,
    keep_zero_sites: bool,
) -> list[dict[str, object]]:
    grouped: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        sites = int(row["sites_compared"])
        if sites == 0 and not keep_zero_sites:
            continue
        if min_sites > 0 and sites < min_sites:
            continue
        grouped[str(row["query_sample"])].append(row)

    if not grouped:
        raise ValueError("No rows remain after filtering. Lower --min-sites or inspect the gtcheck files.")

    ranked_rows: list[dict[str, object]] = []
    for query_sample in grouped:
        ordered = sorted(
            grouped[query_sample],
            key=lambda row: (
                float("inf") if row["discordance_per_site"] is None else row["discordance_per_site"],
                row["discordance"],
                -float(row["avg_neg_log_hwe"]),
                -int(row["matching_genotypes"]),
            ),
        )
        for rank, row in enumerate(ordered[:top_n], start=1):
            ranked = dict(row)
            ranked["rank"] = rank
            ranked_rows.append(ranked)
    return ranked_rows


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def qc_from_vcf(path: Path) -> dict[str, object]:
    total = called = missing = het = hom_ref = hom_alt = partial_or_non_diploid = 0
    with open_text(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.rstrip("\n").split("\t")
                if len(header) < 10:
                    raise ValueError(f"VCF has no sample column: {path}")
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            total += 1
            fmt_keys = fields[8].split(":")
            sample_values = fields[9].split(":")
            fmt_values = dict(zip(fmt_keys, sample_values))
            gt = fmt_values.get("GT", ".").replace("|", "/")
            if gt in {".", "./."}:
                missing += 1
                continue
            alleles = gt.split("/")
            if len(alleles) != 2 or "." in alleles:
                partial_or_non_diploid += 1
                missing += 1
                continue
            called += 1
            if alleles[0] != alleles[1]:
                het += 1
            elif alleles[0] == "0":
                hom_ref += 1
            else:
                hom_alt += 1

    return {
        "query_vcf": str(path),
        "vcf_total_sites": total,
        "vcf_called_sites": called,
        "vcf_missing_sites": missing,
        "vcf_partial_or_non_diploid_sites": partial_or_non_diploid,
        "vcf_call_rate": called / total if total else None,
        "vcf_missing_rate": missing / total if total else None,
        "vcf_het_sites": het,
        "vcf_hom_ref_sites": hom_ref,
        "vcf_hom_alt_sites": hom_alt,
        "vcf_het_rate_among_called": het / called if called else None,
    }


def build_sample_summary(
    ranked_rows: list[dict[str, object]],
    panel_size: Optional[int],
    query_vcf_dir: Optional[Path],
    query_vcf_suffix: str,
) -> list[dict[str, object]]:
    grouped: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in ranked_rows:
        grouped[str(row["query_sample"])].append(row)

    summaries: list[dict[str, object]] = []
    for query_sample in grouped:
        rows = sorted(grouped[query_sample], key=lambda row: int(row["rank"]))
        best = rows[0]
        second = rows[1] if len(rows) > 1 else None
        max_sites = max(int(row["sites_compared"]) for row in rows)
        summary: dict[str, object] = {
            "query_sample": query_sample,
            "query_sample_basename": best["query_sample_basename"],
            "top_panel_sample": best["panel_sample"],
            "top_discordance": best["discordance"],
            "top_sites_compared": best["sites_compared"],
            "top_match_fraction": best["match_fraction"],
            "top_confidence_score": best["confidence_score"],
            "second_panel_sample": second["panel_sample"] if second else None,
            "second_match_fraction": second["match_fraction"] if second else None,
            "match_fraction_gap_rank1_rank2": (
                float(best["match_fraction"]) - float(second["match_fraction"])
                if second and best["match_fraction"] is not None and second["match_fraction"] is not None
                else None
            ),
            "max_sites_compared_any_hit": max_sites,
            "top_sites_compared_fraction_of_panel": (
                int(best["sites_compared"]) / panel_size if panel_size else None
            ),
            "max_sites_compared_fraction_of_panel": max_sites / panel_size if panel_size else None,
        }

        if query_vcf_dir:
            candidate = query_vcf_dir / f"{best['query_sample_basename']}{query_vcf_suffix}"
            if candidate.exists():
                try:
                    summary.update(qc_from_vcf(candidate))
                except Exception as exc:
                    summary["query_vcf"] = str(candidate)
                    summary["query_vcf_qc_error"] = str(exc)
            else:
                summary["query_vcf"] = None
        summaries.append(summary)
    return summaries


def write_tsv(path: Path, rows: list[dict[str, object]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            delimiter="\t",
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({column: fmt(row.get(column)) for column in columns})


def main() -> int:
    args = parse_args()
    inputs = expand_inputs(args.input)
    if not inputs:
        print("ERROR: no gtcheck input files matched", file=sys.stderr)
        return 1

    rows: list[dict[str, object]] = []
    for path in inputs:
        if not path.exists():
            print(f"WARNING: missing file: {path}", file=sys.stderr)
            continue
        try:
            rows.extend(read_gtcheck(path))
        except Exception as exc:
            print(f"WARNING: {exc}", file=sys.stderr)

    if not rows:
        print("ERROR: no valid gtcheck files were parsed", file=sys.stderr)
        return 1

    try:
        ranked_rows = rank_hits(rows, args.top_n, args.min_sites, args.keep_zero_sites)
    except ValueError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    output = Path(args.output)
    write_tsv(output, ranked_rows, TOP_HIT_COLUMNS)
    print(f"Wrote: {output}")

    if args.sample_summary != "none":
        summary_path = (
            Path(args.sample_summary)
            if args.sample_summary
            else output.with_name(f"{output.stem}.sample_summary.tsv")
        )
        query_vcf_dir = Path(args.query_vcf_dir) if args.query_vcf_dir else None
        summaries = build_sample_summary(ranked_rows, args.panel_size, query_vcf_dir, args.query_vcf_suffix)
        write_tsv(summary_path, summaries, SAMPLE_SUMMARY_COLUMNS)
        print(f"Wrote: {summary_path}")

    current_query = None
    for row in ranked_rows:
        if row["query_sample"] != current_query:
            current_query = row["query_sample"]
            print(f"\n=== {current_query} ===")
            print("rank\tpanel_sample\tdiscordance_per_site\tsites_compared\tmatch_fraction\tconfidence_score")
        print(
            "\t".join(
                [
                    fmt(row["rank"]),
                    fmt(row["panel_sample"]),
                    fmt(row["discordance_per_site"]),
                    fmt(row["sites_compared"]),
                    fmt(row["match_fraction"]),
                    fmt(row["confidence_score"]),
                ]
            )
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
