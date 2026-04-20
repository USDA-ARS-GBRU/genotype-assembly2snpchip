#!/usr/bin/env python3
"""
summarize_gtcheck_top_hits.py

Parse one or more bcftools gtcheck output TSV files and report the top N best hits
for each query sample.

Expected bcftools gtcheck output format includes lines like:
DCv2    <query_sample>    <genotyped_sample>    <discordance>    <avg_neg_log_hwe>    <sites_compared>    <matching_genotypes>

Ranking:
  1. discordance_per_site ascending (smaller is better)
  2. discordance ascending
  3. avg_neg_log_hwe descending (bigger is better)
  4. matching_genotypes descending

Notes:
- Rows with sites_compared == 0 are excluded from ranking by default because
  discordance_per_site and match_fraction are undefined for those rows.
- Use --keep-zero-sites if you want such rows retained in the output table.
- In GT-vs-GT mode with bcftools gtcheck -u GT,GT -E 0, discordance is
  interpretable as the number of mismatching genotypes.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import List

import pandas as pd


EXPECTED_COLUMNS = [
    "record_type",
    "query_sample",
    "genotyped_sample",
    "discordance",
    "avg_neg_log_hwe",
    "sites_compared",
    "matching_genotypes",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Summarize top matches from bcftools gtcheck output files.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "-i", "--input",
        nargs="+",
        required=True,
        help="One or more bcftools gtcheck output TSV files.",
    )

    parser.add_argument(
        "-n", "--top-n",
        type=int,
        default=10,
        help="Number of best hits to report per query sample.",
    )

    parser.add_argument(
        "-o", "--output",
        default="gtcheck_top_hits.tsv",
        help="Output TSV filename.",
    )

    parser.add_argument(
        "--min-sites",
        type=int,
        default=0,
        help="Minimum number of compared sites required to keep a hit.",
    )

    parser.add_argument(
        "--keep-zero-sites",
        action="store_true",
        help=(
            "Keep rows with sites_compared == 0 in the output. By default these are "
            "dropped before ranking because per-site metrics are undefined."
        ),
    )

    return parser.parse_args()


def read_gtcheck_file(path: Path) -> pd.DataFrame:
    rows: List[List[str]] = []

    with path.open("r") as f:
        for line in f:
            if line.startswith("#") or line.startswith("INFO") or not line.strip():
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 7 or fields[0] != "DCv2":
                continue

            rows.append(fields[:7])

    if not rows:
        raise ValueError(f"No DCv2 records found in {path}")

    df = pd.DataFrame(rows, columns=EXPECTED_COLUMNS)

    for col in ["discordance", "avg_neg_log_hwe", "sites_compared", "matching_genotypes"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    # Drop rows with malformed numeric content.
    before = len(df)
    df = df.dropna(subset=["discordance", "avg_neg_log_hwe", "sites_compared", "matching_genotypes"]).copy()
    dropped = before - len(df)
    if dropped:
        print(f"WARNING: dropped {dropped} malformed row(s) from {path}", file=sys.stderr)

    # Compute safe per-site metrics.
    denom = df["sites_compared"].replace(0, pd.NA)
    df["discordance_per_site"] = df["discordance"] / denom
    df["match_fraction"] = df["matching_genotypes"] / denom
    df["source_file"] = str(path)
    df["query_sample_basename"] = df["query_sample"].map(lambda x: Path(str(x)).name)

    return df


def rank_hits(df: pd.DataFrame, top_n: int, min_sites: int, keep_zero_sites: bool) -> pd.DataFrame:
    work = df.copy()

    if not keep_zero_sites:
        work = work[work["sites_compared"] > 0].copy()

    if min_sites > 0:
        work = work[work["sites_compared"] >= min_sites].copy()

    if work.empty:
        raise ValueError(
            "No rows remain after filtering. Try lowering --min-sites or use "
            "--keep-zero-sites if you explicitly want zero-site rows retained."
        )

    ranked_frames = []

    for query_sample, sub in work.groupby("query_sample", sort=False):
        # For safety, sort any NA per-site values to the end.
        ranked = sub.sort_values(
            by=[
                "discordance_per_site",
                "discordance",
                "avg_neg_log_hwe",
                "matching_genotypes",
            ],
            ascending=[True, True, False, False],
            na_position="last",
        ).copy()

        ranked["rank"] = range(1, len(ranked) + 1)
        ranked_frames.append(ranked.head(top_n))

    if not ranked_frames:
        raise ValueError("No per-sample ranking tables could be generated.")

    out = pd.concat(ranked_frames, ignore_index=True)

    return out[
        [
            "query_sample",
            "query_sample_basename",
            "genotyped_sample",
            "rank",
            "discordance",
            "discordance_per_site",
            "avg_neg_log_hwe",
            "sites_compared",
            "matching_genotypes",
            "match_fraction",
            "source_file",
        ]
    ]


def main() -> int:
    args = parse_args()

    frames = []
    for f in args.input:
        path = Path(f)
        if not path.exists():
            print(f"WARNING: missing file: {path}", file=sys.stderr)
            continue
        try:
            frames.append(read_gtcheck_file(path))
        except Exception as e:
            print(f"WARNING: {e}", file=sys.stderr)

    if not frames:
        print("ERROR: no valid input files", file=sys.stderr)
        return 1

    combined = pd.concat(frames, ignore_index=True)

    try:
        top_hits = rank_hits(
            combined,
            top_n=args.top_n,
            min_sites=args.min_sites,
            keep_zero_sites=args.keep_zero_sites,
        )
    except ValueError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1

    top_hits.to_csv(args.output, sep="\t", index=False)

    print(f"\nWrote: {args.output}\n")

    for q, sub in top_hits.groupby("query_sample", sort=False):
        print(f"=== {q} ===")
        print(
            sub[
                [
                    "rank",
                    "genotyped_sample",
                    "discordance_per_site",
                    "avg_neg_log_hwe",
                    "sites_compared",
                    "match_fraction",
                ]
            ].to_string(index=False)
        )
        print()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
