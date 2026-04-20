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

Example Usage:
    python summarize_gtcheck_top_hits.py \
      -i benning.tsv clark.tsv jack.tsv \
      -o all_top10.tsv

Wildcard:
    python summarize_gtcheck_top_hits.py \
      -i gtcheck_outputs/*.tsv \
      -o all_top10.tsv

Change number of hits:
    python summarize_gtcheck_top_hits.py \
      -i *.tsv \
      -n 20 \
      -o top20.tsv
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
        help="One or more bcftools gtcheck output TSV files"
    )

    parser.add_argument(
        "-n", "--top-n",
        type=int,
        default=10,
        help="Number of best hits to report per query sample"
    )

    parser.add_argument(
        "-o", "--output",
        default="gtcheck_top_hits.tsv",
        help="Output TSV filename"
    )

    parser.add_argument(
        "--min-sites",
        type=int,
        default=0,
        help="Minimum number of compared sites required to keep a hit"
    )

    return parser.parse_args()


def read_gtcheck_file(path: Path) -> pd.DataFrame:
    rows: List[List[str]] = []

    with path.open("r") as f:
        for line in f:
            if line.startswith("#") or line.startswith("INFO") or not line.strip():
                continue

            fields = line.strip().split("\t")
            if len(fields) < 7 or fields[0] != "DCv2":
                continue

            rows.append(fields[:7])

    if not rows:
        raise ValueError(f"No DCv2 records found in {path}")

    df = pd.DataFrame(rows, columns=EXPECTED_COLUMNS)

    df["discordance"] = pd.to_numeric(df["discordance"], errors="coerce")
    df["avg_neg_log_hwe"] = pd.to_numeric(df["avg_neg_log_hwe"], errors="coerce")
    df["sites_compared"] = pd.to_numeric(df["sites_compared"], errors="coerce")
    df["matching_genotypes"] = pd.to_numeric(df["matching_genotypes"], errors="coerce")

    df["discordance_per_site"] = df["discordance"] / df["sites_compared"]
    df["match_fraction"] = df["matching_genotypes"] / df["sites_compared"]
    df["source_file"] = str(path)
    df["query_sample_basename"] = df["query_sample"].map(lambda x: Path(str(x)).name)

    return df


def rank_hits(df: pd.DataFrame, top_n: int, min_sites: int) -> pd.DataFrame:
    if min_sites > 0:
        df = df[df["sites_compared"] >= min_sites]

    ranked_frames = []

    for query_sample, sub in df.groupby("query_sample", sort=False):
        ranked = sub.sort_values(
            by=[
                "discordance_per_site",
                "discordance",
                "avg_neg_log_hwe",
                "matching_genotypes",
            ],
            ascending=[True, True, False, False],
        ).copy()

        ranked["rank"] = range(1, len(ranked) + 1)
        ranked_frames.append(ranked.head(top_n))

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
    top_hits = rank_hits(combined, args.top_n, args.min_sites)

    top_hits.to_csv(args.output, sep="\t", index=False)

    print(f"\nWrote: {args.output}\n")

    for q, sub in top_hits.groupby("query_sample"):
        print(f"=== {q} ===")
        print(sub[[
            "rank",
            "genotyped_sample",
            "discordance_per_site",
            "avg_neg_log_hwe",
            "sites_compared"
        ]].to_string(index=False))
        print()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())