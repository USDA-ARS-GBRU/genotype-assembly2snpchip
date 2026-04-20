#!/usr/bin/env python3
"""
Summarize one or more bcftools gtcheck output TSV files and optionally attach
single-sample VCF QC metrics for each query sample.

Expected bcftools gtcheck output rows:
DCv2    <query_sample>    <genotyped_sample>    <discordance>
        <avg_neg_log_hwe>    <sites_compared>    <matching_genotypes>

Primary ranking:
  1. discordance_per_site ascending
  2. discordance ascending
  3. avg_neg_log_hwe descending
  4. matching_genotypes descending

Optional extras:
- confidence_score = match_fraction * log10(sites_compared)
- per-sample QC summary from the single-sample filtered VCFs used as gtcheck queries
"""
from __future__ import annotations

import argparse
import gzip
import math
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

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
        description=(
            "Summarize top matches from bcftools gtcheck output files and optionally "
            "attach per-sample query VCF QC metrics."
        ),
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
        help="Output TSV filename for ranked top hits.",
    )
    parser.add_argument(
        "--sample-summary-output",
        default=None,
        help=(
            "Optional TSV filename for one-row-per-query summary. If omitted, a file "
            "named '<output stem>.sample_summary.tsv' is written next to --output."
        ),
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
    parser.add_argument(
        "--panel-size",
        type=int,
        default=None,
        help=(
            "Total number of panel markers. When provided, report each query sample's "
            "maximum sites_compared as a fraction of the panel."
        ),
    )
    parser.add_argument(
        "--vcf-dir",
        default=None,
        help=(
            "Directory containing single-sample filtered query VCFs. If provided, the "
            "script will try to find one VCF per query sample and compute call rate, "
            "missing rate, and heterozygosity metrics."
        ),
    )
    parser.add_argument(
        "--vcf-suffix",
        default=".soy50k.filtered.vcf.gz",
        help="Suffix used to match query sample basename to its VCF in --vcf-dir.",
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

    before = len(df)
    df = df.dropna(subset=["discordance", "avg_neg_log_hwe", "sites_compared", "matching_genotypes"]).copy()
    dropped = before - len(df)
    if dropped:
        print(f"WARNING: dropped {dropped} malformed row(s) from {path}", file=sys.stderr)

    denom = df["sites_compared"].replace(0, pd.NA)
    df["discordance_per_site"] = df["discordance"] / denom
    df["match_fraction"] = df["matching_genotypes"] / denom
    # Slightly penalize tiny site counts while remaining simple to explain.
    df["confidence_score"] = df["match_fraction"] * df["sites_compared"].map(
        lambda x: math.log10(x) if pd.notna(x) and x and x > 0 else float("nan")
    )
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

    ranked_frames: List[pd.DataFrame] = []
    for query_sample, sub in work.groupby("query_sample", sort=False):
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
            "confidence_score",
            "source_file",
        ]
    ]


def open_textmaybe_gzip(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open("r")


def normalize_gt(gt: str) -> Tuple[str, ...]:
    gt = gt.replace("|", "/")
    parts = tuple(gt.split("/"))
    return parts


def qc_from_single_sample_vcf(path: Path) -> Dict[str, object]:
    total = called = missing = het = hom_ref = hom_alt = non_diploid = 0
    with open_textmaybe_gzip(path) as f:
        sample_col_idx: Optional[int] = None
        for line in f:
            if not line.strip():
                continue
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.rstrip("\n").split("\t")
                if len(header) < 10:
                    raise ValueError(f"VCF appears to have no sample column: {path}")
                sample_col_idx = 9
                continue
            if sample_col_idx is None:
                continue
            total += 1
            fields = line.rstrip("\n").split("\t")
            fmt_keys = fields[8].split(":")
            sample_vals = fields[sample_col_idx].split(":")
            fmt = dict(zip(fmt_keys, sample_vals))
            gt = fmt.get("GT", ".")
            if gt in {".", "./.", ".|."}:
                missing += 1
                continue
            alleles = normalize_gt(gt)
            if len(alleles) != 2 or any(a == "." for a in alleles):
                non_diploid += 1
                missing += 1
                continue
            called += 1
            a, b = alleles
            if a != b:
                het += 1
            elif a == "0":
                hom_ref += 1
            else:
                hom_alt += 1

    call_rate = called / total if total else float("nan")
    missing_rate = missing / total if total else float("nan")
    het_rate_called = het / called if called else float("nan")
    return {
        "query_vcf": str(path),
        "vcf_total_sites": total,
        "vcf_called_sites": called,
        "vcf_missing_sites": missing,
        "vcf_non_diploid_or_partial_missing_sites": non_diploid,
        "vcf_call_rate": call_rate,
        "vcf_missing_rate": missing_rate,
        "vcf_het_sites": het,
        "vcf_hom_ref_sites": hom_ref,
        "vcf_hom_alt_sites": hom_alt,
        "vcf_het_rate_among_called": het_rate_called,
    }


def find_query_vcf(vcf_dir: Path, query_sample_basename: str, vcf_suffix: str) -> Optional[Path]:
    candidate = vcf_dir / f"{query_sample_basename}{vcf_suffix}"
    if candidate.exists():
        return candidate
    return None


def build_sample_summary(
    top_hits: pd.DataFrame,
    panel_size: Optional[int],
    vcf_dir: Optional[Path],
    vcf_suffix: str,
) -> pd.DataFrame:
    records = []
    for query_sample, sub in top_hits.groupby("query_sample", sort=False):
        best = sub.sort_values("rank").iloc[0]
        max_sites_compared = int(sub["sites_compared"].max()) if len(sub) else 0
        rec: Dict[str, object] = {
            "query_sample": query_sample,
            "query_sample_basename": best["query_sample_basename"],
            "top_hit": best["genotyped_sample"],
            "top_rank": int(best["rank"]),
            "top_discordance": int(best["discordance"]),
            "top_sites_compared": int(best["sites_compared"]),
            "top_matching_genotypes": int(best["matching_genotypes"]),
            "top_match_fraction": float(best["match_fraction"]),
            "top_confidence_score": float(best["confidence_score"]),
            "max_sites_compared_any_hit": max_sites_compared,
        }
        if panel_size:
            rec["max_sites_compared_fraction_of_panel"] = max_sites_compared / panel_size
            rec["top_sites_compared_fraction_of_panel"] = int(best["sites_compared"]) / panel_size
        if vcf_dir is not None:
            vcf_path = find_query_vcf(vcf_dir, str(best["query_sample_basename"]), vcf_suffix)
            if vcf_path is None:
                rec["query_vcf"] = pd.NA
            else:
                try:
                    rec.update(qc_from_single_sample_vcf(vcf_path))
                except Exception as e:
                    rec["query_vcf"] = str(vcf_path)
                    rec["vcf_qc_error"] = str(e)
        records.append(rec)
    return pd.DataFrame(records)


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

    summary_out = (
        Path(args.sample_summary_output)
        if args.sample_summary_output
        else Path(args.output).with_name(f"{Path(args.output).stem}.sample_summary.tsv")
    )
    vcf_dir = Path(args.vcf_dir) if args.vcf_dir else None
    sample_summary = build_sample_summary(top_hits, args.panel_size, vcf_dir, args.vcf_suffix)
    sample_summary.to_csv(summary_out, sep="\t", index=False)

    print(f"\nWrote: {args.output}")
    print(f"Wrote: {summary_out}\n")

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
                    "confidence_score",
                ]
            ].to_string(index=False)
        )
        print()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
