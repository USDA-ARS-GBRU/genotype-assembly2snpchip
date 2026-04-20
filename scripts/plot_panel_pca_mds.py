#!/usr/bin/env python3
"""
Create PCA and/or MDS figures from a SNP-chip panel VCF plus query VCFs.

The script reads VCF/VCF.gz files directly in Python, converts diploid GT calls
to alternate-allele dosage values (0, 1, 2), imputes remaining missing values by
marker mean, and plots assemblies in the context of the SNP-chip panel.
"""
from __future__ import annotations

import argparse
import gzip
import math
import random
from pathlib import Path
from typing import Iterable, Optional

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.manifold import MDS
from sklearn.preprocessing import StandardScaler


PLOT_STYLE = {
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.color": "#e6e6e6",
    "grid.linewidth": 0.8,
    "axes.edgecolor": "#444444",
    "axes.labelcolor": "#222222",
    "font.size": 10,
    "axes.titlesize": 13,
    "axes.titleweight": "bold",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot PCA/MDS from a SNP-chip panel VCF and assembly-derived query VCFs.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--panel-vcf", required=True, help="Multi-sample SNP-chip panel VCF/VCF.gz.")
    parser.add_argument(
        "--query-vcfs",
        nargs="*",
        default=[],
        help="One or more single-sample assembly-derived query VCF/VCF.gz files.",
    )
    parser.add_argument(
        "--panel-samples",
        default=None,
        help="Optional file with panel sample names to include, one per line.",
    )
    parser.add_argument(
        "--max-panel-samples",
        type=int,
        default=None,
        help="Randomly downsample panel samples to this count before PCA/MDS.",
    )
    parser.add_argument(
        "--max-sites",
        type=int,
        default=None,
        help="Randomly downsample markers to this count after filtering.",
    )
    parser.add_argument(
        "--min-site-call-rate",
        type=float,
        default=0.8,
        help="Minimum nonmissing genotype fraction required to keep a marker.",
    )
    parser.add_argument("--min-maf", type=float, default=0.01, help="Minimum minor allele frequency.")
    parser.add_argument(
        "--metadata",
        default=None,
        help=(
            "Optional TSV with at least a sample column. Optional columns such as group, label, "
            "and source are used for plotting."
        ),
    )
    parser.add_argument("--sample-column", default="sample", help="Sample-name column in --metadata.")
    parser.add_argument("--group-column", default="group", help="Metadata column for point color.")
    parser.add_argument("--label-column", default="label", help="Metadata column for text labels.")
    parser.add_argument("--out-dir", default="figures", help="Directory for output files.")
    parser.add_argument("--prefix", default="panel_pca", help="Output filename prefix.")
    parser.add_argument(
        "--method",
        choices=["pca", "mds", "both"],
        default="pca",
        help="Dimensional reduction method to run.",
    )
    parser.add_argument(
        "--max-mds-samples",
        type=int,
        default=500,
        help="Refuse MDS above this sample count unless you raise the limit.",
    )
    parser.add_argument(
        "--formats",
        nargs="+",
        default=["png", "svg"],
        choices=["png", "svg", "pdf"],
        help="Figure formats to write.",
    )
    parser.add_argument("--dpi", type=int, default=300, help="DPI for PNG output.")
    parser.add_argument("--seed", type=int, default=13, help="Random seed for downsampling.")
    return parser.parse_args()


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def site_key(fields: list[str]) -> str:
    return f"{fields[0]}:{fields[1]}:{fields[3]}:{fields[4]}"


def gt_to_dosage(sample_field: str, fmt_keys: list[str]) -> float:
    values = sample_field.split(":")
    fmt = dict(zip(fmt_keys, values))
    gt = fmt.get("GT", ".")
    gt = gt.replace("|", "/")
    if gt in {".", "./.", ".|."}:
        return np.nan
    alleles = gt.split("/")
    if any(allele == "." for allele in alleles):
        return np.nan
    try:
        allele_values = [int(allele) for allele in alleles]
    except ValueError:
        return np.nan
    if len(allele_values) == 1:
        allele_values = allele_values * 2
    if len(allele_values) != 2:
        return np.nan
    return float(sum(1 for allele in allele_values if allele > 0))


def read_sample_list(path: Optional[Path]) -> Optional[set[str]]:
    if path is None:
        return None
    with path.open("r", encoding="utf-8") as handle:
        return {line.strip() for line in handle if line.strip() and not line.startswith("#")}


def choose_indices(names: list[str], keep: Optional[set[str]], max_count: Optional[int], seed: int) -> list[int]:
    indices = [idx for idx, name in enumerate(names) if keep is None or name in keep]
    if max_count is not None and len(indices) > max_count:
        rng = random.Random(seed)
        indices = sorted(rng.sample(indices, max_count))
    return indices


def read_panel_vcf(
    path: Path,
    keep_samples: Optional[set[str]],
    max_panel_samples: Optional[int],
    seed: int,
) -> tuple[list[str], list[str], np.ndarray]:
    sample_names: list[str] = []
    sample_indices: list[int] = []
    site_keys: list[str] = []
    columns: list[list[float]] = []

    with open_text(path) as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.rstrip("\n").split("\t")
                sample_names_all = header[9:]
                selected = choose_indices(sample_names_all, keep_samples, max_panel_samples, seed)
                sample_indices = selected
                sample_names = [sample_names_all[idx] for idx in selected]
                continue
            if not line.strip() or not sample_names:
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            if "," in fields[4]:
                continue
            fmt_keys = fields[8].split(":")
            dosage = [gt_to_dosage(fields[9 + idx], fmt_keys) for idx in sample_indices]
            site_keys.append(site_key(fields))
            columns.append(dosage)

    if not columns:
        raise ValueError(f"No biallelic genotype records parsed from {path}")
    matrix = np.asarray(columns, dtype=np.float32).T
    return sample_names, site_keys, matrix


def query_sample_name(path: Path, header_samples: list[str]) -> str:
    if len(header_samples) == 1:
        return header_samples[0]
    name = path.name
    for suffix in [".vcf.gz", ".vcf", ".bcf"]:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return path.stem


def read_query_vcf(path: Path, panel_keys: list[str]) -> tuple[str, np.ndarray]:
    panel_index = {key: idx for idx, key in enumerate(panel_keys)}
    vector = np.full(len(panel_keys), np.nan, dtype=np.float32)
    sample_name: Optional[str] = None

    with open_text(path) as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.rstrip("\n").split("\t")
                sample_name = query_sample_name(path, header[9:])
                continue
            if not line.strip() or sample_name is None:
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10 or "," in fields[4]:
                continue
            key = site_key(fields)
            idx = panel_index.get(key)
            if idx is None:
                continue
            vector[idx] = gt_to_dosage(fields[9], fields[8].split(":"))

    if sample_name is None:
        raise ValueError(f"No VCF header with sample found in {path}")
    return sample_name, vector


def marker_filter(matrix: np.ndarray, min_call_rate: float, min_maf: float) -> np.ndarray:
    called = np.isfinite(matrix)
    call_rate = called.mean(axis=0)
    means = np.nanmean(matrix, axis=0)
    allele_freq = means / 2.0
    maf = np.minimum(allele_freq, 1.0 - allele_freq)
    keep = (call_rate >= min_call_rate) & np.isfinite(maf) & (maf >= min_maf)
    if not np.any(keep):
        raise ValueError("No markers remain after call-rate and MAF filtering")
    return keep


def impute_and_standardize(matrix: np.ndarray) -> np.ndarray:
    work = matrix.astype(np.float32, copy=True)
    means = np.nanmean(work, axis=0)
    inds = np.where(~np.isfinite(work))
    work[inds] = np.take(means, inds[1])
    scaler = StandardScaler(with_mean=True, with_std=True)
    scaled = scaler.fit_transform(work)
    scaled = np.nan_to_num(scaled, nan=0.0, posinf=0.0, neginf=0.0)
    return scaled


def read_metadata(path: Optional[Path], sample_column: str) -> pd.DataFrame:
    if path is None:
        return pd.DataFrame(columns=[sample_column])
    meta = pd.read_csv(path, sep="\t")
    if sample_column not in meta.columns:
        raise ValueError(f"Metadata file must contain sample column: {sample_column}")
    return meta


def attach_metadata(
    coords: pd.DataFrame,
    metadata: pd.DataFrame,
    sample_column: str,
    group_column: str,
    label_column: str,
) -> pd.DataFrame:
    out = coords.copy()
    if not metadata.empty:
        out = out.merge(metadata, left_on="sample", right_on=sample_column, how="left")
    if group_column not in out.columns:
        out[group_column] = out["source"]
    out[group_column] = out[group_column].fillna(out["source"])
    if label_column not in out.columns:
        out[label_column] = ""
    out[label_column] = out[label_column].fillna("")
    return out


def save_figure(fig: plt.Figure, out_dir: Path, prefix: str, stem: str, formats: Iterable[str], dpi: int) -> list[Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = []
    for fmt in formats:
        path = out_dir / f"{prefix}_{stem}.{fmt}"
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
        paths.append(path)
    plt.close(fig)
    return paths


def plot_coordinates(
    coords: pd.DataFrame,
    x_col: str,
    y_col: str,
    title: str,
    group_column: str,
    label_column: str,
    x_label: str,
    y_label: str,
) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(9.5, 7.2))
    sns.scatterplot(
        data=coords,
        x=x_col,
        y=y_col,
        hue=group_column,
        style="source",
        s=70,
        linewidth=0.5,
        edgecolor="white",
        alpha=0.9,
        ax=ax,
    )
    query_rows = coords[coords["source"] == "query"]
    for _, row in query_rows.iterrows():
        label = str(row.get(label_column) or row["sample"])
        ax.annotate(label, xy=(row[x_col], row[y_col]), xytext=(5, 5), textcoords="offset points", fontsize=8)
    ax.set_title(title)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False)
    return fig


def run_pca(matrix: np.ndarray, samples: list[str], sources: list[str]) -> tuple[pd.DataFrame, np.ndarray]:
    n_components = min(3, matrix.shape[0], matrix.shape[1])
    if n_components < 2:
        raise ValueError("PCA requires at least two samples and two markers")
    model = PCA(n_components=n_components, random_state=13)
    coords = model.fit_transform(matrix)
    data = {
        "sample": samples,
        "source": sources,
        "PC1": coords[:, 0],
        "PC2": coords[:, 1],
    }
    if n_components >= 3:
        data["PC3"] = coords[:, 2]
    return pd.DataFrame(data), model.explained_variance_ratio_


def run_mds(matrix: np.ndarray, samples: list[str], sources: list[str], seed: int) -> pd.DataFrame:
    model = MDS(n_components=2, metric=True, normalized_stress="auto", random_state=seed, n_init=4)
    coords = model.fit_transform(matrix)
    return pd.DataFrame({"sample": samples, "source": sources, "MDS1": coords[:, 0], "MDS2": coords[:, 1]})


def main() -> int:
    args = parse_args()
    sns.set_theme(style="whitegrid", rc=PLOT_STYLE)

    keep_samples = read_sample_list(Path(args.panel_samples) if args.panel_samples else None)
    panel_samples, keys, panel_matrix = read_panel_vcf(
        Path(args.panel_vcf),
        keep_samples=keep_samples,
        max_panel_samples=args.max_panel_samples,
        seed=args.seed,
    )

    query_names: list[str] = []
    query_vectors: list[np.ndarray] = []
    for query_vcf in args.query_vcfs:
        name, vector = read_query_vcf(Path(query_vcf), keys)
        query_names.append(name)
        query_vectors.append(vector)

    if query_vectors:
        matrix = np.vstack([panel_matrix, np.vstack(query_vectors)])
    else:
        matrix = panel_matrix
    samples = panel_samples + query_names
    sources = ["panel"] * len(panel_samples) + ["query"] * len(query_names)

    keep_markers = marker_filter(matrix, args.min_site_call_rate, args.min_maf)
    if args.max_sites is not None and int(keep_markers.sum()) > args.max_sites:
        rng = np.random.default_rng(args.seed)
        kept_indices = np.flatnonzero(keep_markers)
        selected = rng.choice(kept_indices, size=args.max_sites, replace=False)
        new_keep = np.zeros_like(keep_markers, dtype=bool)
        new_keep[selected] = True
        keep_markers = new_keep

    filtered = matrix[:, keep_markers]
    scaled = impute_and_standardize(filtered)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    metadata = read_metadata(Path(args.metadata) if args.metadata else None, args.sample_column)
    outputs: list[Path] = []

    if args.method in {"pca", "both"}:
        pca_coords, variance = run_pca(scaled, samples, sources)
        pca_coords = attach_metadata(pca_coords, metadata, args.sample_column, args.group_column, args.label_column)
        pca_path = out_dir / f"{args.prefix}_pca_coordinates.tsv"
        pca_coords.to_csv(pca_path, sep="\t", index=False)
        outputs.append(pca_path)
        fig = plot_coordinates(
            pca_coords,
            "PC1",
            "PC2",
            "Panel PCA with assembly-derived query samples",
            args.group_column,
            args.label_column,
            f"PC1 ({variance[0] * 100:.1f}% variance)",
            f"PC2 ({variance[1] * 100:.1f}% variance)",
        )
        outputs.extend(save_figure(fig, out_dir, args.prefix, "pca_pc1_pc2", args.formats, args.dpi))
        if "PC3" in pca_coords.columns:
            fig = plot_coordinates(
                pca_coords,
                "PC2",
                "PC3",
                "Panel PCA: PC2 vs PC3",
                args.group_column,
                args.label_column,
                f"PC2 ({variance[1] * 100:.1f}% variance)",
                f"PC3 ({variance[2] * 100:.1f}% variance)",
            )
            outputs.extend(save_figure(fig, out_dir, args.prefix, "pca_pc2_pc3", args.formats, args.dpi))

    if args.method in {"mds", "both"}:
        if scaled.shape[0] > args.max_mds_samples:
            raise ValueError(
                f"MDS requested for {scaled.shape[0]} samples, above --max-mds-samples {args.max_mds_samples}. "
                "Use --panel-samples or --max-panel-samples, or raise the limit if you have enough memory/time."
            )
        mds_coords = run_mds(scaled, samples, sources, args.seed)
        mds_coords = attach_metadata(mds_coords, metadata, args.sample_column, args.group_column, args.label_column)
        mds_path = out_dir / f"{args.prefix}_mds_coordinates.tsv"
        mds_coords.to_csv(mds_path, sep="\t", index=False)
        outputs.append(mds_path)
        fig = plot_coordinates(
            mds_coords,
            "MDS1",
            "MDS2",
            "MDS of panel and assembly-derived query samples",
            args.group_column,
            args.label_column,
            "MDS1",
            "MDS2",
        )
        outputs.extend(save_figure(fig, out_dir, args.prefix, "mds1_mds2", args.formats, args.dpi))

    summary_path = out_dir / f"{args.prefix}_matrix_summary.tsv"
    pd.DataFrame(
        [
            {
                "panel_samples": len(panel_samples),
                "query_samples": len(query_names),
                "markers_before_filtering": len(keys),
                "markers_after_filtering": int(keep_markers.sum()),
                "min_site_call_rate": args.min_site_call_rate,
                "min_maf": args.min_maf,
            }
        ]
    ).to_csv(summary_path, sep="\t", index=False)
    outputs.append(summary_path)

    for path in outputs:
        print(f"Wrote: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
