#!/usr/bin/env python3
"""
Create publication-friendly figures from gtcheck top-hit summary tables.

Inputs are the TSV files written by scripts/summarize_gtcheck_top_hits.py.
Outputs are PNG/SVG/PDF figures suitable for README files, reports, and slide
decks.
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable, Optional

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
import pandas as pd
import seaborn as sns


PLOT_STYLE = {
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.color": "#e6e6e6",
    "grid.linewidth": 0.8,
    "axes.edgecolor": "#444444",
    "axes.labelcolor": "#222222",
    "xtick.color": "#333333",
    "ytick.color": "#333333",
    "font.size": 10,
    "axes.titlesize": 13,
    "axes.titleweight": "bold",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot gtcheck top-hit and sample-QC summary tables.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--top-hits",
        required=True,
        help="TSV produced by scripts/summarize_gtcheck_top_hits.py.",
    )
    parser.add_argument(
        "--sample-summary",
        default=None,
        help=(
            "Optional one-row-per-query TSV from summarize_gtcheck_top_hits.py. "
            "If omitted, rank-gap metrics are derived from --top-hits when possible."
        ),
    )
    parser.add_argument("--out-dir", default="figures", help="Directory for output figures.")
    parser.add_argument("--prefix", default="gtcheck", help="Filename prefix for output figures.")
    parser.add_argument("--max-rank", type=int, default=10, help="Maximum rank to include.")
    parser.add_argument(
        "--max-heatmap-columns",
        type=int,
        default=30,
        help="Maximum number of panel samples to show as heatmap columns.",
    )
    parser.add_argument(
        "--formats",
        nargs="+",
        default=["png", "svg"],
        choices=["png", "svg", "pdf"],
        help="Figure formats to write.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="DPI for raster outputs such as PNG.",
    )
    return parser.parse_args()


def read_top_hits(path: Path, max_rank: int) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    required = {
        "query_sample",
        "panel_sample",
        "rank",
        "match_fraction",
        "sites_compared",
        "discordance_per_site",
    }
    missing = sorted(required - set(df.columns))
    if missing:
        raise ValueError(f"{path} is missing required columns: {', '.join(missing)}")
    df = df.copy()
    df["rank"] = pd.to_numeric(df["rank"], errors="coerce")
    for col in ["match_fraction", "sites_compared", "discordance_per_site", "confidence_score"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df[df["rank"].between(1, max_rank)].copy()
    if "query_sample_basename" not in df.columns:
        df["query_sample_basename"] = df["query_sample"].astype(str).map(lambda value: Path(value).name)
    return df.sort_values(["query_sample", "rank"])


def derive_sample_summary(top_hits: pd.DataFrame) -> pd.DataFrame:
    records = []
    for query, sub in top_hits.sort_values("rank").groupby("query_sample", sort=False):
        sub = sub.sort_values("rank")
        best = sub.iloc[0]
        second = sub.iloc[1] if len(sub) > 1 else None
        second_match = second["match_fraction"] if second is not None else pd.NA
        gap = best["match_fraction"] - second_match if second is not None else pd.NA
        records.append(
            {
                "query_sample": query,
                "query_sample_basename": best.get("query_sample_basename", Path(str(query)).name),
                "top_panel_sample": best["panel_sample"],
                "top_match_fraction": best["match_fraction"],
                "top_sites_compared": best["sites_compared"],
                "second_panel_sample": second["panel_sample"] if second is not None else pd.NA,
                "second_match_fraction": second_match,
                "match_fraction_gap_rank1_rank2": gap,
                "max_sites_compared_any_hit": sub["sites_compared"].max(),
            }
        )
    return pd.DataFrame(records)


def read_sample_summary(path: Optional[Path], top_hits: pd.DataFrame) -> pd.DataFrame:
    if path is None:
        summary = derive_sample_summary(top_hits)
    else:
        summary = pd.read_csv(path, sep="\t")
    for col in [
        "top_match_fraction",
        "top_sites_compared",
        "second_match_fraction",
        "match_fraction_gap_rank1_rank2",
        "max_sites_compared_any_hit",
        "vcf_call_rate",
        "vcf_missing_rate",
        "vcf_het_rate_among_called",
        "top_sites_compared_fraction_of_panel",
        "max_sites_compared_fraction_of_panel",
    ]:
        if col in summary.columns:
            summary[col] = pd.to_numeric(summary[col], errors="coerce")
    if "query_sample_basename" not in summary.columns:
        summary["query_sample_basename"] = summary["query_sample"].astype(str).map(lambda value: Path(value).name)
    return summary.sort_values(
        ["top_match_fraction", "top_sites_compared", "query_sample_basename"],
        ascending=[False, False, True],
    )


def save_figure(fig: plt.Figure, out_dir: Path, prefix: str, stem: str, formats: Iterable[str], dpi: int) -> list[Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = []
    for fmt in formats:
        path = out_dir / f"{prefix}_{stem}.{fmt}"
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
        paths.append(path)
    plt.close(fig)
    return paths


def ordered_queries(summary: pd.DataFrame) -> list[str]:
    return summary["query_sample"].astype(str).tolist()


def plot_top_hits_lollipop(top_hits: pd.DataFrame, summary: pd.DataFrame) -> plt.Figure:
    order = ordered_queries(summary)
    df = top_hits[top_hits["query_sample"].isin(order)].copy()
    df["query_sample"] = pd.Categorical(df["query_sample"], categories=order[::-1], ordered=True)
    height = max(4.5, 0.38 * len(order) + 1.8)
    fig, ax = plt.subplots(figsize=(12, height))

    for query, sub in df.groupby("query_sample", observed=True):
        y = list(df["query_sample"].cat.categories).index(query)
        ax.hlines(y=y, xmin=sub["match_fraction"].min(), xmax=sub["match_fraction"].max(), color="#bbbbbb", lw=1.0, zorder=1)

    rank_min = float(df["rank"].min())
    rank_max = float(df["rank"].max())
    if rank_max == rank_min:
        df["marker_size"] = 120.0
    else:
        df["marker_size"] = 70.0 + (df["rank"] - rank_min) * (220.0 - 70.0) / (rank_max - rank_min)

    color_norm = Normalize(vmin=float(df["sites_compared"].min()), vmax=float(df["sites_compared"].max()))
    cmap = plt.get_cmap("viridis")

    ax.scatter(
        df["match_fraction"],
        df["query_sample"],
        c=df["sites_compared"],
        s=df["marker_size"],
        cmap=cmap,
        norm=color_norm,
        edgecolors="white",
        linewidths=0.6,
        zorder=3,
    )
    top = df[df["rank"] == 1]
    for _, row in top.iterrows():
        ax.annotate(
            str(row["panel_sample"]),
            xy=(row["match_fraction"], row["query_sample"]),
            xytext=(5, 5),
            textcoords="offset points",
            fontsize=8,
            color="#333333",
        )
    ax.set_title("Top gtcheck hits by assembly")
    ax.set_xlabel("Match fraction")
    ax.set_ylabel("Assembly")
    ax.set_xlim(max(0, df["match_fraction"].min() - 0.02), 1.005)

    colorbar = fig.colorbar(ScalarMappable(norm=color_norm, cmap=cmap), ax=ax, pad=0.02)
    colorbar.set_label("Sites compared")

    rank_values = sorted(df["rank"].dropna().astype(int).unique())
    if len(rank_values) > 4:
        rank_values = rank_values[:3] + [rank_values[-1]]
    legend_handles = []
    for rank in rank_values:
        size = float(df.loc[df["rank"] == rank, "marker_size"].iloc[0])
        legend_handles.append(
            Line2D(
                [0],
                [0],
                marker="o",
                color="none",
                label=f"Rank {rank}",
                markerfacecolor="#666666",
                markeredgecolor="white",
                markeredgewidth=0.6,
                markersize=(size ** 0.5) / 1.3,
            )
        )
    ax.legend(
        handles=legend_handles,
        loc="center left",
        bbox_to_anchor=(1.16, 0.5),
        frameon=False,
        title="Marker size",
    )
    return fig


def plot_rank_gap(summary: pd.DataFrame) -> plt.Figure:
    df = summary.copy()
    order = ordered_queries(df)
    df["query_sample"] = pd.Categorical(df["query_sample"], categories=order[::-1], ordered=True)
    height = max(4.5, 0.36 * len(order) + 1.8)
    fig, ax = plt.subplots(figsize=(11, height))
    sns.scatterplot(
        data=df,
        x="match_fraction_gap_rank1_rank2",
        y="query_sample",
        hue="top_match_fraction",
        size="top_sites_compared",
        sizes=(45, 190),
        palette="mako",
        edgecolor="white",
        linewidth=0.7,
        ax=ax,
    )
    for _, row in df.iterrows():
        if pd.notna(row.get("match_fraction_gap_rank1_rank2")):
            ax.annotate(
                str(row.get("top_panel_sample", "")),
                xy=(row["match_fraction_gap_rank1_rank2"], row["query_sample"]),
                xytext=(5, 4),
                textcoords="offset points",
                fontsize=8,
            )
    ax.axvline(0.01, color="#d95f02", ls="--", lw=1, alpha=0.7)
    ax.set_title("Separation between best and second-best hits")
    ax.set_xlabel("Rank 1 match fraction minus rank 2 match fraction")
    ax.set_ylabel("Assembly")
    ax.set_xlim(left=0)
    ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False, title="Top match / sites")
    return fig


def plot_match_vs_sites(summary: pd.DataFrame) -> plt.Figure:
    df = summary.copy()
    fig, ax = plt.subplots(figsize=(10, 7))
    sns.scatterplot(
        data=df,
        x="top_sites_compared",
        y="top_match_fraction",
        hue="match_fraction_gap_rank1_rank2",
        size="top_sites_compared",
        sizes=(45, 220),
        palette="rocket_r",
        edgecolor="white",
        linewidth=0.7,
        ax=ax,
    )
    for _, row in df.iterrows():
        ax.annotate(
            str(row["query_sample_basename"]),
            xy=(row["top_sites_compared"], row["top_match_fraction"]),
            xytext=(5, 5),
            textcoords="offset points",
            fontsize=8,
        )
    ax.set_title("Best-hit match fraction vs sites compared")
    ax.set_xlabel("Top-hit sites compared")
    ax.set_ylabel("Top-hit match fraction")
    ax.set_ylim(max(0, df["top_match_fraction"].min() - 0.03), 1.005)
    ax.set_xlim(left=0)
    ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False, title="Rank gap / sites")
    return fig


def plot_top_hit_heatmap(top_hits: pd.DataFrame, summary: pd.DataFrame, max_columns: int) -> plt.Figure:
    order = ordered_queries(summary)
    panel_order: list[str] = []
    for query in order:
        sub = top_hits[top_hits["query_sample"] == query].sort_values("rank")
        for panel in sub["panel_sample"].astype(str):
            if panel not in panel_order:
                panel_order.append(panel)
            if len(panel_order) >= max_columns:
                break
        if len(panel_order) >= max_columns:
            break

    matrix = (
        top_hits[top_hits["panel_sample"].astype(str).isin(panel_order)]
        .pivot_table(index="query_sample", columns="panel_sample", values="match_fraction", aggfunc="max")
        .reindex(index=order, columns=panel_order)
    )
    height = max(4, 0.38 * len(order) + 2.0)
    width = max(8, 0.34 * len(panel_order) + 4.5)
    fig, ax = plt.subplots(figsize=(width, height))
    sns.heatmap(
        matrix,
        cmap="YlGn",
        vmin=max(0, float(top_hits["match_fraction"].min()) - 0.01),
        vmax=1,
        linewidths=0.35,
        linecolor="white",
        cbar_kws={"label": "Match fraction"},
        ax=ax,
    )
    ax.set_title("Top-hit match fraction heatmap")
    ax.set_xlabel("Panel sample")
    ax.set_ylabel("Assembly")
    ax.tick_params(axis="x", labelrotation=55)
    return fig


def plot_query_qc_bars(summary: pd.DataFrame) -> Optional[plt.Figure]:
    metric_labels = {
        "vcf_call_rate": "Call rate",
        "vcf_missing_rate": "Missing rate",
        "vcf_het_rate_among_called": "Het rate",
        "top_sites_compared_fraction_of_panel": "Top sites / panel",
        "max_sites_compared_fraction_of_panel": "Max sites / panel",
    }
    available = [col for col in metric_labels if col in summary.columns and summary[col].notna().any()]
    if not available:
        return None
    df = summary[["query_sample", "query_sample_basename", *available]].copy()
    long = df.melt(
        id_vars=["query_sample", "query_sample_basename"],
        value_vars=available,
        var_name="metric",
        value_name="fraction",
    )
    long["metric"] = long["metric"].map(metric_labels)
    order = ordered_queries(summary)
    long["query_sample"] = pd.Categorical(long["query_sample"], categories=order[::-1], ordered=True)

    height = max(4.5, 0.4 * len(order) + 1.8)
    fig, ax = plt.subplots(figsize=(11, height))
    sns.barplot(
        data=long,
        x="fraction",
        y="query_sample",
        hue="metric",
        orient="h",
        palette="Set2",
        ax=ax,
    )
    ax.set_title("Query VCF QC metrics")
    ax.set_xlabel("Fraction")
    ax.set_ylabel("Assembly")
    ax.set_xlim(0, 1)
    ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False, title="Metric")
    return fig


def main() -> int:
    args = parse_args()
    sns.set_theme(style="whitegrid", rc=PLOT_STYLE)

    top_hits = read_top_hits(Path(args.top_hits), args.max_rank)
    if top_hits.empty:
        raise SystemExit("ERROR: no top-hit rows remain after rank filtering")
    summary = read_sample_summary(Path(args.sample_summary) if args.sample_summary else None, top_hits)

    out_dir = Path(args.out_dir)
    outputs: list[Path] = []
    figures = [
        ("top_hits_lollipop", plot_top_hits_lollipop(top_hits, summary)),
        ("rank1_rank2_gap", plot_rank_gap(summary)),
        ("match_fraction_vs_sites", plot_match_vs_sites(summary)),
        ("top_hit_heatmap", plot_top_hit_heatmap(top_hits, summary, args.max_heatmap_columns)),
    ]
    qc_fig = plot_query_qc_bars(summary)
    if qc_fig is not None:
        figures.append(("query_qc_bars", qc_fig))

    for stem, fig in figures:
        outputs.extend(save_figure(fig, out_dir, args.prefix, stem, args.formats, args.dpi))

    for path in outputs:
        print(f"Wrote: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
