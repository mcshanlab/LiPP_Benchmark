"""
Read all *_lipsimscore.csv files in the current directory and plot the
distribution of "Closest_tanimoto" values for each file/category.
"""

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

plt.rcParams.update(
    {
        "font.size": 14,
        "axes.titlesize": 16,
        "axes.labelsize": 14,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
        "legend.fontsize": 12,
        "legend.title_fontsize": 13,
    }
)


def plot_box_with_points(
    categories: list[str],
    values_by_category: list[list[float]],
    fig_height: float,
    all_values: list[float],
) -> None:
    fig2, ax2 = plt.subplots(figsize=(10, fig_height))
    y_positions = np.arange(1, len(categories) + 1)
    ax2.boxplot(
        values_by_category,
        positions=y_positions,
        vert=False,
        patch_artist=True,
        showfliers=False,
        boxprops={"facecolor": "#4c72b0", "alpha": 0.5},
        medianprops={"color": "black"},
        whiskerprops={"color": "black"},
        capprops={"color": "black"},
    )

    rng = np.random.default_rng(0)
    for idx, vals in enumerate(values_by_category, start=1):
        if not vals:
            continue
        jitter = (rng.random(len(vals)) - 0.5) * 0.4
        ax2.scatter(
            vals,
            np.full(len(vals), idx) + jitter,
            s=10,
            color="black",
            alpha=0.5,
            linewidths=0,
        )

    ax2.set_yticks(y_positions)
    ax2.set_yticklabels(categories, fontsize=15)
    ax2.set_xlabel("Tanimoto similarity of the closest matching lipid", fontsize=18)
    ax2.set_ylabel("Protein Leiden Cluster", fontsize=18)
    ax2.tick_params(axis="x", labelsize=15)
    ax2.tick_params(axis="y", labelsize=15)
    # ax2.set_title(f"{rmsd_column}")
    ax2.grid(axis="x", linestyle="--", alpha=0.4)
    if all(0.0 <= v <= 1.0 for v in all_values):
        ax2.set_xlim(0.0, 1.2)

    out_file2 = Path("closest_tanimoto_bar_points.png")
    fig2.tight_layout()
    fig2.savefig(out_file2, dpi=300)
    print(f"Saved plot to {out_file2.resolve()}")


def plot_box_with_colored_points(
    categories: list[str],
    values_by_category: list[list[float]],
    protids_by_category: list[list[str]],
    fig_height: float,
    rmsd_column: str = "AF_lipid_RMSD_spy",
    all_scores_file: Path = Path("All_scores.csv"),
) -> None:
    rmsd_lookup_df = pd.read_csv(all_scores_file)
    if "protid" not in rmsd_lookup_df.columns or rmsd_column not in rmsd_lookup_df.columns:
        raise KeyError(f"All_scores.csv must contain columns: protid, {rmsd_column}")

    rmsd_lookup = (
        rmsd_lookup_df[["protid", rmsd_column]]
        .dropna(subset=["protid"])
        .drop_duplicates(subset=["protid"], keep="first")
        .set_index("protid")[rmsd_column]
    )
    rmsd_lookup = pd.to_numeric(rmsd_lookup, errors="coerce")

    def _rmsd_color(v: float) -> str:
        if pd.isna(v):
            return "#9e9e9e"  # missing
        v = float(v)
        if v < 2:
            return "#1f77b4"  # 0-2
        if v < 5:
            return "#2ca02c"  # 2-5
        if v < 10:
            return "#ff7f0e"  # 5-10
        elif v > 10:
            return "#d62728"  # >10

    fig3, ax3 = plt.subplots(figsize=(10, fig_height))
    y_positions = np.arange(1, len(categories) + 1)
    ax3.boxplot(
        values_by_category,
        positions=y_positions,
        vert=False,
        patch_artist=True,
        showfliers=False,
        boxprops={"facecolor": "#4c72b0", "alpha": 0.25},
        medianprops={"color": "black"},
        whiskerprops={"color": "black"},
        capprops={"color": "black"},
    )

    rng = np.random.default_rng(0)
    for idx, (vals, protids) in enumerate(zip(values_by_category, protids_by_category), start=1):
        if not vals:
            continue
        jitter = (rng.random(len(vals)) - 0.5) * 0.4
        rmsd_vals = [rmsd_lookup.get(p, np.nan) for p in protids]
        colors = [_rmsd_color(v) for v in rmsd_vals]
        ax3.scatter(
            vals,
            np.full(len(vals), idx) + jitter,
            s=16,
            c=colors,
            alpha=0.9,
            linewidths=0,
        )

    ax3.set_yticks(y_positions)
    ax3.set_yticklabels(categories)
    ax3.set_xlabel("Tanimoto similarity of the closest matching lipid")
    ax3.set_ylabel("Protein Leiden Cluster")
    #ax3.set_title(f"{rmsd_column}")
    ax3.grid(axis="x", linestyle="--", alpha=0.4)

    legend_elements = [
        Line2D([0], [0], marker="o", color="w", label=">10", markerfacecolor="#d62728", markersize=6),
        Line2D([0], [0], marker="o", color="w", label="5-10", markerfacecolor="#ff7f0e", markersize=6),
        Line2D([0], [0], marker="o", color="w", label="2-5", markerfacecolor="#2ca02c", markersize=6),
        Line2D([0], [0], marker="o", color="w", label="0-2", markerfacecolor="#1f77b4", markersize=6),
    ]
    ax3.legend(
    handles=legend_elements,
    title=rmsd_column,
    loc="upper left",           # legend corner to use as anchor
    bbox_to_anchor=(1.02, 1.0), # just outside right of axes
    borderaxespad=0.0,
    frameon=True)

    out_file3 = Path(f"closest_tanimoto_bar_points_colored_by_{rmsd_column}.png")
    fig3.tight_layout()
    fig3.savefig(out_file3, dpi=300)
    print(f"Saved plot to {out_file3.resolve()}")


def main() -> None:
    csv_files = sorted(Path(".").glob("*_lipsimscore.csv"))
    if not csv_files:
        raise FileNotFoundError("No *_lipsimscore.csv files found in current directory.")
    all_scores_file = Path("All_scores.csv")

    categories: list[str] = []
    values_by_category: list[list[float]] = []
    protids_by_category: list[list[str]] = []

    for csv_file in csv_files:
        category = csv_file.stem.replace("_lipsimscore", "")
        df = pd.read_csv(csv_file)

        if "Closest_tanimoto" not in df.columns:
            print(f"Skipping {csv_file.name}: no 'Closest_tanimoto' column.")
            continue

        vals = pd.to_numeric(df["Closest_tanimoto"], errors="coerce")
        valid_mask = vals.notna()
        vals = vals[valid_mask]
        if vals.empty:
            print(f"Skipping {csv_file.name}: no non-NaN Closest_tanimoto values.")
            continue
        if "protid" not in df.columns:
            print(f"Skipping {csv_file.name}: no 'protid' column.")
            continue

        categories.append(category)
        values_by_category.append(vals.tolist())
        protids_by_category.append(df.loc[valid_mask, "protid"].astype(str).tolist())

    if not categories:
        raise ValueError("No categories had non-NaN values in 'Closest_tanimoto'.")

    fig_height = max(4, 0.45 * len(categories) + 1.5)
    fig, ax = plt.subplots(figsize=(10, fig_height))

    vp = ax.violinplot(
        values_by_category,
        positions=range(1, len(categories) + 1),
        vert=False,
        showmeans=False,
        showmedians=True,
        showextrema=False,
    )

    for body in vp["bodies"]:
        body.set_alpha(0.6)

    ax.set_yticks(range(1, len(categories) + 1))
    ax.set_yticklabels(categories)
    ax.set_xlabel("Closest_tanimoto")
    ax.set_ylabel("Category")
    ax.set_title("Distribution of Closest_tanimoto by category")
    ax.grid(axis="x", linestyle="--", alpha=0.4)

    all_values = [v for sublist in values_by_category for v in sublist]
    if all(0.0 <= v <= 1.0 for v in all_values):
        ax.set_xlim(0.0, 1.0)

    out_file = Path("closest_tanimoto_distribution.png")
    fig.tight_layout()
    fig.savefig(out_file, dpi=300)
    print(f"Saved plot to {out_file.resolve()}")

    plot_box_with_points(categories, values_by_category, fig_height, all_values)
    all_scores_cols = pd.read_csv(all_scores_file, nrows=0).columns.tolist()
    rmsd_spy_columns = [c for c in all_scores_cols if c.endswith("_lipid_RMSD_spy")]
    if not rmsd_spy_columns:
        raise ValueError("No columns ending with '_lipid_RMSD_spy' found in All_scores.csv.")

    for rmsd_col in rmsd_spy_columns:
        plot_box_with_colored_points(
            categories,
            values_by_category,
            protids_by_category,
            fig_height,
            rmsd_column=rmsd_col,
            all_scores_file=all_scores_file,
        )


if __name__ == "__main__":
    main()
