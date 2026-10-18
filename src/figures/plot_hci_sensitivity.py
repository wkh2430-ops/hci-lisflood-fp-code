"""Plot fixed-RSE HCI sensitivity summaries."""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=Path("data/processed/sensitivity/HCI_sensitivity_summary.csv"))
    parser.add_argument("--output", type=Path, default=Path("figures/reproduced"))
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    data = pd.read_csv(args.input)
    colors = {"HYK": "#1f77b4", "JHT": "#ff7f0e", "GC": "#2ca02c"}
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.8), layout="constrained")
    for ax, parameter, title in zip(axes, ("Qbank", "beta"), ("Bankfull-discharge sensitivity", "Activation-steepness sensitivity")):
        subset = data[data.parameter == parameter]
        scenarios = list(dict.fromkeys(subset.scenario))
        x = np.arange(len(scenarios))
        for code in ("HYK", "JHT", "GC"):
            values = subset[subset.station_code == code].set_index("scenario").loc[scenarios, "HCI_mean_change_pct_vs_exact_baseline"]
            ax.plot(x, values, "o-", label=code, color=colors[code], lw=1.5, ms=4)
        ax.axhline(0, color="0.45", lw=0.8)
        ax.set_xticks(x, scenarios, rotation=30, ha="right")
        ax.set_ylabel("Change in mean HCI (%)")
        ax.set_title(title)
        ax.grid(axis="y", ls="--", alpha=0.25)
    axes[0].legend(frameon=False)
    for extension in ("png", "pdf", "svg"):
        fig.savefig(args.output / f"Fig_S_HCI_Sensitivity.{extension}", dpi=300, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
