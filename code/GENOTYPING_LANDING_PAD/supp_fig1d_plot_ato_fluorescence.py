# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 10:39:46 2026

@author: ig8
"""

#!/usr/bin/env python3
"""Plot the ATO fluorescence CSV as a labelled bar chart.

Requires matplotlib. No pandas or seaborn dependency.

Usage:
    python plot_ato_fluorescence.py "ATO_fluorescence (1).csv"
    python plot_ato_fluorescence.py input.csv --output-dir figures --color "#C44E70"

Writes ato_fluorescence.png (600 dpi), .pdf, and .svg.
Bars show individual observations; no averaging or error bars are applied.
"""

import argparse
import csv
import math
from pathlib import Path
import textwrap


def read_data(csv_path):
    """Read and validate one fluorescence/condition/population series."""
    required = {
        "time_ATO_wk", "percent_fluorescent_positive", "fluorescence",
        "constitutive_inducible", "FACS_population",
    }
    with Path(csv_path).open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError("Missing CSV columns: " + ", ".join(sorted(missing)))
        rows = list(reader)
    if not rows:
        raise ValueError("The CSV contains no measurements.")

    metadata = {}
    for key in ("fluorescence", "constitutive_inducible", "FACS_population"):
        values = {(row.get(key) or "").strip() for row in rows}
        if "" in values or len(values) != 1:
            raise ValueError("Expected one non-empty value for %s; filter the CSV to one series." % key)
        metadata[key] = values.pop()

    points = []
    for row_number, row in enumerate(rows, start=2):
        try:
            week = float(row["time_ATO_wk"])
            percent = float(row["percent_fluorescent_positive"])
        except (TypeError, ValueError) as exc:
            raise ValueError("Invalid time or percentage on CSV row %d." % row_number) from exc
        if not math.isfinite(week) or not math.isfinite(percent) or not 0 <= percent <= 100:
            raise ValueError("CSV row %d needs a finite time and a percentage from 0 to 100." % row_number)
        points.append((week, percent))
    points.sort()
    if len({week for week, _ in points}) != len(points):
        raise ValueError("Repeated time points found. Plot replicates explicitly before summarising them.")
    return points, metadata


def plot_data(points, metadata, output_dir, color="#C44E70"):
    """Create a figure with a true numeric time axis and export three formats."""
    try:
        import matplotlib
        matplotlib.use("Agg")  # Save figures without opening a desktop window.
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise RuntimeError("Matplotlib is required. Install it with: python -m pip install matplotlib") from exc

    weeks, percentages = zip(*points)
    # Keep the bars separated even if future measurements are closer together.
    gaps = [right - left for left, right in zip(weeks, weeks[1:])]
    width = 0.65 * min(gaps) if gaps else 0.65
    reporter = metadata["fluorescence"]
    text_color = "#253245"

    with plt.rc_context({
        "font.family": "DejaVu Sans", "font.size": 11,
        "text.color": text_color, "axes.labelcolor": text_color,
        "xtick.color": text_color, "ytick.color": text_color,
        "svg.fonttype": "none", "pdf.fonttype": 42,
    }):
        fig, ax = plt.subplots(figsize=(18, 5), facecolor="white")
        fig.subplots_adjust(left=0.13, right=0.97, bottom=0.22, top=0.79)

        fig.text(0.13, 0.94, "%s positivity over time" % reporter,
                 fontsize=19, fontweight="bold", ha="left", va="top")
        fig.text(0.13, 0.885, "%s expression" % metadata["constitutive_inducible"].capitalize(),
                 fontsize=14, color="#64748B", ha="left", va="top")

        ax.set_axisbelow(True)
        ax.grid(axis="y", color="#E5E9EE", linewidth=0.0)
        ax.bar(weeks, percentages, width=width, color=color, edgecolor="none", zorder=3)
        for week, percent in points:
            ax.annotate("%g%%" % percent, (week, percent), xytext=(0, 6),
                        textcoords="offset points", ha="center", va="bottom",
                        fontsize=14, fontweight="semibold")

        ax.set_xlabel("ATO time (weeks)", labelpad=15)
        ax.set_ylabel("%s-positive cells (%%)" % reporter, labelpad=15)
        ax.set_xticks(weeks)
        ax.set_xticklabels(["%g" % week for week in weeks])
        ax.set_xlim(min(weeks) - width, max(weeks) + width)
        ax.set_ylim(0, 110)  # Headroom for the label above a 100% bar.
        ax.set_yticks(range(0, 101, 20))
        ax.tick_params(axis="both", length=0, pad=8)
        for side in ("top", "right", "left"):
            ax.spines[side].set_visible(False)
        ax.spines["bottom"].set_color("#AAB4C0")

        gate = textwrap.fill("FACS population: " + metadata["FACS_population"], width=90)
        fig.text(0.13, 0.10, gate, fontsize=13, color="#64748B", va="top")
        fig.text(0.13, 0.045, "One recorded measurement per time point.",
                 fontsize=12, color="#64748B", va="top")

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        paths = []
        for extension in ("png", "pdf", "svg"):
            path = output_dir / ("ato_fluorescence." + extension)
            fig.savefig(path, dpi=600, bbox_inches="tight", facecolor="white")
            paths.append(path.resolve())
        plt.close(fig)
    return paths


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", type=Path, help="Path to the ATO fluorescence CSV")
    parser.add_argument("--output-dir", type=Path, default=Path("."))
    parser.add_argument("--color", default="#C44E70", help="Bar colour as a name or hex code")
    args = parser.parse_args()
    try:
        points, metadata = read_data(args.csv)
        paths = plot_data(points, metadata, args.output_dir, args.color)
    except (OSError, ValueError, RuntimeError) as exc:
        parser.exit(1, "Error: %s\n" % exc)
    for path in paths:
        print("Saved:", path)


if __name__ == "__main__":
    main()
