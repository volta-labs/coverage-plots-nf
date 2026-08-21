#!/usr/bin/env python3
"""
Coverage heterogeneity plot.

Auto-discovers samples from qc_*/hcs_coverage_raw.bed in the target directory.
Heterogeneity = fraction of regions with norm. coverage < LOW or > HIGH.

Usage:
  python3 plot_coverage_heterogeneity.py \
      --dir ~/Downloads/ILMN148_bam_bai_sg_hcs \
      --run-name ILMN148 \
      --groups "cfDNA:Callisto-cfDNA,SGULP13:SOPHiA-gDNA,A3:Callisto-gDNA"

  python3 plot_coverage_heterogeneity.py \
      --dir ~/Downloads/ILMN153_bam_bai_sg_hcs \
      --run-name ILMN153 \
      --groups "SGULP11:SGHC25 Lane 1,SGULP13:SGHC25 Lane 4,cfDNA:SGHC25 Lane 4"
"""

import argparse, glob, os, re, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# ── thresholds ────────────────────────────────────────────────────────────────
LOW_THRESH  = 0.2
HIGH_THRESH = 5.0

# ── color palette (auto-assigned per group in order of first appearance) ──────
PALETTE = ["#4472C4", "#70AD47", "#C0504D", "#ED7D31", "#7030A0", "#00B0F0"]

CHROMS = [f"chr{i}" for i in range(1, 23)]

# ── helpers ───────────────────────────────────────────────────────────────────
def clean_sample_name(folder_name):
    """Strip qc_ prefix, _aligned_sorted and _S## suffixes."""
    name = re.sub(r"^qc_", "", folder_name)
    name = re.sub(r"_aligned_sorted$", "", name)
    name = re.sub(r"_S\d+$", "", name)
    return name

def assign_group(sample_name, patterns, default):
    for pattern, group in patterns:
        if re.search(pattern, sample_name, re.IGNORECASE):
            return group
    return default

def parse_groups(groups_str):
    """Parse 'pattern:Group Name,pattern2:Group 2' into list of (pattern, group)."""
    if not groups_str:
        return []
    pairs = []
    for item in groups_str.split(","):
        if ":" in item:
            pattern, group = item.split(":", 1)
            pairs.append((pattern.strip(), group.strip()))
    return pairs

def load_sample_sheet(tsv_path):
    """Load a TSV with columns 'sample' and 'group'. Returns {sample: group} dict."""
    df = pd.read_csv(tsv_path, sep="\t")
    df.columns = [c.strip().lower() for c in df.columns]
    if "sample" not in df.columns or "group" not in df.columns:
        sys.exit(f"ERROR: sample sheet must have 'sample' and 'group' columns, got: {list(df.columns)}")
    return dict(zip(df["sample"].str.strip(), df["group"].str.strip()))

def discover_samples(directory, group_patterns, default_group, sample_map=None):
    beds = sorted(glob.glob(os.path.join(directory, "qc_*", "hcs_coverage_raw.bed")))
    if not beds:
        sys.exit(f"ERROR: no qc_*/hcs_coverage_raw.bed found in {directory}")
    samples = []
    for bed in beds:
        folder = os.path.basename(os.path.dirname(bed))
        label  = clean_sample_name(folder)
        if sample_map is not None:
            group = sample_map.get(label, default_group)
        else:
            group = assign_group(label, group_patterns, default_group)
        samples.append((label, bed, group))
    return samples

def compute_heterogeneity(samples):
    records = []
    for label, bed_path, group in samples:
        df = pd.read_csv(bed_path, sep="\t", header=None,
                         names=["chrom","start","end","coverage"])
        df = df[df["chrom"].isin(CHROMS)].copy()
        if df.empty or df["coverage"].sum() == 0:
            print(f"  WARNING: {label} — no coverage, skipping")
            continue
        mean_cov   = df["coverage"].mean()
        df["norm"] = df["coverage"] / mean_cov
        total      = len(df)
        outside    = ((df["norm"] < LOW_THRESH) | (df["norm"] > HIGH_THRESH)).sum()
        het_pct    = outside / total * 100
        records.append({"label": label, "group": group,
                         "het": het_pct, "total": total, "outside": int(outside)})
        print(f"  {label:35s}  {het_pct:.3f}%  ({outside}/{total})")
    return records

def make_plot(records, title, group_order, colors, out_file):
    results = pd.DataFrame(records)
    results["group_idx"] = results["group"].apply(
        lambda g: group_order.index(g) if g in group_order else len(group_order)
    )
    results = results.sort_values(["group_idx","label"], ascending=[True, True])
    results = results.iloc[::-1].reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(8, max(4, len(results) * 0.45 + 1.5)))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    bar_colors = [colors.get(g, "#999999") for g in results["group"]]
    y_pos = range(len(results))
    bars  = ax.barh(list(y_pos), results["het"].values,
                    color=bar_colors, edgecolor="white", linewidth=0.4, height=0.65)

    for bar, val in zip(bars, results["het"].values):
        ax.text(bar.get_width() + 0.01, bar.get_y() + bar.get_height() / 2,
                f"{val:.2f}%", va="center", ha="left", fontsize=7.5, color="#444444")

    ax.set_yticks(list(y_pos))
    ax.set_yticklabels(results["label"].values, fontsize=8.5)
    ax.set_xlabel("Coverage heterogeneity (%)", fontsize=10)
    ax.set_title(title, fontsize=11, pad=10)
    ax.set_xlim(0, max(results["het"].max() * 1.25, 1.0))
    ax.spines[["top","right"]].set_visible(False)
    ax.axvline(0, color="#cccccc", lw=0.6)
    ax.xaxis.grid(True, linestyle="--", alpha=0.4, color="#aaaaaa")
    ax.set_axisbelow(True)

    legend_elements = [
        Patch(facecolor=colors[g], label=g)
        for g in group_order if g in results["group"].values
    ]
    if legend_elements:
        ax.legend(handles=legend_elements, title="Group", fontsize=8,
                  title_fontsize=8.5, loc="lower right",
                  framealpha=0.85, edgecolor="#cccccc")

    ax.text(0.01, -0.07,
            f"Heterogeneity = fraction of regions with norm. coverage "
            f"< {LOW_THRESH} or > {HIGH_THRESH}",
            transform=ax.transAxes, fontsize=7, color="#777777")

    plt.tight_layout(pad=1.0)
    plt.savefig(out_file, dpi=180, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out_file}")

# ── main ──────────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser()
parser.add_argument("--dir",          required=True, help="Directory containing qc_* subfolders")
parser.add_argument("--run-name",     required=True, help="Run name for plot title and output filename")
parser.add_argument("--groups",       default="",
                    help='Group patterns, e.g. "cfDNA:Callisto-cfDNA,SGULP13:SOPHiA-gDNA"')
parser.add_argument("--sample-sheet", default=None,
                    help="TSV file with columns 'sample' and 'group' — overrides --groups")
args = parser.parse_args()

directory     = os.path.expanduser(args.dir)
run_name      = args.run_name
default_group = run_name
group_patterns = parse_groups(args.groups)

sample_map = None
if args.sample_sheet:
    sample_map = load_sample_sheet(os.path.expanduser(args.sample_sheet))
    print(f"  Loaded sample sheet: {len(sample_map)} entries")

print(f"\n=== {run_name} ===")
samples = discover_samples(directory, group_patterns, default_group, sample_map=sample_map)
print(f"Found {len(samples)} samples\n")

records = compute_heterogeneity(samples)
if not records:
    sys.exit("No data loaded.")

# assign colors in sorted group order (natural sort: 32-1, 32-2, ..., 33-1, ...)
import re as _re
def _natural_key(s):
    return [int(c) if c.isdigit() else c.lower() for c in _re.split(r'(\d+)', s)]

group_order = sorted(set(g for _, _, g in samples), key=_natural_key)
colors = {g: PALETTE[i % len(PALETTE)] for i, g in enumerate(group_order)}

make_plot(records,
          title    = f"Coverage heterogeneity per sample — {run_name}",
          group_order = group_order,
          colors   = colors,
          out_file = f"coverage_heterogeneity_{run_name}.png")
