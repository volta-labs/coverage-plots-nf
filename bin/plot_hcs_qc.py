#!/usr/bin/env python3
"""
SOPHiA-style QC plots — chromosome coverage + GC bias.

Usage:
  python3 plot_hcs_qc.py --targets FILE --gc FILE --outdir DIR --sample NAME
"""

import argparse, os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
from scipy.stats import gaussian_kde

# ── parse arguments ───────────────────────────────────────────────────────────
parser = argparse.ArgumentParser()
parser.add_argument("--targets", default="qc_hcs/hcs_coverage_raw.bed")
parser.add_argument("--gc",      default="qc_hcs/target_gc.csv")
parser.add_argument("--outdir",  default="qc_hcs")
parser.add_argument("--sample",  default="HCS_v2.0_inferred")
args = parser.parse_args()

TARGET_COV_BED = args.targets
GC_CACHE       = args.gc
OUT_DIR        = args.outdir
SAMPLE_NAME    = args.sample

CHROMS = [f"chr{i}" for i in range(1, 23)]

# SOPHiA color palette
COL_BLUE   = "#4472C4"   # High GC
COL_ORANGE = "#E87722"   # High AT / dropout
COL_GREY   = "#9E9E9E"   # Neutral

if not os.path.exists(TARGET_COV_BED):
    sys.exit(f"ERROR: missing {TARGET_COV_BED}  —  run bedtools coverage first")

# ── load coverage data ────────────────────────────────────────────────────────
print("Loading coverage files...")
targets = pd.read_csv(TARGET_COV_BED, sep="\t", header=None,
                      names=["chrom","start","end","coverage"])

# Normalise chrom names: add 'chr' prefix if absent (handles BEDs with '1','2'... or 'chr1','chr2'...)
targets["chrom"] = targets["chrom"].astype(str).apply(
    lambda c: c if c.startswith("chr") else f"chr{c}"
)

# filter to standard chroms only (excludes alt/unplaced contigs)
targets = targets[targets["chrom"].isin(CHROMS)].copy()

target_mean = targets["coverage"].mean()
targets["norm"] = (targets["coverage"] / target_mean).clip(0.02, 30)
total_bp = int((targets["end"] - targets["start"]).sum())

print(f"  Targets   : {len(targets):,}  ({total_bp/1000:.1f} kbp)")
print(f"  Mean depth: {target_mean:.0f}x")

# ── load GC content ───────────────────────────────────────────────────────────
if not os.path.exists(GC_CACHE):
    sys.exit("ERROR: target_gc.csv not found. Run fetch_all_gc.py first.")

gc_df = pd.read_csv(GC_CACHE).dropna(subset=["gc"])
print(f"  GC cache  : {len(gc_df)} regions")

# ── compute per-sample norm_cov by joining GC cache with this sample's coverage
gc_df = gc_df[["chrom","start","end","gc"]].merge(
    targets[["chrom","start","end","norm"]],
    on=["chrom","start","end"], how="inner"
).rename(columns={"norm": "norm_cov"})
gc_df["norm_cov"] = gc_df["norm_cov"].clip(0.02, 30)
gc_df = gc_df.dropna(subset=["gc","norm_cov"])

# ── colour each target by GC ──────────────────────────────────────────────────
targets_plot = targets.merge(
    gc_df[["chrom","start","end","gc"]],
    on=["chrom","start","end"], how="left"
)
targets_plot["chrom"] = pd.Categorical(
    targets_plot["chrom"], categories=CHROMS, ordered=True
)
targets_plot = targets_plot.sort_values(["chrom","start"]).reset_index(drop=True)

def target_color(gc):
    if pd.isna(gc):  return COL_GREY
    if gc > 55:      return COL_BLUE
    if gc < 38:      return COL_ORANGE
    return COL_GREY

targets_plot["color"] = targets_plot["gc"].apply(target_color)

# ── PLOT 1: chromosome coverage (per-target, SOPHiA style) ───────────────────
print("\nGenerating chromosome coverage plot...")

fig, ax = plt.subplots(figsize=(20, 3.8))
fig.patch.set_facecolor("white")
ax.set_facecolor("white")

x_pos = 0
xticks, xlabels = [], []
bg_colors = ["#F5F5F5", "#EBEBEB"]

for ci, chrom in enumerate(CHROMS):
    cdata = targets_plot[targets_plot["chrom"] == chrom]
    if cdata.empty:
        continue
    n  = len(cdata)
    xs = np.arange(x_pos, x_pos + n)
    ax.axvspan(x_pos, x_pos + n, alpha=1.0, color=bg_colors[ci % 2], zorder=0)
    ax.scatter(xs, cdata["norm"].values, c=cdata["color"].values,
               s=14, alpha=0.85, zorder=2, linewidths=0, rasterized=True)
    xticks.append(x_pos + n // 2)
    xlabels.append(chrom.replace("chr",""))
    x_pos += n + 3

# reference lines
ax.axhline(1.0, color="#555555", lw=0.8, ls="--", alpha=0.7, zorder=3)
ax.axhline(0.2, color="#999999", lw=0.5, ls="--", alpha=0.5, zorder=3)

ax.set_yscale("log")
ax.set_ylim(0.03, 10)
ax.set_yticks([0.04, 0.2, 1, 5])
ax.set_yticklabels(["0.04","0.2","1","5"], fontsize=8)
ax.set_xticks(xticks)
ax.set_xticklabels(xlabels, fontsize=7.5)
ax.set_xlim(0, x_pos)
ax.set_ylabel("Normalized Coverage", fontsize=9)
ax.set_xlabel("Chromosome", fontsize=9)
ax.spines[["top","right"]].set_visible(False)

ax.set_title(
    f"{SAMPLE_NAME} — {total_bp/1000:.2f}kbp; "
    f"mean_depth — {target_mean:.0f}x; "
    f"regions — {len(targets)}",
    fontsize=8.5, loc="left", pad=6
)

legend_elements = [
    Line2D([0],[0], marker="o", color="w", markerfacecolor=COL_BLUE,
           markersize=7, label="High GC"),
    Line2D([0],[0], marker="o", color="w", markerfacecolor=COL_GREY,
           markersize=7, label="Neutral"),
    Line2D([0],[0], marker="o", color="w", markerfacecolor=COL_ORANGE,
           markersize=7, label="High AT"),
]
ax.legend(handles=legend_elements, fontsize=7.5, loc="upper right",
          framealpha=0.85, edgecolor="#cccccc")

plt.tight_layout(pad=0.5)
out1 = f"{OUT_DIR}/{SAMPLE_NAME}_chromosome_coverage.png"
plt.savefig(out1, dpi=180, bbox_inches="tight")
plt.close()
print(f"  Saved: {out1}")

# ── PLOT 2: GC bias ───────────────────────────────────────────────────────────
print("Generating GC bias plot...")

fig, ax = plt.subplots(figsize=(4.2, 4.2))
fig.patch.set_facecolor("white")

x = gc_df["gc"].values
y = gc_df["norm_cov"].values
mask = (x > 15) & (x < 85) & (y > 0.03) & (y < 20)
x, y = x[mask], y[mask]
ly = np.log10(y)

# KDE density for scatter colouring
if len(x) > 10:
    try:
        kde = gaussian_kde(np.vstack([x, ly]))
        density = kde(np.vstack([x, ly]))
        density = (density - density.min()) / (density.max() - density.min() + 1e-9)
    except Exception:
        density = np.ones(len(x)) * 0.5
else:
    density = np.ones(len(x)) * 0.5

# compute severity BEFORE reordering (x and y still aligned)
low_gc_med  = np.median(y[x < 40]) if np.any(x < 40) else 1
high_gc_med = np.median(y[x > 60]) if np.any(x > 60) else 1
gc_ratio    = high_gc_med / max(low_gc_med, 0.01)
severity    = "Low" if gc_ratio < 1.5 else ("Moderate" if gc_ratio < 3.5 else "High")

# sort so dense points render on top
order = np.argsort(density)
x, ly_plot, density = x[order], ly[order], density[order]

cmap = plt.get_cmap("YlOrRd")
sc = ax.scatter(x, 10**ly_plot, c=density, cmap=cmap,
                s=2, alpha=0.75, linewidths=0, zorder=2,
                norm=mcolors.PowerNorm(gamma=0.5, vmin=0, vmax=1))

# reference lines matching SOPHiA grid
for yref in [0.2, 1.0, 5.0]:
    ax.axhline(yref, color="#888888", lw=0.7, ls="--", alpha=0.55, zorder=3)

ax.set_yscale("log")
ax.set_ylim(0.03, 25)
ax.set_yticks([0.04, 0.2, 1, 5, 25])
ax.set_yticklabels(["0.04","0.2","1","5","25"], fontsize=9)
ax.set_xlim(15, 85)
ax.set_xticks([25, 50, 75])
ax.set_xticklabels(["25","50","75"], fontsize=9)
ax.set_xlabel("GC", fontsize=10)
ax.set_ylabel("Normalized Coverage", fontsize=10)
ax.set_title(f"{SAMPLE_NAME} — GC Bias", fontsize=10, fontweight="normal")
ax.spines[["top","right"]].set_visible(False)

ax.text(0.05, 0.96, f"GC Bias: {severity}",
        transform=ax.transAxes, fontsize=9, va="top",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                  edgecolor="#cccccc", alpha=0.9))

plt.tight_layout(pad=0.8)
out2 = f"{OUT_DIR}/{SAMPLE_NAME}_gc_bias.png"
plt.savefig(out2, dpi=180, bbox_inches="tight")
plt.close()
print(f"  Saved: {out2}")

print(f"\nDone!\n  open {out1}\n  open {out2}")
