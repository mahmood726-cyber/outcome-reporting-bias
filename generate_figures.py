#!/usr/bin/env python
"""Generate publication-quality figures for the ORB manuscript.

Produces three figures:
  1. Classification bar chart (Low/Moderate/High risk)
  2. Excess significance distribution histogram
  3. I-squared by ORB risk class (box plot)

Saves 300 dpi PNG + PDF to figures/ directory.
"""

import sys
import io
import os
import csv
import json

# Windows cp1252 safety
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# ── Paths ──────────────────────────────────────────────────────────────
BASE = os.path.dirname(os.path.abspath(__file__))
DATA_CSV = os.path.join(BASE, "data", "output", "orb_results.csv")
DATA_JSON = os.path.join(BASE, "data", "output", "orb_summary.json")
FIG_DIR = os.path.join(BASE, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

# ── Load data ──────────────────────────────────────────────────────────
with open(DATA_CSV, encoding="utf-8") as f:
    rows = list(csv.DictReader(f))

with open(DATA_JSON, encoding="utf-8") as f:
    summary = json.load(f)

print(f"Loaded {len(rows)} reviews from {DATA_CSV}")

# ── Colour palette (colour-blind safe) ─────────────────────────────────
C_LOW = "#4CAF50"       # green
C_MOD = "#FF9800"       # amber
C_HIGH = "#F44336"      # red
C_HIST = "#5C6BC0"      # indigo
PALETTE = [C_LOW, C_MOD, C_HIGH]
CLASS_LABELS = ["Low risk", "Moderate risk", "High risk"]
CLASS_KEYS = ["Low_Risk", "Moderate_Risk", "High_Risk"]

# ── Shared style ───────────────────────────────────────────────────────
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 11,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.dpi": 150,
})


def save(fig, name):
    """Save figure as 300 dpi PNG and PDF."""
    for ext in ("png", "pdf"):
        path = os.path.join(FIG_DIR, f"{name}.{ext}")
        fig.savefig(path, dpi=300, bbox_inches="tight", facecolor="white")
        print(f"  Saved {path}")
    plt.close(fig)


# ═══════════════════════════════════════════════════════════════════════
# FIGURE 1 — Classification bar chart
# ═══════════════════════════════════════════════════════════════════════
print("\nFigure 1: Classification bar chart")
counts = [summary["classification"][k] for k in CLASS_KEYS]
total = sum(counts)
pcts = [100 * c / total for c in counts]

fig1, ax1 = plt.subplots(figsize=(6, 4.5))
bars = ax1.bar(CLASS_LABELS, counts, color=PALETTE, edgecolor="white", width=0.6)

# Annotate bars with count and percentage
for bar, count, pct in zip(bars, counts, pcts):
    ax1.text(
        bar.get_x() + bar.get_width() / 2,
        bar.get_height() + 5,
        f"{count}\n({pct:.1f}%)",
        ha="center", va="bottom", fontsize=11, fontweight="bold",
    )

ax1.set_ylabel("Number of meta-analyses", fontsize=12)
ax1.set_title(
    "Figure 1. ORB risk classification of 403 Cochrane meta-analyses",
    fontsize=12, fontweight="bold", pad=12,
)
ax1.set_ylim(0, max(counts) * 1.25)
ax1.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
fig1.tight_layout()
save(fig1, "figure1_classification")


# ═══════════════════════════════════════════════════════════════════════
# FIGURE 2 — Excess significance histogram
# ═══════════════════════════════════════════════════════════════════════
print("\nFigure 2: Excess significance histogram")
es_vals = [float(r["excess_significance"]) for r in rows]

fig2, ax2 = plt.subplots(figsize=(7, 4.5))
n_bins = 40
ax2.hist(es_vals, bins=n_bins, color=C_HIST, edgecolor="white", alpha=0.85)
ax2.axvline(0, color="black", linestyle="--", linewidth=1.2, label="Zero (no excess)")
ax2.axvline(1, color=C_MOD, linestyle=":", linewidth=1.2, label="Excess = 1")
ax2.axvline(2, color=C_HIGH, linestyle=":", linewidth=1.2, label="Excess = 2")

# Annotation
mean_es = summary["excess_significance"]["mean"]
ax2.annotate(
    f"Mean = {mean_es:.2f}",
    xy=(mean_es, 0),
    xytext=(mean_es + 2.5, ax2.get_ylim()[1] * 0.7),
    arrowprops=dict(arrowstyle="->", color="black"),
    fontsize=10, fontstyle="italic",
)

ax2.set_xlabel("Excess significance (observed - expected significant studies)", fontsize=11)
ax2.set_ylabel("Number of meta-analyses", fontsize=11)
ax2.set_title(
    "Figure 2. Distribution of excess significance across 403 reviews",
    fontsize=12, fontweight="bold", pad=12,
)
ax2.legend(loc="upper right", fontsize=9, framealpha=0.9)
fig2.tight_layout()
save(fig2, "figure2_excess_significance")


# ═══════════════════════════════════════════════════════════════════════
# FIGURE 3 — I-squared by ORB risk class (box plot)
# ═══════════════════════════════════════════════════════════════════════
print("\nFigure 3: I-squared by ORB risk class")
i2_by_class = []
for key in CLASS_KEYS:
    vals = [float(r["I2"]) for r in rows if r["orb_class"] == key]
    i2_by_class.append(vals)

fig3, ax3 = plt.subplots(figsize=(6, 5))
bp = ax3.boxplot(
    i2_by_class,
    tick_labels=CLASS_LABELS,
    patch_artist=True,
    widths=0.5,
    showmeans=True,
    meanprops=dict(marker="D", markerfacecolor="white", markeredgecolor="black", markersize=6),
    medianprops=dict(color="black", linewidth=1.5),
    flierprops=dict(marker="o", markerfacecolor="grey", markersize=4, alpha=0.5),
)

for patch, colour in zip(bp["boxes"], PALETTE):
    patch.set_facecolor(colour)
    patch.set_alpha(0.7)

# Reference lines for conventional heterogeneity thresholds
for threshold, label in [(25, "Low"), (50, "Moderate"), (75, "Substantial")]:
    ax3.axhline(threshold, color="grey", linestyle="--", linewidth=0.8, alpha=0.6)
    ax3.text(3.55, threshold + 1, f"{threshold}%", fontsize=8, color="grey", va="bottom")

ax3.set_ylabel("I-squared (%)", fontsize=12)
ax3.set_ylim(-5, 105)
ax3.set_title(
    "Figure 3. Heterogeneity (I-squared) by ORB risk class",
    fontsize=12, fontweight="bold", pad=12,
)

# Add sample size annotations
for i, (key, n) in enumerate(zip(CLASS_KEYS, [len(v) for v in i2_by_class]), start=1):
    ax3.text(i, -3, f"n={n}", ha="center", fontsize=9, color="grey")

fig3.tight_layout()
save(fig3, "figure3_i2_by_class")


# ═══════════════════════════════════════════════════════════════════════
print(f"\nAll figures saved to {FIG_DIR}")
print("Done.")
