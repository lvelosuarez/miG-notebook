"""
Generate updated visual abstract for mSystems submission.

Key changes vs original:
- Result 3: "Database Selection Matters"  → "Classifiers Behave Differently"
- Result 4: "Cross-Cohort Integration Robust" → "meteor Converges; sylph Retains Asymmetry"
- Methods: added classifier row (meteor / sylph)
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
import numpy as np

# ── Canvas ──────────────────────────────────────────────────────────────────
FIG_W, FIG_H = 5.3, 8.5
fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
ax.set_xlim(0, FIG_W)
ax.set_ylim(0.75, FIG_H)
ax.axis("off")
fig.patch.set_facecolor("white")

# ── Palette ─────────────────────────────────────────────────────────────────
TAN       = "#F5E6C8"
D_BROWN   = "#4A3B2A"
BLUE      = "#4D7FAF"
ORANGE    = "#C4925A"
LT_BLUE   = "#E8F2FA"
LT_ORANGE = "#FBF3E8"
EDGE_TAN  = "#C8B090"
EDGE_BLUE = "#8AADCC"
EDGE_ORA  = "#C8A070"
WHITE     = "#FFFFFF"
GREY      = "#777777"


def rbox(ax, x, y, w, h, fc=WHITE, ec=EDGE_TAN, lw=1.5, r=0.18):
    """Draw a rounded rectangle in data coordinates."""
    patch = FancyBboxPatch(
        (x, y), w, h,
        boxstyle=f"round,pad=0,rounding_size={r}",
        facecolor=fc, edgecolor=ec, linewidth=lw,
        transform=ax.transData, clip_on=False,
    )
    ax.add_patch(patch)


def stick_person(ax, cx, cy, color, scale=0.22):
    """Draw a simple stick figure centred at (cx, cy)."""
    r = scale * 0.38          # head radius
    ax.add_patch(mpatches.Circle((cx, cy + scale * 0.55), r,
                                  facecolor=color, edgecolor=color, linewidth=0,
                                  transform=ax.transData, zorder=3))
    # body
    ax.plot([cx, cx], [cy + scale * 0.15, cy + scale * 0.45],
            color=color, lw=1.8, transform=ax.transData, zorder=3)
    # arms
    ax.plot([cx - scale * 0.32, cx + scale * 0.32],
            [cy + scale * 0.32, cy + scale * 0.32],
            color=color, lw=1.8, transform=ax.transData, zorder=3)
    # legs
    ax.plot([cx, cx - scale * 0.26], [cy + scale * 0.15, cy - scale * 0.10],
            color=color, lw=1.8, transform=ax.transData, zorder=3)
    ax.plot([cx, cx + scale * 0.26], [cy + scale * 0.15, cy - scale * 0.10],
            color=color, lw=1.8, transform=ax.transData, zorder=3)


# ═══════════════════════════════════════════════════════════════════════════
# SECTION 1 – BACKGROUND  (y = 6.5 → 8.5)
# ═══════════════════════════════════════════════════════════════════════════
SEC1_TOP = 8.35

# Header box
rbox(ax, 0.22, SEC1_TOP - 0.52, 2.00, 0.48, fc=TAN, ec=EDGE_TAN, lw=1.5)
ax.text(1.22, SEC1_TOP - 0.28, "Background",
        ha="center", va="center", fontsize=14, fontweight="bold", color=D_BROWN)

# Bullet text
bg = (
    "- Saliva commonly used\n"
    "  to sequence host DNA\n"
    "  in large-scale biobanks\n"
    "\n"
    "- Existing datasets also\n"
    "  contain microbial DNA,\n"
    "  but extracted using\n"
    "  host-focused protocols"
)
ax.text(0.28, SEC1_TOP - 0.70, bg,
        ha="left", va="top", fontsize=9.5, color=D_BROWN, linespacing=1.45)

# ── Mouth icon (right side) ──────────────────────────────────────────────
mx, my = 3.85, 7.75
# outer lips
ax.add_patch(mpatches.Ellipse((mx, my), 1.35, 0.68,
                               facecolor="#D93040", edgecolor="#A01020", lw=2,
                               transform=ax.transData, zorder=5))
# dark interior
ax.add_patch(mpatches.Ellipse((mx, my - 0.06), 1.00, 0.36,
                               facecolor="#111111",
                               transform=ax.transData, zorder=6))
# upper lip (covers top of dark area)
ax.add_patch(mpatches.Ellipse((mx, my + 0.14), 1.10, 0.38,
                               facecolor="#D93040",
                               transform=ax.transData, zorder=7))
# upper lip centre dip
ax.add_patch(mpatches.Ellipse((mx, my + 0.28), 0.38, 0.24,
                               facecolor="#D93040",
                               transform=ax.transData, zorder=8))
# tongue
ax.add_patch(mpatches.Ellipse((mx + 0.12, my - 0.10), 0.52, 0.30,
                               facecolor="#D04858",
                               transform=ax.transData, zorder=7))
# teeth hint
ax.add_patch(mpatches.Ellipse((mx, my + 0.08), 0.80, 0.15,
                               facecolor="#FFFFFF", alpha=0.6,
                               transform=ax.transData, zorder=7))
# saliva drop
ax.plot([mx + 0.12, mx + 0.12], [my - 0.26, my - 0.55],
        color="#90C0E8", lw=2.5, transform=ax.transData, zorder=8)
ax.add_patch(mpatches.Ellipse((mx + 0.12, my - 0.62), 0.14, 0.18,
                               facecolor="#90C0E8",
                               transform=ax.transData, zorder=8))

# ── Monitor / chart icon ─────────────────────────────────────────────────
mnx, mny = 3.85, 6.82
# screen
rbox(ax, mnx - 0.65, mny - 0.30, 1.30, 0.72, fc=WHITE, ec=BLUE, lw=2.0, r=0.08)
# chart lines
xv = np.linspace(mnx - 0.52, mnx + 0.52, 30)
ax.plot(xv, mny + 0.12 + 0.10 * np.sin(np.linspace(0, 3 * np.pi, 30)),
        color=BLUE, lw=2, transform=ax.transData)
ax.plot(xv, mny - 0.02 + 0.07 * np.sin(np.linspace(1.2, 4.2 * np.pi, 30)),
        color=ORANGE, lw=2, transform=ax.transData)
ax.plot(xv, mny + 0.20 + 0.05 * np.cos(np.linspace(0, 2 * np.pi, 30)),
        color="#80A878", lw=1.5, transform=ax.transData)
# stand
ax.plot([mnx, mnx], [mny - 0.30, mny - 0.52], color=BLUE, lw=2.5, transform=ax.transData)
ax.plot([mnx - 0.28, mnx + 0.28], [mny - 0.52, mny - 0.52], color=BLUE, lw=2.5, transform=ax.transData)

# Divider
ax.axhline(6.35, xmin=0.04, xmax=0.96, color=EDGE_TAN, lw=1.2)

# ═══════════════════════════════════════════════════════════════════════════
# SECTION 2 – METHODS  (y = 3.85 → 6.35)
# ═══════════════════════════════════════════════════════════════════════════
# Header
rbox(ax, 0.22, 5.85, 1.80, 0.48, fc=TAN, ec=EDGE_TAN, lw=1.5)
ax.text(1.12, 6.09, "Methods",
        ha="center", va="center", fontsize=14, fontweight="bold", color=D_BROWN)

# ── miG box ──────────────────────────────────────────────────────────────
rbox(ax, 0.22, 4.70, 2.30, 1.00, fc=WHITE, ec=EDGE_BLUE, lw=1.5, r=0.15)
for i in range(5):
    stick_person(ax, 0.52 + i * 0.38, 5.30, BLUE, scale=0.22)
ax.text(0.38, 4.88, "human genome sequencing", fontsize=7.5, color=D_BROWN)
ax.text(0.38, 4.76, "miG", fontsize=10, fontweight="bold", color=D_BROWN)
ax.text(2.38, 4.76, "n = 39", fontsize=9, color=D_BROWN, ha="right")

# ── ASAL box ─────────────────────────────────────────────────────────────
rbox(ax, 2.80, 4.70, 2.30, 1.00, fc=WHITE, ec=EDGE_ORA, lw=1.5, r=0.15)
for i in range(2):
    stick_person(ax, 3.08 + i * 0.38, 5.30, ORANGE, scale=0.22)
for i in range(2):
    stick_person(ax, 3.84 + i * 0.38, 5.30, BLUE, scale=0.22)
ax.text(2.96, 4.88, "microbiome-specific extraction", fontsize=7.5, color=D_BROWN)
ax.text(2.96, 4.76, "ASAL", fontsize=10, fontweight="bold", color=D_BROWN)
ax.text(4.96, 4.76, "n = 14", fontsize=9, color=D_BROWN, ha="right")

# ── Classifier row (NEW) ─────────────────────────────────────────────────
rbox(ax, 0.22, 3.90, 4.88, 0.72, fc=LT_BLUE, ec=EDGE_BLUE, lw=1.2, r=0.15)
ax.text(2.66, 4.58, "─── Classifiers ───",
        ha="center", va="center", fontsize=8, color=GREY)

# meteor
ax.add_patch(mpatches.Circle((0.68, 4.22), 0.14,
                               facecolor=BLUE, edgecolor=BLUE,
                               transform=ax.transData, zorder=5))
ax.text(0.68, 4.22, "m", ha="center", va="center",
        fontsize=9, fontweight="bold", color=WHITE, zorder=6)
ax.text(1.00, 4.31, "meteor", fontsize=10, fontweight="bold", color=BLUE, va="center")
ax.text(1.00, 4.14, "saliva MSP DB (853 genomes)", fontsize=8, color=GREY, va="center")

# vs divider
ax.text(2.66, 4.22, "vs", ha="center", va="center", fontsize=11, color=GREY, style="italic")

# sylph
ax.add_patch(mpatches.Circle((3.30, 4.22), 0.14,
                               facecolor=D_BROWN, edgecolor=D_BROWN,
                               transform=ax.transData, zorder=5))
ax.text(3.30, 4.22, "s", ha="center", va="center",
        fontsize=9, fontweight="bold", color=WHITE, zorder=6)
ax.text(3.62, 4.31, "sylph", fontsize=10, fontweight="bold", color=D_BROWN, va="center")
ax.text(3.62, 4.14, "GTDB R220 (k-mer)", fontsize=8, color=GREY, va="center", clip_on=False)

# Divider
ax.axhline(3.78, xmin=0.04, xmax=0.96, color=EDGE_TAN, lw=1.2)

# ═══════════════════════════════════════════════════════════════════════════
# SECTION 3 – RESULTS  (y = 0.12 → 3.78)
# ═══════════════════════════════════════════════════════════════════════════
# Header
rbox(ax, 0.22, 3.28, 1.72, 0.48, fc=TAN, ec=EDGE_TAN, lw=1.5)
ax.text(1.08, 3.52, "Results",
        ha="center", va="center", fontsize=14, fontweight="bold", color=D_BROWN)

# ── Result 1: Depth Drives Stability ─────────────────────────────────────
ry1 = 2.95
ax.annotate("", xy=(0.86, ry1 + 0.22), xytext=(0.86, ry1 - 0.06),
            arrowprops=dict(arrowstyle="<->", color=D_BROWN, lw=2.2))
ax.text(1.25, ry1 + 0.08, "Depth Drives Stability",
        ha="left", va="center", fontsize=11, color=D_BROWN)

# ── Result 2: Minimal Extraction Bias ────────────────────────────────────
ry2 = 2.42
# Syringe-like icon: a circle with a cross
ax.add_patch(mpatches.Circle((0.86, ry2), 0.18,
                               facecolor=WHITE, edgecolor=BLUE, lw=2.0,
                               transform=ax.transData, zorder=4))
ax.plot([0.86, 0.86], [ry2 - 0.18, ry2 - 0.38], color=BLUE, lw=2.2, transform=ax.transData)
ax.add_patch(mpatches.Circle((0.86, ry2 - 0.46), 0.10,
                               facecolor=BLUE, edgecolor=BLUE, lw=0,
                               transform=ax.transData, zorder=4))
ax.plot([0.68, 1.04], [ry2, ry2], color=BLUE, lw=2.2, transform=ax.transData)
ax.text(1.25, ry2, "Minimal Extraction Bias",
        ha="left", va="center", fontsize=11, color=D_BROWN)

# ── Result 3: Classifiers Behave Differently (NEW) ───────────────────────
ry3 = 1.82
# Balance / scale icon: two arms and a pivot
ax.plot([0.54, 1.18], [ry3, ry3], color=D_BROWN, lw=2.2, transform=ax.transData)  # beam
ax.plot([0.86, 0.86], [ry3 + 0.20, ry3 - 0.12], color=D_BROWN, lw=2.2, transform=ax.transData)  # pivot
ax.add_patch(mpatches.Circle((0.86, ry3 + 0.22), 0.06,
                               facecolor=D_BROWN, transform=ax.transData, zorder=5))
# left pan (meteor – higher / converging)
ax.plot([0.54, 0.54], [ry3, ry3 - 0.16], color=D_BROWN, lw=1.5, transform=ax.transData)
ax.plot([0.42, 0.66], [ry3 - 0.16, ry3 - 0.16], color=BLUE, lw=2.0, transform=ax.transData)
ax.text(0.54, ry3 - 0.30, "m", ha="center", fontsize=8, fontweight="bold", color=BLUE)
# right pan (sylph – lower / retaining asymmetry)
ax.plot([1.18, 1.18], [ry3, ry3 - 0.26], color=D_BROWN, lw=1.5, transform=ax.transData)
ax.plot([1.06, 1.30], [ry3 - 0.26, ry3 - 0.26], color=D_BROWN, lw=2.0, transform=ax.transData)
ax.text(1.18, ry3 - 0.40, "s", ha="center", fontsize=8, fontweight="bold", color=D_BROWN)
ax.text(1.25, ry3, "Classifiers Behave",
        ha="left", va="center", fontsize=11, color=D_BROWN)
ax.text(1.25, ry3 - 0.28, "Differently",
        ha="left", va="center", fontsize=11, color=D_BROWN)

# ── Result 4: meteor Converges; sylph Retains Asymmetry (NEW) ────────────
ry4 = 1.08
# Split arrow icon: two paths from same start
# meteor → convergence (arrow pointing to centre dot)
ax.annotate("", xy=(0.86, ry4 + 0.14), xytext=(0.58, ry4 + 0.28),
            arrowprops=dict(arrowstyle="-|>", color=BLUE, lw=2.0))
ax.annotate("", xy=(0.86, ry4 + 0.14), xytext=(0.58, ry4 + 0.00),
            arrowprops=dict(arrowstyle="-|>", color=BLUE, lw=2.0))
ax.add_patch(mpatches.Circle((0.86, ry4 + 0.14), 0.08,
                               facecolor=BLUE, transform=ax.transData, zorder=5))
# sylph → divergence (arrows spreading out)
ax.annotate("", xy=(1.18, ry4 + 0.28), xytext=(1.10, ry4 + 0.14),
            arrowprops=dict(arrowstyle="-|>", color=D_BROWN, lw=2.0))
ax.annotate("", xy=(1.18, ry4 + 0.00), xytext=(1.10, ry4 + 0.14),
            arrowprops=dict(arrowstyle="-|>", color=D_BROWN, lw=2.0))
ax.add_patch(mpatches.Circle((1.10, ry4 + 0.14), 0.06,
                               facecolor=D_BROWN, transform=ax.transData, zorder=5))

ax.text(1.35, ry4 + 0.22, "meteor Converges;",
        ha="left", va="center", fontsize=11, color=BLUE)
ax.text(1.35, ry4 - 0.02, "sylph Retains Asymmetry",
        ha="left", va="center", fontsize=11, color=D_BROWN)

# ═══════════════════════════════════════════════════════════════════════════
# Save
# ═══════════════════════════════════════════════════════════════════════════
out_base = (
    "/Users/lourdes/Dropbox/Curro/Projects/open/GazelADN/miG-notebook/visual_abstract"
)
plt.tight_layout(pad=0.1)
plt.savefig(out_base + ".png", dpi=150, bbox_inches="tight", facecolor="white")
plt.savefig(out_base + ".svg", bbox_inches="tight", facecolor="white")
print(f"Saved {out_base}.png and {out_base}.svg")
