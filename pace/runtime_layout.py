"""
Reproduce the PACE runtime processor layout plot.

The chart shows which MPI processor ranks are assigned to each component
of a coupled climate model (e.g. E3SM/CESM).  Each component is drawn as
a filled rectangle whose y-extent gives the processor-rank range for that
component.  CPL (the coupler) spans all ranks; ATM shares the same rank
range as CPL but is drawn with a narrower width to make it visible.

Processor layout for PACE run 224501
-------------------------------------
Component  | color        | ranks (y)       | x-width
-----------+--------------+-----------------+-----------
ICE        | cyan         |    0 –  356     | full
LND        | lime green   |  356 –  375     | full
OCN        | blue-purple  | 1374 – 2312     | full
CPL        | orange       | 2312 – 2419     | full
ATM        | light blue   | 2312 – 2419     | 107 (narrow, drawn over CPL)
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ---------------------------------------------------------------------------
# Data: (label, face_color, x_start, x_end, y_start, y_end)
# x / y values are MPI processor ranks (0-based).
# CPL and ATM share the same rank range; ATM is drawn on top with a narrower
# x-extent so both are visible simultaneously.
# ---------------------------------------------------------------------------
TOTAL_PROCS = 2419

components = [
    # name    color          x0     x1           y0    y1
    ("ICE",  "cyan",          0, TOTAL_PROCS,    0,   356),
    ("LND",  "#00CC00",       0, TOTAL_PROCS,  356,   375),
    ("OCN",  "#8888FF",       0, TOTAL_PROCS, 1374,  2312),
    ("CPL",  "orange",        0, TOTAL_PROCS, 2312,  2419),
    # ATM shares the CPL rank range; its narrower x-width distinguishes it
    ("ATM",  "#87CEEB",       0,         107, 2312,  2419),
]

# ---------------------------------------------------------------------------
# Build figure
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(10, 7))

for name, color, x0, x1, y0, y1 in components:
    rect = mpatches.Rectangle(
        (x0, y0), x1 - x0, y1 - y0,
        linewidth=0.8,
        edgecolor="black",
        facecolor=color,
        zorder=2,
    )
    ax.add_patch(rect)

    # Place the label at the visual centre of the rectangle
    cx = (x0 + x1) / 2
    cy = (y0 + y1) / 2
    ax.text(
        cx, cy, name,
        ha="center", va="center",
        fontsize=14, fontweight="bold",
        zorder=3,
    )

# ---------------------------------------------------------------------------
# Axes formatting
# ---------------------------------------------------------------------------
ax.set_xlim(0, TOTAL_PROCS)
ax.set_ylim(0, TOTAL_PROCS)

# Y-axis ticks match the component boundaries visible in the original plot
yticks = [0, 356, 375, 1374, 1536, 2312, 2419]
ax.set_yticks(yticks)
ax.set_yticklabels([f"{v:.1f}" for v in yticks])

ax.set_xlabel("Processors#", fontsize=12)

ax.tick_params(axis="both", labelsize=10)

plt.tight_layout()
plt.savefig("runtime_layout.png", dpi=100, bbox_inches="tight")
plt.show()
