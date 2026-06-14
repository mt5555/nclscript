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
TOTAL_PROCS = 4096
CORES_PER_GPU = 64
TOTAL_GPUS = TOTAL_PROCS // CORES_PER_GPU

t_ice1=0
t_ice2=72.2

t_lnd1=t_ice2
t_lnd2=t_lnd1 + 3.6

t_ocn1=t_lnd2
t_ocn2=t_ocn1 + 199.7

t_atm1=t_lnd2
t_atm2=t_atm1+187.6

t_cpl1 = max(t_atm2,t_ocn2)
t_cpl2 = t_cpl1 + 106.3

total = 307.2
if total > t_cpl2:
  print("Adding other overheads to CPL: original=",t_cpl2," new=",total)
  t_cpl2=total
  

c_ocn1=0
c_ocn2=TOTAL_PROCS-256
c_atm1=c_ocn2
c_atm2=c_atm1+256


  

components = [
    # name    color          x0     x1           y0    y1
    ("ICE",  "cyan",          0, TOTAL_PROCS, t_ice1,  t_ice2),
    ("LND",  "#00CC00",       0, TOTAL_PROCS, t_lnd1,  t_lnd2),
    ("OCN",  "#8888FF",       c_ocn1, c_ocn2, t_ocn1,  t_ocn2),
    ("CPL",  "orange",        0, TOTAL_PROCS, t_cpl1,  t_cpl2),
    # ATM shares the CPL rank range; its narrower x-width distinguishes it
    ("ATM",  "#87CEEB",       c_atm1,  c_atm2, t_atm1,  t_atm2),
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
ax.set_ylim(0, t_cpl2)

# Y-axis ticks match the component boundaries visible in the original plot
yticks = [t_ice1, t_ocn1, t_cpl1, t_cpl2]
ax.set_yticks(yticks)
ax.set_yticklabels([f"{v:.1f}" for v in yticks])

xticks = [0, 1024, 2048, 3072, 4096]
ax.set_xticks(xticks)
ax.set_xticklabels([f"{v:.0f}" for v in xticks])

ax.set_xlabel("Processor Cores", fontsize=12)
ax.set_ylabel("Seconds per Model Day", fontsize=12)

# Top x-axis in equivalent GPU counts.
ax_top = ax.secondary_xaxis(
    "top",
    functions=(
        lambda cores: cores / CORES_PER_GPU,
        lambda gpus: gpus * CORES_PER_GPU,
    ),
)
gpu_ticks = [0, 16, 32, 48, TOTAL_GPUS]
ax_top.set_xticks(gpu_ticks)
ax_top.set_xticklabels([f"{v:.0f}" for v in gpu_ticks])
ax_top.set_xlabel("GPU Count", fontsize=12)

ax.tick_params(axis="both", labelsize=10)

plt.tight_layout()
plt.savefig("runtime_layout.png", dpi=100, bbox_inches="tight")
plt.show()
