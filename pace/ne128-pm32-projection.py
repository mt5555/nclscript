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
import math

# ---------------------------------------------------------------------------
# Data: (label, face_color, x_start, x_end, y_start, y_end)
# x / y values are MPI processor ranks (0-based).
# CPL and ATM share the same rank range; ATM is drawn on top with a narrower
# x-extent so both are visible simultaneously.
# ---------------------------------------------------------------------------
TOTAL_NODES=32
TOTAL_PROCS = TOTAL_NODES*64    # GPU nodes only have 64 cores per node
TOTAL_GPUS = TOTAL_NODES*4      # PM-GPU, do we use 8 per node or 4? 
GPU_SCALE = 24  # GPU size = size of 24 CPUs

ATM_GPUS=16*4
OCN_GPUS=16*4

t_ice1=0
t_ice2=72.2

t_lnd1=t_ice2
t_lnd2=t_lnd1 + 3.6

t_ocn1=t_lnd2
t_ocn2=t_ocn1 + 199.7 / 2  # assuming 2x faster

t_atm1=t_lnd2
t_atm2=t_atm1+187.6 /2    # assuming 2x faster

t_cpl1 = max(t_atm2,t_ocn2)
t_cpl2 = t_cpl1 + 21.5

total = 307.2
#if total > t_cpl2:
#  print("Adding other overheads to CPL: original=",t_cpl2," new=",total)
#  t_cpl2=total


# atm is using 16 nodes = 64 GPUs, and 64 cores
# ocn is using 8 nodes = 32 GPUs

c_atm1=TOTAL_PROCS
c_atm2=c_atm1 + ATM_GPUS*GPU_SCALE
c_ocn1=c_atm2
c_ocn2=c_ocn1 + OCN_GPUS*GPU_SCALE


  

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


# atm cpus:  64 1 cpu wide boxes:
for x0 in range(0, TOTAL_PROCS, round(TOTAL_PROCS/ATM_GPUS)):
  x1=x0+1
  y0=t_atm1
  y1=t_atm2
  color="#87CEEB"
  rect = mpatches.Rectangle(
    (x0, y0), x1 - x0, y1 - y0,
    linewidth=0.0,
    edgecolor=color,
    facecolor=color,
    zorder=2,
  )
  ax.add_patch(rect)

for x00 in range(0, TOTAL_PROCS, round(TOTAL_PROCS/OCN_GPUS)):
  x0=x00+2
  x1=x0+1
  y0=t_ocn1
  y1=t_ocn2
  color="#8888FF"
  rect = mpatches.Rectangle(
    (x0, y0), x1 - x0, y1 - y0,
    linewidth=0.0,
    edgecolor=color,
    facecolor=color,
    zorder=2,
  )
  ax.add_patch(rect)


  



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
ax.set_xlim(0, TOTAL_PROCS + TOTAL_GPUS*GPU_SCALE)
ax.set_ylim(0, t_cpl2)

# Y-axis ticks match the component boundaries visible in the original plot
yticks = [t_ice1, t_ocn1, t_cpl1, t_cpl2]
ax.set_yticks(yticks)
ax.set_yticklabels([f"{v:.1f}" for v in yticks])

xticks = [0, 1024, 2048, 3072, 4096]
ax.set_xticks(xticks)
ax.set_xticklabels([f"{v:.0f}" for v in xticks])

ax.set_xlabel("Processor Cores (64 per node)", fontsize=12)
ax.xaxis.set_label_coords(0.2,-.05)
ax.set_ylabel("Seconds per Model Day", fontsize=12)

# Top x-axis: (cores - midpoint) / GPU_SCALE, so the centre of the plot is 0.
_mid = TOTAL_PROCS
ax_top = ax.secondary_xaxis(
    "top",
    functions=(
        lambda cores: (cores - _mid) / GPU_SCALE,
        lambda gpu_units: gpu_units * GPU_SCALE + _mid,
    ),
)
top_ticks = [0, 128, 256]
ax_top.set_xticks(top_ticks)
ax_top.set_xticklabels([f"{v:.0f}" for v in top_ticks])
ax_top.set_xlabel(f"                                                                   GPUs (4 per node)", fontsize=12)

ax.tick_params(axis="both", labelsize=10)

plt.tight_layout()
plt.savefig("ne128-pm32.png", dpi=300, bbox_inches="tight")
plt.show()
