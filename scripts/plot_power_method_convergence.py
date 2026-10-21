"""Show normalized power iterates aligning with the dominant eigenvector."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Arc


plt.rcParams.update({"font.size": 15, "svg.fonttype": "path"})
blue, teal = "#2563eb", "#0f766e"
ink, gray = "#334155", "#cbd5e1"

A = np.diag([1.0, 0.7])
initial_angle = np.deg2rad(70.0)
q = np.array([np.cos(initial_angle), np.sin(initial_angle)])
steps = [0, 2, 4, 6, 10]
iterates = [q.copy()]
for _ in range(steps[-1]):
    q = A @ q
    q /= np.linalg.norm(q)
    iterates.append(q.copy())

fig, axes = plt.subplots(1, len(steps), figsize=(12, 3.4))


def arrow(ax, endpoint, color):
    ax.annotate(
        "", xy=endpoint, xytext=(0, 0),
        arrowprops={"arrowstyle": "-|>", "color": color, "lw": 2.8,
                    "mutation_scale": 17, "shrinkA": 0, "shrinkB": 0},
        zorder=3,
    )


for ax, k in zip(axes, steps):
    q = iterates[k]
    angle = np.arctan2(q[1], q[0])
    # For this diagonal example the tangent contracts by exactly 0.7 per step.
    assert np.isclose(np.linalg.norm(q), 1.0)
    assert np.isclose(q[1] / q[0], 0.7**k * np.tan(initial_angle))

    unit_arc = np.linspace(0, np.pi / 2, 150)
    ax.plot(np.cos(unit_arc), np.sin(unit_arc), color=gray,
            linewidth=1.2, linestyle=(0, (3, 3)), zorder=0)
    ax.plot([0, 0], [-0.07, 1.06], color=gray, linewidth=0.9, zorder=0)

    wedge = np.linspace(0, angle, 80)
    ax.fill(np.r_[0, np.cos(wedge), 0],
            np.r_[0, np.sin(wedge), 0], color=blue, alpha=0.055, zorder=0)
    ax.add_patch(Arc((0, 0), 0.7, 0.7, theta1=0,
                     theta2=np.rad2deg(angle), color=blue, linewidth=1.2))
    arrow(ax, (1, 0), teal)
    arrow(ax, q, blue)
    ax.scatter([0], [0], s=15, color=ink, zorder=4)

    ax.text(*(1.13 * q + np.array([0.0, 0.03])),
            rf"$\mathbf{{q}}_{{{k}}}$", color=blue, fontsize=19,
            ha="center", va="bottom")
    ax.text(1.03, -0.12, r"$\mathbf{x}_1$", color=teal,
            fontsize=18, ha="center", va="top")
    ax.text(-0.07, -0.1, "0", color=ink, fontsize=12, ha="right")
    ax.set_title("Start" if k == 0 else f"{k} steps",
                 fontsize=17, color=ink, pad=14)
    ax.text(0.55, -0.43, rf"Angle: ${np.rad2deg(angle):.1f}^\circ$",
            color=ink, fontsize=15, ha="center")
    ax.set_xlim(-0.15, 1.35)
    ax.set_ylim(-0.55, 1.35)
    ax.set_aspect("equal")
    ax.axis("off")

fig.tight_layout(pad=0.6, w_pad=1.5)
destination = Path(__file__).resolve().parents[1] / "_static"
fig.savefig(destination / "power_method_convergence.svg", bbox_inches="tight")
fig.savefig("/tmp/nla-power-method-convergence.png", dpi=180, bbox_inches="tight")
plt.close(fig)
