"""Illustrate projection removal and normalization in Gram-Schmidt."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


plt.rcParams.update({"font.size": 15, "svg.fonttype": "path"})
origin = np.zeros(2)
q1 = np.array([1., 0.])
a2 = np.array([2.6, 2.2])
r12 = q1 @ a2
projection = r12 * q1
remainder = a2 - projection
q2 = remainder / np.linalg.norm(remainder)

assert np.isclose(q1 @ remainder, 0.)
assert np.isclose(q1 @ q2, 0.)
assert np.isclose(np.linalg.norm(q1), 1.)
assert np.isclose(np.linalg.norm(q2), 1.)
assert np.allclose(a2, projection + np.linalg.norm(remainder) * q2)

blue, purple = "#2563eb", "#7c3aed"
orange, teal = "#c26a13", "#0f766e"
ink, gray = "#334155", "#94a3b8"
fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.8))


def arrow(ax, start, end, color, linestyle="-", linewidth=2.5, alpha=1.):
    ax.annotate(
        "", xy=end, xytext=start,
        arrowprops={"arrowstyle": "-|>", "color": color,
                    "lw": linewidth, "mutation_scale": 19,
                    "shrinkA": 0, "shrinkB": 0,
                    "linestyle": linestyle, "alpha": alpha},
        zorder=3,
    )


left, right = axes
left.set_title(r"Remove the component along $q_1$", fontsize=16, pad=17, color=ink)
left.plot([-0.35, 3.7], [0, 0], color=gray, linewidth=1.2, zorder=0)
left.fill([0, projection[0], a2[0]], [0, 0, a2[1]], color="#f1f5f9", zorder=0)
arrow(left, origin, projection, orange)
arrow(left, origin, a2, blue)
# Translating the remainder to the projection's tip shows vector subtraction.
arrow(left, projection, a2, purple)
arrow(left, origin, q1, teal, linewidth=3.)
left.text(2.62, 2.43, r"$a_2$", color=blue, fontsize=21, ha="center")
left.text(2.82, 1.05, r"$a_2^{(1)}$", color=purple, fontsize=21)
left.text(1.9, -0.36, r"$r_{12}q_1$", color=orange, fontsize=18, ha="center")
left.text(0.63, -0.36, r"$q_1$", color=teal, fontsize=20, ha="center")
left.text(-0.16, -0.31, r"$0$", color=ink, fontsize=14)
size = 0.2
left.plot([projection[0] - size, projection[0] - size, projection[0]],
          [0, size, size], color=ink, linewidth=1.4)
left.text(1.7, -0.85, r"$a_2^{(1)}=a_2-(q_1^Ta_2)q_1$",
          ha="center", color=ink, fontsize=17)

right.set_title("Normalize the perpendicular remainder", fontsize=16, pad=17, color=ink)
right.plot([-0.45, 1.6], [0, 0], color=gray, linewidth=1.2, zorder=0)
arrow(right, origin, remainder, purple, linestyle="--", linewidth=1.8, alpha=0.5)
arrow(right, origin, q1, teal, linewidth=3.)
arrow(right, origin, q2, purple, linewidth=3.)
right.text(0.18, 2.21, r"$a_2^{(1)}$", color=purple, fontsize=21, va="center")
right.text(0.22, 0.77, r"$q_2$", color=purple, fontsize=21)
right.text(1.15, 0.04, r"$q_1$", color=teal, fontsize=21)
right.text(1., -0.3, "1", color=teal, fontsize=14, ha="center")
right.text(-0.26, 1., "1", color=purple, fontsize=14, ha="center", va="center")
right.text(-0.17, -0.31, r"$0$", color=ink, fontsize=14)
right.plot([size, size, 0], [0, size, size], color=ink, linewidth=1.4)
right.text(0.6, -0.85, r"$q_2=a_2^{(1)}/\|a_2^{(1)}\|_2$",
           ha="center", color=ink, fontsize=17)

for ax in axes:
    ax.scatter([0], [0], s=22, color=ink, zorder=4)
    ax.set_aspect("equal")
    ax.set_ylim(-1.05, 2.65)
    ax.axis("off")
left.set_xlim(-0.45, 4.05)
right.set_xlim(-1.65, 2.85)
fig.tight_layout(pad=0.8, w_pad=2.5)

destination = Path(__file__).resolve().parents[1] / "_static"
fig.savefig(destination / "gram_schmidt_geometry.svg", bbox_inches="tight")
fig.savefig("/tmp/nla-gram-schmidt-geometry.png", dpi=180, bbox_inches="tight")
plt.close(fig)
