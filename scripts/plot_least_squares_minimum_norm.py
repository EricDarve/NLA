"""Illustrate the minimum-norm point in an affine least-squares solution set."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import numpy as np


plt.rcParams.update({"font.size": 15, "svg.fonttype": "path"})
origin = np.zeros(2)
minimum = np.array([1.6, 0.])
z = np.array([0., 2.])
other = minimum + z

# This rank-one example matches the output-space projection used in the book.
A = np.diag([2., 0.])
b = np.array([3.2, 3.2])
assert np.allclose(A @ z, 0.)
assert np.allclose(A @ minimum, A @ other)
assert np.allclose(minimum, np.linalg.pinv(A) @ b)
assert np.isclose(minimum @ z, 0.)
assert np.isclose(np.linalg.norm(other) ** 2,
                  np.linalg.norm(minimum) ** 2 + np.linalg.norm(z) ** 2)

blue, purple, orange = "#2563eb", "#7c3aed", "#c26a13"
ink, gray = "#334155", "#94a3b8"
fig, ax = plt.subplots(figsize=(9, 6))

ax.plot([-1.9, 4.35], [0, 0], color=gray, linewidth=1.2, zorder=0)
ax.plot([0, 0], [-1.8, 3.25], color=gray, linewidth=1.6, zorder=0)
ax.plot([minimum[0], minimum[0]], [-1.25, 3.25],
        color="#64748b", linewidth=1.7, linestyle="--", zorder=0)
ax.add_patch(Circle(origin, np.linalg.norm(minimum),
                    fill=False, edgecolor=purple, linewidth=1.6,
                    linestyle=":", alpha=0.55, zorder=0))
ax.fill([0, minimum[0], other[0]], [0, 0, other[1]],
        color="#f1f5f9", zorder=0)


def arrow(start, end, color):
    ax.annotate("", xy=end, xytext=start,
                arrowprops={"arrowstyle": "-|>", "color": color,
                            "lw": 2.5, "mutation_scale": 19,
                            "shrinkA": 0, "shrinkB": 0}, zorder=3)


arrow(origin, other, blue)
arrow(origin, minimum, purple)
arrow(minimum, other, orange)
ax.scatter([0, minimum[0]], [0, 0], color=ink, s=23, zorder=4)

ax.text(1.83, 2.07, r"$x=x^*+z$", color=blue, fontsize=20)
ax.text(1.82, 0.92, r"$z$", color=orange, fontsize=21)
ax.text(1.05, -0.34, r"$x^*$", color=purple, fontsize=21, ha="center")
ax.text(-0.18, -0.35, r"$0$", color=ink, fontsize=15)
ax.text(-0.18, 3.1, r"$N(A)$", color=ink, fontsize=18, ha="right")
ax.text(3.3, -0.35, r"$\operatorname{range}(A^T)$", color=ink, fontsize=17, ha="center")
ax.annotate("All least-squares\nminimizers", xy=(1.6, 2.95), xytext=(2.15, 2.95),
            ha="left", va="center", fontsize=17, color=ink,
            arrowprops={"arrowstyle": "-", "color": gray, "lw": 1.2})
ax.text(-0.25, -1.94, r"Radius $\|x^*\|_2$", color=purple, fontsize=17, ha="center")

# The first marker identifies orthogonality to the null space at the origin.
# The second shows the Pythagorean decomposition of any other minimizer.
size = 0.18
ax.plot([size, size, 0], [0, -size, -size], color=ink, linewidth=1.4)
ax.plot([minimum[0] - size, minimum[0] - size, minimum[0]],
        [0, size, size], color=ink, linewidth=1.4)

ax.set_aspect("equal")
ax.set_xlim(-2., 4.5)
ax.set_ylim(-2.2, 3.5)
ax.axis("off")
fig.tight_layout(pad=0.5)

destination = Path(__file__).resolve().parents[1] / "_static"
fig.savefig(destination / "least_squares_minimum_norm.svg", bbox_inches="tight")
fig.savefig("/tmp/nla-least-squares-minimum-norm.png", dpi=180, bbox_inches="tight")
plt.close(fig)
