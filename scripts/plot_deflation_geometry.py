"""Illustrate eigenvalue preservation under orthogonal Schur deflation."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


plt.rcParams.update({"font.size": 15, "svg.fonttype": "path"})
blue, purple = "#2563eb", "#7c3aed"
orange, teal = "#c26a13", "#0f766e"
ink, gray = "#334155", "#94a3b8"

A = np.array([[4., 1.], [0., 2.]])
origin = np.zeros(2)
q1, q2 = np.eye(2).T
x2 = np.array([-1., 2.])
P = np.eye(2) - np.outer(q1, q1)
projected = P @ x2
image = A @ q2

assert np.allclose(A @ q1, 4 * q1)
assert np.allclose(A @ x2, 2 * x2)
assert np.allclose(projected / np.linalg.norm(projected), q2)
assert not np.allclose(image, 2 * q2)
assert np.allclose(P @ image, 2 * q2)
assert np.allclose(image, q1 + 2 * q2)


def arrow(ax, start, end, color, linewidth=2.5, linestyle="-", alpha=1.):
    ax.annotate(
        "", xy=end, xytext=start,
        arrowprops={"arrowstyle": "-|>", "color": color,
                    "lw": linewidth, "mutation_scale": 19,
                    "shrinkA": 0, "shrinkB": 0,
                    "linestyle": linestyle, "alpha": alpha},
        zorder=3,
    )


fig, (left, right) = plt.subplots(1, 2, figsize=(11.8, 5.5))
for ax in (left, right):
    ax.plot([-1.55, 1.6], [0, 0], color=gray, linewidth=1.1, zorder=0)
    ax.plot([0, 0], [-0.15, 2.65], color=gray, linewidth=1.1,
            linestyle=(0, (3, 3)), zorder=0)
    arrow(ax, origin, q1, teal, linewidth=3.)
    ax.text(0.75, -0.27, r"$q_{\star,1}$", color=teal,
            fontsize=19, ha="center")
    ax.text(0.08, 2.55, r"$\mathcal{S}_1^\perp$", color=ink, fontsize=16)
    ax.scatter([0], [0], s=22, color=ink, zorder=5)
    ax.text(-0.12, -0.21, "0", color=ink, fontsize=13, ha="right")
    # The known direction and its orthogonal complement meet at a right angle.
    ax.plot([0.16, 0.16, 0], [0, 0.16, 0.16], color=ink, linewidth=1.2)
    ax.set_aspect("equal")
    ax.set_xlim(-1.65, 1.85)
    ax.set_ylim(-0.6, 2.85)
    ax.axis("off")

left.set_title("Project the eigenvector, then normalize", fontsize=16,
               color=ink, pad=18)
left.fill([0, x2[0], 0], [0, x2[1], projected[1]], color=blue,
          alpha=0.045, zorder=0)
arrow(left, origin, x2, blue)
arrow(left, origin, projected, purple, linewidth=2., linestyle="--", alpha=0.65)
arrow(left, origin, q2, purple, linewidth=3.)
arrow(left, x2, projected, orange, linewidth=2.)
left.text(-1.12, 2.02, r"$x_2$", color=blue, fontsize=21,
          ha="right", va="center")
left.text(0.16, 1.95, r"$P_1x_2$", color=purple, fontsize=19)
left.text(0.16, 0.85, r"$q_{\star,2}$", color=purple, fontsize=19)
left.text(-0.5, 2.2, "project", color=orange, fontsize=14, ha="center")
left.text(0.1, -0.54, r"$P_1x_2=2q_{\star,2}$",
          color=ink, fontsize=18, ha="center")

right.set_title("Multiply the new direction, then project", fontsize=16,
                color=ink, pad=18)
right.fill([0, image[0], 0], [0, image[1], 2], color=blue,
           alpha=0.045, zorder=0)
right.plot([1, 1], [0, 2], color=gray, linewidth=1.1,
           linestyle=(0, (3, 3)), zorder=0)
arrow(right, origin, image, blue)
arrow(right, origin, 2 * q2, purple, linewidth=2.5)
arrow(right, image, 2 * q2, orange, linewidth=2.)
right.text(1.13, 2.02, r"$Aq_{\star,2}$", color=blue,
           fontsize=20, va="center")
right.text(-0.16, 1.95, r"$2q_{\star,2}$", color=purple,
           fontsize=20, ha="right")
right.text(0.5, 2.2, "project", color=orange, fontsize=14, ha="center")
right.text(0.1, -0.54, r"$P_1Aq_{\star,2}=2q_{\star,2}$",
           color=ink, fontsize=18, ha="center")

fig.tight_layout(pad=0.8, w_pad=3.)
destination = Path(__file__).resolve().parents[1] / "_static"
fig.savefig(destination / "deflation_geometry.svg", bbox_inches="tight")
fig.savefig("/tmp/nla-deflation-geometry.png", dpi=180, bbox_inches="tight")
plt.close(fig)
