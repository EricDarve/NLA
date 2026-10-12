"""Generate the projection diagram used in the normal-equations notes."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


plt.rcParams.update({"font.size": 15, "svg.fonttype": "path"})
fig, ax = plt.subplots(figsize=(9, 5))

origin = np.zeros(2)
p = np.array([3.2, 0.0])
b = np.array([3.2, 3.2])
candidate = np.array([6.0, 0.0])
residual = b - p
displacement = p - candidate

# These coordinates depict a one-dimensional column space in R^2.
assert residual @ displacement == 0
assert np.isclose(
    np.linalg.norm(b - candidate) ** 2,
    np.linalg.norm(residual) ** 2 + np.linalg.norm(displacement) ** 2,
)

blue = "#2563eb"
purple = "#7c3aed"
orange = "#c26a13"
teal = "#0f766e"
gray = "#64748b"

ax.fill(
    [p[0], candidate[0], b[0]],
    [p[1], candidate[1], b[1]],
    color="#f1f5f9",
    zorder=0,
)
ax.plot([-0.7, 7.5], [0, 0], color="#94a3b8", linewidth=1.5, zorder=1)


def arrow(start, end, color, linestyle="-"):
    ax.annotate(
        "",
        xy=end,
        xytext=start,
        arrowprops={
            "arrowstyle": "-|>",
            "color": color,
            "lw": 2.5,
            "mutation_scale": 19,
            "shrinkA": 0,
            "shrinkB": 0,
            "linestyle": linestyle,
        },
        zorder=3,
    )


arrow(origin, b, blue)
arrow(origin, p, purple)
arrow(p, b, orange)
arrow(candidate, b, gray, "--")
arrow(candidate, p, teal)

ax.text(3.2, 3.5, r"$b$", color=blue, fontsize=21, ha="center")
ax.text(3.4, 0.95, r"$r^*=b-p$", color=orange, fontsize=18)
ax.text(5.15, 1.8, r"$b-Ax$", color=gray, fontsize=18)
ax.text(3.2, -0.42, r"$p=Ax^*$", color=purple, fontsize=18, ha="center")
ax.text(6.0, -0.42, r"$Ax$", color=gray, fontsize=18, ha="center")
ax.text(4.65, -0.85, r"$p-Ax$", color=teal, fontsize=17, ha="center")
ax.text(6.35, 0.23, r"$\operatorname{range}(A)$", color="#334155", fontsize=17)
ax.text(-0.2, -0.42, r"$0$", color="#334155", fontsize=16)
ax.scatter([0, p[0], candidate[0]], [0, 0, 0], s=25, color="#334155", zorder=4)

# Mark the perpendicular residual and in-subspace displacement.
size = 0.25
ax.plot(
    [p[0], p[0] + size, p[0] + size],
    [size, size, 0],
    color="#334155",
    linewidth=1.5,
    zorder=4,
)

ax.set_aspect("equal")
ax.set_xlim(-0.8, 8.05)
ax.set_ylim(-1.05, 3.95)
ax.axis("off")
fig.tight_layout(pad=0.4)

destination = Path(__file__).resolve().parents[1] / "_static"
fig.savefig(destination / "normal_equations_projection.svg", bbox_inches="tight")
fig.savefig("/tmp/nla-normal-equations-projection.png", dpi=170, bbox_inches="tight")
plt.close(fig)
