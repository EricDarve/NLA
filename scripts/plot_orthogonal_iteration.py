"""Show both columns evolving in the deflation example."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


plt.rcParams.update({"font.size": 15, "svg.fonttype": "path"})
blue, purple = "#2563eb", "#7c3aed"
ink, gray = "#334155", "#cbd5e1"
A = np.array([[4., 1.], [0., 2.]])
Q0 = np.array([[1., -1.], [1., 1.]]) / np.sqrt(2.)
Q = Q0.copy()
iterates = [Q.copy()]
for _ in range(6):
    Q, R = np.linalg.qr(A @ Q)
    phase = np.sign(np.diag(R))
    Q = Q * phase
    iterates.append(Q.copy())

assert np.allclose(iterates[1], np.array([[5., -2.], [2., 5.]])/np.sqrt(29.))
steps = [0, 1, 3, 6]
fig, axes = plt.subplots(1, len(steps), figsize=(12, 3.9))


def arrow(ax, vector, color, linewidth=2.8, linestyle="-"):
    ax.annotate(
        "", xy=vector, xytext=(0, 0),
        arrowprops={"arrowstyle": "-|>", "color": color,
                    "lw": linewidth, "mutation_scale": 18,
                    "shrinkA": 0, "shrinkB": 0, "linestyle": linestyle},
        zorder=3,
    )


for ax, k in zip(axes, steps):
    Q = iterates[k]
    assert np.allclose(Q.T @ Q, np.eye(2))
    direction = np.linalg.matrix_power(A, k) @ Q0[:, 0]
    direction /= np.linalg.norm(direction)
    assert np.allclose(Q[:, 0], direction)
    circle = np.linspace(0, np.pi, 160)
    ax.plot(np.cos(circle), np.sin(circle), color=gray, linewidth=1.1,
            linestyle=(0, (3, 3)), zorder=0)
    arrow(ax, (1, 0), gray, linewidth=1.3, linestyle="--")
    arrow(ax, (0, 1), gray, linewidth=1.3, linestyle="--")
    for j, color in enumerate((blue, purple)):
        arrow(ax, Q[:, j], color)
        label = 1.17 * Q[:, j] + np.array([0, 0.04])
        ax.text(*label, rf"$q_{{{k},{j+1}}}$", color=color,
                fontsize=18, ha="center", va="bottom")
    corner = np.array([0.14*Q[:, 0], 0.14*(Q[:, 0]+Q[:, 1]), 0.14*Q[:, 1]])
    ax.plot(corner[:, 0], corner[:, 1], color=ink, linewidth=1.1)
    ax.scatter([0], [0], s=17, color=ink, zorder=4)
    ax.text(-0.07, -0.22, "0", color=ink, fontsize=12)
    ax.set_title("Start" if k == 0 else f"{k} step" + ("s" if k > 1 else ""),
                 color=ink, fontsize=17, pad=17)
    ax.set_aspect("equal")
    ax.set_xlim(-1.35, 1.4)
    ax.set_ylim(-0.35, 1.55)
    ax.axis("off")

fig.tight_layout(pad=0.6, w_pad=2.)
destination = Path(__file__).resolve().parents[1] / "_static"
fig.savefig(destination / "orthogonal_iteration.svg", bbox_inches="tight")
fig.savefig("/tmp/nla-orthogonal-iteration.png", dpi=180, bbox_inches="tight")
plt.close(fig)
