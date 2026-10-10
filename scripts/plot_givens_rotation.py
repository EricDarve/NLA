"""Generate the geometric illustration used in the Givens rotation notes."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


plt.rcParams.update({"font.size": 14, "svg.fonttype": "path"})
fig, ax = plt.subplots(figsize=(7, 5.5))
x = np.array([3., 4.])
r = np.linalg.norm(x)
theta = np.arctan2(x[1], x[0])
c, s = x / r
G = np.array([[c, s], [-s, c]])
assert np.allclose(G.T @ G, np.eye(2))
assert np.allclose(G @ x, [r, 0])

ax.axhline(0, color="#cbd5e1", linewidth=1)
ax.axvline(0, color="#cbd5e1", linewidth=1)

def arrow(start, end, color):
    ax.annotate("", xy=end, xytext=start,
                arrowprops={"arrowstyle": "-|>", "color": color,
                            "lw": 2.5, "mutation_scale": 19,
                            "shrinkA": 0, "shrinkB": 0})

angles = np.linspace(theta, 0, 120)
arc = r * np.column_stack([np.cos(angles), np.sin(angles)])
ax.plot(arc[:, 0], arc[:, 1], color="#c26a13", linewidth=2)
arrow(arc[-6], arc[-1], "#c26a13")
arrow(np.zeros(2), x, "#2563eb")
arrow(np.zeros(2), np.array([r, 0.]), "#7c3aed")

# The small arc marks the angle; the outer arc shows the clockwise motion.
angle_arc = 1.35 * np.column_stack([np.cos(angles), np.sin(angles)])
ax.plot(angle_arc[:, 0], angle_arc[:, 1], color="#64748b", linewidth=1.4)
ax.text(1.45, .55, r"$\theta$", color="#64748b", fontsize=17)
ax.text(2.35, 4.35, r"$x=(a,b)^T$", color="#2563eb", fontsize=18)
ax.text(3.35, -.65, r"$Gx=r e_1$", color="#7c3aed", fontsize=18)
ax.text(4.7, 2.85, r"$G$", color="#c26a13", fontsize=19)
ax.text(5.6, -.25, r"$e_1$", color="#64748b")
ax.text(-.6, 4.9, r"$e_2$", color="#64748b")
ax.scatter([0], [0], s=24, color="#334155", zorder=5)
ax.text(-.4, -.45, r"$0$", color="#334155")

ax.set_aspect("equal")
ax.set_xlim(-.8, 6.)
ax.set_ylim(-.85, 5.4)
ax.axis("off")
fig.tight_layout(pad=.5)
destination = Path(__file__).resolve().parents[1] / "_static"
fig.savefig(destination / "givens_rotation.svg", bbox_inches="tight")
fig.savefig("/tmp/nla-givens-rotation.png", dpi=160, bbox_inches="tight")
plt.close(fig)
