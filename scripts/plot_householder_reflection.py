"""Generate the geometric illustration used in the Householder notes."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


plt.rcParams.update({"font.size": 14, "svg.fonttype": "path"})
fig, ax = plt.subplots(figsize=(8, 6))
x = np.array([3., 4.])
target = np.array([-5., 0.])
v = x - target
midpoint = (x + target) / 2
normal = v / np.linalg.norm(v)
tangent = np.array([-normal[1], normal[0]])

# The reflecting line passes through the origin and the midpoint.
assert np.isclose(np.linalg.norm(x), np.linalg.norm(target))
assert np.isclose(v @ midpoint, 0)
P = np.eye(2) - 2 * np.outer(v, v) / (v @ v)
assert np.allclose(P @ x, target)

ax.axhline(0, color="#cbd5e1", linewidth=1, zorder=0)
ax.axvline(0, color="#cbd5e1", linewidth=1, zorder=0)
ax.plot([-2.7, 1.15], [5.4, -2.3], color="#334155", linewidth=2)
ax.text(-2.65, 5.65, "Reflection line", color="#334155")

def arrow(start, end, color):
    ax.annotate("", xy=end, xytext=start,
                arrowprops={"arrowstyle": "-|>", "color": color,
                            "lw": 2.5, "mutation_scale": 19,
                            "shrinkA": 0, "shrinkB": 0})

arrow(np.zeros(2), x, "#2563eb")
arrow(np.zeros(2), target, "#7c3aed")
# Translating the difference vector makes the reflection geometry visible.
arrow(target, x, "#c26a13")
ax.text(3.25, 4.05, r"$x$", color="#2563eb", fontsize=19)
ax.text(-5.25, -.65, r"$\sigma e_1$", color="#7c3aed", fontsize=19)
ax.text(-3.7, 2., r"$v=x-\sigma e_1$", color="#c26a13", fontsize=17)
ax.text(4.1, -.4, r"$e_1$", color="#64748b")
ax.scatter([0], [0], s=24, color="#334155", zorder=5)
ax.text(-.4, -.48, r"$0$", color="#334155")

# A right-angle marker shows that v is normal to the reflecting line.
size = .3
corner = np.array([midpoint + size * tangent,
                   midpoint + size * (tangent + normal),
                   midpoint + size * normal])
ax.plot(corner[:, 0], corner[:, 1], color="#334155", linewidth=1.4)
ax.scatter(*midpoint, s=22, color="#334155", zorder=5)

ax.set_aspect("equal")
ax.set_xlim(-5.8, 4.6)
ax.set_ylim(-2.5, 6.1)
ax.axis("off")
fig.tight_layout(pad=.5)
destination = Path(__file__).resolve().parents[1] / "_static"
fig.savefig(destination / "householder_reflection.svg", bbox_inches="tight")
fig.savefig("/tmp/nla-householder-reflection.png", dpi=160, bbox_inches="tight")
plt.close(fig)
