"""Show a two-dimensional Schur subspace converging in three dimensions."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


plt.rcParams.update({"font.size": 14, "svg.fonttype": "path"})
blue, teal = "#2563eb", "#0f766e"
orange, ink, gray = "#c26a13", "#334155", "#b6c2d1"
A = np.diag([1.05, 1., 0.65])
G = np.array([[1., 0.], [0., 1.], [0.9, 1.2]])
Q, _ = np.linalg.qr(G, mode="reduced")
steps = [0, 1, 3, 7]
iterates = [Q.copy()]
for _ in range(steps[-1]):
    Q, _ = np.linalg.qr(A @ Q, mode="reduced")
    iterates.append(Q.copy())

Pi = np.diag([1., 1., 0.])
angle = np.linspace(0, 2*np.pi, 160)
circle = np.array([np.cos(angle), np.sin(angle)])
target = np.vstack((circle, np.zeros_like(angle)))
fig = plt.figure(figsize=(11.2, 8.4))


def line(ax, points, **kwargs):
    ax.plot(points[0], points[1], points[2], **kwargs)


for panel, k in enumerate(steps, start=1):
    Q = iterates[k]
    d = np.linalg.norm((np.eye(3)-Pi) @ Q, 2)
    # The highest unit vector realizes the maximum distance to the target plane.
    highest = Q @ (Q[2, :] / d)
    foot = Pi @ highest
    assert np.isclose(np.linalg.norm(highest), 1.)
    assert np.isclose(np.linalg.norm(highest-foot), d)
    assert np.isclose(highest[2], d)
    expected_slopes = np.array([.9*(.65/1.05)**k, 1.2*.65**k])
    assert np.allclose(Q[2, :], expected_slopes @ Q[:2, :])
    assert np.allclose(Q.T @ Q, np.eye(2))
    assert np.isclose(d, np.linalg.norm(expected_slopes)
                      / np.sqrt(1+np.linalg.norm(expected_slopes)**2))
    assert d <= np.linalg.norm([.9, 1.2])*.65**k + 1e-14

    # Plot each plane independently of the signs and rotation of its QR basis.
    normal = np.cross(Q[:, 0], Q[:, 1])
    if normal[2] < 0:
        normal *= -1
    hinge = np.cross(np.array([0., 0., 1.]), normal)
    hinge /= np.linalg.norm(hinge)
    basis = np.column_stack((hinge, highest))
    assert np.allclose(basis.T @ basis, np.eye(2))
    current = basis @ circle

    ax = fig.add_subplot(2, 2, panel, projection="3d")
    ax.set_proj_type("ortho")
    ax.view_init(elev=26, azim=143)
    ax.set_box_aspect((1, 1, 1), zoom=1.4)
    ax.set_xlim(-1.15, 1.15)
    ax.set_ylim(-1.15, 1.15)
    ax.set_zlim(-1.15, 1.15)
    ax.set_axis_off()

    for j in range(3):
        end = 1.17*np.eye(3)[:, j]
        ax.quiver(0, 0, 0, *end, color=gray, linewidth=1.,
                  arrow_length_ratio=.06, zorder=0)
        ax.text(*(1.15*end), rf"$e_{{{j+1}}}$", color=ink, fontsize=13,
                ha="right" if j < 2 else "left",
                va="center" if j < 2 else "bottom")

    # Transparent, equal-radius patches reveal the orientation of both planes.
    for vertices, color, alpha in ((target, teal, .12), (current, blue, .23)):
        ax.add_collection3d(Poly3DCollection(
            [vertices.T], facecolor=color, edgecolor="none", alpha=alpha))
    line(ax, target, color=teal, linewidth=1.5, linestyle=(0, (4, 3)))
    line(ax, current, color=blue, linewidth=1.8)
    for j in range(2):
        line(ax, np.column_stack((-basis[:, j], basis[:, j])),
             color=blue, linewidth=.8, alpha=.4)
        axis = np.eye(3)[:, j]
        line(ax, np.column_stack((-axis, axis)), color=teal,
             linewidth=.7, alpha=.4)

    line(ax, np.column_stack((foot, highest)), color=orange,
         linewidth=2.4, linestyle=(0, (3, 2)))
    ax.scatter(*highest, color=blue, s=22, depthshade=False)
    ax.scatter(*foot, color=teal, s=20, depthshade=False)
    ax.scatter(0, 0, 0, color=ink, s=13, depthshade=False)
    title = "Start" if k == 0 else f"{k} step" + ("s" if k > 1 else "")
    ax.set_title(title, color=ink, fontsize=17, pad=0)
    ax.text2D(.5, .04, rf"$d_{{{k}}}={d:.3f}$", transform=ax.transAxes,
              ha="center", color=orange, fontsize=16)

legend = [Patch(facecolor=blue, alpha=.3, label="Current plane"),
          Patch(facecolor=teal, alpha=.2, label="True Schur plane"),
          Line2D([0], [0], color=orange, linewidth=2., linestyle="--",
                 label=r"Subspace error $d_k$")]
fig.legend(handles=legend, loc="upper center", bbox_to_anchor=(.5, .99),
           frameon=False, ncol=3, fontsize=14)
fig.subplots_adjust(left=.02, right=.98, bottom=.02, top=.9,
                    wspace=.02, hspace=.07)
destination = Path(__file__).resolve().parents[1] / "_static"
fig.savefig(destination / "subspace_convergence.svg", bbox_inches="tight")
fig.savefig("/tmp/nla-subspace-convergence.png", dpi=180, bbox_inches="tight")
plt.close(fig)
