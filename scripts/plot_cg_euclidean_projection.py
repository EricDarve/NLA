"""Illustrate the Euclidean trial construction preceding the CG derivation.

The plotted projections are calculated from an SPD system. They require the
exact solution and are not CG iterates. Requires NumPy and Matplotlib.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import numpy as np

plt.rcParams.update({"font.size": 15, "svg.fonttype": "path",
                     "svg.hashsalt": "cg-euclidean-projection"})
BLUE, PURPLE, ORANGE = "#2563eb", "#7c3aed", "#c26a13"
TEAL, INK, GRAY = "#0f766e", "#334155", "#cbd5e1"


def calculate():
    a = np.array([[2., -1.], [-1., 2.]])
    b = np.array([3., 0.])
    x = np.linalg.solve(a, b)
    p1 = b/np.linalg.norm(b)
    new = a @ p1 - p1*(p1 @ a @ p1)
    p2 = -new/np.linalg.norm(new)  # Choose the upward sign for readability.
    p = np.column_stack((p1, p2))
    mu = p.T @ x
    x1 = mu[0]*p1
    x2 = x1+mu[1]*p2
    np.testing.assert_allclose(p.T @ p, np.eye(2))
    np.testing.assert_allclose(x1, [2., 0.])
    np.testing.assert_allclose(x2, x)
    np.testing.assert_allclose(p1 @ (x-x1), 0.)
    assert np.linalg.matrix_rank(np.column_stack((b, a @ b))) == 2
    # Verify the perpendicular projection is closer than other points on K_1.
    for t in np.linspace(-1, 3, 30):
        np.testing.assert_allclose(np.linalg.norm(x-t*p1)**2,
                                   np.linalg.norm(x-x1)**2+(t-mu[0])**2)
    return x, x1, x2, mu


def arrow(ax, start, end, color, width=2.4, style="->"):
    ax.annotate("", xy=end, xytext=start,
                arrowprops={"arrowstyle": style, "lw": width, "color": color,
                            "shrinkA": 0, "shrinkB": 0}, zorder=4)


def plot():
    x, x1, x2, mu = calculate()
    fig, axes = plt.subplots(1, 2, figsize=(11.8, 5.5))
    for ax in axes:
        ax.set_aspect("equal")
        ax.set_xlim(-.35, 3.15)
        ax.set_ylim(-.55, 2.18)
        ax.axis("off")
        ax.plot([-.25, 3.05], [0, 0], color=TEAL, lw=1.4, zorder=1)
        ax.scatter(0, 0, color=INK, s=28, zorder=5)
        ax.text(-.08, -.23, r"$0=x_0$", ha="center", color=INK, fontsize=13)
        ax.scatter(*x, color=INK, s=45, zorder=6)

    ax = axes[0]
    ax.add_patch(Circle(x, 1., fill=False, edgecolor=GRAY, lw=1.4))
    arrow(ax, (0, 0), x, INK, width=1.7)
    arrow(ax, (0, 0), x1, BLUE)
    ax.plot([2, 2], [0, 1], color=ORANGE, lw=1.9, linestyle=(0, (4, 3)))
    ax.plot([1.86, 1.86, 2], [0, .14, .14], color=ORANGE, lw=1.2)
    ax.scatter(*x1, color=BLUE, s=40, zorder=6)
    ax.text(1.03, -.27, r"$\mathcal{K}_1=\operatorname{span}\{p_1\}$", color=TEAL,
            ha="center", fontsize=13)
    ax.text(2.07, -.3, r"$x_1=\mu_1p_1$", color=BLUE, fontsize=13)
    ax.text(2.12, .45, r"$x-x_1$", color=ORANGE, fontsize=14)
    ax.text(2.1, 1.13, r"$x$ (unknown)", color=INK, fontsize=14)
    ax.text(.85, .55, r"$x$", color=INK, fontsize=15)
    ax.set_title("First space: the closest point on a line", color=INK, fontsize=16, pad=16)

    ax = axes[1]
    ax.fill_between([-.25, 3.05], -.4, 1.7, color=BLUE, alpha=.035, zorder=0)
    ax.plot([0, 0], [-.2, 1.65], color=PURPLE, lw=1.3, zorder=1)
    arrow(ax, (0, 0), x, GRAY, width=1.5)
    arrow(ax, (0, 0), x1, BLUE)
    arrow(ax, x1, x2, ORANGE)
    ax.plot([1.85, 1.85, 2], [0, .15, .15], color=GRAY, lw=1.2)
    ax.scatter(*x1, color=BLUE, s=40, zorder=6)
    ax.text(1., -.27, r"$\mu_1p_1$", color=BLUE, ha="center", fontsize=14)
    ax.text(2.15, .38, r"$\mu_2p_2$", color=ORANGE, fontsize=14)
    ax.text(2.13, 1.13, r"$x_2=x$", color=INK, fontsize=14)
    ax.text(-.12, 1.5, r"$p_2$", color=PURPLE, ha="right", fontsize=14)
    ax.text(2.88, -.25, r"$p_1$", color=TEAL, fontsize=14)
    ax.text(.27, 1.53, r"$\mathcal{K}_2=\operatorname{span}\{p_1,p_2\}$",
            color=INK, fontsize=13)
    ax.set_title("Second space: add the perpendicular part", color=INK, fontsize=16, pad=16)
    fig.text(.5, .97, r"Euclidean best approximations in nested Krylov spaces",
             ha="center", color=INK, fontsize=17)
    fig.text(.5, .09, r"$x_1=2p_1,\qquad x_2=2p_1+p_2=x$",
             ha="center", color=INK, fontsize=16)
    fig.text(.5, .025, r"The geometry is simple, but the coefficients $\mu_j=p_j^Tx$ require the unknown solution.",
             ha="center", color=INK, fontsize=13)
    fig.subplots_adjust(left=.025, right=.985, top=.83, bottom=.2, wspace=.12)
    root = Path(__file__).resolve().parents[1]
    fig.savefig(root/"_static/cg_euclidean_projection.svg", bbox_inches="tight",
                metadata={"Date": None})
    for suffix, dpi in (("", 180), ("-book-width", 65)):
        fig.savefig(f"/tmp/nla-cg-euclidean{suffix}.png", dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"Verified SPD example, orthonormal directions, projections, and coefficients {mu}.")


if __name__ == "__main__":
    plot()
