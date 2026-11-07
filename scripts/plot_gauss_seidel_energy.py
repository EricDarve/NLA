"""Plot exact coordinate minimization for a two-dimensional SPD quadratic."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


plt.rcParams.update({
    "font.size": 15,
    "svg.fonttype": "path",
    "svg.hashsalt": "gauss-seidel-energy",
})
BLUE, ORANGE = "#2563eb", "#c26a13"
INK, GRAY, TEAL = "#334155", "#94a3b8", "#0f766e"
A = np.array([[1., -.8], [-.8, 1.]])


def trajectory(sweeps=6):
    points = [np.array([2.4, 1.8])]
    for _ in range(sweeps):
        for i in (0, 1):
            previous = points[-1]
            current = previous.copy()
            current[i] = .8 * current[1-i]
            # The derivative along the updated coordinate vanishes.
            assert abs((A @ current)[i]) < 1e-14
            drop = .5*(previous @ A @ previous - current @ A @ current)
            np.testing.assert_allclose(drop, .5*A[i, i]*(current[i]-previous[i])**2)
            assert drop > 0
            points.append(current)
    np.testing.assert_allclose(points[2], [1.44, 1.152])
    return np.array(points)


def plot():
    points = trajectory()
    eigenvalues, eigenvectors = np.linalg.eigh(A)
    assert np.all(eigenvalues > 0)
    angle = np.linspace(0, 2*np.pi, 700)
    circle = np.array([np.cos(angle), np.sin(angle)])

    fig, ax = plt.subplots(figsize=(8.8, 7.8))
    for j, point in enumerate(points):
        level = .5*point @ A @ point
        ellipse = eigenvectors @ (np.sqrt(2*level/eigenvalues)[:, None]*circle)
        np.testing.assert_allclose(.5*np.sum(ellipse*(A @ ellipse), axis=0), level)
        color = GRAY if j == 0 else (BLUE if j % 2 else ORANGE)
        ax.plot(*ellipse, color=color, linewidth=1., alpha=.55, zorder=1)

    for j, (start, end) in enumerate(zip(points[:-1], points[1:])):
        color = BLUE if j % 2 == 0 else ORANGE
        ax.annotate("", xy=end, xytext=start, zorder=4,
                    arrowprops={"arrowstyle": "->", "color": color,
                                "linewidth": 2., "shrinkA": 0, "shrinkB": 0,
                                "mutation_scale": 12 if j < 6 else 8})
        ax.scatter(*end, s=25 if j < 6 else 12, facecolor="white",
                   edgecolor=color, linewidth=1.1, zorder=5)
    ax.scatter(*points[0], color=INK, s=46, zorder=6)
    ax.scatter(0, 0, color=TEAL, edgecolor="white", s=65, zorder=7)

    ax.plot([1.05, 1.85], [1.8, 1.8], color=BLUE, linestyle=(0, (3, 3)), lw=1)
    ax.plot([1.44, 1.44], [.82, 1.48], color=ORANGE, linestyle=(0, (3, 3)), lw=1)
    ax.annotate("Horizontal minimum\n" + r"$\partial\Phi/\partial z_1=0$",
                xy=points[1], xytext=(.35, 2.55), color=BLUE, fontsize=14,
                arrowprops={"arrowstyle": "-", "color": BLUE, "lw": .9})
    ax.annotate("Vertical minimum\n" + r"$\partial\Phi/\partial z_2=0$",
                xy=points[2], xytext=(2., 1.12), color=ORANGE, fontsize=14,
                arrowprops={"arrowstyle": "-", "color": ORANGE, "lw": .9})
    ax.text(2.58, 1.94, r"$x^{(0)}=(2.4,\,1.8)$", color=INK, fontsize=14)
    ax.annotate(r"$x^{(1)}=(1.44,\,1.152)$" + "\nOne complete sweep",
                xy=points[2], xytext=(2., .42), color=INK, fontsize=13,
                arrowprops={"arrowstyle": "-", "color": GRAY, "lw": .9})
    ax.annotate("Solution " + r"$x=0$", xy=(0, 0), xytext=(-1.3, -.65),
                color=TEAL, fontsize=14,
                bbox={"facecolor": "white", "edgecolor": "none", "pad": 1.5},
                arrowprops={"arrowstyle": "-", "color": TEAL, "lw": .9})

    ax.set_aspect("equal")
    ax.set_xlim(-2.9, 4.1)
    ax.set_ylim(-2.7, 3.)
    ax.spines[["top", "right"]].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_position("zero")
        ax.spines[spine].set_color(GRAY)
    ax.set_xticks([-2, -1, 1, 2, 3])
    ax.set_yticks([-2, -1, 1, 2])
    ax.tick_params(colors=GRAY, labelsize=12)
    ax.text(3.95, -.22, r"$z_1$", fontsize=16, color=INK)
    ax.text(-.35, 2.84, r"$z_2$", fontsize=16, color=INK)

    fig.text(.09, .94, "Gauss-Seidel: minimize along one coordinate",
             color=INK, fontsize=19, weight="bold")
    fig.text(.09, .895, r"$\Phi(z)=\frac{1}{2}(z_1^2-1.6z_1z_2+z_2^2)$",
             color=INK, fontsize=17)
    handles = [Line2D([], [], color=BLUE, lw=2, label=r"Update $z_1\leftarrow0.8z_2$"),
               Line2D([], [], color=ORANGE, lw=2, label=r"Update $z_2\leftarrow0.8z_1$")]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(.5, .055),
               ncol=2, frameon=False, fontsize=14)
    fig.text(.5, .025, "Each contour shows constant energy; arrows show coordinate updates.",
             ha="center", fontsize=12, color=INK)
    fig.subplots_adjust(left=.08, right=.98, bottom=.14, top=.86)
    root = Path(__file__).resolve().parents[1]
    fig.savefig(root / "_static/gauss_seidel_energy.svg", bbox_inches="tight",
                metadata={"Date": None})
    for suffix, dpi in (("", 160), ("-book-width", 85)):
        fig.savefig(f"/tmp/nla-gauss-seidel-energy{suffix}.png", dpi=dpi,
                    bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    plot()
    print("Verified coordinate minima, energy decrease, contour levels, and first sweep.")
