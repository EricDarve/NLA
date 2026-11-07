"""Compute and compare actual 2D Jacobi and Gauss-Seidel iterates."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({"font.size": 15, "svg.fonttype": "path",
                     "svg.hashsalt": "jacobi-gauss-seidel"})
BLUE, ORANGE, PURPLE = "#2563eb", "#c26a13", "#7c3aed"
INK, GRAY, TEAL = "#334155", "#cbd5e1", "#0f766e"
A = np.array([[1., -.8], [-.8, 1.]])
START = np.array([2.4, 1.8])
SWEEPS = 4


def compute():
    jacobi, gs, midpoints = [START.copy()], [START.copy()], []
    for _ in range(SWEEPS):
        old = jacobi[-1]
        new = old - (A @ old)/np.diag(A)
        jacobi.append(new)
        current = gs[-1].copy()
        for i in range(2):
            current[i] -= (A[i] @ current)/A[i, i]
            assert abs((A @ current)[i]) < 1e-14
            if i == 0:
                midpoints.append(current.copy())
        gs.append(current.copy())
    jacobi, gs, midpoints = map(np.array, (jacobi, gs, midpoints))
    # Independent checks using the splitting matrices.
    d = np.diag(np.diag(A))
    lower, upper = -np.tril(A, -1), -np.triu(A, 1)
    rj = np.linalg.solve(d, lower + upper)
    rg = np.linalg.solve(d-lower, upper)
    for k in range(SWEEPS + 1):
        np.testing.assert_allclose(jacobi[k], np.linalg.matrix_power(rj, k) @ START)
        np.testing.assert_allclose(gs[k], np.linalg.matrix_power(rg, k) @ START)
    np.testing.assert_allclose(jacobi[1], [1.44, 1.92])
    np.testing.assert_allclose(gs[1], [1.44, 1.152])
    full_path = np.array([START] + [p for pair in zip(midpoints, gs[1:]) for p in pair])
    assert np.all(np.diff(np.einsum('ij,ij->i', full_path @ A, full_path)) < 0)
    return jacobi, gs, midpoints


def arrow(ax, start, end, color):
    ax.annotate("", xy=end, xytext=start, zorder=3,
                arrowprops={"arrowstyle": "->", "color": color,
                            "lw": 2., "shrinkA": 0, "shrinkB": 0,
                            "mutation_scale": 12})


def plot():
    jacobi, gs, midpoints = compute()
    fig, axes = plt.subplots(1, 2, figsize=(11.8, 5.9))
    w, v = np.linalg.eigh(A)
    t = np.linspace(0, 2*np.pi, 800)
    circle = np.array([np.cos(t), np.sin(t)])
    for ax in axes:
        for level in (.01, .03, .08, .2, .5, .5*START @ A @ START):
            ellipse = v @ (np.sqrt(2*level/w)[:, None]*circle)
            ax.plot(*ellipse, color=GRAY, lw=1, zorder=0)
        ax.set_aspect("equal")
        ax.set_xlim(-.25, 2.8)
        ax.set_ylim(-.25, 2.3)
        ax.spines[["top", "right"]].set_visible(False)
        for spine in ("left", "bottom"):
            ax.spines[spine].set_position("zero")
            ax.spines[spine].set_color(GRAY)
        ax.tick_params(colors=INK, labelsize=11)
        ax.set_xticks([1, 2])
        ax.set_yticks([1, 2])
        ax.text(2.65, -.18, r"$z_1$", color=INK)
        ax.text(-.2, 2.17, r"$z_2$", color=INK)
        ax.scatter(0, 0, s=60, color=TEAL, edgecolor="white", zorder=6)
        ax.text(.08, -.19, "Solution", color=TEAL, fontsize=12)

    for start, end in zip(jacobi[:-1], jacobi[1:]):
        arrow(axes[0], start, end, PURPLE)
    for start, middle, end in zip(gs[:-1], midpoints, gs[1:]):
        arrow(axes[1], start, middle, BLUE)
        arrow(axes[1], middle, end, ORANGE)
        axes[1].scatter(*middle, s=23, facecolor="white", edgecolor=BLUE, zorder=4)

    for ax, points, color, offsets in (
        (axes[0], jacobi, PURPLE, [(8, 8), (-4, 12), (10, -7), (-22, 4), (9, -14)]),
        (axes[1], gs, INK, [(8, 8), (10, -4), (10, -4), (10, -4), (8, -13)]),
    ):
        ax.scatter(points[:, 0], points[:, 1], s=36, color=color, zorder=5)
        for k, (point, offset) in enumerate(zip(points, offsets)):
            ax.annotate(str(k), xy=point, xytext=offset, textcoords="offset points",
                        fontsize=13, weight="bold", color=color,
                        bbox={"facecolor": "white", "edgecolor": "none", "pad": 1.})

    axes[0].set_title("Jacobi: both use old values", color=PURPLE, fontsize=17, pad=18)
    axes[1].set_title("Gauss-Seidel: reuse the new value", color=INK, fontsize=17, pad=18)
    axes[0].text(.5, -.13, "Each arrow is one complete sweep.",
                 transform=axes[0].transAxes, ha="center", color=PURPLE, fontsize=13)
    axes[1].text(.5, -.13, "Horizontal + vertical = one sweep.",
                 transform=axes[1].transAxes, ha="center", color=INK, fontsize=13)
    fig.text(.5, .96, r"Same energy, same start: $x^{(0)}=(2.4,\,1.8)^T$",
             ha="center", color=INK, fontsize=17)
    fig.text(.5, .035, "Numbers mark completed sweeps; both panels use the same contours and scale.",
             ha="center", color=INK, fontsize=13)
    fig.subplots_adjust(left=.06, right=.98, bottom=.18, top=.82, wspace=.17)
    root = Path(__file__).resolve().parents[1]
    fig.savefig(root / "_static/jacobi_gauss_seidel.svg", bbox_inches="tight",
                metadata={"Date": None})
    for suffix, dpi in (("", 170), ("-book-width", 65)):
        fig.savefig(f"/tmp/nla-jacobi-gauss-seidel{suffix}.png", dpi=dpi,
                    bbox_inches="tight")
    plt.close(fig)
    print("Completed sweeps: k, Jacobi, Gauss-Seidel")
    for k in range(SWEEPS + 1):
        print(k, jacobi[k], gs[k])


if __name__ == "__main__":
    plot()
