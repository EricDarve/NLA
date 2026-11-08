"""Compare actual Jacobi, Gauss-Seidel, and optimal SOR paths in 2D.

Uses the same SPD quadratic and starting point as plot_jacobi_gauss_seidel.py.
"""

from pathlib import Path

from plot_jacobi_gauss_seidel import (
    A, START, BLUE, ORANGE, PURPLE, INK, GRAY, TEAL, arrow, compute,
)
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams["svg.hashsalt"] = "sor-geometry"
SWEEPS = 3


def calculate():
    jacobi, gs, gs_midpoints = compute()
    rho_j = abs(A[0, 1])/np.sqrt(A[0, 0]*A[1, 1])
    omega = 2/(1+np.sqrt(1-rho_j**2))
    np.testing.assert_allclose(omega, 1.25)
    sor, midpoints = [START.copy()], []
    for _ in range(SWEEPS):
        current = sor[-1].copy()
        for i in range(2):
            before = current.copy()
            correction = -(A[i] @ current)/A[i, i]
            current[i] += omega*correction
            drop = .5*(before @ A @ before-current @ A @ current)
            np.testing.assert_allclose(drop, .5*omega*(2-omega)*A[i, i]*correction**2)
            assert drop > 0
            if i == 0:
                midpoints.append(current.copy())
        sor.append(current.copy())
    sor, midpoints = np.array(sor), np.array(midpoints)
    np.testing.assert_allclose(sor[1], [1.2, .75])
    d, lower, upper = np.diag(np.diag(A)), -np.tril(A, -1), -np.triu(A, 1)

    def iteration(w):
        return np.linalg.solve(d-w*lower, (1-w)*d+w*upper)

    g = iteration(omega)
    for k, point in enumerate(sor):
        np.testing.assert_allclose(point, np.linalg.matrix_power(g, k) @ START)
    # At the optimum both eigenvalues equal omega-1; the double root is
    # sensitive to floating-point perturbations.
    np.testing.assert_allclose(np.linalg.eigvals(g), [.25, .25], atol=1e-8)
    for w in np.linspace(.001, 1.999, 2001):
        assert max(abs(np.linalg.eigvals(iteration(w)))) >= .25-1e-8
    return jacobi[:SWEEPS+1], gs[:SWEEPS+1], gs_midpoints[:SWEEPS], sor, midpoints


def plot():
    jacobi, gs, gs_midpoints, sor, sor_midpoints = calculate()
    fig, axes = plt.subplots(1, 3, figsize=(11.8, 5.1))
    w, v = np.linalg.eigh(A)
    t = np.linspace(0, 2*np.pi, 800)
    circle = np.array([np.cos(t), np.sin(t)])
    for ax in axes:
        for level in (.003, .01, .03, .08, .2, .5, .5*START @ A @ START):
            ellipse = v @ (np.sqrt(2*level/w)[:, None]*circle)
            ax.plot(*ellipse, color=GRAY, lw=1, zorder=0)
        ax.set_aspect("equal")
        ax.set_xlim(-.28, 2.8)
        ax.set_ylim(-.28, 2.3)
        ax.spines[["top", "right"]].set_visible(False)
        for spine in ("left", "bottom"):
            ax.spines[spine].set_position("zero")
            ax.spines[spine].set_color(GRAY)
        ax.tick_params(colors=INK, labelsize=11)
        ax.set_xticks([1, 2])
        ax.set_yticks([1, 2])
        ax.text(2.65, -.23, r"$z_1$", color=INK, fontsize=13)
        ax.text(-.25, 2.16, r"$z_2$", color=INK, fontsize=13)
        ax.scatter(0, 0, s=40, color=TEAL, edgecolor="white", zorder=6)

    for start, end in zip(jacobi[:-1], jacobi[1:]):
        arrow(axes[0], start, end, PURPLE)
    for ax, points, midpoints in ((axes[1], gs, gs_midpoints),
                                  (axes[2], sor, sor_midpoints)):
        for start, middle, end in zip(points[:-1], midpoints, points[1:]):
            arrow(ax, start, middle, BLUE)
            arrow(ax, middle, end, ORANGE)
            ax.scatter(*middle, s=20, facecolor="white", edgecolor=BLUE, zorder=4)
    for ax, points, color, offsets in (
        (axes[0], jacobi, PURPLE, [(7, 7), (-3, 10), (8, -8), (-17, 5)]),
        (axes[1], gs, INK, [(7, 7), (8, -4), (8, -4), (8, -6)]),
        (axes[2], sor, INK, [(7, 7), (8, -4), (8, -4), (10, -4)]),
    ):
        ax.scatter(points[:, 0], points[:, 1], s=28, color=color, zorder=5)
        for k, (point, offset) in enumerate(zip(points, offsets)):
            ax.annotate(str(k), xy=point, xytext=offset, textcoords="offset points",
                        fontsize=13, weight="bold", color=color,
                        bbox={"facecolor": "white", "edgecolor": "none", "pad": 1.})

    for ax, title, subtitle in zip(axes,
            ("Jacobi", "Gauss-Seidel", r"SOR: $\omega=1.25$"),
            ("Use old values", "Minimize along each coordinate", "Move beyond each coordinate minimum")):
        ax.set_title(title, color=INK, fontsize=17, pad=15)
        ax.text(.5, -.18, subtitle, transform=ax.transAxes, ha="center",
                color=INK, fontsize=11)
    fig.text(.5, .96, r"Three computed sweeps from $x^{(0)}=(2.4,\,1.8)^T$",
             ha="center", color=INK, fontsize=17)
    fig.text(.5, .085, "Numbers mark completed sweeps; identical contours and axes in all panels.",
             ha="center", color=INK, fontsize=12)
    fig.text(.5, .035, "Blue: first-coordinate update. Orange: second-coordinate update. Teal: solution.",
             ha="center", color=INK, fontsize=12)
    fig.subplots_adjust(left=.04, right=.98, bottom=.25, top=.83, wspace=.16)
    root = Path(__file__).resolve().parents[1]
    fig.savefig(root / "_static/sor_geometry.svg", bbox_inches="tight", metadata={"Date": None})
    for suffix, dpi in (("", 180), ("-book-width", 65)):
        fig.savefig(f"/tmp/nla-sor-geometry{suffix}.png", dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print("omega_opt=1.25; spectral radii: Jacobi 0.8, GS 0.64, SOR 0.25")
    for k in range(SWEEPS+1):
        print(f"k={k}: Jacobi {jacobi[k]}, GS {gs[k]}, SOR {sor[k]}")


if __name__ == "__main__":
    plot()
