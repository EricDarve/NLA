"""Compare asymptotic decay factors for the five-point Poisson problem.

Run with Python, NumPy, and Matplotlib. These are spectral-radius predictions,
not measured errors or finite-iteration error bounds. GS and SOR use red-black
ordering. The small-matrix checks independently verify the spectral formulas.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


plt.rcParams.update({
    "font.size": 15,
    "svg.fonttype": "path",
    "svg.hashsalt": "poisson-convergence",
})
COLORS = ("#2563eb", "#7c3aed", "#c26a13")
LABELS = ("Jacobi", "Gauss-Seidel", "Optimal SOR")
INK, GRAY, TEAL = "#334155", "#cbd5e1", "#0f766e"
TARGET = 1e-6


def radii(n):
    angle = np.pi / (n + 1)
    return np.array([
        np.cos(angle),
        np.cos(angle)**2,
        (1 - np.sin(angle)) / (1 + np.sin(angle)),
    ])


def sweep_counts(n):
    rho = radii(n)
    counts = np.ceil(np.log(TARGET) / np.log(rho)).astype(int)
    assert np.all(rho**counts <= TARGET)
    assert np.all(rho**(counts - 1) > TARGET)
    return counts


def verify():
    for n in (3, 5, 8):
        t = (2*np.eye(n) - np.diag(np.ones(n-1), 1)
             - np.diag(np.ones(n-1), -1))
        a = np.kron(np.eye(n), t) + np.kron(t, np.eye(n))
        order = sorted(range(n*n), key=lambda j: (j//n + j % n) % 2)
        a = a[np.ix_(order, order)]
        diagonal = np.diag(np.diag(a))
        lower, upper = -np.tril(a, -1), -np.triu(a, 1)
        omega = 2 / (1 + np.sin(np.pi / (n + 1)))
        matrices = (
            np.linalg.solve(diagonal, lower + upper),
            np.linalg.solve(diagonal - lower, upper),
            np.linalg.solve(diagonal - omega*lower,
                            (1-omega)*diagonal + omega*upper),
        )
        observed = [max(abs(np.linalg.eigvals(r))) for r in matrices]
        # The repeated eigenvalues at optimal omega are sensitive to roundoff.
        np.testing.assert_allclose(observed, radii(n), atol=1e-7, rtol=0)
    np.testing.assert_array_equal(sweep_counts(63), [11463, 5732, 141])


def plot():
    fig, axes = plt.subplots(1, 2, figsize=(11.8, 5.1))
    for ax in axes:
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.grid(axis="y", which="major", color=GRAY, alpha=.6, linewidth=.8)
        ax.set_axisbelow(True)
        ax.spines[["top", "right"]].set_visible(False)
        ax.spines[["left", "bottom"]].set_color(GRAY)
        ax.tick_params(colors=INK, labelsize=12)
        ax.tick_params(which="minor", bottom=False, left=False)

    ax = axes[0]
    steps = np.geomspace(1, 16000, 700)
    counts = sweep_counts(63)
    for rho, count, color in zip(radii(63), counts, COLORS):
        # Stop just below the visible range to avoid underflow far off the plot.
        visible = steps <= np.log(1e-8) / np.log(rho)
        ax.plot(steps[visible], rho**steps[visible], color=color, linewidth=2.2)
        ax.scatter(count, rho**count, color=color, s=38, zorder=4)
    ax.axhline(TARGET, color=TEAL, linewidth=1.2, linestyle=(0, (4, 3)))
    ax.text(2, TARGET*1.6, r"Target: $10^{-6}$", color=TEAL, fontsize=12)
    ax.set_title(r"Decay model: $63\times63$ grid", color=INK, fontsize=17, pad=15)
    ax.set_ylabel(r"Asymptotic decay factor $\rho^k$", color=INK, labelpad=9)
    ax.set_xlabel(r"Sweep count $k$ (log scale)", color=INK, labelpad=9)
    ax.set_xlim(1, 16000)
    ax.set_ylim(1e-7, 2)
    ax.set_xticks([1, 10, 100, 1000, 10000], ["1", "10", "100", "1,000", "10,000"])
    ax.set_yticks([1, 1e-2, 1e-4, 1e-6])

    ax = axes[1]
    grids = np.array([7, 15, 31, 63, 127])
    counts = np.array([sweep_counts(int(n)) for n in grids])
    for j, color in enumerate(COLORS):
        ax.plot(grids, counts[:, j], color=color, linewidth=2.2,
                marker="o", markersize=5)
    ax.set_title(r"Sweeps to reach $\rho^k\leq10^{-6}$", color=INK, fontsize=17, pad=15)
    ax.set_ylabel("Predicted sweep count", color=INK, labelpad=9)
    ax.set_xlabel(r"Grid points per side $N$ (log scale)", color=INK, labelpad=9)
    ax.set_xlim(6, 150)
    ax.set_ylim(10, 1e5)
    ax.set_xticks(grids, [str(n) for n in grids])
    ax.set_yticks([10, 100, 1000, 10000, 100000])

    handles = [Line2D([], [], color=color, linewidth=2.2, label=label)
               for color, label in zip(COLORS, LABELS)]
    fig.legend(handles=handles, loc="upper center", ncol=3,
               bbox_to_anchor=(.53, 1.01), frameon=False, fontsize=14)
    fig.subplots_adjust(left=.09, right=.98, bottom=.17, top=.83, wspace=.36)
    root = Path(__file__).resolve().parents[1]
    fig.savefig(root / "_static/poisson_convergence.svg", bbox_inches="tight",
                metadata={"Date": None})
    for suffix, dpi in (("", 180), ("-book-width", 65)):
        fig.savefig(f"/tmp/nla-poisson-convergence{suffix}.png", dpi=dpi,
                    bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    verify()
    plot()
    for n in (7, 15, 31, 63, 127):
        print(f"N={n}: Jacobi, Gauss-Seidel, optimal SOR: {sweep_counts(n)}")
