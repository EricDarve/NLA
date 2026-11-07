"""Solve two Poisson grids and plot measured error after every complete sweep.

All three methods start from zero. The exact solution u=x(1-x)y(1-y)
satisfies both the PDE and its five-point discretization exactly.
No spectral-radius formula is used to generate the error curves.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

plt.rcParams.update({"font.size": 15, "svg.fonttype": "path",
                     "svg.hashsalt": "poisson-measured-convergence"})
COLORS = ("#2563eb", "#7c3aed", "#c26a13")
METHODS = ("Jacobi", "Gauss-Seidel", "Optimal SOR")
INK, GRAY, TEAL = "#334155", "#cbd5e1", "#0f766e"
TARGET = 1e-6


def neighbors(u):
    return u[:-2, 1:-1] + u[2:, 1:-1] + u[1:-1, :-2] + u[1:-1, 2:]


def problem(n):
    h = 1/(n+1)
    x, y = np.meshgrid(np.arange(1, n+1)*h, np.arange(1, n+1)*h)
    exact = x*(1-x)*y*(1-y)
    rhs = h*h*2*(x*(1-x) + y*(1-y))
    padded = np.pad(exact, 1)
    np.testing.assert_allclose(4*exact-neighbors(padded), rhs, atol=2e-16)
    red = np.indices((n, n)).sum(axis=0) % 2 == 0
    return exact, rhs, red


def sweep(u, rhs, red, method, omega):
    interior = u[1:-1, 1:-1]
    if method == "Jacobi":
        # Evaluate all right-hand sides before overwriting any unknown.
        interior[:] = (rhs + neighbors(u))/4
    else:
        for mask in (red, ~red):
            # Recompute after the first color, so the second sees new values.
            provisional = (rhs + neighbors(u))/4
            interior[mask] += omega*(provisional[mask]-interior[mask])


def verify_sweeps():
    """Compare the grid implementation with independently formed splittings."""
    n = 4
    _, rhs, red = problem(n)
    t = 2*np.eye(n)-np.diag(np.ones(n-1), 1)-np.diag(np.ones(n-1), -1)
    a = np.kron(np.eye(n), t) + np.kron(t, np.eye(n))
    perm = np.r_[np.flatnonzero(red.ravel()), np.flatnonzero(~red.ravel())]
    a = a[np.ix_(perm, perm)]
    d, l, upper = np.diag(np.diag(a)), -np.tril(a, -1), -np.triu(a, 1)
    start = np.random.default_rng(42).normal(size=(n, n))
    for method, omega in zip(METHODS, (1., 1., 2/(1+np.sin(np.pi/(n+1))))):
        u = np.pad(start, 1)
        state = start.ravel()[perm].copy()
        for _ in range(3):
            if method == "Jacobi":
                state = np.linalg.solve(d, (l+upper)@state + rhs.ravel()[perm])
            else:
                state = np.linalg.solve(d-omega*l,
                        ((1-omega)*d+omega*upper)@state + omega*rhs.ravel()[perm])
            sweep(u, rhs, red, method, omega)
            np.testing.assert_allclose(u[1:-1, 1:-1].ravel()[perm], state, atol=1e-14)


def solve(n, method):
    exact, rhs, red = problem(n)
    omega = 2/(1+np.sin(np.pi/(n+1))) if method == "Optimal SOR" else 1.
    u = np.zeros((n+2, n+2))
    norm_exact = np.linalg.norm(exact.ravel())
    errors = [1.]
    for k in range(1, 25001):
        sweep(u, rhs, red, method, omega)
        error = np.linalg.norm((u[1:-1, 1:-1]-exact).ravel())/norm_exact
        assert np.isfinite(error)
        errors.append(error)
        if error <= TARGET:
            break
    else:
        raise RuntimeError(f"{method} did not converge for N={n}")
    assert errors[-2] > TARGET >= errors[-1]
    assert np.all(u[[0, -1], :] == 0) and np.all(u[:, [0, -1]] == 0)
    residual = rhs - (4*u[1:-1, 1:-1]-neighbors(u))
    assert np.linalg.norm(residual.ravel())/np.linalg.norm(rhs.ravel()) < 1e-4
    print(f"N={n}, {method}: {k} sweeps; relative error {error:.8e}")
    return np.array(errors)


def plot():
    verify_sweeps()
    fig, axes = plt.subplots(1, 2, figsize=(11.8, 5.3))
    all_counts = []
    for ax, n in zip(axes, (31, 63)):
        counts = []
        for method, color in zip(METHODS, COLORS):
            errors = solve(n, method)
            steps = np.arange(1, len(errors))
            ax.plot(steps, errors[1:], color=color, linewidth=2.2)
            ax.scatter(steps[-1], errors[-1], color=color, s=38, zorder=4)
            counts.append(int(steps[-1]))
        all_counts.append(counts)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(1, 16000)
        ax.set_ylim(4e-7, 2)
        ax.set_xticks([1, 10, 100, 1000, 10000], ["1", "10", "100", "1,000", "10,000"])
        ax.set_yticks([1, 1e-2, 1e-4, 1e-6])
        ax.grid(axis="y", which="major", color=GRAY, alpha=.6, linewidth=.8)
        ax.set_axisbelow(True)
        ax.spines[["top", "right"]].set_visible(False)
        ax.spines[["left", "bottom"]].set_color(GRAY)
        ax.tick_params(colors=INK, labelsize=12)
        ax.tick_params(which="minor", bottom=False, left=False)
        ax.axhline(TARGET, color=TEAL, lw=1.2, linestyle=(0, (4, 3)))
        ax.set_title(rf"Measured error: ${n}\times{n}$ grid", color=INK, fontsize=17, pad=15)
        ax.set_xlabel(r"Sweep count $k$ (log scale)", color=INK, labelpad=9)
    axes[0].set_ylabel("Relative solution error", color=INK, labelpad=9)
    handles = [Line2D([], [], color=c, lw=2.2, label=m) for c, m in zip(COLORS, METHODS)]
    fig.legend(handles=handles, loc="upper center", ncol=3,
               bbox_to_anchor=(.53, 1.01), frameon=False, fontsize=14)
    fig.text(.53, .02, r"Same zero initial guess; dashed line marks relative error $10^{-6}$.",
             ha="center", color=INK, fontsize=12)
    fig.subplots_adjust(left=.09, right=.98, bottom=.24, top=.83, wspace=.26)
    for ax, counts in zip(axes, all_counts):
        center = (ax.get_position().x0 + ax.get_position().x1)/2
        fig.text(center, .085,
                 f"Sweeps: J {counts[0]:,} | GS {counts[1]:,} | SOR {counts[2]:,}",
                 ha="center", color=INK, fontsize=11)
    root = Path(__file__).resolve().parents[1]
    fig.savefig(root / "_static/poisson_measured_convergence.svg", bbox_inches="tight",
                metadata={"Date": None})
    for suffix, dpi in (("", 180), ("-book-width", 65)):
        fig.savefig(f"/tmp/nla-poisson-measured{suffix}.png", dpi=dpi, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    plot()
