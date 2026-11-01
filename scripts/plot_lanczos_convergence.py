"""Compare Chebyshev trial filters, the Lanczos error, and its upper bound.

Run from any directory with Python, NumPy, and Matplotlib installed.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from numpy.polynomial.chebyshev import chebval


plt.rcParams.update({
    "font.size": 15,
    "svg.fonttype": "path",
    "svg.hashsalt": "lanczos-convergence",
})
blue, purple, orange = "#2563eb", "#7c3aed", "#c26a13"
teal, ink, gray = "#0f766e", "#334155", "#cbd5e1"


def lanczos(A, q1, steps):
    """Three-term recurrence with full reorthogonalization for this example."""
    n = len(q1)
    Q = np.zeros((n, steps + 1))
    Q[:, 0] = q1 / np.linalg.norm(q1)
    alpha, beta = [], []

    for j in range(steps):
        w = A @ Q[:, j]
        if j > 0:
            w -= beta[j - 1] * Q[:, j - 1]
        a = Q[:, j] @ w
        w -= a * Q[:, j]
        # Remove small components reintroduced by rounding.
        for _ in range(2):
            for i in range(j + 1):
                w -= (Q[:, i] @ w) * Q[:, i]
        b = np.linalg.norm(w)
        assert b > 100 * np.finfo(float).eps * np.linalg.norm(A)
        alpha.append(a)
        beta.append(b)
        Q[:, j + 1] = w / b

    T = np.diag(alpha) + np.diag(beta[:-1], 1) + np.diag(beta[:-1], -1)
    return Q, T, np.asarray(beta)


def chebyshev_filter(values, k, eigenvalues):
    """The normalized degree-(k-1) trial polynomial from the theorem."""
    largest, second, smallest = eigenvalues[0], eigenvalues[1], eigenvalues[-1]
    mapped = (2 * np.asarray(values) - second - smallest) / (second - smallest)
    rho = (largest - second) / (second - smallest)
    coefficients = [0.] * (k - 1) + [1.]
    denominator = chebval(1 + 2 * rho, coefficients)
    return chebval(mapped, coefficients) / denominator, denominator


def example_data():
    """Check the projection, trial filters, and eigenvalue error bounds."""
    n, steps = 40, 12
    eigenvalues = np.r_[2.5, np.linspace(2., 1., n - 1)]
    A = np.diag(eigenvalues)
    q1 = np.ones(n) / np.sqrt(n)
    spread = eigenvalues[0] - eigenvalues[-1]
    tan_squared = np.sum(q1[1:] ** 2) / q1[0] ** 2
    norm_A = np.linalg.norm(A, 2)
    eps = np.finfo(float).eps
    Q, T, beta = lanczos(A, q1, steps)

    np.testing.assert_allclose(Q.T @ Q, np.eye(steps + 1), atol=1e-14)
    extended_T = np.vstack((T, np.zeros(steps)))
    extended_T[-1, -1] = beta[-1]
    np.testing.assert_allclose(A @ Q[:, :steps], Q @ extended_T, atol=1e-13)

    records = []
    old_theta = -np.inf
    for k in range(1, steps + 1):
        Qk, Tk = Q[:, :k], T[:k, :k]
        np.testing.assert_allclose(Tk, Qk.T @ A @ Qk, atol=1e-13)
        theta, Y = np.linalg.eigh(Tk)
        largest_ritz = theta[-1]
        error = eigenvalues[0] - largest_ritz
        assert largest_ritz >= old_theta - 20 * eps * norm_A
        old_theta = largest_ritz
        # Stop the plotted range before roundoff obscures the positive error.
        assert error > 100 * eps * norm_A

        p_values, denominator = chebyshev_filter(eigenvalues, k, eigenvalues)
        np.testing.assert_allclose(p_values[0], 1., atol=1e-14)
        assert np.max(np.abs(p_values[1:])) <= 1 / denominator + 20 * eps
        bound = spread * tan_squared / denominator ** 2
        np.testing.assert_allclose(
            bound, spread * tan_squared * np.max(np.abs(p_values[1:])) ** 2,
            rtol=1e-13,
        )

        # The trial vector belongs to the Krylov space and supplies the bound.
        trial = q1 * p_values
        distance = np.linalg.norm(trial - Qk @ (Qk.T @ trial))
        assert distance < 1e-12
        trial_error = eigenvalues[0] - (trial @ A @ trial) / (trial @ trial)
        assert error <= trial_error + 50 * eps * norm_A
        assert trial_error <= bound + 50 * eps * norm_A

        u = Qk @ Y[:, -1]
        residual = np.linalg.norm(A @ u - largest_ritz * u)
        np.testing.assert_allclose(
            residual, beta[k - 1] * abs(Y[-1, -1]),
            rtol=1e-9, atol=50 * eps * norm_A,
        )
        records.append({
            "k": k,
            "ritz_value": largest_ritz,
            "error": error,
            "bound": bound,
            "filter_values": np.abs(p_values[1:]),
        })

    return eigenvalues, records


def plot(eigenvalues, records):
    snapshots = {3: (blue, "o"), 7: (purple, "s"), 11: (orange, "^")}
    fig, (left, right) = plt.subplots(1, 2, figsize=(11.8, 5.1))

    for ax in (left, right):
        ax.set_yscale("log")
        ax.grid(axis="y", which="major", color=gray, alpha=.6, linewidth=.8)
        ax.set_axisbelow(True)
        ax.spines[["top", "right"]].set_visible(False)
        ax.spines[["left", "bottom"]].set_color(gray)
        ax.tick_params(colors=ink, labelsize=12)
        ax.tick_params(axis="y", which="minor", left=False)

    left.set_title("Suppressing unwanted components", fontsize=17, color=ink, pad=15)
    for record in records:
        if record["k"] not in snapshots:
            continue
        color, marker = snapshots[record["k"]]
        values = record["filter_values"]
        assert np.all(values > 0)
        left.scatter(eigenvalues[1:], values, color=color, marker=marker,
                     s=25, alpha=.85, edgecolors="none", zorder=3)
    left.axhline(1., color=teal, linewidth=1.6, linestyle=(0, (4, 3)))
    left.text(1.5, 1.6, r"Wanted component: $p(2.5)=1$", ha="center",
              color=teal, fontsize=13)
    smallest = min(record["filter_values"].min()
                   for record in records if record["k"] in snapshots)
    left.set_ylim(10. ** np.floor(np.log10(smallest)), 5.)
    left.set_xlim(.97, 2.03)
    left.set_xticks([1., 1.25, 1.5, 1.75, 2.])
    left.set_xlabel(r"Unwanted eigenvalues $\lambda_i$", fontsize=15, color=ink, labelpad=9)
    left.set_ylabel(r"Filter magnitude $|p(\lambda_i)|$", fontsize=15, color=ink, labelpad=9)

    right.set_title("Convergence of the largest Ritz value", fontsize=17, color=ink, pad=15)
    steps = [record["k"] for record in records]
    right.plot(steps, [record["error"] for record in records],
               color=ink, linewidth=1.8, marker="o", markersize=4.5,
               label="Ritz-value error", zorder=2)
    right.plot(steps, [record["bound"] for record in records],
               color=teal, linewidth=1.6, linestyle=(0, (4, 3)),
               label="Chebyshev bound", zorder=1)
    for record in records:
        if record["k"] not in snapshots:
            continue
        color, marker = snapshots[record["k"]]
        right.scatter(record["k"], record["error"], color=color,
                      marker=marker, s=72, zorder=4)
    right.set_xlim(.7, 12.3)
    right.set_ylim(1e-12, 2e2)
    right.set_xticks([1, 3, 5, 7, 9, 11])
    right.set_xlabel(r"Krylov subspace size $k$", fontsize=15, color=ink, labelpad=9)
    right.set_ylabel(r"Eigenvalue error $\lambda_1-\theta_1^{(k)}$",
                     fontsize=15, color=ink, labelpad=9)
    right.legend(loc="upper right", fontsize=12, frameon=False,
                 handlelength=1.8, handletextpad=.5)

    handles = [Line2D([], [], color=color, marker=marker, linestyle="none",
                      markersize=7, label=rf"$k={k}$")
               for k, (color, marker) in snapshots.items()]
    fig.legend(handles=handles, loc="upper center", ncol=3, frameon=False,
               bbox_to_anchor=(.53, 1.01), fontsize=14, handletextpad=.4,
               columnspacing=2.2)
    fig.subplots_adjust(left=.09, right=.98, bottom=.16, top=.83, wspace=.42)

    destination = Path(__file__).resolve().parents[1] / "_static"
    fig.savefig(destination / "lanczos_convergence.svg", bbox_inches="tight",
                metadata={"Date": None})
    fig.savefig("/tmp/nla-lanczos-convergence.png", dpi=180, bbox_inches="tight")
    fig.savefig("/tmp/nla-lanczos-convergence-book-width.png", dpi=65,
                bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    eigenvalues, records = example_data()
    plot(eigenvalues, records)
    print("Verified Lanczos projections, Chebyshev filters, and error bounds.")
    for record in records:
        if record["k"] in (3, 7, 11):
            print(f"k={record['k']}: Ritz-value error {record['error']:.2e}; "
                  f"Chebyshev bound {record['bound']:.2e}")
