"""Illustrate Arnoldi convergence with computed Ritz-root filters and residuals."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


plt.rcParams.update({
    "font.size": 15,
    "svg.fonttype": "path",
    "svg.hashsalt": "arnoldi-convergence",
})
blue, purple, orange = "#2563eb", "#7c3aed", "#c26a13"
teal, ink, gray = "#0f766e", "#334155", "#cbd5e1"


def arnoldi(A, q1, steps):
    """Two-pass modified Gram-Schmidt, as in the Arnoldi implementation."""
    Q = np.zeros((len(q1), steps + 1))
    H = np.zeros((steps + 1, steps))
    Q[:, 0] = q1 / np.linalg.norm(q1)
    for j in range(steps):
        w = A @ Q[:, j]
        for _ in range(2):
            for i in range(j + 1):
                h = Q[:, i] @ w
                H[i, j] += h
                w -= h * Q[:, i]
        H[j + 1, j] = np.linalg.norm(w)
        assert H[j + 1, j] > 100 * np.finfo(float).eps * np.linalg.norm(A)
        Q[:, j + 1] = w / H[j + 1, j]
    return Q, H


def example_data():
    """Construct the example and check the identities used in the figure."""
    n, steps = 40, 10
    eigenvalues = np.r_[10., np.linspace(2., 1., n - 1)]
    X = np.eye(n) + .2 * np.diag(np.ones(n - 1), 1)
    X /= np.linalg.norm(X, axis=0)
    # Solve for A rather than forming X^{-1} explicitly.
    A = np.linalg.solve(X.T, (X * eigenvalues).T).T
    q1 = X @ np.ones(n)
    q1 /= np.linalg.norm(q1)
    c = np.linalg.solve(X, q1)
    norm_A = np.linalg.norm(A, 2)
    eps = np.finfo(float).eps

    assert np.linalg.norm(A - A.T) > 1e-2
    assert np.linalg.cond(X) < 1.6
    np.testing.assert_allclose(A @ X, X * eigenvalues, atol=1e-13)
    np.testing.assert_allclose(c, np.full(n, c[0]), atol=1e-14)

    Q, H = arnoldi(A, q1, steps)
    np.testing.assert_allclose(Q.T @ Q, np.eye(steps + 1), atol=1e-14)
    np.testing.assert_allclose(A @ Q[:, :steps], Q @ H, atol=1e-13)

    records = []
    for k in range(1, steps + 1):
        Qk, Hk = Q[:, :k], H[:k, :k]
        theta, Y = np.linalg.eig(Hk)
        wanted = np.argmin(np.abs(theta - eigenvalues[0]))
        y = Y[:, wanted] / np.linalg.norm(Y[:, wanted])
        u = Qk @ y
        np.testing.assert_allclose(np.linalg.norm(u), 1., atol=1e-14)

        residual = np.linalg.norm(A @ u - theta[wanted] * u)
        predicted = abs(H[k, k - 1] * y[-1])
        np.testing.assert_allclose(
            residual, predicted, rtol=1e-9, atol=50 * eps * norm_A)

        # p(z) = product_{j != wanted} (z - theta_j) / (10 - theta_j).
        # The true eigenvalue is used only to normalize this analysis filter.
        roots = np.delete(theta, wanted)
        p_values = np.prod(
            (eigenvalues[:, None] - roots) / (eigenvalues[0] - roots),
            axis=1,
        )
        np.testing.assert_allclose(p_values[0], 1., atol=1e-14)
        filtered = X @ (c * p_values)

        # Independently evaluate the factored polynomial on the matrix.
        direct = q1.astype(complex)
        for root in roots:
            direct = (A @ direct - root * direct) / (eigenvalues[0] - root)
        np.testing.assert_allclose(direct, filtered, atol=1e-13)

        # The normalized filter gives the same direction as the Ritz vector.
        parallel_error = np.linalg.norm(filtered - u * np.vdot(u, filtered))
        assert parallel_error / np.linalg.norm(filtered) < 1e-11
        assert residual > 0
        records.append({
            "k": k,
            "theta": theta[wanted],
            "relative_residual": residual / norm_A,
            "filter_values": np.abs(p_values[1:]),
        })

    # Keep the initial increase: the example does not imply monotone residuals.
    assert records[1]["relative_residual"] > records[0]["relative_residual"]
    return eigenvalues, records


def plot(eigenvalues, records):
    snapshots = {2: (blue, "o"), 4: (purple, "s"), 6: (orange, "^")}
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
        k = record["k"]
        if k not in snapshots:
            continue
        color, marker = snapshots[k]
        values = record["filter_values"]
        assert np.all(values > 0)
        left.scatter(eigenvalues[1:], values, color=color, marker=marker,
                     s=25, alpha=.85, edgecolors="none", zorder=3)

    left.axhline(1., color=teal, linewidth=1.6, linestyle=(0, (4, 3)))
    left.text(1.5, 1.6, r"Wanted component: $p(10)=1$", ha="center",
              color=teal, fontsize=13)
    smallest = min(record["filter_values"].min()
                   for record in records if record["k"] in snapshots)
    left.set_ylim(10. ** np.floor(np.log10(smallest)), 5.)
    left.set_xlim(.97, 2.03)
    left.set_xticks([1., 1.25, 1.5, 1.75, 2.])
    left.set_xlabel(r"Unwanted eigenvalues $\lambda_i$", fontsize=15, color=ink, labelpad=9)
    left.set_ylabel(r"Filter magnitude $|p(\lambda_i)|$", fontsize=15, color=ink, labelpad=9)

    right.set_title("Convergence of the Ritz pair", fontsize=17, color=ink, pad=15)
    steps = [record["k"] for record in records]
    residuals = [record["relative_residual"] for record in records]
    right.plot(steps, residuals, color=ink, linewidth=1.8,
               marker="o", markersize=4.5, zorder=2)
    for record in records:
        if record["k"] not in snapshots:
            continue
        color, marker = snapshots[record["k"]]
        right.scatter(record["k"], record["relative_residual"],
                      color=color, marker=marker, s=72, zorder=4)
    right.set_xlim(.7, 10.3)
    right.set_ylim(1e-14, 1.)
    right.set_xticks(steps)
    right.set_xlabel(r"Krylov subspace size $k$", fontsize=15, color=ink, labelpad=9)
    right.set_ylabel(r"Relative Ritz residual $\|Au-\theta u\|_2/\|A\|_2$",
                     fontsize=15, color=ink, labelpad=9)

    handles = [Line2D([], [], color=color, marker=marker, linestyle="none",
                      markersize=7, label=rf"$k={k}$")
               for k, (color, marker) in snapshots.items()]
    fig.legend(handles=handles, loc="upper center", ncol=3, frameon=False,
               bbox_to_anchor=(.53, 1.01), fontsize=14, handletextpad=.4,
               columnspacing=2.2)
    fig.subplots_adjust(left=.09, right=.98, bottom=.16, top=.83, wspace=.42)

    destination = Path(__file__).resolve().parents[1] / "_static"
    fig.savefig(destination / "arnoldi_convergence.svg", bbox_inches="tight",
                metadata={"Date": None})
    fig.savefig("/tmp/nla-arnoldi-convergence.png", dpi=180, bbox_inches="tight")
    # Check the figure at approximately the book's normal content width.
    fig.savefig("/tmp/nla-arnoldi-convergence-book-width.png", dpi=65,
                bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    eigenvalues, records = example_data()
    plot(eigenvalues, records)
    print("Verified filter normalization, Ritz directions, and residual identities.")
    for record in records:
        if record["k"] in (2, 4, 6):
            print(f"k={record['k']}: largest unwanted filter value "
                  f"{record['filter_values'].max():.2e}; "
                  f"relative Ritz residual {record['relative_residual']:.2e}")
