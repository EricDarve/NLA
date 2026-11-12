"""Compare GMRES for matrices with identical spectra and different eigenbases.

Run with Python, NumPy, and Matplotlib. The upper-triangular matrices are
constructed directly, without a numerically unstable similarity transform.
"""

from fractions import Fraction
from math import factorial
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


plt.rcParams.update({
    "font.size": 15,
    "svg.fonttype": "path",
    "svg.hashsalt": "gmres-nonnormal-convergence",
})
blue, purple, orange = "#2563eb", "#7c3aed", "#c26a13"
teal, ink, gray = "#0f766e", "#334155", "#cbd5e1"


def eigenbasis(A, superdiagonal):
    """Return unit-column eigenvectors and their 2-norm condition number.

    The diagonal spacing h=1/64 is exact in binary arithmetic. With S the
    unit superdiagonal shift, unnormalized eigenvectors form V=exp(a*S),
    a=superdiagonal/h, and V^{-1}=exp(-a*S). The exponential sums terminate.
    Using the explicit inverse avoids estimating a tiny singular value or
    numerically inverting an extremely ill-conditioned eigenvector matrix.
    """
    n = len(A)
    a = 64 * superdiagonal
    assert a == int(a)  # Exact rational check below uses the chosen integer a.
    coefficients = np.ones(n)
    for j in range(1, n):
        coefficients[j] = coefficients[j - 1] * a / j
    V, Vinv = np.zeros_like(A), np.zeros_like(A)
    for i in range(n):
        for j in range(i, n):
            V[i, j] = coefficients[j - i]
            Vinv[i, j] = (-1.) ** (j - i) * coefficients[j - i]

    # Verify the inverse formula exactly, rather than testing a product
    # whose cancellation would lose many digits in floating-point arithmetic.
    exact = [Fraction(int(a) ** j, factorial(j)) for j in range(n)]
    for degree in range(n):
        convolution = sum(exact[j] * (-1) ** (degree - j) * exact[degree - j]
                          for j in range(degree + 1))
        assert convolution == (1 if degree == 0 else 0)

    lengths = np.linalg.norm(V, axis=0)
    X = V / lengths
    Xinv = lengths[:, None] * Vinv
    np.testing.assert_allclose(np.linalg.norm(X, axis=0), 1., atol=1e-15)
    np.testing.assert_allclose(A @ X, X * np.diag(A), atol=2e-14)
    condition = np.linalg.norm(X, 2) * np.linalg.norm(Xinv, 2)
    return X, condition


def gmres_residuals(A, b):
    """Unrestarted GMRES, with two orthogonalization passes and QR solves.

    Each recorded residual is computed explicitly. Full reorthogonalization
    keeps the displayed comparison separate from loss of orthogonality.
    """
    n = len(b)
    beta = np.linalg.norm(b)
    Q = np.zeros((n, n + 1))
    Q[:, 0] = b / beta
    H = np.zeros((n + 1, n))
    residuals = [1.]
    condition_A = np.linalg.cond(A)
    for j in range(n):
        w = A @ Q[:, j]
        for _ in range(2):
            for i in range(j + 1):
                coefficient = Q[:, i] @ w
                H[i, j] += coefficient
                w -= coefficient * Q[:, i]
        H[j + 1, j] = np.linalg.norm(w)
        k = j + 1
        Qk, Hk = Q[:, :k], H[:k + 1, :k]
        rhs = np.zeros(k + 1)
        rhs[0] = beta
        Z, R = np.linalg.qr(Hk, mode="complete")
        g = Z.T @ rhs
        y = np.linalg.solve(R[:k, :], g[:k])
        residual = np.linalg.norm(b - A @ (Qk @ y)) / beta

        np.testing.assert_allclose(residual, abs(g[-1]) / beta, atol=3e-14)
        np.testing.assert_allclose(Qk.T @ Qk, np.eye(k), atol=2e-14)
        np.testing.assert_allclose(
            A @ Qk, Qk @ H[:k, :k] + np.outer(w, np.eye(k)[-1]),
            atol=2e-14,
        )
        # The small least-squares problems remain well conditioned even
        # though the eigenvector matrices are extremely ill conditioned.
        assert np.linalg.cond(Hk) <= condition_A * (1 + 1e-11)
        assert residual <= residuals[-1] + 3e-14
        residuals.append(residual)
        if residual < 1e-12:
            break
        assert H[j + 1, j] > 1e-14
        Q[:, j + 1] = w / H[j + 1, j]
    return np.asarray(residuals)


def example_data():
    n = 33
    eigenvalues = .75 + np.arange(n) / 64
    b = np.random.default_rng(5).normal(size=n)
    b /= np.linalg.norm(b)
    records = []
    for superdiagonal, color in ((0., blue), (.25, purple), (1., orange)):
        A = np.diag(eigenvalues) + superdiagonal * np.diag(np.ones(n - 1), 1)
        np.testing.assert_array_equal(np.linalg.eigvals(A), eigenvalues)
        _, condition = eigenbasis(A, superdiagonal)
        residuals = gmres_residuals(A, b)
        target_step = int(np.flatnonzero(residuals <= 1e-8)[0])
        if superdiagonal == 0:
            # The simple eigenvalue polynomial (1-z)^k supplies this bound.
            assert np.all(residuals <= .25 ** np.arange(len(residuals)) + 3e-14)
        records.append({
            "superdiagonal": superdiagonal, "color": color,
            "condition": condition, "residuals": residuals,
            "condition_A": np.linalg.cond(A),
            "target_step": target_step,
        })
    assert [r["target_step"] for r in records] == [9, 16, 31]
    return eigenvalues, records


def plot(eigenvalues, records):
    fig, (left, right) = plt.subplots(1, 2, figsize=(11.8, 5.1))
    for ax in (left, right):
        ax.spines[["top", "right"]].set_visible(False)
        ax.spines[["left", "bottom"]].set_color(gray)
        ax.tick_params(colors=ink, labelsize=12)
        ax.set_axisbelow(True)

    left.set_title("Identical eigenvalue clusters", fontsize=17, color=ink, pad=15)
    left.axvspan(.75, 1.25, color=gray, alpha=.18)
    left.axvline(0, color=gray, linewidth=1.2, linestyle=(0, (4, 3)))
    for row, record in zip((2, 1, 0), records):
        left.scatter(eigenvalues, np.full(len(eigenvalues), row),
                     color=record["color"], s=16, zorder=3)
    left.set_xlim(-.05, 1.4)
    left.set_ylim(-.5, 2.5)
    left.set_xticks([0, .5, .75, 1., 1.25])
    left.set_yticks([0, 1, 2], ["1", "0.25", "0"])
    left.set_xlabel(r"Eigenvalue $\lambda$", color=ink, labelpad=9)
    left.set_ylabel("First superdiagonal entry", color=ink, labelpad=9)

    right.set_title("Different GMRES convergence", fontsize=17, color=ink, pad=15)
    right.set_yscale("log")
    right.grid(axis="y", which="major", color=gray, alpha=.6, linewidth=.8)
    right.tick_params(axis="y", which="minor", left=False)
    for record in records:
        values = record["residuals"]
        right.plot(np.arange(len(values)), values, color=record["color"],
                   linewidth=2., marker="o", markersize=3., markevery=3)
    right.axhline(1e-8, color=teal, linewidth=1.1, linestyle=(0, (4, 3)))
    right.text(.70, 1.6e-8, r"Target: $10^{-8}$",
               transform=right.get_yaxis_transform(), ha="center", va="bottom",
               color=teal, fontsize=12)
    right.set_xlim(0, 33)
    right.set_ylim(1e-12, 2.)
    right.set_xticks([0, 10, 20, 30])
    right.set_yticks([1., 1e-3, 1e-6, 1e-9, 1e-12])
    right.set_xlabel(r"Iteration count $k$", color=ink, labelpad=9)
    right.set_ylabel(r"Relative residual $\|r_k\|_2/\|r_0\|_2$",
                     color=ink, labelpad=9)

    labels = [r"$\kappa_2(X)=1$", r"$\kappa_2(X)\approx3.5\times10^{13}$",
              r"$\kappa_2(X)\approx2.0\times10^{32}$"]
    handles = [Line2D([], [], color=r["color"], linewidth=2., label=label)
               for r, label in zip(records, labels)]
    fig.legend(handles=handles, loc="upper center", ncol=3,
               bbox_to_anchor=(.53, 1.08), frameon=False, fontsize=13,
               title=r"Eigenvector basis conditioning: $\kappa_2(X)$, not $\kappa_2(A)$",
               title_fontsize=13,
               handlelength=1.8, columnspacing=1.8)
    fig.subplots_adjust(left=.08, right=.98, bottom=.17, top=.82, wspace=.42)

    destination = Path(__file__).resolve().parents[1] / "_static"
    fig.savefig(destination / "gmres_nonnormal_convergence.svg",
                bbox_inches="tight", metadata={"Date": None})
    fig.savefig("/tmp/nla-gmres-nonnormal.png", dpi=180, bbox_inches="tight")
    fig.savefig("/tmp/nla-gmres-nonnormal-book-width.png", dpi=65,
                bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    eigenvalues, records = example_data()
    plot(eigenvalues, records)
    print("Verified spectra, eigenvectors, explicit inverse formula, Arnoldi "
          "relations, and QR residual norms.")
    for record in records:
        print(f"superdiagonal={record['superdiagonal']:g}: "
              f"kappa_2(X)={record['condition']:.3e}; "
              f"kappa_2(A)={record['condition_A']:.3e}; "
              f"target reached in {record['target_step']} steps")
