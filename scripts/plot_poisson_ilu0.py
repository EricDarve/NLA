"""Compare the complete Poisson and strict ILU(0)-preconditioned spectra.

Run with Python, NumPy, and Matplotlib. No SciPy or generic ILU package is
needed: the natural-order five-point pattern is factored without any fill,
reordering, pivoting, or diagonal modification.
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({
    "font.size": 15,
    "svg.fonttype": "path",
    "svg.hashsalt": "poisson-ilu0",
})
BLUE, PURPLE, TEAL = "#2563eb", "#7c3aed", "#0f766e"
INK, GRAY = "#334155", "#cbd5e1"


def stencil(v, m):
    """Apply h^2 A on one vector or a matrix of column vectors."""
    grid = np.asarray(v).reshape((m, m) + v.shape[1:])
    out = 4 * grid.copy()
    out[1:] -= grid[:-1]
    out[:-1] -= grid[1:]
    out[:, 1:] -= grid[:, :-1]
    out[:, :-1] -= grid[:, 1:]
    return out.reshape(v.shape)


class PoissonILU0:
    """Natural-order strict ILU(0) for h^2 A, stored as M0 = C C^T.

    If M0 = L U with L unit lower triangular, its upper factor is
    U = D L^T. Only the west/south entries of L and east/north entries of
    U are retained. Their off-diagonal U entries remain -1 because the
    five-point adjacency graph contains no triangles. The diagonal obeys
    d_i = 4 - sum(1/d_j) over the earlier neighbors j.
    """
    def __init__(self, m):
        self.m = m
        self.n = m * m
        d = np.empty(self.n)
        for i in range(self.n):
            d[i] = 4.
            if i % m:
                d[i] -= 1 / d[i - 1]
            if i >= m:
                d[i] -= 1 / d[i - m]
        assert np.all(d > 0)
        self.d = d
        self.root = np.sqrt(d)

    def forward(self, v):
        """Solve C w = v, including for multiple right-hand sides."""
        w = np.array(v, dtype=float, copy=True)
        for i in range(self.n):
            if i % self.m:
                w[i] += w[i - 1] / self.root[i - 1]
            if i >= self.m:
                w[i] += w[i - self.m] / self.root[i - self.m]
            w[i] /= self.root[i]
        return w

    def backward(self, v):
        """Solve C^T w = v."""
        w = np.array(v, dtype=float, copy=True)
        for i in range(self.n - 1, -1, -1):
            if i % self.m < self.m - 1:
                w[i] += w[i + 1] / self.root[i]
            if i < self.n - self.m:
                w[i] += w[i + self.m] / self.root[i]
            w[i] /= self.root[i]
        return w

    def solve(self, v):
        return self.backward(self.forward(v))

    def matrices(self):
        """Explicit L, U for checks, not for iterative applications."""
        L = np.eye(self.n)
        for i in range(self.n):
            if i % self.m:
                L[i, i - 1] = -1 / self.d[i - 1]
            if i >= self.m:
                L[i, i - self.m] = -1 / self.d[i - self.m]
        U = self.d[:, None] * L.T
        return L, U


def generic_ilu0(A):
    """Independent pattern-restricted Gaussian elimination for small checks."""
    packed = A.copy()
    pattern = A != 0
    for i in range(len(A)):
        for j in range(i):
            if not pattern[i, j]:
                continue
            packed[i, j] /= packed[j, j]
            for k in range(j + 1, len(A)):
                if pattern[i, k]:
                    packed[i, k] -= packed[i, j] * packed[j, k]
    return np.tril(packed, -1) + np.eye(len(A)), np.triu(packed)


def check_zero_fill():
    for m in (2, 3, 5):
        A = stencil(np.eye(m * m), m)
        ilu = PoissonILU0(m)
        L, U = ilu.matrices()
        reference_L, reference_U = generic_ilu0(A)
        np.testing.assert_allclose(L, reference_L, atol=1e-15)
        np.testing.assert_allclose(U, reference_U, atol=1e-15)
        assert not np.any((L + U)[A == 0])
        # The product can have fill; the no-fill condition concerns the factors.
        M = L @ U
        np.testing.assert_allclose(M, M.T, atol=1e-15)
        assert np.linalg.eigvalsh(M)[0] > 0
        rhs = np.arange(1., m * m + 1)
        np.testing.assert_allclose(M @ ilu.solve(rhs), rhs, atol=3e-14)


def gmres_history(matvec, b, solve=None, tolerance=1e-8, limit=160):
    """Unrestarted (optionally right-preconditioned) GMRES.

    Two orthogonalization passes preserve the basis; incremental Givens
    rotations supply the residual estimate. Every plotted norm is also
    checked explicitly in the original system. Additional solves used to
    form these intermediate approximations are for verification/plotting.
    """
    n = len(b)
    Q = np.zeros((n, limit + 1))
    H = np.zeros((limit + 1, limit))
    R = np.zeros_like(H)
    beta = np.linalg.norm(b)
    Q[:, 0] = b / beta
    g = np.zeros(limit + 1)
    g[0] = beta
    rotations, history = [], [1.]
    for j in range(limit):
        direction = Q[:, j] if solve is None else solve(Q[:, j])
        w = matvec(direction)
        for _ in range(2):
            coefficients = Q[:, :j + 1].T @ w
            H[:j + 1, j] += coefficients
            w -= Q[:, :j + 1] @ coefficients
        H[j + 1, j] = np.linalg.norm(w)
        assert H[j + 1, j] > 1e-14
        Q[:, j + 1] = w / H[j + 1, j]
        R[:, j] = H[:, j]
        for i, (c, s) in enumerate(rotations):
            a, h = R[i:i + 2, j]
            R[i:i + 2, j] = [c * a + s * h, -s * a + c * h]
        a, h = R[j:j + 2, j]
        magnitude = np.hypot(a, h)
        c, s = a / magnitude, h / magnitude
        rotations.append((c, s))
        R[j:j + 2, j] = [magnitude, 0.]
        a, h = g[j:j + 2]
        g[j:j + 2] = [c * a + s * h, -s * a + c * h]
        k = j + 1
        y = np.linalg.solve(R[:k, :k], g[:k])
        x = Q[:, :k] @ y
        if solve is not None:
            x = solve(x)
        residual = np.linalg.norm(b - matvec(x)) / beta
        np.testing.assert_allclose(residual, abs(g[k]) / beta, atol=3e-13)
        assert residual <= history[-1] + 1e-13
        history.append(residual)
        if residual <= tolerance:
            np.testing.assert_allclose(Q[:, :k + 1].T @ Q[:, :k + 1],
                                       np.eye(k + 1), atol=2e-14)
            return np.asarray(history)
    raise RuntimeError("GMRES failed to reach the requested tolerance.")


def example_data(m=30):
    check_zero_fill()
    n, h = m * m, 1 / (m + 1)
    ilu = PoissonILU0(m)
    A0 = stencil(np.eye(n), m)
    # Work with h^2 A and h^2 M internally. Their preconditioned operator
    # equals M^{-1} A, independent of this common scaling.
    left = ilu.forward(A0)
    transformed = ilu.forward(left.T).T
    np.testing.assert_allclose(transformed, transformed.T, atol=1e-15)
    values, vectors = np.linalg.eigh(transformed)
    assert np.all(values > 0)
    # Verify all generalized eigenpairs A0 u = lambda M0 u.
    u = ilu.backward(vectors)
    np.testing.assert_allclose(ilu.solve(stencil(u, m)), u * values,
                               atol=3e-13)
    p = np.arange(1, m + 1)
    one_dimensional = 4 * np.sin(p * np.pi * h / 2) ** 2
    analytic = np.sort((one_dimensional[:, None] + one_dimensional).ravel())
    np.testing.assert_allclose(np.linalg.eigvalsh(A0), analytic, atol=5e-14)
    # A sine-product eigenvector independently checks the indexing and scale.
    i, j = np.meshgrid(p, p, indexing="ij")
    for a, b in ((1, 1), (2, 7), (m, m)):
        v = (np.sin(a * np.pi * i * h) * np.sin(b * np.pi * j * h)).ravel()
        eigenvalue = one_dimensional[a - 1] + one_dimensional[b - 1]
        np.testing.assert_allclose(stencil(v, m), eigenvalue * v, atol=4e-14)

    rhs = np.random.default_rng(302).normal(size=n)
    rhs /= np.linalg.norm(rhs)
    raw = gmres_history(lambda v: stencil(v, m) / h**2, rhs)
    preconditioned = gmres_history(lambda v: stencil(v, m) / h**2, rhs,
                                  solve=lambda r: h**2 * ilu.solve(r))
    return analytic / h**2, values, raw, preconditioned


def plot(raw_values, preconditioned_values, raw_history, ilu_history):
    fig = plt.figure(figsize=(11.8, 7.6))
    grid = fig.add_gridspec(2, 2, height_ratios=(1, 1.2), hspace=.70, wspace=.30)
    left = fig.add_subplot(grid[0, 0])
    right = fig.add_subplot(grid[0, 1])
    bottom = fig.add_subplot(grid[1, :])
    for ax in (left, right, bottom):
        ax.spines[["top", "right"]].set_visible(False)
        ax.spines[["left", "bottom"]].set_color(GRAY)
        ax.tick_params(colors=INK, labelsize=12)
        ax.set_axisbelow(True)
        ax.grid(axis="y", color=GRAY, alpha=.5, linewidth=.8)

    left.hist(raw_values, bins=np.linspace(0, 8000, 25),
              color=BLUE, alpha=.85, edgecolor="white", linewidth=.6)
    left.set_title("Before preconditioning", fontsize=17, color=INK, pad=14)
    left.set_xlabel(r"Eigenvalues of $A$", color=INK, labelpad=8)
    left.set_ylabel("Number of eigenvalues", color=INK, labelpad=8)
    left.set_xlim(0, 8000)
    left.set_xticks([0, 2000, 4000, 6000, 8000])

    right.axvspan(.75, 1.25, color=TEAL, alpha=.08)
    right.hist(preconditioned_values, bins=np.linspace(0, 1.25, 31),
               color=PURPLE, alpha=.85, edgecolor="white", linewidth=.6)
    right.set_title("After ILU(0)", fontsize=17, color=INK, pad=14)
    right.set_xlabel(r"Eigenvalues of $M^{-1}A$", color=INK, labelpad=8)
    right.set_ylabel("Number of eigenvalues", color=INK, labelpad=8)
    right.set_xlim(0, 1.25)
    right.set_xticks([0, .25, .5, .75, 1., 1.25])
    fraction = np.mean((preconditioned_values >= .75) & (preconditioned_values <= 1.25))
    right.text(.04, .92, f"{100 * fraction:.0f}% lie in [0.75, 1.25]",
               transform=right.transAxes, color=TEAL, fontsize=13, va="top")

    bottom.set_yscale("log")
    bottom.set_title("Effect on GMRES convergence", fontsize=17, color=INK, pad=14)
    for history, color, label in ((raw_history, BLUE, "Without preconditioning"),
                                  (ilu_history, PURPLE, "With ILU(0)")):
        steps = np.arange(len(history))
        bottom.plot(steps, history, color=color, linewidth=2.1, label=label)
        bottom.scatter(steps[-1], history[-1], color=color, s=35, zorder=3)
        bottom.annotate(f"{steps[-1]} steps", xy=(steps[-1], history[-1]),
                        xytext=(steps[-1] + 3, 3e-7), fontsize=13, color=color,
                        arrowprops={"arrowstyle": "->", "color": color},
                        va="bottom")
    bottom.axhline(1e-8, color=TEAL, linewidth=1.2, linestyle=(0, (4, 3)))
    bottom.text(57, 1.5e-8, r"Target: $10^{-8}$", color=TEAL, fontsize=12)
    bottom.set_xlim(0, 105)
    bottom.set_ylim(1e-9, 2.)
    bottom.set_yticks([1., 1e-2, 1e-4, 1e-6, 1e-8])
    bottom.set_xlabel(r"GMRES iteration $k$", color=INK, labelpad=8)
    bottom.set_ylabel(r"Relative residual $\|b-Ax_k\|_2/\|b\|_2$",
                      color=INK, labelpad=8)
    bottom.legend(loc="upper right", frameon=False, fontsize=13)
    fig.subplots_adjust(left=.09, right=.97, bottom=.09, top=.93)
    destination = Path(__file__).resolve().parents[1] / "_static"
    fig.savefig(destination / "poisson_ilu0_spectrum.svg", bbox_inches="tight",
                metadata={"Date": None})
    fig.savefig("/tmp/nla-poisson-ilu0-spectrum.png", dpi=180, bbox_inches="tight")
    fig.savefig("/tmp/nla-poisson-ilu0-spectrum-book-width.png", dpi=65,
                bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    raw, preconditioned, raw_history, ilu_history = example_data()
    plot(raw, preconditioned, raw_history, ilu_history)
    for name, values in (("A", raw), ("M^-1 A", preconditioned)):
        print(f"{name}: min {values.min():.9g}; max {values.max():.9g}; "
              f"eigenvalue ratio {values.max() / values.min():.6g}")
    print(f"GMRES steps: {len(raw_history) - 1} without preconditioning; "
          f"{len(ilu_history) - 1} with strict ILU(0).")
    print("Verified strict zero fill, sine eigenvalues, all preconditioned eigenpairs, "
          "and explicit GMRES residuals.")
