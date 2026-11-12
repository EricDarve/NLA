"""Compare computed CG energy errors for spread and clustered SPD spectra.

Run from the repository with Python, NumPy, and Matplotlib installed.
The diagonal matrices express the experiment in an eigenvector basis.
No convergence curve is prescribed: every plotted error comes from CG iterates.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

BLUE, PURPLE, TEAL = "#2563eb", "#7c3aed", "#0f766e"
INK, GRAY = "#334155", "#cbd5e1"
plt.rcParams.update({
    "font.size": 15, "text.color": INK, "axes.labelcolor": INK,
    "xtick.color": INK, "ytick.color": INK, "svg.fonttype": "path",
    "svg.hashsalt": "cg-spectrum-convergence",
})


def cg_errors(eigenvalues, exact, tolerance=1e-12, max_steps=120):
    """Run ordinary double-precision CG from zero; measure the true error."""
    b = eigenvalues * exact
    x = np.zeros_like(b)
    r = b.copy()
    p = r.copy()
    initial_energy = np.sqrt(exact @ (eigenvalues * exact))
    errors = [1.]
    for k in range(1, max_steps + 1):
        ap = eigenvalues * p
        rr = r @ r
        alpha = rr / (p @ ap)
        x += alpha * p
        r -= alpha * ap
        error = exact - x
        energy = np.sqrt(error @ (eigenvalues * error)) / initial_energy
        errors.append(energy)

        # Check the recursively updated residual against the original system.
        np.testing.assert_allclose(
            r, b - eigenvalues * x, rtol=0.,
            atol=5e-13 * np.linalg.norm(b),
        )
        if k <= 6:
            # Independent minimization in the Krylov space:
            # minimize ||A^(1/2)(exact - B y)||_2 by dense least squares.
            scaled = eigenvalues / eigenvalues.max()
            basis = np.column_stack([b * scaled**j for j in range(k)])
            root = np.sqrt(eigenvalues)
            coeff = np.linalg.lstsq(root[:, None] * basis, root * exact,
                                    rcond=None)[0]
            np.testing.assert_allclose(x, basis @ coeff, rtol=2e-8, atol=2e-10)
        if energy <= tolerance:
            break
        p = r + (r @ r) / rr * p
    else:
        raise AssertionError("The example did not reach the requested accuracy.")
    errors = np.asarray(errors)
    assert np.all(np.diff(errors) <= 5e-15), "Energy error must decrease."
    rho = (np.sqrt(eigenvalues.max() / eigenvalues.min()) - 1) / (
        np.sqrt(eigenvalues.max() / eigenvalues.min()) + 1)
    assert np.all(errors <= 2 * rho**np.arange(len(errors)) + 5e-15)
    return errors


def plot():
    n = 120
    spectra = [
        np.linspace(1., 100., n),
        np.concatenate([
            np.linspace(1., 1.1, 40),
            np.linspace(19.95, 20.05, 40),
            np.linspace(99.9, 100., 40),
        ]),
    ]
    # Identical eigenvector coefficients in the initial error for both systems.
    exact = np.ones(n) / np.sqrt(n)
    histories = []
    for eigenvalues in spectra:
        assert len(np.unique(eigenvalues)) == n
        assert eigenvalues.min() == 1. and eigenvalues.max() == 100.
        histories.append(cg_errors(eigenvalues, exact))

    fig, (left, right) = plt.subplots(1, 2, figsize=(11.8, 5.8),
                                     gridspec_kw={"width_ratios": [1, 1.15]})
    for ax in (left, right):
        ax.spines[["top", "right"]].set_visible(False)
        ax.spines[["left", "bottom"]].set_color(GRAY)
        ax.tick_params(labelsize=12)
    left.set_title("Same endpoints, different spectra", fontsize=17, pad=20)
    for values, row, color in zip(spectra, (1., 0.), (BLUE, PURPLE)):
        left.hlines(row, 1, 100, color=GRAY, lw=1, zorder=0)
        left.scatter(values, np.full(n, row), marker="|", s=240,
                     linewidths=1.2, color=color, zorder=3)
    left.set_yticks([1, 0], ["Spread", "Clustered"])
    left.tick_params(axis="y", length=0)
    left.set_ylim(-.8, 1.7)
    left.set_xlim(-4, 105)
    left.set_xticks([1, 20, 40, 60, 80, 100])
    left.set_xlabel(r"Eigenvalue $\lambda$")
    left.text(50, 1.4, "120 distinct eigenvalues in each system",
              fontsize=12, ha="center", color=INK)
    for center in (1.05, 20., 99.95):
        left.annotate("40", xy=(center, -.12), xytext=(center, -.48),
                      ha="center", fontsize=12, color=PURPLE,
                      arrowprops={"arrowstyle": "-", "color": PURPLE, "lw": .9})
    left.text(.5, -.2, r"Both have $\lambda_{\min}=1$, $\lambda_{\max}=100$",
              ha="center", transform=left.transAxes, fontsize=13)

    right.set_title("Measured CG energy error", fontsize=17, pad=20)
    max_k = max(len(h) - 1 for h in histories)
    steps = np.arange(max_k + 1)
    rho = 9. / 11.
    right.semilogy(steps, 2 * rho**steps, color=TEAL, ls="--", lw=2,
                   label=r"Common bound $2(9/11)^k$")
    for history, color, name in zip(histories, (BLUE, PURPLE), ("Spread", "Clustered")):
        right.semilogy(np.arange(len(history)), history, color=color, lw=2.4,
                       label=name)
        right.plot(len(history) - 1, history[-1], "o", color=color, ms=5)
    right.axhline(1e-12, color=GRAY, ls=":", lw=1.2)
    right.text(44, 1.6e-12, r"Stopping level: $10^{-12}$", fontsize=11, color=INK)
    right.set_xlim(0, 76)
    right.set_ylim(1e-13, 3)
    right.set_yticks([1., 1e-3, 1e-6, 1e-9, 1e-12])
    right.set_xlabel(r"CG iteration $k$")
    right.set_ylabel(r"$\|x-x_k\|_A\,/\,\|x-x_0\|_A$", fontsize=15)
    right.grid(axis="y", color=GRAY, alpha=.6, lw=.7)
    right.legend(loc="upper right", frameon=False, fontsize=11)
    fig.text(.5, .98, r"Same condition number $\kappa_2(A)=100$; faster convergence with clustering",
             ha="center", fontsize=18, color=INK)
    fig.subplots_adjust(left=.09, right=.98, bottom=.23, top=.81, wspace=.43)
    root = Path(__file__).resolve().parents[1]
    fig.savefig(root / "_static/cg_convergence.svg", bbox_inches="tight",
                metadata={"Date": None})
    for suffix, dpi in (("", 180), ("-book-width", 65)):
        fig.savefig(f"/tmp/nla-cg-convergence{suffix}.png", dpi=dpi,
                    bbox_inches="tight")
    plt.close(fig)
    for name, history in zip(("Spread", "Clustered"), histories):
        print(f"{name}: {len(history)-1} steps, relative energy error {history[-1]:.3e}")
    print("Verified SPD spectra, CG residuals, independent Krylov minimization, "
          "energy decrease, and the common condition-number bound.")


if __name__ == "__main__":
    plot()
