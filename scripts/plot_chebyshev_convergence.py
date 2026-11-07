"""Plot the interval error bounds for Chebyshev and Richardson iteration.

Run from any directory with Python, NumPy, and Matplotlib installed.
The curves are theoretical bounds, not errors measured for a chosen system.
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
    "svg.hashsalt": "chebyshev-convergence",
})
blue, purple = "#2563eb", "#7c3aed"
orange, ink, gray = "#c26a13", "#334155", "#cbd5e1"
TOLERANCE = 1e-6


def comparison(kappa):
    """Compute exact interval bounds and verify the polynomial interpretation."""
    q = (np.sqrt(kappa) - 1) / (np.sqrt(kappa) + 1)
    richardson_factor = (kappa - 1) / (kappa + 1)
    last = int(np.ceil(np.log(TOLERANCE) / np.log(richardson_factor)))
    tick_spacing = 20 if kappa == 10 else 200
    plot_limit = int(np.ceil(1.08 * last / tick_spacing)) * tick_spacing
    steps = np.arange(plot_limit + 1)
    chebyshev = 2 * q ** steps / (1 + q ** (2 * steps))
    richardson = richardson_factor ** steps
    first_chebyshev = int(np.flatnonzero(chebyshev <= TOLERANCE)[0])
    first_richardson = int(np.flatnonzero(richardson <= TOLERANCE)[0])

    np.testing.assert_allclose(chebyshev[0], 1.)
    np.testing.assert_allclose(richardson[0], 1.)
    np.testing.assert_allclose(chebyshev[1], richardson[1])
    assert np.all(chebyshev <= richardson + 5 * np.finfo(float).eps)
    assert np.all(np.diff(chebyshev) < 0)
    assert np.all(np.diff(richardson) < 0)

    # On [1, kappa], both interval maxima are attained at the endpoints.
    # Check the formulas against independently evaluated polynomials.
    alpha, beta = 1., float(kappa)
    c, d = (alpha + beta) / 2, (beta - alpha) / 2
    eigenvalues = np.linspace(alpha, beta, 1001)
    for k in (0, 1, 5, 20, first_chebyshev):
        coefficients = [0.] * k + [1.]
        p = chebval((c - eigenvalues) / d, coefficients)
        p /= chebval(c / d, coefficients)
        np.testing.assert_allclose(np.max(np.abs(p)), chebyshev[k], rtol=2e-13)
        richardson_p = (1 - eigenvalues / c) ** k
        np.testing.assert_allclose(
            np.max(np.abs(richardson_p)), richardson[k], rtol=2e-13,
        )
    for bound, first in ((chebyshev, first_chebyshev),
                         (richardson, first_richardson)):
        assert bound[first] <= TOLERANCE < bound[first - 1]

    return steps, chebyshev, richardson, first_chebyshev, first_richardson


def plot():
    fig, axes = plt.subplots(1, 2, figsize=(11.8, 5.1))
    for ax, kappa in zip(axes, (10, 100)):
        steps, chebyshev, richardson, kc, kr = comparison(kappa)
        ax.set_yscale("log")
        ax.grid(axis="y", which="major", color=gray, alpha=.6, linewidth=.8)
        ax.set_axisbelow(True)
        ax.spines[["top", "right"]].set_visible(False)
        ax.spines[["left", "bottom"]].set_color(gray)
        ax.tick_params(colors=ink, labelsize=12)
        ax.tick_params(axis="y", which="minor", left=False)

        ax.plot(steps, richardson, color=blue, linewidth=2.2)
        ax.plot(steps, chebyshev, color=purple, linewidth=2.2)
        ax.axhline(TOLERANCE, color=orange, linewidth=1.3,
                   linestyle=(0, (4, 3)))
        ax.text(.55, TOLERANCE * 1.4, r"Target: $10^{-6}$",
                transform=ax.get_yaxis_transform(), ha="center", va="bottom",
                color=orange, fontsize=12)

        for count, values, color, label_position in (
            (kc, chebyshev, purple, (.38, .39)),
            (kr, richardson, blue, (.80, .39)),
        ):
            ax.scatter(count, values[count], color=color, s=38, zorder=4)
            ax.annotate(
                f"{count} steps", xy=(count, values[count]),
                xytext=label_position, textcoords="axes fraction",
                color=color, fontsize=13,
                arrowprops={"arrowstyle": "->", "color": color,
                            "linewidth": 1.1},
                ha="center", va="bottom",
            )

        ax.set_title(rf"Spectral ratio $\kappa={kappa}$",
                     fontsize=17, color=ink, pad=15)
        ax.set_xlim(0, steps[-1])
        ax.set_ylim(1e-7, 2.)
        ax.set_yticks([1., 1e-2, 1e-4, 1e-6])
        ax.set_xticks(np.arange(0, 81, 20) if kappa == 10
                      else np.arange(0, 801, 200))
        ax.set_xlabel(r"Iteration count $k$", color=ink, labelpad=9)
    axes[0].set_ylabel("Relative error bound", color=ink, labelpad=9)

    handles = [Line2D([], [], color=purple, linewidth=2.2,
                      label="Chebyshev bound"),
               Line2D([], [], color=blue, linewidth=2.2,
                      label="Richardson bound")]
    fig.legend(handles=handles, loc="upper center", ncol=2,
               bbox_to_anchor=(.53, 1.01), frameon=False, fontsize=14,
               columnspacing=2.2)
    fig.subplots_adjust(left=.09, right=.98, bottom=.16, top=.83, wspace=.26)

    destination = Path(__file__).resolve().parents[1] / "_static"
    fig.savefig(destination / "chebyshev_convergence.svg", bbox_inches="tight",
                metadata={"Date": None})
    fig.savefig("/tmp/nla-chebyshev-convergence.png", dpi=180,
                bbox_inches="tight")
    fig.savefig("/tmp/nla-chebyshev-convergence-book-width.png", dpi=65,
                bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    plot()
    print("Verified interval maxima and tolerance crossings for both methods.")
    for kappa in (10, 100):
        *_, kc, kr = comparison(kappa)
        print(f"kappa={kappa}: Chebyshev {kc} steps; Richardson {kr} steps.")
