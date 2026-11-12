"""Show an A-inner product as an ordinary dot product after coordinate scaling."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Arc
import numpy as np

plt.rcParams.update({"font.size": 15, "svg.fonttype": "path",
                     "svg.hashsalt": "cg-inner-product"})
BLUE, PURPLE, TEAL = "#2563eb", "#7c3aed", "#0f766e"
INK, GRAY = "#334155", "#cbd5e1"


def plot():
    a = np.diag([4., 1.])
    scale = np.diag(np.sqrt(np.diag(a)))
    y, z = np.array([.5, 1.]), np.array([.5, -1.])
    ty, tz = scale @ y, scale @ z
    np.testing.assert_allclose(scale.T @ scale, a)
    np.testing.assert_allclose(y @ z, -.75)
    np.testing.assert_allclose(y @ a @ z, ty @ tz)
    np.testing.assert_allclose(ty @ tz, 0.)
    for v in (y, z):
        np.testing.assert_allclose(v @ a @ v, np.linalg.norm(scale @ v)**2)

    fig, axes = plt.subplots(1, 2, figsize=(11.8, 5.9))
    for j, ax in enumerate(axes):
        ax.set_aspect("equal")
        ax.set_xlim(-.7, 1.9)
        ax.set_ylim(-1.55, 1.55)
        # The same grid is mapped along with the vectors: horizontal
        # distances double, while vertical distances stay unchanged.
        for v in np.arange(-1.5, 2.1, .5):
            ax.axvline(v*(2 if j else 1), color=GRAY, lw=.7, alpha=.55, zorder=0)
            ax.axhline(v, color=GRAY, lw=.7, alpha=.55, zorder=0)
        ax.spines[["top", "right"]].set_visible(False)
        for spine in ("left", "bottom"):
            ax.spines[spine].set_position("zero")
            ax.spines[spine].set_color(INK)
        ax.set_xticks([.5, 1., 1.5])
        ax.set_yticks([-1., -.5, .5, 1.])
        ax.tick_params(colors=INK, labelsize=11)
        ax.text(1.78, -.2, r"$v_1$" if j == 0 else r"$\widetilde v_1$", color=INK)
        ax.text(-.24, 1.4, r"$v_2$" if j == 0 else r"$\widetilde v_2$", color=INK)
        ax.scatter(0, 0, color=INK, s=22, zorder=6)

    for ax, vectors, labels in (
        (axes[0], (y, z), (r"$y=(\frac{1}{2},1)^T$", r"$z=(\frac{1}{2},-1)^T$")),
        (axes[1], (ty, tz), (r"$\widetilde y=(1,1)^T$", r"$\widetilde z=(1,-1)^T$")),
    ):
        for v, label, color in zip(vectors, labels, (BLUE, PURPLE)):
            ax.annotate("", xy=v, xytext=(0, 0), zorder=4,
                        arrowprops={"arrowstyle": "->", "lw": 2.5, "color": color,
                                    "shrinkA": 0, "shrinkB": 0})
            ax.text(v[0]+.08, v[1]+(.07 if v[1]>0 else -.17), label,
                    color=color, fontsize=14,
                    bbox={"facecolor": "white", "edgecolor": "none", "pad": 1.})
    angle = np.degrees(np.arctan2(y[1], y[0]))
    axes[0].add_patch(Arc((0, 0), .65, .65, theta1=-angle, theta2=angle,
                          color=TEAL, lw=1.5, zorder=3))
    corner = np.array([[.14, .14], [.28, 0.], [.14, -.14]])
    axes[1].plot(corner[:, 0], corner[:, 1], color=TEAL, lw=1.5, zorder=5)
    axes[0].set_title("Original coordinates", color=INK, fontsize=18, pad=18)
    axes[1].set_title("First coordinate multiplied by 2", color=INK, fontsize=18, pad=18)
    for ax, formula in zip(axes, (r"$y^Tz=-\frac{3}{4}$", r"$\widetilde y^T\widetilde z=y^TAz=0$")):
        ax.text(.5, -.15, formula, transform=ax.transAxes, ha="center",
                color=TEAL, fontsize=17)
    fig.text(.5, .97, r"$A=\operatorname{diag}(4,1):\quad \widetilde v=A^{1/2}v=(2v_1,v_2)^T$",
             ha="center", color=INK, fontsize=18)
    fig.text(.5, .025, "An ordinary dot product in rescaled coordinates is the A-inner product in the original coordinates.",
             ha="center", color=INK, fontsize=12)
    fig.subplots_adjust(left=.055, right=.97, bottom=.2, top=.82, wspace=.3)
    root = Path(__file__).resolve().parents[1]
    fig.savefig(root/"_static/cg_inner_product.svg", bbox_inches="tight", metadata={"Date": None})
    for suffix, dpi in (("", 180), ("-book-width", 65)):
        fig.savefig(f"/tmp/nla-cg-inner-product{suffix}.png", dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print("Verified weighted dot product, transformed dot product, and energy norms.")


if __name__ == "__main__":
    plot()
