"""Shared layout and colors for orthogonal triangularization figures."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle
import numpy as np


def plot_sequence(states, active_rows, active_columns, eliminated_groups,
                  titles, actions, row_labels, footer, stem, active_masks=None):
    """Highlight the next update, optionally using masks for two-sided actions."""
    plt.rcParams.update({"font.size": 13, "svg.fonttype": "path"})
    nrows, ncols = states[0].shape
    cell = 0.55
    pitch = (ncols - 1) * cell + 1.6
    left_edges = [stage * pitch for stage in range(len(states))]
    right = left_edges[-1] + (ncols - 1) * cell + 0.8
    width = right + 1.85
    fig, ax = plt.subplots(figsize=(width * 12.5 / 11.85, 5.1))

    blue_fill, blue_edge = "#dbeafe", "#60a5fa"
    fill_in_color = "#2563eb"
    green, orange = "#15803d", "#c26a13"
    ink, muted = "#334155", "#94a3b8"
    lo, hi = -(nrows - 1) * cell - 0.32, 0.32
    action_y = lo - 0.5

    for stage, (left, state) in enumerate(zip(left_edges, states)):
        center = left + (ncols - 1) * cell / 2
        ax.text(center, 0.72, titles[stage], ha="center", va="center", fontsize=15, color=ink)
        ax.text(center, action_y, actions[stage], ha="center", va="top", fontsize=12.5, color=ink)
        created_entries = {entry for group in eliminated_groups[:stage] for entry in group}

        for row in range(nrows):
            for col in range(ncols):
                x, y = left + col * cell, -row * cell
                if active_masks is not None and stage < len(active_masks):
                    active = active_masks[stage][row, col]
                else:
                    active = (stage < len(active_rows) and row in active_rows[stage]
                              and col in active_columns[stage])
                if active:
                    ax.add_patch(Rectangle(
                        (x - 0.25, y - 0.24), 0.5, 0.48,
                        facecolor=blue_fill, edgecolor=blue_edge, linewidth=1., zorder=0,
                    ))
                zero = np.isclose(state[row, col], 0.)
                created = (row, col) in created_entries
                fill_in = np.isclose(states[0][row, col], 0.) and not zero
                color = (green if created else muted if zero
                         else fill_in_color if fill_in else ink)
                ax.text(
                    x, y, "0" if zero else r"$\bullet$",
                    ha="center", va="center", fontsize=21 if zero else 19,
                    color=color, fontweight="bold" if created else "normal",
                )

        for x, direction in [(left - 0.38, 1), (left + (ncols - 1) * cell + 0.38, -1)]:
            ax.plot(
                [x + direction * 0.12, x, x, x + direction * 0.12],
                [hi, hi, lo, lo], color=ink, linewidth=1.5,
            )
        if stage < len(active_rows):
            for row, col in eliminated_groups[stage]:
                ax.add_patch(Circle((left + col * cell, -row * cell), 0.22,
                                    fill=False, edgecolor=orange, linewidth=2.2))
            ax.annotate("", xy=(left_edges[stage + 1] - 0.63, -(nrows - 1) * cell / 2),
                        xytext=(left + (ncols - 1) * cell + 0.57, -(nrows - 1) * cell / 2),
                        arrowprops={"arrowstyle": "-|>", "color": ink,
                                    "lw": 1.7, "mutation_scale": 17})

    for row, label in enumerate(row_labels):
        ax.text(-0.6, -row * cell, label, ha="right", va="center", fontsize=12, color=ink)

    ax.text((left_edges[-1] + (ncols - 1) * cell) / 2, action_y - 0.78, footer,
            ha="center", va="center", fontsize=15, color=ink)
    ax.set_aspect("equal")
    ax.set_xlim(-1.85, right)
    ax.set_ylim(action_y - 1.13, 1.1)
    ax.axis("off")
    fig.tight_layout(pad=0.5)

    destination = Path(__file__).resolve().parents[1] / "_static"
    fig.savefig(destination / f"{stem}.svg", bbox_inches="tight")
    fig.savefig(f"/tmp/nla-{stem.replace('_', '-')}.png", dpi=180, bbox_inches="tight")
    plt.close(fig)
