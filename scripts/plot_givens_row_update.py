"""Illustrate the left-to-right Givens update after adding a matrix row."""

import numpy as np
from givens_matrix_diagrams import plot_sequence

# Compute the states rather than merely sketch their zero patterns.
initial = np.array([[4., 2., -1.], [0., 3., 2.], [0., 0., 5.], [3., 1., 2.]])
current = initial.copy()
states = [current.copy()]
for k in range(3):
    a, b = current[k, k], current[3, k]
    length = np.hypot(a, b)
    c, s = a / length, b / length
    rotation = np.array([[c, s], [-s, c]])
    current[[k, 3], k:] = rotation @ current[[k, 3], k:]
    current[3, k] = 0.
    assert np.allclose(np.tril(current[:3], -1), 0.)
    assert np.allclose(current[3, :k + 1], 0.)
    assert np.allclose(current.T @ current, initial.T @ initial)
    states.append(current.copy())

plot_sequence(
    states,
    active_rows=[(k, 3) for k in range(3)],
    active_columns=[range(k, 3) for k in range(3)],
    eliminated_groups=[[(3, k)] for k in range(3)],
    titles=["Start", "One entry cleared", "Two entries cleared", "Triangular again"],
    actions=[
        "Rotate row 1\nwith the added row",
        "Rotate row 2\nwith the added row",
        "Rotate row 3\nwith the added row",
        "Updated factor above;\nadded row is zero",
    ],
    row_labels=["row 1", "row 2", "row 3", "added row"],
    footer="Clear the added row from left to right.",
    stem="givens_row_update",
)
