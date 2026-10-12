"""Illustrate adjacent-row Givens QR, working upward within each column."""

import numpy as np
from givens_matrix_diagrams import plot_sequence


initial = np.array([[4., 2., 1.], [3., 5., 2.], [1., 2., 6.]])
current = initial.copy()
transform = np.eye(3)
row_pairs = [(1, 2), (0, 1), (1, 2)]
columns = [0, 0, 1]
targets = [(2, 0), (1, 0), (2, 1)]
states = [current.copy()]

for step, ((upper, lower), col) in enumerate(zip(row_pairs, columns)):
    a, b = current[upper, col], current[lower, col]
    length = np.hypot(a, b)
    c, s = a / length, b / length
    rotation = np.array([[c, s], [-s, c]])
    current[[upper, lower], col:] = rotation @ current[[upper, lower], col:]
    transform[[upper, lower], :] = rotation @ transform[[upper, lower], :]
    current[lower, col] = 0.
    for row, column in targets[:step + 1]:
        assert np.isclose(current[row, column], 0.)
    assert np.allclose(current.T @ current, initial.T @ initial)
    states.append(current.copy())

assert np.allclose(np.tril(current, -1), 0.)
assert np.allclose(transform @ initial, current)
assert np.allclose(transform.T @ transform, np.eye(3))

plot_sequence(
    states,
    active_rows=row_pairs,
    active_columns=[range(col, 3) for col in columns],
    eliminated_groups=[[entry] for entry in targets],
    titles=["Start", "Bottom entry cleared", "First column finished", "Upper triangular"],
    actions=[
        "Column 1: rotate\nrows 2 and 3",
        "Column 1: rotate\nrows 1 and 2",
        "Column 2: rotate\nrows 2 and 3",
        "Upper triangular\nfactor",
    ],
    row_labels=["row 1", "row 2", "row 3"],
    footer="Work upward within each column, then move to the next column.",
    stem="givens_dense_qr",
)
