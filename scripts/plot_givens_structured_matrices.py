"""Generate matching Hessenberg and tridiagonal Givens QR diagrams."""

import numpy as np
from givens_matrix_diagrams import plot_sequence


def elimination_states(initial, tridiagonal=False):
    current = initial.copy()
    states = [current.copy()]
    for k in range(3):
        a, b = current[k, k], current[k + 1, k]
        length = np.hypot(a, b)
        c, s = a / length, b / length
        rotation = np.array([[c, s], [-s, c]])
        current[[k, k + 1], k:] = rotation @ current[[k, k + 1], k:]
        current[k + 1, k] = 0.
        assert np.allclose(np.tril(current, -2), 0.)
        assert np.allclose(np.tril(current[:, :k + 1], -1), 0.)
        assert np.allclose(current.T @ current, initial.T @ initial)
        if tridiagonal:
            assert np.allclose(np.triu(current, 3), 0.)
        states.append(current.copy())
    assert np.allclose(np.tril(current, -1), 0.)
    return states


hessenberg = np.array([
    [4., 2., 3., 1.],
    [1., 5., 2., 3.],
    [0., 1., 6., 2.],
    [0., 0., 1., 7.],
])
tridiagonal = np.diag([4., 5., 6., 7.]) + np.diag(np.ones(3), 1) + np.diag(np.ones(3), -1)

for initial, stem, footer, banded in [
    (hessenberg, "givens_hessenberg", "Clear the subdiagonal from top to bottom.", False),
    (tridiagonal, "givens_tridiagonal", "Each rotation touches at most three columns.", True),
]:
    states = elimination_states(initial, tridiagonal=banded)
    if banded:
        # Only these two positions become nonzero outside the original band.
        fill_in = np.isclose(initial, 0.) & ~np.isclose(states[-1], 0.)
        assert set(map(tuple, np.argwhere(fill_in))) == {(0, 2), (1, 3)}
    plot_sequence(
        states,
        active_rows=[(k, k + 1) for k in range(3)],
        active_columns=[range(k, min(k + 3, 4) if banded else 4) for k in range(3)],
        eliminated_groups=[[(k + 1, k)] for k in range(3)],
        titles=["Start", "One entry cleared", "Two entries cleared", "Upper triangular"],
        actions=["Rotate rows\n1 and 2", "Rotate rows\n2 and 3", "Rotate rows\n3 and 4", "Upper triangular\nfactor"],
        row_labels=["row 1", "row 2", "row 3", "row 4"],
        footer=footer,
        stem=stem,
    )
