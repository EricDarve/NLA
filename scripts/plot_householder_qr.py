"""Show Householder QR clearing a whole column segment at each step."""

import numpy as np
from givens_matrix_diagrams import plot_sequence


initial = np.array([[4., 2., 1.], [3., 5., 2.], [1., 2., 6.], [2., 1., 4.]])
current = initial.copy()
transform = np.eye(4)
states = [current.copy()]

for k in range(3):
    x = current[k:, k].copy()
    sigma = -np.copysign(np.linalg.norm(x), x[0])
    v = x.copy()
    v[0] -= sigma
    reflection = np.eye(len(x)) - 2 * np.outer(v, v) / (v @ v)
    target = np.zeros(len(x))
    target[0] = sigma
    assert np.allclose(reflection @ x, target)
    assert np.allclose(reflection.T @ reflection, np.eye(len(x)))

    current[k:, k:] = reflection @ current[k:, k:]
    transform[k:, :] = reflection @ transform[k:, :]
    current[k, k] = sigma
    current[k + 1:, k] = 0.
    assert np.allclose(np.tril(current[:, :k + 1], -1), 0.)
    assert np.allclose(current.T @ current, initial.T @ initial)
    states.append(current.copy())

assert np.allclose(np.tril(current, -1), 0.)
assert np.allclose(transform @ initial, current)
assert np.allclose(transform.T @ transform, np.eye(4))

plot_sequence(
    states,
    active_rows=[range(k, 4) for k in range(3)],
    active_columns=[range(k, 3) for k in range(3)],
    eliminated_groups=[[(row, k) for row in range(k + 1, 4)] for k in range(3)],
    titles=["Start", "First column finished", "Two columns finished", "Upper trapezoidal"],
    actions=[
        "Reflect rows\n1 through 4",
        "Reflect rows\n2 through 4",
        "Reflect rows\n3 and 4",
        "Triangular block above;\nzero row below",
    ],
    row_labels=["row 1", "row 2", "row 3", "row 4"],
    footer="One reflection clears the whole column segment below the diagonal.",
    stem="householder_qr",
)
