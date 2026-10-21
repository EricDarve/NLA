"""Show two-sided Householder reduction in the existing matrix-diagram style."""

import numpy as np
from givens_matrix_diagrams import plot_sequence


initial = np.array([
    [4., 2., 1., 3.],
    [3., 5., 2., 1.],
    [1., 2., 6., 2.],
    [2., 1., 4., 7.],
])
current = initial.copy()
transform = np.eye(4)
states = [current.copy()]
active_masks = []
targets = []

for k in range(2):
    x = current[k + 1:, k].copy()
    sigma = -np.copysign(np.linalg.norm(x), x[0])
    v = x.copy()
    v[0] -= sigma
    reflection = np.eye(len(x)) - 2 * np.outer(v, v) / (v @ v)
    embedded = np.eye(4)
    embedded[k + 1:, k + 1:] = reflection
    assert np.allclose(embedded.T @ embedded, np.eye(4))

    left_update = embedded @ current
    assert np.allclose(left_update[k + 2:, k], 0.)
    assert np.allclose(np.tril(left_update[:, :k], -2), 0.)
    current = left_update @ embedded
    # The right update fixes every column already finished, including column k.
    assert np.allclose(current[:, :k + 1], left_update[:, :k + 1])
    current[k + 2:, k] = 0.
    transform = transform @ embedded
    assert np.allclose(current, transform.T @ initial @ transform)
    assert np.allclose(np.tril(current[:, :k + 1], -2), 0.)
    assert np.isclose(np.linalg.norm(current), np.linalg.norm(initial))
    states.append(current.copy())

    # Left updates act on the current column and trailing block; right updates
    # also change the top rows of the later columns. Their union is L-shaped.
    active = np.zeros((4, 4), dtype=bool)
    active[:, k + 1:] = True
    active[k + 1:, k] = True
    active_masks.append(active)
    targets.append([(row, k) for row in range(k + 2, 4)])

assert np.allclose(np.tril(current, -2), 0.)
assert np.all(np.abs(np.diag(current, -1)) > 1e-8)
assert np.array_equal(np.isclose(current, 0.), np.tril(np.ones((4, 4)), -2).astype(bool))

plot_sequence(
    states,
    active_rows=[range(1, 4), range(2, 4)],
    active_columns=[range(4), range(1, 4)],
    active_masks=active_masks,
    eliminated_groups=targets,
    titles=["Start", "First column finished", "Upper Hessenberg"],
    actions=[
        "Reflect rows 2–4,\nthen columns 2–4",
        "Reflect rows 3–4,\nthen columns 3–4",
        "Keep the first\nsubdiagonal",
    ],
    row_labels=["row 1", "row 2", "row 3", "row 4"],
    footer="Reflect from both sides; leave the first subdiagonal.",
    stem="householder_hessenberg",
)
