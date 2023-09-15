Not all matrices are diagonalizable. But all square matrices have a Schur decomposition:
$$
A = QTQ^H
$$
- $H$: transpose conjugate
- $T$: triangular with **eigenvalues on the diagonal**
- $Q$: complex unitary; $Q^H Q =I$.

Schur decompositions can be accurately computed because they rely on unitary matrices. This is one of the best decompositions to compute eigenvalues.

[[Eigenvalues]], [[Orthogonal matrix and projector]]