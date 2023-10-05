Not all matrices are diagonalizable. But all square matrices have a Schur decomposition:
$$
A = QTQ^H
$$
- $T$: triangular with **eigenvalues on the diagonal**
- $Q$: square complex unitary; $Q^H Q =I$.

Schur decompositions can be accurately computed because they rely on unitary matrices. This is one of the best decompositions to compute eigenvalues.

**Proof of the existence of the Schur decomposition**

We prove the result by induction on the size of the matrix. For $n=1$, the result is true.

Assume it is true for matrices of size less than $n$.

We know that all matrices have at least one eigenvalue. Denote $\lambda$ and $x$ the eigenvalue and eigenvector with $\|x\|_2 = 1$. We can define a unitary basis
$$
Q_1 = \big[ x, x_1, \dots, x_{n-1} \big]
$$
We have
$$
Q_1^H A Q_1 = \begin{pmatrix}
\lambda & T_{12} \\
0 & T_{22}
\end{pmatrix}
$$
By induction, we can find a unitary matrix $Q_2$ of size $n-1$ such that
$$
Q_2^H T_{22} Q_2 = T^*_{22}
$$
where $T^*_{22}$ is upper triangular. Define the unitary matrix
$$
Q_3 = Q_1 \begin{pmatrix}
1 & \\
 & Q_2
\end{pmatrix}
$$
Then
$$
Q_3^H A Q_3 =
\begin{pmatrix}
\lambda & * \\
0 & T^*_{22}
\end{pmatrix}
$$
is upper triangular. $\square$

[[Eigenvalues]], [[Orthogonal matrix and projector]], [[Hermitian and symmetric matrices]]