For any square matrix $A$ there exists at least a scalar $\lambda \in \mathbb C$ and $x \in \mathbb C^n$ such that $Ax = \lambda x$. They are called an eigenvalue and eigenvector.

A matrix $A$ is diagonalizable if there exists a basis $x_1$, ..., $x_n$ of eigenvectors. In that case, we can write:
$$
A = X \Lambda X^{-1}
$$
where column $i$ of $X$ is $x_i$ and $\Lambda$ is a diagonal matrix with $\lambda_i$ on the diagonal:
$$
A x_i = \lambda_i \, x_i
$$
A key result is that $A^k$ has a very simple expression:
$$
A^k = X \Lambda^k X^{-1}
$$
where $\Lambda^k$ is diagonal with $\lambda_i^k$ on the diagonal.

[[Matrix-vector and matrix-matrix product]], [[Invertible matrix]]