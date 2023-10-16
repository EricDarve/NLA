- The previous sections have shown that the QR factorization exists.
- Let's assume that $A$ is full column rank, then we will prove that its QR factorization is unique if we require that $r_{ii} > 0$.

Proof: consider $A^T A$. Since $A$ is full column rank, this matrix is [[Symmetric Positive Definite Matrices|SPD]]. From $A=QR$, we get:
$$
A^T A = R^T Q^T Q R = R^T R
$$
Denote by $L = R^T$. $L$ is a lower triangular matrix with positive entries on the diagonal. We have $A^T A = LL^T$. So $LL^T$ is the Cholesky factorization of $A$. This factorization is [[Existence of the Cholesky factorization|unique]]. So the factor $R$ is unique. But we also have:
$$
Q = AR^{-1}
$$
so the factor $Q$ is also unique. $\square$

[[Symmetric Positive Definite Matrices]], [[Existence of the Cholesky factorization]], [[QR factorization]], [[QR using Householder transformations]], [[QR using Givens transformations]], [[Gram-Schmidt]]