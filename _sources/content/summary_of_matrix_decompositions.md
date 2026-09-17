# Summary of Matrix Decompositions

The decompositions in this chapter describe a matrix in coordinates that reveal its structure. Eigenvalues describe **time evolution** along individual modes; singular values describe **size and scale in space**.

Here $H$ denotes conjugate transpose. A unitary matrix satisfies $Q^HQ=I$; for real matrices, the corresponding condition is orthogonality, $Q^TQ=I$.

## The Main Decompositions

| Decomposition | When it exists | Structure of the factors |
| :--- | :--- | :--- |
| **Complex Schur**<br>$A=QTQ^H$ | Every square complex matrix, including real matrices viewed over $\mathbb{C}$. | $Q$ is unitary; $T$ is upper triangular. Its diagonal entries are the eigenvalues. The columns of $Q$ need not be eigenvectors. |
| **Real Schur**<br>$A=QSQ^T$ | Every square real matrix. | $Q$ is real orthogonal; $S$ is real block upper triangular. Its $1\times1$ diagonal blocks contain real eigenvalues; its $2\times2$ diagonal blocks have nonreal conjugate eigenvalue pairs. |
| **Eigendecomposition**<br>$A=X\Lambda X^{-1}$ | Exactly when a square matrix has a basis of eigenvectors over the chosen field. | $X$ is invertible and contains the eigenvectors; $\Lambda$ is diagonal and contains their eigenvalues. The eigenvectors need not be orthogonal. |
| **Unitary diagonalization**<br>$A=Q\Lambda Q^H$ | Exactly when a complex square matrix is normal: $A^HA=AA^H$. | $Q$ contains an orthonormal eigenvector basis; $\Lambda$ contains the eigenvalues. For real symmetric matrices, $Q$ and $\Lambda$ can both be real. |
| **Singular value decomposition**<br>$A=U\Sigma V^H$ | Every $m\times n$ real or complex matrix. | In the full SVD, $U$ is $m\times m$ and $V$ is $n\times n$, both unitary. The $m\times n$ matrix $\Sigma$ is rectangular diagonal with nonnegative entries. For real $A$, the factors can be real and the formula is $A=U\Sigma V^T$. |

## How the Decompositions Relate

- **Schur form always exists; an eigenvector basis may not.** A diagonal Schur form is a unitary diagonalization and exists exactly for normal matrices. A general eigendecomposition allows a nonorthogonal basis, so it is not a special case of Schur form.
- **Real and complex diagonalization differ.** A real matrix may need complex eigenvectors. In particular, real orthogonal diagonalization is possible exactly for real symmetric matrices; real normality alone is not enough.
- **The SVD uses separate input and output bases.** It pairs directions through $Av_i=\sigma_i u_i$. Its positive singular values are the square roots of the positive eigenvalues of both $A^HA$ and $AA^H$; the numbers of zero eigenvalues can differ when $A$ is rectangular. For normal matrices, the singular values are the eigenvalue magnitudes, arranged in decreasing order.

## Key Consequences

**Time evolution.** If $A=X\Lambda X^{-1}$, then

$$
A^k=X\Lambda^kX^{-1}.
$$

Each mode is multiplied by $\lambda_i^k$. The geometry of the eigenvectors matters too: nonorthogonal modes can combine to produce large amplification even when each eigenvalue has magnitude less than one.

**Size, scale, and rank.** With $p=\min(m,n)$ and $\sigma_1\geq\cdots\geq\sigma_p\geq0$,

$$
\|A\|_2=\sigma_1,
\qquad
\|A\|_F^2=\sum_{i=1}^p\sigma_i^2.
$$

The rank is the number $r$ of positive singular values. The first $r$ left and right singular vectors span $R(A)$ and $R(A^H)$; the remaining vectors in the full SVD span $N(A^H)$ and $N(A)$, respectively.

**Low-rank approximation.** The compact SVD keeps all $r$ positive singular values and is exact. Keeping only the largest $k<r$ gives

$$
A_k=\sum_{i=1}^k\sigma_i u_i v_i^H,
$$

the best approximation among matrices of rank at most $k$ in both the spectral and Frobenius norms. Its errors are

$$
\|A-A_k\|_2=\sigma_{k+1},
\qquad
\|A-A_k\|_F^2=\sum_{i=k+1}^p\sigma_i^2.
$$

The preceding sections give the proofs and examples for [eigendecomposition and Schur forms](eigendecomposition.md), [normal matrices](normal_matrices.md), and the [SVD](singular_value_decomposition.md). Later chapters develop the numerical methods for computing and using these decompositions.
