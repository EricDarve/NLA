We now use [[Key idea of iterative methods for eigenvalue computation|the algorithm]] to compute $Q_k$ and $H_k$ to approximate the eigenvectors and eigenvalues of $A$.

We approximate the eigenvalues of $A$ using the eigenvalues of $H_k$:
$$
\lambda(A) \approx \lambda(H_k).
$$
We briefly outline why this is a reasonable idea. Consider the eigenvectors and eigenvalues of $H_k$:
$$
H_k = X_k \Lambda_k X_k^{-1}
$$
Use the [[Key idea of iterative methods for eigenvalue computation#^d63511|previous equation]] with $Q_k$ and $h_{k+1,k}$:
$$
Q_k (X_k \Lambda_k X_k^{-1}) + h_{k+1,k} \, q_{k+1} e_k^T = A Q_k
$$
Multiply by $X_k$ to the right:
$$
A (Q_k X_k) = (Q_k X_k) \Lambda_k + h_{k+1,k} \, q_{k+1} \, x_{k,}
$$
We see that we have approximate eigenvectors and eigenvalues of $A$ assuming that $h_{k+1,k}$ is small:
$$
A (Q_k X_k) \approx (Q_k X_k) \Lambda_k
$$

- $\Lambda_k$ are called the [[Convergence of the orthogonal iteration|Ritz eigenvalues.]]
- $Q_k X_k$ are the approximate eigenvectors.

**Arnoldi algorithm:**

- Step 1: compute $H_k$ using the iteration starting from $q_1$.
- Step 2: Compute the eigenvalue of $H_k$.

$A q_j$ and $q_k^T A q_j$ only required!

Low flop count if $A$ is sparse.