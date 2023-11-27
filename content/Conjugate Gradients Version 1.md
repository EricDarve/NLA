Let's use a different norm compared to the [[Krylov methods for sparse systems|previous derivation.]] This will only apply to matrices $A$ that are [[Symmetric Positive Definite Matrices|symmetric positive definite]].

We will consider the $A$-norm:
$$
\begin{gather}
\| x - x^{(k)} \|_A \\[.5em]
\| z \|_A = \sqrt{z^T A z}
\end{gather}
$$
If $A$ is [[Symmetric Positive Definite Matrices|symmetric positive definite,]] then this is a [[Operator and matrix norms|norm.]]

If we look for a solution that [[Krylov methods for sparse systems|minimizes that norm]] in the [[Krylov subspace|Krylov subspace]] we obtain a simple equation:
$$
\| x - x^{(k)} \|_A^2 = (x - Q_k y)^T A (x - Q_k y)
$$
The solution minimizes 
$$
y^T Q_k^T A Q_k y  - 2 b^T Q_k y + b^T x
$$
We can prove that the solution is
$$
(Q_k^T A Q_k) y_k = Q_k^T b
$$

We recognize $Q_k^T A Q_k$ from [[Lanczos process|Lanczos.]] This is the same matrix. We project $A$ using $Q_k$ and solve the resulting $k \times k$ system.

**Proof.** The function to minimize is
$$
\begin{align}
\| x - x^{(k)} \|_A^2 & = (x - Q_k y)^T A (x - Q_k y) \\[.2em]
& = (A^{1/2} x - A^{1/2} Q_k y)^T (A^{1/2} x - A^{1/2} Q_k y) \\[.2em]
& = \| A^{1/2} x - A^{1/2} Q_k y \|_2^2
\end{align}
$$
Using our previous reasoning for [[Method of normal equation|least-squares]], we find that
$$
A^{1/2} x - A^{1/2} Q_k y_k  \; \perp A^{1/2} Q_k
$$

$$Q_k^T A^{1/2} ( A^{1/2} x - A^{1/2} Q_k y_k) = 0$$

$$Q_k^T A x - Q^T A Q_k y_k = 0$$

$$(Q_k^T A Q_k) y_k = Q_k^T b$$
$\square$

The final equation is:
$$
(Q_k^T A Q_k) y_k = Q_k^T b
$$
Because we are using the $\| \; \|_A$ norm, we are able to get rid of $x$ and have $b$ instead.

This leads to:

### Conjugate gradients: version 1

1. Compute $Q_k$ using Lanczos.
2. $(Q_k^T A Q_k) \; y_k = Q_k^T b = \|b\|_2 \, e_1.$
3. $x^{(k)} = Q_k \, y_k.$

- Space cost: $O(kn).$ This is required because of the last step $Q_k y_k$ which requires saving $Q_k$.
- Time cost: $O(k \, {\rm nnz} + kn).$

This process can be made much more efficient. In particular, the space cost can be reduced to $O(n)$. The time cost can be reduced as well but will remain $O(k \, {\rm nnz} + kn).$

But it takes some insight to make it happen.