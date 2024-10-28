[[Orthogonal iteration|We saw]] that $Q_k$ converges to $Q$, and $T_k = Q_k^H A Q_k$ converges to an upper triangular matrix.

The rate of convergence is not straightforward to derive.

Here is a sketch of a proof. Start with
$$
A^k = X \Lambda^k X^{-1}.
$$
Span of $Q_k[1:i]$ converges to span of $X[1:i]$ = span of $Q[1:i].$ The convergence rate is given by
$$
\Big| \frac{\lambda_{i+1}}{\lambda_i} \Big|
$$
**Proposition 1: Convergence of [[Orthogonal iteration|orthogonal iteration]]**

The norm of the block $T_k[i+1:n,1:i]$ decays like 
$$
\Big| \frac{\lambda_{i+1}}{\lambda_i} \Big|^k.
$$
$\square$

**Proposition 2: Convergence of the Ritz eigenvalues**

Assume that we start with a random $Q_0$ with $r$ columns. In that case, $T_k$ has dimension $r \times r$. The eigenvalues of $T_k$ are called Ritz eigenvalues. The $i$th Ritz eigenvalue converges to $\lambda_i$ with rate
$$
\Big| \frac{\lambda_{r+1}}{\lambda_i} \Big|^k
$$
$\square$

Let's go back to the case of $Q_0 = I$ and $T_k$ of size $n \times n$. Consider a case where $\lambda_j \gg \lambda_{j+1}$. Then we have the following block structure for $T_k$.

$$T_k =
\begin{pmatrix}
* & * \\
\epsilon & *
\end{pmatrix}
$$
We converge quickly to a $2 \times 2$ block upper triangular matrix.

More generally if we have sufficient separation between eigenvalues, i.e., $|\lambda_i| \gg |\lambda_{i+1}|,$ then convergence to an [[Schur decomposition|upper triangular matrix]] is very fast.