# Orthogonal Iteration and Eigenvalues

:::{admonition} Key Idea
:class: tip

Orthogonal (subspace) iteration produces a sequence of orthonormal bases $Q_k$. The projected matrix,

$$
T_k = Q_k^H A Q_k
$$

is the Rayleigh–Ritz compression of $A$ onto the subspace $\operatorname{range}(Q_k)$. As this subspace converges to an invariant subspace of $A$, the matrix $T_k$ converges to the (upper triangular) Schur form of $A$ restricted to that subspace. Consequently, the eigenvalues of $T_k$ (and ultimately its diagonal entries) become increasingly accurate approximations of the corresponding eigenvalues of $A$.
:::

## From Subspaces to Eigenvalue Approximations

Recall that one step of orthogonal iteration (using the thin QR factorization) is:

$$
Z_{k+1} = A Q_k, \qquad Z_{k+1} = Q_{k+1} R_{k+1}.
$$

We define the Rayleigh–Ritz matrix as:

$$
T_k = Q_k^H A Q_k \in \mathbb{C}^{p\times p}.
$$

By construction, $\operatorname{range}(Q_k)$ provides an orthonormal basis for a $p$-dimensional subspace, and $T_k$ represents the matrix $A$ projected onto that basis.

## Convergence of $T_k$ to the Schur Form

Assume $A$ has the Schur decomposition $A = Q T Q^H$, where the eigenvalues on the diagonal of $T$ are ordered by magnitude, $|\lambda_1| > |\lambda_2| > \dots > |\lambda_n| > 0$.

If orthogonal iteration is initialized with a matrix $Q_0$ whose columns are linearly independent and have a nonzero component in the desired invariant subspace, then $\operatorname{range}(Q_k)$ will converge. Under standard spectral separation conditions, it converges to the dominant $p$-dimensional $A$-invariant subspace, which is spanned by the first $p$ Schur vectors.

Let $Q_p = Q[:,1:p]$ be the $n \times p$ matrix of the first $p$ Schur vectors, and let $T_p = T[1:p,1:p]$ be the corresponding $p \times p$ principal submatrix of $T$. As the subspaces converge, so do the projected matrices:

$$
\operatorname{range}(Q_k) \to \operatorname{range}(Q_p) \quad\Longrightarrow\quad T_k = Q_k^H A Q_k \to T_p
$$

This means $T_k$ asymptotically becomes upper triangular, and its eigenvalues (which are its diagonal entries) converge to the $p$ dominant eigenvalues of $A$, $\{\lambda_1, \dots, \lambda_p\}$.

## Shifts: Accelerating Eigenvalue Convergence

Shifting modifies the iteration to dramatically accelerate convergence by emphasizing specific eigenvalues.

### The Problem with Slow Convergence

Without shifting, the convergence of $\mathrm{range}(Q_k)$ to the dominant invariant subspace is governed by the eigenvalue separation. The convergence rate typically depends on the ratio of the largest "unwanted" eigenvalue ($\lambda_{p+1}$) to the smallest "wanted" one ($\lambda_p$):

$$
\frac{|\lambda_{p+1}|}{|\lambda_p|}
$$

If this ratio is close to 1 (i.e., the spectral gap is small), convergence can be impractically slow.

### The Shifted Iteration

To accelerate this, we introduce a scalar **shift** $\mu$ and apply the iteration to the matrix $A - \mu I$ instead:

$$
Z_k = (A - \mu I) Q_k, \qquad Z_k = Q_{k+1} R_{k+1}.
$$

This works because $A$ and $A - \mu I$ share the same invariant subspaces (and Schur vectors), so $\mathrm{range}(Q_{k+1})$ still converges to an invariant subspace. However, the convergence rate is now determined by the eigenvalues of the *shifted* matrix, which are $\lambda_i - \mu$.

### The Full QR Algorithm ($p=n$)

This strategy becomes exceptionally powerful when we set $p=n$, which transforms orthogonal iteration into the **QR algorithm**. In this special case, $Q_k$ and $T_k$ are square $n \times n$ matrices.

We can analyze the convergence by partitioning $T_k$:

$$
T_k = \begin{bmatrix}
T_{k,11} & T_{k,12} \\
T_{k,21} & T_{k,22}
\end{bmatrix}
$$

Here, $T_{k,22}$ is a scalar (the bottom-right entry) and $T_{k,21}$ is a row vector. For the $p=n$ case, the analysis shows that the sub-diagonal block $T_{k,21}$ converges to zero at a rate governed by the ratio of the smallest-magnitude eigenvalues:

$$
|T_{k,21}| \to 0 \quad \text{with rate} \quad \frac{|\lambda_n|}{|\lambda_{n-1}|}
$$

As $T_{k,21} \to 0$, $T_k$ becomes block upper triangular, and the scalar $T_{k,22}$ converges to the eigenvalue $\lambda_n$.

When we apply a shift $\mu$, the convergence rate of $T_{k,21}$ to zero becomes:

$$
\frac{|\lambda_n - \mu|}{|\lambda_{n-1} - \mu|}
$$

:::{admonition} Key Idea
:class: tip

This new ratio is the key to acceleration. To find the smallest eigenvalue $\lambda_n$, we can choose a shift $\mu$ that is a good approximation of $\lambda_n$. This makes the numerator $|\lambda_n - \mu|$ very small, forcing rapid convergence of $T_{k,21}$ to 0.
:::

The diagonal entries of $T_k = Q_k^H A Q_k$ provide our best available approximations to the eigenvalues. This enables a powerful feedback loop:

* Set $p=n$, so that $Q_k$ is a full $n \times n$ unitary matrix at each step
* At step $k$, compute the projected matrix $T_k = Q_k^H A Q_k$.
* Use the bottom-right entry of this matrix, $\mu_k = (T_k)_{nn}$, as the shift for the next iteration. This value, known as the Rayleigh quotient shift, is an excellent, adaptively improving approximation of $\lambda_n$.

This strategy of adaptively choosing a shift based on the last entry of $T_k$ is the core idea that transforms subspace iteration into the highly efficient QR iteration with shifts.