# Orthogonal Iteration and Eigenvalues

(sec:oi-eigs)=
:::{admonition} Key Idea
:class: tip

Orthogonal (subspace) iteration produces a sequence of orthonormal bases $Q_k$. The **projected matrix**,

$$
T_k = Q_k^H A Q_k
$$

is the Rayleigh–Ritz compression of $A$ onto the subspace $\operatorname{range}(Q_k)$. As $Q_k$ converges to the Schur vectors of $A$, the matrix $T_k$ approaches the (upper triangular) Schur form $T$. Consequently, the diagonal entries of $T_k$ become increasingly accurate approximations of the eigenvalues of $A$.
:::

## From Subspaces to Eigenvalue Approximations

Recall that one step of orthogonal iteration (using the thin QR factorization) is:

$$
Z_{k+1} = A Q_k, \qquad Z_{k+1} = Q_{k+1} R_{k+1}.
$$

We define the Rayleigh–Ritz matrix as:

$$
T_k = Q_k^H A Q_k \;\in\; \mathbb{C}^{p\times p}.
$$

By construction, $\operatorname{range}(Q_k)$ provides an orthonormal basis for a $p$-dimensional subspace, and $T_k$ represents the matrix $A$ projected onto that basis.

## Convergence of $T_k$ to the Schur Form

Assume $A$ has the Schur decomposition:

$$
A = Q T Q^H.
$$

If orthogonal iteration is initialized with a matrix $Q_0$ of full column rank (whose columns have a nonzero component in the desired invariant subspace), then $\operatorname{range}(Q_k)$ will converge. Under standard spectral separation conditions, this subspace converges to the corresponding $p$-dimensional $A$-invariant subspace spanned by the appropriate columns of $Q$.

Consequently:

$$
Q_k \;\to\; Q[:,1:p] \quad\Longrightarrow\quad
T_k = Q_k^H A Q_k \;\to\; T[1:p,1:p],
$$

This means $T_k$ asymptotically becomes upper triangular, with the target eigenvalues on its diagonal.

## Shifts: Accelerating Eigenvalue Convergence

Shifting modifies the iteration to dramatically accelerate convergence by emphasizing specific eigenvalues.

Without shifting, the convergence of the $p$-dimensional subspace $\mathrm{range}(Q_k)$ to the dominant invariant subspace (spanned by the first $p$ Schur vectors) is governed by the eigenvalue separation. The convergence rate typically depends on the ratio of the largest "unwanted" eigenvalue ($\lambda_{p+1}$) to the smallest "wanted" one ($\lambda_p$):

$$
\frac{|\lambda_{p+1}|}{|\lambda_p|}
$$

If this ratio is close to 1 (i.e., the spectral gap is small), convergence can be impractically slow.

To accelerate this, we introduce a scalar **shift** $\mu$ and apply the iteration to the matrix $A - \mu I$ instead:

$$
Z_k = (A - \mu I)\,Q_k, \qquad Z_k = Q_{k+1}\,R_{k+1}.
$$

This works because $A$ and $A - \mu I$ share the same invariant subspaces (and Schur vectors), so $\mathrm{range}(Q_{k+1})$ still converges to the same desired subspace. However, the convergence rate is now determined by the eigenvalues of the *shifted* matrix:

$$
\frac{|\lambda_{p+1}(A - \mu I)|}{|\lambda_p(A - \mu I)|} = \frac{|\lambda_{p+1} - \mu|}{|\lambda_p - \mu|}
$$

```{admonition} Accelerating the Rate
:class: tip

This new ratio is the key to acceleration. Convergence becomes very fast if this ratio is much less than 1 ($\ll 1$). This can be achieved by choosing a shift $\mu$ such that:

$$
|\lambda_{p+1} - \mu| \ll |\lambda_p - \mu|
$$

By picking $\mu$ to be a good estimate of an "unwanted" eigenvalue (like $\lambda_{p+1}$), we can make the numerator extremely small, forcing rapid convergence.
```

The focus now shifts to the projected matrix $T_k = Q_k^H A Q_k$, as its eigenvalues, $\text{eig}(T_k)$, serve as our best available approximations to the true eigenvalues of $A$.

This enables a powerful feedback loop:

1.  Run the iteration to compute the basis $Q_k$.
2.  Compute the projected matrix $T_k = Q_k^H A Q_k$ and find its eigenvalues (known as "Ritz values").
3.  Use one of these Ritz values as the shift $\mu$ for the *next* iteration.

This strategy of adaptively choosing shifts based on the current Rayleigh-Ritz approximations is the core idea that transforms subspace iteration into the highly efficient **QR iteration with shifts**.