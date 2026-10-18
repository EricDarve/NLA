# Orthogonal Iteration and Eigenvalues

:::{admonition} Key Idea
:class: tip

Orthogonal (subspace) iteration produces a sequence of orthonormal bases $Q_k$. The projected matrix,

$$
T_k = Q_k^H A Q_k
$$

is the Rayleigh–Ritz compression of $A$ onto the subspace $\operatorname{range}(Q_k)$. As this subspace converges to the desired invariant subspace, the eigenvalues of $T_k$ approach the wanted eigenvalues of $A$. To obtain Schur directions and read eigenvalues from the diagonal, we need convergence of the nested leading column spaces, as established in [orthogonal iteration](orthogonal_iteration.md).
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

The columns of $Q_k$ form an orthonormal basis for a $p$-dimensional subspace, and $T_k$ represents the action of $A$ projected onto that subspace and expressed in this basis.

## Subspace Convergence and Schur Form

The preceding section distinguishes two conclusions. Under its diagonalizability, boundary-gap, and invertible starting-coordinate assumptions, the $p$-dimensional subspace converges to $\operatorname{range}(Q_{\star,p})$. This is enough for the eigenvalues of $T_k$ to approach $\lambda_1,\ldots,\lambda_p$.

The basis may still change within that subspace. Changing an orthonormal basis within a fixed invariant subspace changes its matrix representation by unitary similarity, which preserves eigenvalues. Thus convergence of the subspace does not require convergence of $T_k$ to a fixed triangular matrix.

With strict ordering $|\lambda_1|>\cdots>|\lambda_n|>0$ and the starting-coordinate condition for each leading set of columns, the nested leading subspaces converge. The individual columns then approach Schur directions, up to signs or complex phases. Consequently,

$$
(T_k)_{jj}\longrightarrow\lambda_j,
\qquad
(T_k)_{j\ell}\longrightarrow0\quad(j>\ell).
$$

Entries above the diagonal may continue to change with the phases of the basis vectors. For full iteration ($p=n$), $T_k$ is unitarily similar to $A$ at every step, so it already has exactly the same eigenvalues in exact arithmetic. The goal is to approach triangular form so those eigenvalues can be read from the diagonal.

## Shifts: Accelerating Eigenvalue Convergence

Shifting changes the eigenvalue ratios that control convergence. An effective choice can reduce the number of iterations.

### The Problem with Slow Convergence

Without shifting, the convergence of $\mathrm{range}(Q_k)$ to the dominant invariant subspace is governed by the eigenvalue separation. The convergence rate typically depends on the ratio of the largest "unwanted" eigenvalue ($\lambda_{p+1}$) to the smallest "wanted" one ($\lambda_p$):

$$
\frac{|\lambda_{p+1}|}{|\lambda_p|}
$$

If this ratio is close to 1 (i.e., the spectral gap is small), convergence can be impractically slow.

### The Shifted Iteration

To accelerate this, we introduce a scalar **shift** $\mu$ and apply the iteration to the matrix $A - \mu I$ instead:

$$
Z_{k+1} = (A - \mu I) Q_k,
\qquad Z_{k+1} = Q_{k+1} R_{k+1}.
$$

The matrices $A$ and $A-\mu I$ have the same invariant subspaces, and an exact Schur basis for $A$ is also a Schur basis for $A-\mu I$. The eigenvalues become $\lambda_i-\mu$. Their magnitudes can have a different ordering: power or partial subspace iteration now favors eigenvalues farthest from $\mu$. A shift therefore changes both the convergence ratios and potentially which subspace is selected.

### The Full QR Algorithm ($p=n$)

For $p=n$, we can express the basis updates directly as updates of $T_k$, giving the **QR algorithm**. The derivation is developed in [QR iteration](qr_iteration.md). In this case, $Q_k$ and $T_k$ are square $n\times n$ matrices.

For $n\ge2$, partition after the first $n-1$ rows and columns:

$$
T_k = \begin{bmatrix}
T_{k,11} & T_{k,12} \\
T_{k,21} & T_{k,22}
\end{bmatrix}
$$

Here $T_{k,22}$ is the bottom-right scalar and $T_{k,21}$ is a row vector. Under the strict ordering and starting-space conditions above, the unshifted iteration satisfies

$$
\|T_{k,21}\|_2
=O\!\left(\left|\frac{\lambda_n}{\lambda_{n-1}}\right|^k\right).
$$

As $T_{k,21} \to 0$, $T_k$ becomes block upper triangular, and the scalar $T_{k,22}$ converges to the eigenvalue $\lambda_n$.

For a fixed shift $\mu$ that is not an eigenvalue, suppose $\lambda_n$ is the unique closest eigenvalue to $\mu$ and the corresponding starting-space condition holds. The shifted boundary ratio is

$$
\rho_\mu=
\frac{|\lambda_n-\mu|}{\min_{j<n}|\lambda_j-\mu|}<1.
$$

:::{admonition} Key Idea
:class: tip

Choosing $\mu$ close to an isolated eigenvalue $\lambda_n$ makes $\rho_\mu$ small. In full iteration, the other $n-1$ directions separate from this last direction more quickly, so the bottom row approaches block upper triangular form faster. For a fixed shift the bound remains geometric. Updating the shift can give faster local convergence; the [section on shifts](qr_iteration_with_shifts.md) develops those results.
:::

Once $T_k$ is close to Schur form, its diagonal entries provide eigenvalue approximations that can guide the next shift:

* Set $p=n$, so that $Q_k$ is a full $n \times n$ unitary matrix at each step
* At step $k$, compute the projected matrix $T_k = Q_k^H A Q_k$.
* A common choice is the bottom-right entry, $\mu_k=(T_k)_{nn}$. This is the Rayleigh quotient of the last basis vector and is called the **Rayleigh quotient shift**.

Later sections explain how shift choices accelerate QR iteration and how to handle real matrices with complex conjugate eigenvalues.
