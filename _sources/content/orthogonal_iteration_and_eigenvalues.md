# Orthogonal Iteration and Eigenvalues

The [previous section](orthogonal_iteration.md) explains how orthogonal iteration finds invariant subspaces and Schur directions. To compute eigenvalues from these results, we study the matrix $T_k=Q_k^HAQ_k$. Following this matrix directly will lead to the QR iteration.

We keep the same notation: $Q_\star$ is an exact Schur basis in $A=Q_\star TQ_\star^H$, with $|\lambda_1|\ge\cdots\ge|\lambda_n|$, and $Q_k\in\mathbb C^{n\times p}$ is the orthonormal basis at iteration $k$. The first $i$ columns are denoted by $Q_{\star,i}$ and $Q_{k,i}$.

## Extracting Eigenvalues from a Subspace

For $p<n$, orthogonal iteration gives a small subspace in which to look for eigenvectors. A vector in this subspace has the form $\boldsymbol u=Q_k\boldsymbol y$. Multiplication by $A$ can take it outside the subspace. Projecting back gives

$$
Q_kQ_k^HA\boldsymbol u
=Q_k\underbrace{(Q_k^HAQ_k)}_{T_k}\boldsymbol y.
$$

Thus the $p\times p$ matrix

$$
\boxed{T_k=Q_k^HAQ_k}
$$

describes the projected action of $A$ in the current basis. For $p=1$, this is the Rayleigh quotient used in the power method. For general $p$, its eigenvalues are called **Ritz values**. To see why they approximate eigenvalues of $A$, let

$$
T_k\boldsymbol y=\theta\boldsymbol y,
\qquad \|\boldsymbol y\|_2=1,
\qquad \boldsymbol u=Q_k\boldsymbol y.
$$

Using the subspace residual $E_k=AQ_k-Q_kT_k$ from the previous section,

$$
\begin{aligned}
A\boldsymbol u-\theta\boldsymbol u
&=AQ_k\boldsymbol y-Q_kT_k\boldsymbol y\\
&=E_k\boldsymbol y,\\
\|A\boldsymbol u-\theta\boldsymbol u\|_2
&\le\|E_k\|_2.
\end{aligned}
$$

The vector $\boldsymbol u$ is called a **Ritz vector**. It has unit norm, and its eigenvalue residual is small whenever $E_k$ is small. The residual is perpendicular to the current subspace, because

$$
Q_k^HE_k=Q_k^HAQ_k-T_k=0,
\qquad
Q_k^H(A\boldsymbol u-\theta\boldsymbol u)=0.
$$

Thus the Ritz pair satisfies the eigenvector equation after projection onto the subspace. If the subspace is exactly invariant, then $E_k=0$ and $(\theta,\boldsymbol u)$ is an exact eigenpair of $A$. For a nonnormal matrix, the eigenvalue error can be much larger than the residual.

## What Subspace Convergence Guarantees

Suppose the assumptions of the preceding convergence theorem hold: $A=X\Lambda X^{-1}$ is diagonalizable, $|\lambda_p|>|\lambda_{p+1}|$, and the first $p$ rows of $X^{-1}Q_0$ form an invertible matrix. The columns of $X$ must follow the eigenvalue ordering above. Then

$$
\operatorname{range}(Q_k)
\longrightarrow\operatorname{range}(Q_{\star,p}).
$$

On this exact invariant subspace, define

$$
T_{\star,p}=Q_{\star,p}^HAQ_{\star,p}.
$$

This is the leading $p\times p$ block of the exact Schur matrix $T$, with eigenvalues $\lambda_1,\ldots,\lambda_p$. The current basis may keep changing within the converging subspace. We can account for these changes by choosing unitary $p\times p$ matrices $W_k$ that align the two bases:

$$
\begin{aligned}
\|Q_k-Q_{\star,p}W_k\|_2&\longrightarrow0,\\
\|W_kT_kW_k^H-T_{\star,p}\|_2&\longrightarrow0.
\end{aligned}
$$

Unitary changes of basis preserve eigenvalues, and eigenvalues depend continuously on matrix entries. Therefore the Ritz values approach $\lambda_1,\ldots,\lambda_p$, up to reordering and counting multiplicities, even if the basis keeps rotating. **Convergence of the subspace is enough to recover its eigenvalues.** It does not by itself make $T_k$ triangular or make its diagonal entries eigenvalue approximations.

## When Eigenvalues Appear on the Diagonal

For full iteration, $p=n$, $Q_k$ is square and unitary. Consequently, $T_k=Q_k^HAQ_k$ already has exactly the same eigenvalues as $A$ at every step in exact arithmetic. Also, $E_k=0$ at every step, so the full-space residual cannot measure progress. The computational goal is to make $T_k$ triangular so that those eigenvalues can be read from its diagonal.

The previous section proves this under the stronger assumptions

$$
|\lambda_1|>\cdots>|\lambda_n|>0,
$$

with the starting-coordinate condition for every leading set of columns: the first $i$ rows of $X^{-1}Q_{0,i}$ form an invertible $i\times i$ matrix. Each leading column space then converges, so each column approaches its Schur direction up to a sign or complex phase.

Let $D_k$ be a diagonal matrix of these unit-modulus phases. Equivalently,

$$
\begin{aligned}
\|Q_k-Q_\star D_k\|_2&\longrightarrow0,\\
\|T_k-D_k^HTD_k\|_2&\longrightarrow0.
\end{aligned}
$$

Because $T$ is upper triangular, $D_k^HTD_k$ is upper triangular with the same diagonal. Hence

$$
\boxed{(T_k)_{jj}\longrightarrow\lambda_j,
\qquad (T_k)_{j\ell}\longrightarrow0\quad(j>\ell).}
$$

Entries above the diagonal may keep changing with the phases. These phase changes do not affect the limiting diagonal.

:::{admonition} Key Result
:class: tip

For partial iteration, convergence of the $p$-dimensional subspace makes the eigenvalues of $T_k$ converge to the wanted eigenvalues of $A$.

For full iteration, the eigenvalues of $T_k$ are already exact. Convergence of the nested leading subspaces makes its entries below the diagonal tend to zero, exposing those eigenvalues on the diagonal.
:::

## From Orthogonal Iteration to QR Iteration

For computing all eigenvalues, the goal can be expressed entirely in terms of $T_k$: preserve the eigenvalues while driving the entries below the diagonal to zero. This suggests computing $T_{k+1}$ directly from $T_k$, without explicitly updating the full basis $Q_k$ and then forming $Q_k^HAQ_k$ at every step.

As $T_k$ approaches triangular form, its diagonal entries provide eigenvalue approximations that can guide acceleration. A **shift** replaces $A$ by $A-\mu I$, preserving its invariant subspaces while changing the eigenvalue ratios that control convergence. A suitable shift can speed up convergence. Shift selection and convergence are developed in the later section on [shifted QR iteration](qr_iteration_with_shifts.md).

The immediate task is to derive the direct update of $T_k$. The [next section](qr_iteration.md) shows how a QR factorization of $T_k$ gives $T_{k+1}$, producing the **QR iteration**.
