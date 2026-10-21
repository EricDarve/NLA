(lls-methods-summary)=
# Summary of Least-Squares Methods

For $A\in\mathbb{R}^{m\times n}$ and $b\in\mathbb{R}^m$, least squares minimizes

$$
\min_{x\in\mathbb{R}^n}\|Ax-b\|_2.
$$

A minimizer always exists, and the fitted vector $Ax$ is unique. The coefficient vector is unique if and only if $A$ has full column rank. Otherwise, the unique minimum-norm solution is the shortest coefficient vector among all minimizers.

All least-squares minimizers have a residual perpendicular to the column space:

$$
A^T(b-Ax)=0.
$$

The minimum-norm minimizer additionally satisfies $x\perp N(A)$. These two geometric conditions connect the methods below.

## Normal Equations

The [normal-equations method](normal_equations.md) forms $C=A^TA$ and $d=A^Tb$, then solves

$$
Cx=d
$$

using Cholesky factorization and two triangular solves. This solver requires full column rank, which makes $C$ positive definite. The normal equations themselves remain valid for rank-deficient matrices, but $C$ is then singular and the Cholesky procedure described here does not apply.

The accuracy risk comes from forming $A^TA$: for full column rank, $\kappa_2(A^TA)=\kappa_2(A)^2$. Rounding can obscure its small eigenvalues, so a stable Cholesky solve alone does not ensure an accurate least-squares answer. If the columns are nearly dependent, the rounded Gram matrix can even lose positive definiteness.

For dense $A$ with $m\ge n$, the leading cost is $mn^2+\frac13n^3$ flops when symmetry is used to form the Gram matrix. For $m\gg n$, this is about half the factorization arithmetic of Householder QR. The saving must be weighed against the accuracy risk.

## QR Factorization

For full column rank, the [QR method](LS_using_QR.md) uses the reduced factorization $A=Q_1R_1$ and solves

$$
R_1x=Q_1^Tb
$$

by back substitution. Householder QR stores $Q$ implicitly as reflections, so the transformed right-hand side can be computed without forming an explicit orthogonal matrix.

Householder QR gives small factorization and orthogonality errors under the arithmetic assumptions discussed earlier. It avoids forming $A^TA$, and its triangular factor satisfies $\kappa_2(R_1)=\kappa_2(A)$. This avoids squaring the condition number of the matrix used in the solve; it does not remove the sensitivity of the original least-squares problem.

The leading dense factorization cost is $2mn^2-\frac23n^3$ flops for $m\ge n$, plus $O(mn)$ work to transform one right-hand side and $O(n^2)$ for back substitution. The [QR summary](qr_summary.md) compares Householder, Givens, and Gram-Schmidt implementations and their different orthogonality guarantees.

QR exists even when columns are dependent. What fails is the ordinary back-substitution solve with a singular $R_1$. QR-based minimum-norm solvers use an additional orthogonal factorization to obtain a basis for the row space; the full-column-rank algorithm presented here does not include this extra step.

## Singular Value Decomposition

The [SVD method](LS_using_SVD.md) works for any matrix shape or rank. With the compact SVD $A=U_r\Sigma_rV_r^T$, where $r=\operatorname{rank}(A)$ and $\Sigma_r$ contains only positive singular values, the minimum-norm solution is

$$
x^*=V_r\Sigma_r^{-1}U_r^Tb=A^\dagger b.
$$

Compute $U_r^Tb$, divide by the positive singular values, and multiply by $V_r$. There is no need to form the pseudoinverse $A^\dagger$ explicitly. If $A=0$, the minimum-norm solution is $x^*=0$.

The SVD identifies the directions that affect the fit and sets the null-space component to zero. It also exposes small singular values, which can amplify errors when inverted. **SVD does not make an ill-conditioned problem insensitive to perturbations.** Treating small singular values as zero requires a tolerance; doing so solves a truncated approximation to the original matrix.

A standard dense SVD computing the singular vectors needed for the solve costs $O(mn\min(m,n))$. For $m\ge n$, this is $O(mn^2)$, the same order as the other methods, but generally with more arithmetic than Householder QR.

## Choosing a Method

For a dense problem with full column rank, Householder QR is a standard choice for reliable computation. Normal equations can be useful when the matrix is well-conditioned and the arithmetic saving matters. Use SVD when rank information or a minimum-norm solution for dependent columns or an underdetermined system is needed. The choice concerns how we compute the answer; the sensitivity of that answer still depends on the problem.
