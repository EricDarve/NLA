# The Minimal Residual Method (MINRES)

GMRES minimizes the residual over a Krylov subspace for general nonsingular matrices, but its work and storage grow with the subspace. When $A$ is symmetric, we can perform the same minimization using short recurrences. This gives the **minimal residual method (MINRES)**, developed by Paige and Saunders.

Unlike CG, MINRES does not require positive definiteness. In particular, it applies to **symmetric indefinite** matrices, whose eigenvalues include both positive and negative values.

Throughout this section, $A\in\mathbb R^{n\times n}$ is symmetric and nonsingular. We use the same notation as in GMRES: $x$ is the exact solution, $x_0$ the initial guess, $r_0=b-Ax_0$, and $\beta=\|r_0\|_2$. If $r_0=0$, we are already done. Otherwise, the search space is $x_0+\mathcal K_k(A,r_0)$. The derivation assumes exact arithmetic.

## From Arnoldi to Lanczos

For a symmetric matrix, the Arnoldi process reduces to the [Lanczos process](lanczos_iteration.md). The projected matrix becomes symmetric tridiagonal, and generating each basis vector requires only the two preceding vectors.

While Lanczos generates a new direction, its relation is

$$
AQ_k=Q_{k+1}\underline T_k,
$$

where $Q_k=[q_1,\ldots,q_k]$ has orthonormal columns and $q_1=r_0/\beta$. The matrix $\underline T_k$ is $(k+1)\times k$: its first $k$ rows form a symmetric tridiagonal matrix $T_k$, and its last row can have a nonzero entry only in column $k$. Thus it is the tridiagonal specialization of GMRES's rectangular Hessenberg matrix $\underline H_k$.

## The MINRES Least-Squares Problem

As in GMRES, write a candidate as $x_0+Q_k y$. Since $r_0=Q_{k+1}(\beta e_1)$, where $e_1$ is the first coordinate vector in $\mathbb R^{k+1}$,

$$
\begin{aligned}
\|b-A(x_0+Q_k y)\|_2
&=\|r_0-AQ_k y\|_2\\
&=\|Q_{k+1}(\beta e_1-\underline T_k y)\|_2\\
&=\|\beta e_1-\underline T_k y\|_2.
\end{aligned}
$$

The last equality uses the orthonormal columns of $Q_{k+1}$. Therefore MINRES solves

$$
\boxed{
y_k=\underset{y\in\mathbb R^k}{\operatorname{argmin}}
\|\beta e_1-\underline T_k y\|_2,
\qquad x_k=x_0+Q_k y_k.
}
$$

**For a symmetric nonsingular matrix, MINRES and unrestarted GMRES give the same iterates in exact arithmetic when started from the same $x_0$.** They search the same space and minimize the same residual norm. MINRES uses symmetry to carry out this calculation more cheaply.

## Efficiency Gains

We solve the small problem with Givens rotations, just as in GMRES:

$$
G_k^T\cdots G_1^T\underline T_k
=\begin{pmatrix}R_k\\0\end{pmatrix},
\qquad
G_k^T\cdots G_1^T(\beta e_1)
=\begin{pmatrix}t_k\\\rho_k\end{pmatrix}.
$$

Here $R_k$ is upper triangular, $t_k\in\mathbb R^k$, and $\rho_k$ is a scalar. As in GMRES, nonsingularity of $A$ ensures that $R_k$ is nonsingular. We have $R_k y_k=t_k$ and $\|r_k\|_2=|\rho_k|$, where $r_k=b-Ax_k$.

Tridiagonal structure gives $R_k$ only two possible nonzero diagonals above its main diagonal. This allows both the factorization and the solution to be updated with short recurrences. To see why we need not store $Q_k$, define solution directions by

$$
V_k=Q_kR_k^{-1}=[v_1,\ldots,v_k],
\qquad Q_k=V_kR_k.
$$

This is a mathematical definition, not an instruction to form an inverse. Column $k$ of the second identity gives

$$
v_k=\frac{q_k-(R_k)_{k-1,k}v_{k-1}-(R_k)_{k-2,k}v_{k-2}}
{(R_k)_{kk}},
$$

with terms having indices below $1$ omitted. Only the last two solution directions are needed. The previous directions and the first $k-1$ entries of the rotated right-hand side remain unchanged, so

$$
x_k=x_0+V_k t_k
=x_{k-1}+(t_k)_k v_k.
$$

Thus we can update the approximation immediately rather than saving the entire Lanczos basis for a final combination.

Each standard MINRES step uses one matrix-vector product and a fixed number of vector operations. Its cost is $O(\operatorname{nnz}(A)+n)$, where $\operatorname{nnz}(A)$ counts the nonzero matrix entries, and its vector storage is $O(n)$, as for CG. This short-recurrence implementation is described on the [Stanford MINRES page](https://web.stanford.edu/group/SOL/software/minres/).

## Convergence

MINRES inherits the residual-minimizing property of GMRES. Since the search spaces grow,

$$
\|r_{k+1}\|_2\le\|r_k\|_2.
$$

The norm may stay unchanged at some steps; minimization does not guarantee rapid convergence. For a nonsingular system, exact termination occurs in at most $n$ steps, including earlier termination when the Lanczos space stops growing, as in GMRES.

Symmetry gives an orthonormal eigenvector basis, so the eigenvector-conditioning factor in the GMRES bound is one. If $A$ is also positive definite, the [Chebyshev bound](chebyshev_acceleration.md) gives

$$
\|r_k\|_2
\le2\|r_0\|_2
\left(\frac{\sqrt{\kappa}-1}{\sqrt{\kappa}+1}\right)^k,
\qquad
\kappa=\kappa_2(A)=\frac{\lambda_{\max}(A)}{\lambda_{\min}(A)}.
$$

For an indefinite matrix, the spectrum lies on both sides of zero, so this positive-interval bound does not apply. The polynomial must be small on both parts of the spectrum while remaining one at zero. MINRES still minimizes the residual, but its rate depends on the eigenvalue distribution.

In floating-point arithmetic, Lanczos vectors can lose orthogonality. The exact convergence statements need not hold literally, and an explicitly computed residual should be checked before accepting the answer. Reorthogonalization can help, but increases the work and storage beyond the basic short-recurrence method.

## CG vs. MINRES

For SPD $A$, both methods search the same Krylov space but use different objectives:

$$
\begin{aligned}
x_k^{\mathrm{CG}}
&=\underset{v\in x_0+\mathcal K_k}{\operatorname{argmin}}\|x-v\|_A,\\
x_k^{\mathrm{MINRES}}
&=\underset{v\in x_0+\mathcal K_k}{\operatorname{argmin}}\|b-Av\|_2.
\end{aligned}
$$

CG minimizes the energy error; MINRES minimizes the Euclidean residual. Their iterates generally differ. CG guarantees nonincreasing energy error, while MINRES guarantees nonincreasing Euclidean residual norms, in exact arithmetic.

CG is the usual choice for SPD systems. MINRES extends residual minimization to symmetric indefinite systems, where the energy norm used by CG is unavailable. Both retain a fixed amount of vector work and storage per iteration.

Even with a short recurrence, many iterations may be needed. The [next section](preconditioning.md) changes the system seen by the iterative method while preserving the solution of $Ax=b$.
