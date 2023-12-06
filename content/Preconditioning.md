The convergence of both [[Conjugate Gradients algorithm|CG]] and [[GMRES algorithm|GMRES]] depends on the [[Conditioning of a linear system|condition number]] of the matrix and the distribution of the eigenvalues. Convergence improves if the matrix is “closer” to the identity or has clustered eigenvalues and its eigenbasis is well-conditioned.

[[Convergence of the Conjugate Gradients]], [[Convergence of GMRES]]

### Goal of preconditioner

Preconditioners should be easy to apply and should improve the convergence of iterative methods.

The eigenvalues of the preconditioned system should be clustered and lead to fewer iterations in CG/GMRES/MINRES.

## 3 types of preconditioners

There are three main ways to precondition a linear system $Ax=b$.
- Left preconditioning: $(M_1 A) x = M_1 b$
- Right preconditioning: $(A M_2) z = b,$ $x = M_2 z$
- Symmetric preconditioning: $(M_1 A M_2) z = M_1 b,$ $x = M_2 z$

Symmetric preconditioning is required for [[Conjugate Gradients algorithm|CG]] because we need to maintain the [[Symmetric Positive Definite Matrices|SPD]] property. This is the case for [[MINRES]] as well. More specifically, we need to use
$$
M_1 A M_2
$$
with $M_1 = C$ and $M_2 = C^T$. The preconditioned system must be constructed in a symmetric way as
$$
C A C^T.
$$

## Examples of preconditioners

- (Block) diagonal matrix = [[Jacobi iteration|Jacobi]]; $M^{-1} = D$
- [[Gauss-Seidel iteration|Gauss-Seidel]]: $M^{-1} = D - L$
- (Symmetric) [[SOR iteration|Successive Over-Relaxation]]
- Incomplete LU and Cholesky
- Multigrid
- Fast solvers with low-rank compression steps