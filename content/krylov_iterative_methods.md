# Krylov Subspace Methods

## Searching for a Better Approximation

Jacobi, Gauss-Seidel, and SOR improve an approximate solution through inexpensive corrections. Their convergence can nevertheless be slow. Krylov subspace methods use information gathered over successive iterations to choose a correction from an expanding space of possibilities.

We consider a real nonsingular system $Ax=b$. Let $x$ denote the exact solution, $x^{(0)}$ an initial guess, and $r^{(0)}=b-Ax^{(0)}$ its residual. The residual identifies which equations are not yet satisfied. Applying $A$ repeatedly to this vector provides additional search directions, giving the **Krylov subspace**

$$
\mathcal K_k(A,r^{(0)})
=\operatorname{span}\{r^{(0)},Ar^{(0)},\ldots,A^{k-1}r^{(0)}\}.
$$

These are the same spaces used in the earlier [Arnoldi](arnoldi_iteration.md) and [Lanczos](lanczos_iteration.md) methods for eigenvalue approximation. Here we use them to construct a correction to a linear-system solution:

$$
\boxed{x^{(k)}=x^{(0)}+z^{(k)},
\qquad z^{(k)}\in\mathcal K_k(A,r^{(0)}).}
$$

Thus each approximation can combine all the directions available so far. If $x^{(0)}=0$, then $r^{(0)}=b$, and we simply search for $x^{(k)}$ in $\mathcal K_k(A,b)$. This is the starting choice used in the next section.

Building these spaces requires matrix-vector products, so we can work with a sparse matrix or a procedure that computes $Av$ without storing $A$. In practice, we construct a suitable basis rather than forming the powers $A^j$ explicitly.

## The Central Challenge: What Is “Optimal”?

Once the search space is chosen, we still need to select an approximation within it. The most direct goal would be to minimize the Euclidean error $\|x-x^{(k)}\|_2$.

If we had an orthonormal basis for the space, the best approximation would be an orthogonal projection of the exact solution onto the search space. But computing that projection requires information about $x$, the solution we are trying to find. The subspace is available; this error-minimizing choice is not directly computable from it.

We therefore need a criterion that uses the known data $A$ and $b$. Two choices lead to CG and GMRES.

## Two Computable Paths to Optimality

**Conjugate gradient (CG)** applies when $A$ is symmetric positive definite (SPD). It chooses the approximation that minimizes the error in the **energy norm**, $\|v\|_A=\sqrt{v^TAv}$, over the current search space.

Although the exact error is unknown, this minimization is computable. As in the [energy interpretation of Gauss-Seidel](gauss_seidel_method.md), minimizing the squared energy error is equivalent to minimizing the quadratic $\Phi(v)=\tfrac12v^TAv-b^Tv$. We can evaluate and minimize this quadratic using $A$ and $b$. CG performs that minimization over the growing Krylov space, using a short recurrence that stores only a few vectors.

**Generalized minimal residual (GMRES)** applies to general nonsingular matrices. It chooses the approximation that minimizes $\|b-Ax^{(k)}\|_2$ over the current search space. The residual is computable, and the minimization reduces to a smaller linear least-squares problem. Arnoldi provides the orthonormal basis and the small matrices needed for this calculation.

In exact arithmetic, both methods choose the best approximation according to their own criterion. CG minimizes the energy norm of the error; GMRES minimizes the Euclidean norm of the residual. These are different objectives and need not select the same approximation.

## The Main Practical Challenges

An expanding space offers more choices, but we must manage its cost. Standard GMRES stores the Krylov basis, and orthogonalizing each new vector requires more work as the basis grows. Restarting limits storage and work, but can slow convergence. CG exploits symmetry and positive definiteness to avoid storing the full basis.

Optimality within a search space also does not guarantee rapid convergence. The space may need many directions before it contains an accurate approximation. **Preconditioning** addresses this difficulty by changing the system into an equivalent form that is easier for the iteration to solve. The preconditioner must improve convergence while remaining inexpensive to apply.

Finally, we must decide when to stop. The residual is available, while the exact error usually is not. As in the earlier analysis of linear solves, a small residual must be interpreted together with the sensitivity of the problem. Rounding errors can also affect the orthogonality and recurrence relations used by the methods.

## The Path through This Chapter

We begin with [CG](conjugate_gradient.md) for SPD matrices, then develop [GMRES](gmres.md) for general systems. The [comparison of CG and GMRES](cg_vs_gmres.md) brings together their optimality properties, costs, and uses.

Next, [MINRES](minres.md) shows how to minimize the residual while retaining a short recurrence for symmetric matrices, including indefinite ones. The final sections develop [preconditioning](preconditioning.md) and [preconditioned CG](pcg.md), connecting the choice of an inexpensive approximate solve to the convergence of the Krylov method.
