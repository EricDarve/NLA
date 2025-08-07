3 main methods:
- CG: conjugate gradients; for [[Symmetric Positive Definite Matrices|symmetric positive definite matrices]]
- MINRES: for symmetric matrices
- GMRES: for general matrices

Compared to [[Classical iterative methods to solve sparse linear systems|classical iterative methods]] based on splitting, Krylov methods get much more accurate estimates by using the [[Krylov subspace|Krylov subspace:]]
$$
{\mathcal K}(A,b,k) = \text{span}\{ b, A b, A^2 b, \ldots, A^{k-1} b \}
$$
The goal of these methods is to find an “optimal” solution in the Krylov subspace.

We search for solutions of the form
$$
x^{(k)} = Q_k \, y_k
$$
where $Q_k$ is an [[Orthogonal matrix and projector|orthogonal]] basis of the Krylov subspace $K_k$. Vector $y_k$ is defined as
$$
y_k = {\rm argmin}_{y} \| x - Q_k y \|
$$
What norm should we use?

The naive choice is
$$
\| x - x^{(k)} \|_2 = \| x - Q_k y_k \|_2
$$
But solving the least-squares problem requires knowing $Q_k^T x$. This is not possible, unfortunately. So other ideas are required.