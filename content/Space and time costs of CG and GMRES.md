
Algorithm | Cost at iteration $k$ | Space total | Time total
--- | :-: | :-: | :-:
CG | $O(\text{nnz} + n)$ | $O(n)$ | $O(k \, (\text{nnz} + n))$
GMRES | $O(\text{nnz} + kn)$ | $O(kn)$ | $O(k \, \text{nnz} + k^2n)$


- GMRES struggles with large $k$. 
- The algorithm is often restarted when $k$ gets too large to reduce the space and time costs. 
- But restarting slows down the convergence and leads to more iterations.

[[Motivation of iterative methods for eigenvalue computation|Sparse-matrix vector products]], [[Krylov subspace]], [[Krylov methods for sparse systems]], [[Conjugate Gradients algorithm|CG]], [[GMRES algorithm|GMRES]]