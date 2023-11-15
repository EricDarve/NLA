For sparse matrices, $A = LU$ can be [[Motivation of iterative methods for eigenvalue computation|expensive]] when $A$ is large.

We will see more powerful methods later on, but a few simple strategies can be used when $A$ is sparse.

- [[Splitting methods]]
	- We outline the key idea of classical iterative methods using a splitting of matrix $A$ into $A=M-N$ where solving with $M$ is assumed to be computationally very cheap.
- [[Convergence of classical iterative methods]]
	- We establish a key condition for classical iterative methods to converge.
	- The condition is that $\rho(M^{-1}N) < 1$.
- [[Jacobi iteration]]
	- A very simple technique.
	- Very easy to implement on parallel computers and multicore processors.
	- Does not always converge fast.
	- This algorithm converges for matrices that are strictly row or column diagonally dominant.
- [[Gauss-Seidel iteration]]
	- This method has better convergence compared to Jacobi.
	- However, it is more difficult to parallelize and may have lower performance, that is lower [Gflops](https://en.wikipedia.org/wiki/FLOPS) on modern multicore processors.
	- Gauss-Seidel converges for matrices that are strictly row or column diagonally dominant, and for all [[Symmetric Positive Definite Matrices|symmetric positive definite]] matrices.
- [[SOR iteration]]
	- In some cases, we may improve the convergence of Gauss-Seidel by using a relaxation parameter $\omega$.
- [[SOR iteration as a splitting method]]
	- We further analyze the SOR method.
	- We show that it is a type of splitting method.
	- Convergence requires that $0 < \omega < 2$.
	- SOR also converges for all [[Symmetric Positive Definite Matrices|symmetric positive definite]] matrices.
- [[Chebyshev iteration]]
	- This iterative method is very fast but it requires to estimate an interval that contains the eigenvalues.
	- Its convergence rate can be near-optimal in some cases.