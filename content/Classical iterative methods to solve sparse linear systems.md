For sparse matrices, $A = LU$ can be [[Motivation of iterative methods for eigenvalue computation|expensive]] when $A$ is large.

We will see more powerful methods later on, but a few simple strategies can be used when $A$ is sparse.

- [[Splitting methods]]
	- We outline the key idea of classical iterative methods using a splitting of matrix $A$ into $A=M-N$ where solving with $M$ is assumed to be computationally very cheap.
- [[Convergence of classical iterative methods]]
	- We establish a key condition for classical iterative methods to converge.
	- The condition is that $\rho(M^{-1}N) < 1$.