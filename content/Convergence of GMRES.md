[[GMRES|Recall]] that we are minimizing the following norm
$$
\|r^{(k)}\|_2 = \| b - A Q_k y \|_2 
= \| (I - A q(A)) b \|_2
$$

where $q$ is a polynomial of degree $k-1$. To derive an error bound, we should search for a [[Connection between Arnoldi and polynomials of A|polynomial]] $p(x)$ of degree $k$ such that $p(0) = 1$ and which [[Connection between Arnoldi and polynomials of A|minimizes]] $\| p(A) b \|_2$.

**Theorem.** Suppose that $A$ has an eigenvalue decomposition $A = X \Lambda X^{-1}$. Then the GMRES error can be upper-bounded by
$$
  \|r^{(k)}\|_2 \; \le \;
  \|b\|_2 \; 
  \kappa(X)
  \min_{ \substack{p \text{ of degree $k$ }\\p(0) = 1} }
  \max_{ \substack{\lambda \text{ eigenvalue}\\ \text{of $A$}} } |p(\lambda)|
$$
where $\kappa(X) = \|X\|_2 \; \|X^{-1}\|_2$ is the [[Conditioning of a linear system|condition number]] of $X$.

The take-home key ideas to understand the convergence of GMRES are:

1. [[Convergence of the Conjugate Gradients|Distribution of eigenvalues:]] are they distributed in a few compact clusters or uniformly distributed around 0? Are they accumulating near 0?
2. Condition number of the eigenvector basis: is it small or large?

## Failure case

A challenging scenario arises when the eigenvalues of a matrix $A$ are uniformly distributed on the unit circle, which is the circle with a radius of 1 centered at the origin.

Consider, for instance, if $A$ is a permutation matrix that represents a permutation $\pi: \{1, \ldots, n\} \to \{1, \ldots, n\}$, and our goal is to solve $Ax = e_1.$ The solution in this case is the standard basis vector $e_r$, where $\pi(r) = 1.$ In this context, $\mathcal K(A, e_1, k)$ represents a coordinate subspace, spanned by $\{e_1, e_{\pi(1)}, e_{\pi^2(1)}, \ldots, e_{\pi^{k-1}(1)}\}.$ Convergence is achieved only when $k$ is such that $e_r$ is within $\mathcal K(A, e_1, k),$ which occurs if and only if $\pi^{k}(1) = 1.$ This might require $k$ to be as large as $n,$ especially in cases where $\pi$ is a cyclic permutation like $\pi(i) - 1 = i + 1$ mod $n.$ We also notice that, in this case, the solution remains orthogonal to the Krylov subspace until $k=n$.