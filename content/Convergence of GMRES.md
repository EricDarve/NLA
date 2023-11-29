[[GMRES|Recall]] that we are minimize the following norm
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

1. [[Convergence of the Conjugate Gradients|Distribution of eigenvalues:]] are they distributed in a few compact clusters or uniformly distributed around 0? What is the condition number of the matrix?
2. Condition number of the eigenvector basis: is it small or large?