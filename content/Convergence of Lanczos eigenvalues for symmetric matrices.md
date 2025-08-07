Consider the largest eigenvalue of a symmetric matrix. It satisfies:
$$
\lambda_1 = \max_{x \neq 0} \frac{\boldsymbol x^T A \boldsymbol x}{\| \boldsymbol x \|^2_2}
$$
The [[Lanczos process|Lanczos]] approximation uses [[Convergence of the orthogonal iteration|Ritz eigenvalues]] (see also [[Arnoldi process]]):
$$
T_k = Q_k^T A Q_k
$$
From
$$
\lambda_1(T_k) = \max_{y \neq 0 } \frac{ \boldsymbol y^T Q_k^T A Q_k \boldsymbol y}{\| \boldsymbol y \|_2^2}
$$
we immediately find that $\lambda_1(T_k) \le \lambda_1$.

We can do a more precise analysis using [[Connection between Arnoldi and polynomials of A|polynomials]] again. Here is a brief outline. We have, using the [[Connection between Arnoldi and polynomials of A|polynomial]] interpretation of [[Lanczos process|Lanczos]] and [[Krylov subspace|Krylov subspaces:]]
$$
\lambda_1(T_k) = \max_p \frac{ \boldsymbol q_1^T p(A)^T A p(A) \boldsymbol q_1 }{ \boldsymbol q_1^T p(A)^2 \boldsymbol q_1 }
$$
where $p$ is a polynomial of degree $k-1$.

By making the special choice of Chebyshev polynomials, we can prove that:
$$
\lambda_1 \ge \lambda_1(T_k) \ge
\lambda_1 - (\lambda_1 - \lambda_n) \;
\Big( \frac{\tan \phi}{T^{\rm cheb}_{k-1}(1 + 2\rho_1)} \Big)^2
$$
where $\phi$ depends on $q_1$. We defined
$$
\rho_1 = \frac{\lambda_1 - \lambda_2}{\lambda_2 - \lambda_n}
$$
Here is a schematic of the distribution of eigenvalues along the real axis for symmetric matrices:

![[2022-11-02-12-16-26.png|350]]

With Chebyshev polynomials, we have
$$
T^{\rm cheb}_{k-1}(1 + 2\rho_1) = O \big[ (4\rho_1)^k \big]
$$
We obtain a rapid convergence when $\lambda_1$ is well separated from the other eigenvalues and $\rho_1 \gg 1$.

![[2022-11-02-12-17-51.png]]

$T^{\rm cheb}_{k-1}(1 + 2\rho_1)$ becomes very large with $k$ quickly.