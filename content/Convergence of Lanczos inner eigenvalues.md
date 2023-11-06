We can extend the [[Convergence of Lanczos eigenvalues for symmetric matrices|previous result]] for $\lambda_1$.

With $A \to -A$, the same result applies to $\lambda_n$.

There is a similar result for “inner” eigenvalues, but the convergence is slower.
$$
\begin{gather}
\lambda_i \ge \lambda_i(T_k) \ge
\lambda_i - (\lambda_1 - \lambda_n) \;
\Big( \frac{\kappa_i \tan \phi}{T^{\rm cheb}_{k-i}(1 + 2\rho_i)} \Big)^2 \\[1em]
\rho_i = \frac{\lambda_i - \lambda_{i+1}}{\lambda_{i+1} - \lambda_n}, \quad
\kappa_i = \prod_{j=1}^{i-1} \frac{\lambda_j(T_k)-\lambda_n}{\lambda_j(T_k)-\lambda_i}
\end{gather}
$$
For $\lambda_i(T_k)$, the convergence slows down with increasing $i$.

[[Lanczos process]], [[Convergence of Lanczos eigenvalues for symmetric matrices]]