This may not be obvious at first sight, but [[SOR iteration|SOR]] is a splitting method:
$$
x^{(k+1)}_{SOR} = x^{(k)}_{SOR} + \omega \bigg[ D^{-1}
\big( b + L x^{(k+1)}_{SOR} + U x^{(k)}_{SOR} \big) - x^{(k)}_{SOR} \bigg]
$$
This implies:
$$
\big (\frac{1}{\omega} D - L \big) \; x^{(k+1)}_{SOR} = b + \Big(
\big( \frac{1}{\omega} - 1 \big) D + U \Big) \; x^{(k)}_{SOR}
$$
So we indeed have a [[Splitting methods|splitting method]] with the following choices:
$$
\begin{gather}
A = M - N, \\
M = \frac{1}{\omega} D - L, \\
N = \Big( \frac{1}{\omega} - 1 \Big) D + U
\end{gather}
$$
The choice $\omega = 1$ corresponds to [[Gauss-Seidel iteration|Gauss-Seidel.]]

## [[Convergence of classical iterative methods|Convergence]]

**Theorem 1.**

We have
$$
G_{SOR} = M^{-1} N
$$
We can prove that
$$
\rho(G_{SOR}) \ge |\omega - 1|.
$$
Therefore $0 < \omega < 2$ is required for convergence.

**Theorem 2.**

- If $A$ is symmetric positive definite, then $\rho(G_{SOR}) < 1$ for all $0 < \omega < 2$. 
- In that case, SOR converges for all $0 < \omega < 2$. 
- With $\omega = 1$, we recover that [[Gauss-Seidel iteration|Gauss-Seidel]] converges for [[Symmetric Positive Definite Matrices|symmetric positive definite matrices.]]