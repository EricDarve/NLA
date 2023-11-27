This may not be obvious at first sight, but [[SOR iteration|SOR]] is a [[Splitting methods|splitting method:]]
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
G_\text{SOR} = M^{-1} N
$$
We can prove that
$$
\rho(G_\text{SOR}) \ge |\omega - 1|.
$$
Therefore $0 < \omega < 2$ is required for convergence.

**Proof**

The proof uses the [[Determinant|determinant]] of $G_\text{SOR}$:
$$
\det G_\text{SOR} = (\det (\omega^{-1} D - L))^{-1} \;
\det(( \omega^{-1} - 1) D + U)
$$
Using the fact that the matrices are triangular, we get
$$
\det G_\text{SOR} = \omega^n \; (\omega^{-1} - 1)^n
= (1 - \omega)^n
$$
If all the eigenvalues are smaller than 1, then the determinant is in the interval $(-1,1)$ and
$$
|1 - \omega| < 1
$$
which implies that $0 < \omega < 2$ for stability. This is a required but not sufficient condition.

**Theorem 2.**

- If $A$ is [[Symmetric Positive Definite Matrices|symmetric positive definite,]] then $\rho(G_{SOR}) < 1$ for all $0 < \omega < 2$. 
- In that case, SOR converges for all $0 < \omega < 2$. 
- With $\omega = 1$, we recover that [[Gauss-Seidel iteration|Gauss-Seidel]] converges for [[Symmetric Positive Definite Matrices|symmetric positive definite matrices.]]