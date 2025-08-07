Assume, for example, that we can quickly solve $Mx = b$ and $M \approx A$. Then:
$$
\begin{gather}
Ax = b \\
(A - M) x + Mx = b \\
Mx = b - (A-M) x
\end{gather}
$$
This suggests a simple iterative strategy based on the idea of fixed-point iteration.
$$
\begin{gather}
A = M - N \\
M x^{(k+1)} = b + N x^{(k)}
\end{gather}
$$
If $x^{(k)}$ converges, it must converge to the solution of $Ax = b$.

This method, when it converges, requires solving with $M$ instead of $A$. In many cases, this allows solving $Ax=b$ approximately very fast.