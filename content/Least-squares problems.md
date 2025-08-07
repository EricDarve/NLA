Solve $x^* = \text{argmin}_x \; \| Ax - b \|_2$ where $A$ has size $m \times n$.

![[2022-10-11-15-34-53.png]]
Because we use the [[Vector norms|2-norm]], the solution is found when $Ax - b \; \perp \; \text{span}(A)$ ([[Pythagorean theorem]]).

This is equivalent to the equation
$$
\begin{gather}
A^T (Ax - b) = 0 \\
(A^T A) \, x = A^T b
\end{gather}
$$
See [[Dot product|dot product]] and orthogonality.

This equation is a linear equation in $x$ with a square matrix $A^T A$.

