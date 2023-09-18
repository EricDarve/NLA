Numerical analysis precept: consider some algorithm $A(X)$. Suppose the algorithm breaks down when $X=X_0$ (e.g., division by 0). Then the numerical error is typically large when $X \sim X_0$.

Consider this case:
$$
A =
\begin{pmatrix}
0 & 1 \\
1 & \pi
\end{pmatrix}
$$
LU breaks down immediately because $a_{11} = 0$.

What happens when we have a small [[LU algorithm#^pivot]]?

Consider now:
$$
\begin{gather}
A =
\begin{pmatrix}
	\epsilon & 1 \\
	1 & \pi
\end{pmatrix}, \\[1em]
L = \begin{pmatrix}
	1 & 0 \\ \epsilon^{-1} & 1
\end{pmatrix}, \quad
U = \begin{pmatrix}
	\epsilon & 1 \\ 0 & \pi - \epsilon^{-1}
\end{pmatrix}
\end{gather}
$$

- Floating point numbers on a computer can only be represented by a finite number of digits.
- Consider $\pi - \epsilon^{-1}$.
- When $\epsilon \ll 1$, $\epsilon^{-1} \gg \pi$.
- For $\epsilon$ sufficiently small, $\pi - \epsilon^{-1} \equiv -\epsilon^{-1}$ on a computer.
- Because of numerical roundoff errors, on a computer, we get
$$
U = \begin{pmatrix}
\epsilon & 1 \\ 0 & {\color{red}\cancel{\pi}} - \epsilon^{-1}
\end{pmatrix}
$$
for $\epsilon$ sufficiently small. $LU$ is no longer equal to $A$!

If we solve the linear system using the numerical approximations to $L$ and $U$, we will get the wrong result.

A computer program cannot store all the digits of $\epsilon^{-1}$. If $\epsilon^{-1} \approx 10^{20}$, $\pi$ appears around digit 20. If you store only 16 digits of $\epsilon^{-1}$, there is no room to store $\pi$.

Consequence: $\epsilon^{-1} \equiv \epsilon^{-1} - \pi$!

[[Solving linear systems using LU]], [[Triangular factorization]], [[LU algorithm]], [[Existence of LU]]