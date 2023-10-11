**Remedy for instability:** perform row pivoting so that $a_{kk}$ is always the largest entry in the column.

This guarantees that $|l_{ij}| \le 1$.

With row pivoting, our factorization looks like
$$
PA = LU
$$
Note that although $L$ is bounded by definition, there is no guarantee for $U$. But in most practical cases, this algorithm is very accurate.

Let's apply row pivoting to our previous case:
$$
\begin{gather}
A =
\begin{pmatrix}
\epsilon & 1 \\
1 & \pi
\end{pmatrix}, \qquad
P A =
\begin{pmatrix}
1 & \pi \\
\epsilon & 1
\end{pmatrix} 
\\[1em]
L = \begin{pmatrix}
1 & 0 \\ \epsilon & 1
\end{pmatrix}, \qquad
U = \begin{pmatrix}
1 & \pi \\ 0 & 1 - \epsilon \pi
\end{pmatrix}, \qquad PA = LU
\end{gather}
$$
This new factorization is **backward stable and accurate.**

[[Stability of the LU factorization]], [[Forward and backward error]], [[Backward error analysis for LU]]