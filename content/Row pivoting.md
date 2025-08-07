- **Remedy for instability:** perform row pivoting so that $a_{kk}$ is always the largest entry in the column.
- Since $l_{ij} = a_{ij} / a_{jj}$, row pivoting guarantees that $|l_{ij}| \le 1$.

With row pivoting, our factorization looks like
$$
PA = LU
$$
where $P$ is a permutation of the rows of $A$. $P$ is key for the existence of the LU factorization.

- Note that although $L$ is bounded by definition, there is no guarantee for $U$. 
- In most practical cases, this algorithm is very accurate.
- This factorization exists for all square matrices.
- This factorization is unique for all matrices such that $\det(A[1: k, 1: k]) \neq 0$ for all $1 \le k \le n-1$.

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