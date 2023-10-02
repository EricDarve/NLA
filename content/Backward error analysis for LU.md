This is the backward error bound for LU:
$$
(A + E) \; \tilde x = b\\[1em]
|E| \le n \, u \; (2|A| + 4 |\tilde{L}| |\tilde{U}|) + O(u^2)
$$

The backward error can become very large when we have large entries in $L$ or $U$. The LU factorization is **not** a backward stable algorithm.

In our previous example:
$$
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
$$
Error is $O(u \epsilon^{-1})$. It can become arbitrarily large.

Where do these large entries come from?

Recall the LU factorization:
- $u_{k,} = a_{k,}$
- $l_{,k} = a_{,k} / a_{kk}$
- $A \leftarrow A - l_{,k} * u_{k,}$

We have a problem when the pivot $a_{kk}$ becomes very small.

[[Stability of the LU factorization]], [[Forward and backward error]]