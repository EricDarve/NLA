This is the backward error bound for LU:
$$
\begin{gather}
(A + E) \; \tilde x = b \\
|E| \le n \, u \; (2|A| + 4 |\tilde{L}| \: |\tilde{U}|) + O(u^2)
\end{gather}
$$
- $|\tilde{L}|$, $|\tilde{U}|$: matrices obtained by taking the absolute value of the entries of $\tilde{L}$ and $\tilde{U}$. 
- We use the same notation for $|E|$ and $|A|$.
- The backward error can become very large when we have large entries in $L$ or $U$. 
- The LU factorization is **not** a backward stable algorithm.

In our previous example:
$$
\begin{gather}
A =
\begin{pmatrix}
\epsilon & 1 \\
1 & \pi
\end{pmatrix} \\[1em]
L = \begin{pmatrix}
1 & 0 \\ \epsilon^{-1} & 1
\end{pmatrix}, \qquad
U = \begin{pmatrix}
\epsilon & 1 \\ 0 & \pi - \epsilon^{-1}
\end{pmatrix}
\end{gather}
$$
- The backward error is $O(u \epsilon^{-1})$, where $u$ is the [[Floating point arithmetic and unit roundoff error|unit roundoff]].
- It can become arbitrarily large regardless of how small $u$ is. 
- We can choose $\epsilon$ as small as we want and make the backward error large.
- The factorization is not backward stable.

**Where do these large entries come from?**

Recall the LU factorization:
- $u_{k,} = a_{k,}$
- $l_{,k} = a_{,k} / a_{kk}$
- $A \leftarrow A - l_{,k} * u_{k,}$

We have a problem when the [[LU algorithm#^pivot|pivot]] $a_{kk}$ becomes very small.

[[LU algorithm]], [[Stability of the LU factorization]], [[Forward and backward error]]