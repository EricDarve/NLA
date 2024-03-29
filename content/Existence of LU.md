Let's look at a case where the LU factorization breaks down. This happens when a [[LU algorithm#^pivot|pivot]] becomes equal to 0.

Consider:
$$
A = \begin{pmatrix}
1 & 6 & 1 & 0 \\
0 & 1 & 9 & 3 \\
1 & 6 & 1 & 1 \\
0 & 0 & 1 & 9
\end{pmatrix}
$$

Why does the LU algorithm break down?

- $A = LU$
- We have $\det(A) = \det(L) \det(U) = \det(U)$ because $L$ is lower triangular with 1 on the diagonal.
- det($A$) = 0 if and only if det($U$) = 0. This is equivalent to saying that one of the diagonal entries of $U$ is 0. [[LU and determinant]].
- Furthermore, because $A = LU$ and the triangular form of $L$ and $U$ we have that
$$
A[1:k,1:k] = L[1:k,1:k] \: U[1:k,1:k]
$$
- Assume that all the pivots are non-zero up to step $k$. Then $A[1:k,1:k]$ is singular if and only if $u_{kk}$ is zero. This can be stated equivalently as $\det(A[1:k,1:k]) \neq 0$ if and only if $u_{kk} \neq 0$. [[LU and determinant]].
- So, if $\det(A[1:k,1:k]) \neq 0$, $1 \le k \le n-1,$ then all the pivots remain non-zero, and the algorithm completes.

In our example, $u_{33} = 0$ at step $k=3$. This is because the top left $3 \times 3$ block of $A$ is singular. So we cannot proceed to step $k=4$.

[[Triangular factorization]], [[LU algorithm]], [[LU and determinant]]

### Existence of the LU factorization

Theorem: the LU factorization exists if
$$
\det(A[1: k, 1: k]) \neq 0.
$$
for all $1 \le k \le n-1$.

This was proved above. A more precise version of this result is as follows.

Theorem: The LU factorization exists if and only if
$$
{\rm rank}(A[1:k,1:k]) = {\rm rank}(A[1:n,1:k])
$$
for all $1 \le k \le n-1$.

### More details

We make a few observations.

A matrix $A$ is non-singular if and only if ${\rm rank}(A[1:n,1:k]) = k$, $1 \le k \le n$. For non-singular matrices, $\det(A[1: k, 1: k]) \neq 0$, $1 \le k \le n-1$, is a necessary and sufficient condition for the existence of an LU factorization.

Assume that $A$ has an LU factorization. If $u_{kk} = 0$, then $a_{,k}$ is a linear combination of the columns $l_{,j}$, $j < k$. If $u_{kk} = 0$ and $u_{ll} \neq 0$ for $l < k$, then $a_{,k}$ is a linear combination of the columns $a_{,j}$, $j < k$.

When $u_{kk} = 0$, column $l_{,k}$ of $L$ is not uniquely defined, and an infinite number of LU factorizations satisfy $A = LU$.

We also have:

Theorem: There exists a unique $L$ and $U$ if and only if $\det(A[1: k, 1: k]) \neq 0$ for all $1 \le k \le n-1$.

Variant:

Theorem: There exists a unique non-singular $L$ and $U$ if and only if $\det(A[1: k, 1: k]) \neq 0$ for all $1 \le k \le n$.