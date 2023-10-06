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
- If det($A$) = 0, then det($U$) = 0. 
- One of the diagonal entries of $U$ must be 0. 
- From the algorithm, this implies that $a_{kk} = 0$ at some point. 
- The algorithm is dividing by 0. One of the [[LU algorithm#^pivot|pivots]] is 0.

In our example, $u_{33} = 0$ at step $k=3$. We cannot proceed to step $k=4$.

[[Triangular factorization]], [[LU algorithm]]

### Existence of the LU factorization

Theorem: the LU factorization exists if
$$
\det(A[1: k, 1: k]) \neq 0.
$$
for all $1 \le k \le n-1$.

A more precise version of this result is as follows.

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