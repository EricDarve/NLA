For any square matrix $A$, there exists at least a scalar $\lambda \in \mathbb C$ and $x \in \mathbb C^n$ such that $Ax = \lambda x$. They are called an eigenvalue and eigenvector of $A$.

### Existence of eigenvalues

We prove the existence of at least one eigenvalue and eigenvector $Ax = \lambda x$ with $\lambda \in \mathbb C$.

Consider the sequence $x,$ $Ax,$ ..., $A^nx.$ We have $n+1$ vectors so they are linearly dependent. So there exists $a_0,$ ..., $a_n$ such that
$$
a_0 x + \dots + a_n A^n x = 0
$$
Using the fundamental theorem of algebra we can find $c$ and $\lambda_i$ such that
$$
a_0 + a_1 z + \dots + a_n z^n = c (z-\lambda_1) \cdots (z-\lambda_m), \qquad c \neq 0
$$
This shows that
$$
(A-\lambda_1 I) \cdots (A-\lambda_m I) x = 0
$$
If $A-\lambda_i I$ is non-singular for all $i$ then
$$
(A-\lambda_1 I) \cdots (A-\lambda_m I)
$$
is non-singular. This is a contradiction.  So there exists an $i$ such that $A-\lambda_i I$ is singular. So $\lambda_i$ is an eigenvalue.

### Finding all the eigenvalues

We show that there exists a non-singular matrix $X$ and upper triangular matrix $T$ such that
$$
A = X T X^{-1}
$$
We will show that all the eigenvalues of $A$ can be found along the diagonal of $T$.

Let's prove the existence of $X$ and $T$.

We use a proof by induction. The result can be verified when $n=1.$

Assume that the result is true for matrices of size less than $n.$ 

We know that there exists $x$ and $\lambda$ such that $A x = \lambda x.$ Denote by
$$
U = \text{range}(A - \lambda I)
$$
Since $A - \lambda I$ is singular $\dim(U) < n$. 

Moreover, $U$ is stable, that is $AU \subset U$. To prove this, take $u \in U$. Then:
$$
Au = (A - \lambda I) u + \lambda u
$$
Since $(A - \lambda I) u \in U$, we have $Au \in U$. So $U$ is stable. 

We can consider the restriction of $A$ to the subspace $U$. Using the induction hypothesis, there is a basis $u_1,$ ..., $u_m$ such that $A$ restricted to $U$ is upper triangular in that basis.

Let's extend this set of independent vectors to get a full basis of $\mathbb R^n$:
$$
X = \big[ u_1, \dots, u_m, v_1, \dots, v_{n-m} \big]
$$
Take a $v_k$ and
$$
A v_k = (A - \lambda I) v_k + \lambda v_k
$$
So $Av_k \in \text{span}(u_1, \dots, u_m, v_k)$. This proves that
$$
X^{-1} A X = T
$$
where $T$ is upper triangular. So $A = X T X^{-1}$.

The eigenvalues of $A$ are equal to the eigenvalues of $T$. Moreover, we can prove that $T - \lambda I$ is singular if and only if this matrix has a 0 on the diagonal. This implies that all the eigenvalues of $T$ are on its diagonal, and there are $n$ of them. Note that these eigenvalues might be repeated.

[[Matrix-vector and matrix-matrix product]], [[Invertible matrix]]