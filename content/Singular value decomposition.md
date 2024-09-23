One of the most important matrix decompositions in data science and machine learning.

Eigenvalues are great for understanding $A^n$. Does it grow? Does it shrink? How can it be easily modeled?

But it says nothing about the size of $Ax$. How is $A$ transforming $x$ and rescaling the vector?

Because $A$ is a matrix, multiple scales are involved when applying $A$. These scales can be represented systematically using the singular value decomposition.
$$
A = U \Sigma V^T
$$
- $A$: $m \times n$.
- $U$: orthogonal, $m \times m$. Left singular vectors
- $V$: orthogonal, $n \times n$. Right singular vectors
- $\Sigma$: $m \times n$, diagonal matrix with real positive entries = pure scaling = singular values.

![[Drawing 2023-10-04 12.13.49.excalidraw.svg]]

![[Drawing 2023-10-04 2.excalidraw.svg]]

- $A v_i = \sigma_i \, u_i$

Consider a ball $B$ in $R^n$. $A$ transforms this ball into an ellipsoid.

- $Ax = U \Sigma V^Tx$:
	- $V^T x$: point on the unit ball
	- $\Sigma (V^T x)$: point on an ellipsoid; the axes are aligned with the coordinate axes.
	- $U(\Sigma V^T x)$: rotate/reflect the ellipsoid.

The lengths of the axes of this ellipsoid are the singular values of $A$.

![[Drawing 2023-10-04 12.26.03.excalidraw.svg]]

As we can expect, the [[Operator and matrix norms|size of a matrix]] can be related to its singular values:
$$
\| A \|_2 = \sigma_1(A), \qquad
\| A \|_F = \sqrt{\sum_{i=1}^p \sigma_i^2}
$$
We can also define a new operator norm using the SVD, the Schatten $p$-norm:
$$
\| A \|_p^\text{{Schatten}} = \| (\sigma_1, \dots, \sigma_r) \|_{p}
$$
where $r$ is the rank and $\| \; \|_{p}$ is the vector $p$-norm.

[[The four fundamental spaces]]. Assume $A$ is $m \times n$. $r$: number of non-zero singular values = rank of the matrix. Then:
$$
\begin{aligned}
N(A) & = \{v_{r+1}, \ldots, v_n\} \\
R(A) & = \{u_1, \ldots, u_r\} \\
N(A^T) & = \{u_{r+1}, \ldots, u_m\} = R(A)^\perp \\
R(A^T) & = \{v_1, \ldots, v_r\} = N(A)^\perp
\end{aligned}
$$
We recover the [[The four fundamental spaces|four fundamental spaces]] and the [[The four fundamental spaces|rank-nullity theorem]].

Connection with [[Eigenvalues|eigenvalues]]. The eigenvalues of $AA^T$ and $A^T A$ are equal to $\sigma_1^2,$ ..., $\sigma_r^2,$ or 0. The eigenvectors of $AA^T$ are given by $U$, and those of $A^T A$ by $V$:
$$
AA^T = U \Sigma^2 U^T, \qquad
A^TA = V \Sigma^2 V^T
$$
The computational cost of computing the singular value decomposition is $O(n^3)$.

**Proof of the existence of the SVD.** Multiple proofs are possible. Let's look at one of them. Define the matrix
$$
B = \begin{bmatrix} 0 & A \\ A^T & 0 \end{bmatrix}
$$
This matrix is [[Hermitian and symmetric matrices|symmetric]] and therefore is [[Unitarily diagonalizable matrices|unitarily diagonalizable]]:
$$
B = Q \Lambda Q^T
$$
We can restrict this factorization such that $\Lambda$ has only non-zero entries on the diagonal. Assume that $[x;y]$ is an eigenvector of $B$ associated with $\lambda \ne 0$. Then
$$
B \begin{bmatrix} x \\ -y \end{bmatrix} = \begin{bmatrix} 0 & A \\ A^T & 0 \end{bmatrix} \begin{bmatrix} x \\ -y \end{bmatrix} = \begin{bmatrix} -Ay \\ A^Tx \end{bmatrix} = -\lambda \begin{bmatrix} x \\ -y \end{bmatrix}
$$
So $-\lambda$ is also an eigenvalue. Since the eigenvectors must be [[Orthogonal matrix and projector|orthogonal]] we have
$$
x^T x - y^T y = 0
$$
Let's normalize our eigenvector such that
$$
x^T x + y^T y = 2
$$
This implies that $\|x\|_2 = \|y\|_2 = 1$.

Note that since $Ay = \lambda x$ and $A^T x = \lambda y$, we have
$$
A^T A y = \lambda^2 y, \quad AA^T x = \lambda^2 x
$$
So $\lambda^2$ is an eigenvalue of $A^T A$ and $AA^T$. 

Denote by $X$ all the eigenvectors of $AA^T$ and by $Y$ those of $A^TA$ (keeping only the non-zero eigenvalues). We have shown that
$$
Q = \frac{1}{\sqrt{2}} \begin{bmatrix} X & X \\ Y & -Y \end{bmatrix}, \quad
\Lambda = \begin{bmatrix} \Sigma & 0 \\ 0 & -\Sigma \end{bmatrix}
$$
Using
$$
\begin{bmatrix} 0 & A \\ A^T & 0 \end{bmatrix}
= Q \Lambda Q^T
$$
we find that
$$
A = X \Sigma Y^T
$$
This is the SVD of $A$.

$\square$

The SVD is unique if and only if all the singular values are distinct (up to a sign change of columns in $U$ and $V$).

[[The four fundamental spaces]], [[Eigenvalues]], [[Operator and matrix norms]], [[Orthogonal matrix and projector]]