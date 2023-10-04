One of the most important matrix decomposition in data science and machine learning.

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
	- $\Sigma (V^ x)$: point on an ellipsoid; the axes are aligned with the coordinate axes.
	- $U(\Sigma V^T x)$: rotate/reflect the ellipsoid.

The lengths of the axes of this ellipsoid are the singular values of $A$.

![[Drawing 2023-10-04 12.26.03.excalidraw.svg]]

As we can expect, the [[Operator and matrix norms|size of a matrix]] can be related to its singular values:
$$
\| A \|_2 = \sigma_1(A), \qquad
\| A \|_F = \sqrt{\sum_{i=1}^p \sigma_i^2}
$$

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
[[The four fundamental spaces]], [[Eigenvalues]], [[Operator and matrix norms]], [[Orthogonal matrix and projector]]