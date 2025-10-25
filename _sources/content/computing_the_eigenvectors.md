# Computing Eigenvectors from the Schur Decomposition

We will focus on coming sections on computing the eigenvalues and accelerating their convergence using shifts. However, once we have computed the Schur decomposition of $A = Q T Q^H$, we may also want to compute its eigenvectors. This can be achieved by first finding the eigenvectors of the upper triangular matrix $T$. Let's do that now.

If $v$ is an eigenvector of $T$ for an eigenvalue $\lambda$ (meaning $Tv = \lambda v$), then $x = Qv$ is the corresponding eigenvector of $A$. This is because:

$$Ax = (QTQ^H)(Qv) = QTv = Q(\lambda v) = \lambda(Qv) = \lambda x$$

Our task, therefore, is to find a non-zero vector $v$ for a given eigenvalue $\lambda$ of $T$ such that $Tv = \lambda v$. This is equivalent to solving the homogeneous linear system:

$$(T - \lambda I)v = 0$$

## Solving for the Eigenvector

The eigenvalues $\lambda$ of $T$ are its diagonal elements. Let's find the eigenvector $v$ corresponding to a specific eigenvalue $\lambda$. We can partition the matrix $T$ and the vector $v$ to isolate this $\lambda$:

$$
T = \begin{pmatrix} T_{11} & T_{12} & T_{13} \\
0 & \lambda & T_{23} \\
0 & 0 & T_{33}
\end{pmatrix}, \quad v = \begin{pmatrix} v_1 \\ v_2 \\ v_3 \end{pmatrix}
$$

Here, $T_{11}$ and $T_{33}$ are upper triangular matrices, and $v_2$ is a scalar component corresponding to $\lambda$.

The system $(T - \lambda I)v = 0$ can now be written in block form:

$$
\begin{pmatrix} T_{11} - \lambda I & T_{12} & T_{13} \\
0 & 0 & T_{23} \\
0 & 0 & T_{33} - \lambda I
\end{pmatrix}
\begin{pmatrix} v_1 \\ v_2 \\ v_3 \end{pmatrix}
=
\begin{pmatrix} 0 \\ 0 \\ 0 \end{pmatrix}
$$

We solve this system using block back-substitution, starting from the bottom row:

1.  **Third Row:** $(T_{33} - \lambda I)v_3 = 0$.
    The eigenvalues of $T_{33}$ are its diagonal elements. Assuming $\lambda$ is distinct from the eigenvalues of $T_{33}$, the matrix $(T_{33} - \lambda I)$ is invertible. Therefore, the only solution is $v_3 = 0$.

2.  **Second Row:** $0 \cdot v_1 + 0 \cdot v_2 + T_{23} v_3 = 0$.
    Substituting $v_3 = 0$, this equation becomes $0 = 0$. This is always satisfied and provides no information for finding $v_1$ or $v_2$.

3.  **First Row:** $(T_{11} - \lambda I)v_1 + T_{12}v_2 + T_{13}v_3 = 0$.
    With $v_3 = 0$, this simplifies to $(T_{11} - \lambda I)v_1 + T_{12}v_2 = 0$.

We are left with one equation: $(T_{11} - \lambda I)v_1 + T_{12}v_2 = 0$. Since $\lambda$ is an eigenvalue, the full system $(T - \lambda I)v = 0$ is singular and must have a non-trivial solution. The singularity (the $0$ on the diagonal) corresponds to the $v_2$ component, which acts as a free variable.

Because eigenvectors are unique only up to a non-zero scalar multiple, we can set $v_2 = 1$ to find one such eigenvector. The first-row equation becomes:

$$(T_{11} - \lambda I)v_1 + T_{12}(1) = 0 \implies (T_{11} - \lambda I)v_1 = -T_{12}$$

Similar to the $T_{33}$ case, $(T_{11} - \lambda I)$ is invertible because $\lambda$ is not an eigenvalue of $T_{11}$. We can now solve for $v_1$:

$$v_1 = -(T_{11} - \lambda I)^{-1} T_{12}$$

## Solution and Cost

We have found all components of the eigenvector $v$ for the triangular matrix $T$:

$$
v = \begin{pmatrix} v_1 \\ v_2 \\ v_3 \end{pmatrix} = \begin{pmatrix}
- (T_{11}-\lambda I)^{-1} T_{12} \\ 1 \\ 0
\end{pmatrix}
$$

The corresponding eigenvector $x$ of the original matrix $A$ is $x = Qv$.

Computationally, finding $v_1$ involves solving a triangular linear system, which costs $O(n^2)$ operations. The final matrix-vector multiplication $x = Qv$ also costs $O(n^2)$. Thus, the computational cost to find each eigenvector from the Schur decomposition is $O(n^2)$.