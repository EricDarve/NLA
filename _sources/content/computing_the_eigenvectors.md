# Computing Eigenvectors from the Schur Decomposition

Once QR iteration has produced a Schur decomposition $A=Q_\star TQ_\star^H$ for $A\in\mathbb C^{n\times n}$, we may also want its eigenvectors. The columns of $Q_\star$ are Schur vectors and generally are not eigenvectors of $A$. To obtain eigenvectors, we first find those of the upper triangular matrix $T$, then transform them using the unitary matrix $Q_\star$.

If $v$ is an eigenvector of $T$ for an eigenvalue $\lambda$ (meaning $Tv = \lambda v$), then $x = Q_\star v$ is the corresponding eigenvector of $A$. This is because:

$$
\begin{aligned}
Ax&=(Q_\star TQ_\star^H)(Q_\star v)\\
&=Q_\star Tv=Q_\star(\lambda v)
=\lambda x.
\end{aligned}
$$

Our task, therefore, is to find a non-zero vector $v$ for a given eigenvalue $\lambda$ of $T$ such that $Tv = \lambda v$. This is equivalent to solving the homogeneous linear system:

$$(T - \lambda I)v = 0$$

## Solving for the Eigenvector

The eigenvalues of $T$ are its diagonal entries. First assume that $\lambda=t_{jj}$ is a **simple eigenvalue**, meaning that it occurs only once on the diagonal of $T$. This ensures that neither of the other diagonal blocks below has $\lambda$ as an eigenvalue. Partition $T$ and $v$ around position $j$:

$$
T = \begin{pmatrix} T_{11} & T_{12} & T_{13} \\
0 & \lambda & T_{23} \\
0 & 0 & T_{33}
\end{pmatrix}, \quad v = \begin{pmatrix} v_1 \\ v_2 \\ v_3 \end{pmatrix}
$$

Here, $T_{11}$ is $(j-1)\times(j-1)$, $T_{33}$ is $(n-j)\times(n-j)$, and $v_2$ is a scalar. The diagonal blocks are upper triangular. For $j=1$ or $j=n$, omit the corresponding empty blocks.

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
    Since $\lambda$ is simple, it is not a diagonal entry of $T_{33}$. Thus $T_{33}-\lambda I$ is invertible, and the only solution is $v_3=0$.

2.  **Second Row:** $0 \cdot v_1 + 0 \cdot v_2 + T_{23} v_3 = 0$.
    Substituting $v_3 = 0$, this equation becomes $0 = 0$. This is always satisfied and provides no information for finding $v_1$ or $v_2$.

3.  **First Row:** $(T_{11} - \lambda I)v_1 + T_{12}v_2 + T_{13}v_3 = 0$.
    With $v_3 = 0$, this simplifies to $(T_{11} - \lambda I)v_1 + T_{12}v_2 = 0$.

The remaining equation is $(T_{11}-\lambda I)v_1=-T_{12}v_2$. Since $\lambda$ is simple, $T_{11}-\lambda I$ is also invertible. Thus each choice of $v_2$ determines exactly one $v_1$. Choosing $v_2=0$ would give $v_1=v_3=0$, so every nonzero eigenvector has $v_2\ne0$. We can therefore scale it so that $v_2=1$. The first-row equation becomes:

$$
\begin{aligned}
(T_{11}-\lambda I)v_1+T_{12}&=0,\\
(T_{11}-\lambda I)v_1&=-T_{12}.
\end{aligned}
$$

In mathematical notation, the solution is:

$$v_1 = -(T_{11} - \lambda I)^{-1} T_{12}$$

## Solution and Cost

We have found all components of the eigenvector $v$ for the triangular matrix $T$:

$$
v = \begin{pmatrix} v_1 \\ v_2 \\ v_3 \end{pmatrix} = \begin{pmatrix}
- (T_{11}-\lambda I)^{-1} T_{12} \\ 1 \\ 0
\end{pmatrix}
$$

The corresponding eigenvector of $A$ is $x=Q_\star v$. To obtain a unit eigenvector, use $x=Q_\star v/\|v\|_2$, since multiplication by $Q_\star$ preserves the Euclidean norm.

In computation, solve $(T_{11}-\lambda I)v_1=-T_{12}$ by back-substitution rather than forming an inverse. This costs $O(j^2)$ operations, and the multiplication by $Q_\star$ costs at most $O(n^2)$. Thus each eigenvector costs $O(n^2)$ operations after the Schur decomposition is available. If all eigenvalues are simple, computing all $n$ eigenvectors costs $O(n^3)$.

## Scope and Numerical Considerations

For a repeated eigenvalue, the shifted diagonal blocks need not be invertible, so the inverse formula above is not generally valid. The corresponding eigenspace of $T$ is the null space of $T-\lambda I$. The dimension of this eigenspace can be smaller than the eigenvalue's algebraic multiplicity; in that case the matrix is defective. Thus we cannot expect one independent eigenvector for every occurrence on the diagonal.

The derivation uses an upper triangular Schur form over $\mathbb C$. A real Schur form can contain $2\times2$ diagonal blocks for complex conjugate eigenvalues; these require block back-substitution or conversion to complex Schur form.

In floating-point arithmetic, small differences $t_{ii}-\lambda$ can produce large intermediate components. Practical routines use scaled triangular solves to avoid overflow, as in [LAPACK's eigenvector routine](https://www.netlib.org/lapack/complex16/ztrevc.f). Closely spaced eigenvalues can also make individual eigenvectors sensitive to perturbations, even when the Schur decomposition is computed accurately.
