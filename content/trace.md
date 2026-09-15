# The Trace of a Matrix

The **trace** of a square matrix is the sum of its diagonal entries. For $A\in\mathbb{R}^{n\times n}$ or $A\in\mathbb{C}^{n\times n}$,

$$\operatorname{tr}(A)=\sum_{i=1}^n a_{ii}.$$

For example,

$$
A=\begin{pmatrix}
3 & 0 & 1 \\
5 & -2 & 9 \\
4 & 6 & 8
\end{pmatrix}
\quad\Longrightarrow\quad
\operatorname{tr}(A)=3+(-2)+8=9.
$$

In particular, $\operatorname{tr}(I_n)=n$. Only diagonal entries contribute to the trace, so it takes $n-1$ additions to compute it.

## Properties of the Trace

### Linearity

For square matrices $A,B$ of the same size and scalars $\alpha,\beta$,

$$
\operatorname{tr}(\alpha A+\beta B)
=\alpha\operatorname{tr}(A)+\beta\operatorname{tr}(B).
$$

````{prf:proof}
The diagonal entries of $\alpha A+\beta B$ are $\alpha a_{ii}+\beta b_{ii}$. Thus

$$
\begin{aligned}
\operatorname{tr}(\alpha A+\beta B)
&=\sum_{i=1}^n(\alpha a_{ii}+\beta b_{ii}) \\
&=\alpha\sum_{i=1}^n a_{ii}+\beta\sum_{i=1}^n b_{ii}.
\end{aligned}
$$
````

### Cyclic Property

If $A$ is $m\times n$ and $B$ is $n\times m$, then

$$\operatorname{tr}(AB)=\operatorname{tr}(BA).$$

The factors may be rectangular: $AB$ is $m\times m$, while $BA$ is $n\times n$.

````{prf:proof}
Expand the diagonal entries of each product and interchange the finite sums:

$$
\begin{aligned}
\operatorname{tr}(AB)
&=\sum_{i=1}^m\sum_{j=1}^n a_{ij}b_{ji} \\
&=\sum_{j=1}^n\sum_{i=1}^m b_{ji}a_{ij} \\
&=\operatorname{tr}(BA).
\end{aligned}
$$
````

For three factors of compatible sizes, grouping them into two factors gives

$$
\operatorname{tr}(ABC)
=\operatorname{tr}(BCA)
=\operatorname{tr}(CAB).
$$

For example, apply the two-factor identity to $A$ and $BC$ to move $A$ from the front to the back. The same argument works for any number of factors: moving the first factor to the end preserves the trace.

An arbitrary reordering need not preserve the trace. For example, let

$$
A=\begin{pmatrix}1&0\\0&0\end{pmatrix},\qquad
B=\begin{pmatrix}0&1\\0&0\end{pmatrix},\qquad
C=\begin{pmatrix}0&0\\1&0\end{pmatrix}.
$$

Then $ABC=A$ and $ACB=0$, so $\operatorname{tr}(ABC)=1$ while $\operatorname{tr}(ACB)=0$.

The trace is also not multiplicative: in general, $\operatorname{tr}(AB)\neq\operatorname{tr}(A)\operatorname{tr}(B)$. Taking $A=B=I_2$ gives $2$ on the left and $4$ on the right.

### Transpose and Conjugate Transpose

For a square matrix $A$,

$$
\operatorname{tr}(A^T)=\operatorname{tr}(A),\qquad
\operatorname{tr}(A^H)=\overline{\operatorname{tr}(A)}.
$$

Here, $A^H$ denotes the conjugate transpose and the overline denotes complex conjugation.

````{prf:proof}
Transposition leaves the diagonal entries unchanged. The conjugate transpose replaces each diagonal entry $a_{ii}$ by $\overline{a_{ii}}$, so its trace is the complex conjugate of their sum.
````

### Invariance Under a Change of Basis

If $S$ is invertible, then

$$\operatorname{tr}(S^{-1}AS)=\operatorname{tr}(A).$$

Thus similar matrices have the same trace. Although their individual diagonal entries may differ, the trace does not depend on the basis used to represent the linear transformation.

````{prf:proof}
Use the cyclic property to move $S^{-1}$ to the end:

$$
\begin{aligned}
\operatorname{tr}(S^{-1}AS)
&=\operatorname{tr}(ASS^{-1}) \\
&=\operatorname{tr}(A).
\end{aligned}
$$
````

The later section on [applications of eigenvalues](applications_of_eigenvalues.md) uses this property to relate the trace to the sum of the eigenvalues.

## Useful Identities

### Outer Products

For column vectors $u,v\in\mathbb{R}^n$,

$$\operatorname{tr}(uv^T)=\sum_{i=1}^n u_iv_i=v^Tu.$$

The diagonal entries of $uv^T$ are $u_iv_i$, which gives the identity directly. For complex vectors, the corresponding identity with conjugation is

$$\operatorname{tr}(uv^H)=\sum_{i=1}^n u_i\overline{v_i}=v^Hu.$$

### The Frobenius Norm

For a possibly rectangular matrix $A\in\mathbb{C}^{m\times n}$,

$$
\begin{aligned}
\operatorname{tr}(A^HA)
&=\sum_{j=1}^n\sum_{i=1}^m |a_{ij}|^2 \\
&=\|A\|_F^2.
\end{aligned}
$$

Indeed, diagonal entry $j$ of $A^HA$ is the squared norm of column $j$ of $A$. Summing over the columns gives the squared Frobenius norm. For a real matrix, this reduces to $\operatorname{tr}(A^TA)=\|A\|_F^2$.

## Interpretation: Divergence of a Linear Vector Field

Let $A\in\mathbb{R}^{n\times n}$ and consider the vector field $F(x)=Ax$. The **divergence** is the sum of the partial derivatives of each component with respect to its matching coordinate:

$$
\nabla\cdot F=\sum_{i=1}^n\frac{\partial F_i}{\partial x_i}.
$$

Since $F_i(x)=\sum_{j=1}^n a_{ij}x_j$, we have $\partial F_i/\partial x_i=a_{ii}$. Therefore,

$$\nabla\cdot F=\operatorname{tr}(A).$$

If particles move with velocity $F(x)$, so that $\dot{x}=Ax$, the divergence gives the instantaneous rate of volume change per unit volume:

- A positive trace means that a volume of moving particles grows.
- A negative trace means that it shrinks.
- A zero trace means that volume is preserved; the flow is **divergence-free**.

This concerns volume, rather than expansion or contraction in every direction. For example,

$$A=\begin{pmatrix}2&0\\0&-1\end{pmatrix}$$

has positive trace, but the motion $\dot{x}=Ax$ expands in the first coordinate and contracts in the second. If the first diagonal entry is changed to $1$, the trace becomes zero: expansion in one coordinate and contraction in the other preserve area while changing shape.

These statements describe motion under the velocity field $F(x)=Ax$. For a single application of the linear map $x\mapsto Ax$, the volume scaling factor is $|\det(A)|$, as discussed in the preceding section.
