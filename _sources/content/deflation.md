# The Method of Deflation

The [power method](power_method.md) finds the dominant direction. Repeating it from a different starting vector will usually find the same direction again. To find another one, we remove the component along the direction already known **after every multiplication by $A$**. This is **deflation**.

We seek a new direction perpendicular to the known one. For a nonsymmetric matrix, eigenvectors need not be perpendicular, so this new direction need not be an eigenvector of $A$. It will instead be an eigenvector of the projected operator, with the next eigenvalue of $A$.

## A Two-Dimensional Example

Consider

$$
A=\begin{pmatrix}4&1\\0&2\end{pmatrix}.
$$

Its eigenvalues are $\lambda_1=4$ and $\lambda_2=2$. Corresponding eigenvectors are $(1,0)^T$ and $x_2=(-1,2)^T$.

Suppose the power method has found the first unit direction exactly:

$$
q_{\star,1}=\begin{pmatrix}1\\0\end{pmatrix}.
$$

The star marks an exact basis direction, to distinguish it from numerical iterates. The orthogonal projection onto the perpendicular direction is

$$
P_1=I-q_{\star,1}q_{\star,1}^H
=\begin{pmatrix}0&0\\0&1\end{pmatrix}.
$$

Projecting $x_2$ and normalizing gives a perpendicular unit vector:

$$
P_1x_2=\begin{pmatrix}0\\2\end{pmatrix},
\qquad
q_{\star,2}
=\frac{P_1x_2}{\|P_1x_2\|_2}
=\begin{pmatrix}0\\1\end{pmatrix}.
$$

This changes the direction: $q_{\star,2}$ is not an eigenvector of $A$. Multiplication by $A$ brings back a component along the first direction:

$$
Aq_{\star,2}
=\begin{pmatrix}1\\2\end{pmatrix}
=q_{\star,1}+2q_{\star,2}.
$$

Projection removes that component:

$$
P_1Aq_{\star,2}=2q_{\star,2}.
$$

Thus $q_{\star,2}$ is an eigenvector of $P_1A$ with eigenvalue $2$, the second eigenvalue of $A$. Starting perpendicular to the first direction was not enough: we had to project again after multiplication.

:::{figure} ../_static/deflation_geometry.svg
:alt: Two diagrams illustrate deflation for the matrix with rows (4,1) and (0,2). On the left, the blue eigenvector x2=(-1,2) is projected horizontally onto the vertical complement of the known direction q-star-1. Its purple projection (0,2) is normalized to q-star-2=(0,1). On the right, A sends q-star-2 to the blue vector (1,2). An orange horizontal projection removes the known component and leaves the purple vector 2 q-star-2=(0,2).
:width: 820px
:name: fig-deflation-geometry

Left: projection changes the eigenvector's direction; normalization gives $q_{\star,2}$. Right: applying $A$ to this new direction creates a component along the known horizontal direction. Removing it leaves $2q_{\star,2}$. The eigenvalue $2$ is preserved even though the eigenvector has changed.
:::

## Why Projection Preserves the Eigenvalue

Suppose we have found $i$ orthonormal directions. Collect them in

$$
Q_{\star,i}=[q_{\star,1},\ldots,q_{\star,i}],
\qquad
\mathcal S_i=\operatorname{range}(Q_{\star,i}).
$$

Assume that $\mathcal S_i$ is **invariant under $A$**: applying $A$ to any vector in $\mathcal S_i$ keeps it in $\mathcal S_i$. This holds for the first direction because it is an eigenvector. The construction below will preserve this property as we add directions.

Let $P_i=I-Q_{\star,i}Q_{\star,i}^H$ be the orthogonal projection onto $\mathcal S_i^\perp$. Take an eigenvector $x$ with eigenvalue $\lambda$ that is not in $\mathcal S_i$. Split it into its known component and its perpendicular remainder:

$$
x=s+v,
\qquad
s=(I-P_i)x\in\mathcal S_i,
\qquad
v=P_ix\ne0.
$$

Because the known subspace is invariant, $As\in\mathcal S_i$, so $P_iAs=0$. Consequently,

$$
\begin{aligned}
P_iAv
&=P_iA(x-s)\\
&=P_iAx-P_iAs\\
&=\lambda P_ix-0\\
&=\lambda v.
\end{aligned}
$$

The perpendicular remainder is therefore an eigenvector of $P_iA$ with the **same eigenvalue**. Although $Av$ can have a component in the known subspace, its component perpendicular to that subspace is exactly $\lambda v$.

Invariance also gives

$$
P_iA(I-P_i)=0,
\qquad
P_iA=P_iAP_i.
$$

The first identity says that $P_iA$ sends every vector in the known subspace $\mathcal S_i$ to zero.

## Why the New Directions Are Schur Vectors

For the sequential construction, order the eigenvalues so that

$$
|\lambda_1|>|\lambda_2|>\cdots>|\lambda_n|>0.
$$

They are distinct, so their eigenvectors $x_1,\ldots,x_n$ are independent. Suppose that $\mathcal S_i$ spans the first $i$ eigenvectors. Projecting the next one gives

$$
q_{\star,i+1}
=\frac{P_ix_{i+1}}{\|P_ix_{i+1}\|_2},
\qquad
P_iAq_{\star,i+1}
=\lambda_{i+1}q_{\star,i+1}.
$$

The denominator is nonzero because $x_{i+1}\notin\mathcal S_i$. This describes the direction we want; the algorithm below finds it without knowing $x_{i+1}$.

Since $P_i$ retains the perpendicular component, the last identity means

$$
Aq_{\star,i+1}
=\lambda_{i+1}q_{\star,i+1}
+s_{i+1},
\qquad
s_{i+1}\in\mathcal S_i.
$$

Therefore, adding this direction gives a larger invariant subspace. It is also the span of $x_1,\ldots,x_{i+1}$, because we have only subtracted a vector in $\mathcal S_i$ from $x_{i+1}$ and normalized.

Starting with the first eigenvector and repeating this construction produces an orthonormal basis. In that basis, each column satisfies

$$
Aq_{\star,j}
=\sum_{\ell=1}^{j-1}t_{\ell j}q_{\star,\ell}
+\lambda_jq_{\star,j}.
$$

Only earlier directions and the current direction occur. Collecting the columns gives

$$
AQ_\star=Q_\star T,
\qquad
A=Q_\star TQ_\star^H,
$$

where $Q_\star=[q_{\star,1},\ldots,q_{\star,n}]$ is unitary and $T$ is upper triangular. This is the **Schur decomposition**.

The off-diagonal entries of $T$ record the components along earlier directions. They explain why a Schur vector need not be an eigenvector of $A$. The diagonal entries still give the eigenvalues of $A$, and orthogonality gives

$$
\lambda_j=q_{\star,j}^HAq_{\star,j}.
$$

For Hermitian matrices, the eigenvectors can be chosen orthonormal, and these Schur directions are eigenvectors of $A$ itself.

## The Spectrum After Deflation

The geometric argument identifies the next direction. We must also check that it is dominant for the projected operator, so that power iteration can find it. The following result gives the full spectrum after projection. It holds for any ordered Schur decomposition, without the strict-gap or nonzero-eigenvalue assumptions above.

:::{prf:theorem} Exact Schur Deflation
:label: thm:exact-schur-deflation

For $1\le i<n$, let $P_i=I-Q_{\star,i}Q_{\star,i}^H$ and $M_i=P_iA$. The eigenvalues of $M_i$ are $i$ zeros followed by $\lambda_{i+1},\ldots,\lambda_n$, counted with multiplicity. Moreover,

$$
M_iq_{\star,i+1}
=\lambda_{i+1}q_{\star,i+1}.
$$
:::

:::{prf:proof} Deflation in Schur Coordinates.
Partition the Schur form after its first $i$ rows and columns:

$$
T=\begin{pmatrix}T_{11}&T_{12}\\0&T_{22}\end{pmatrix}.
$$

In Schur coordinates, $P_i$ sets the first $i$ components to zero:

$$
\begin{aligned}
Q_\star^HP_iQ_\star
&=\begin{pmatrix}0&0\\0&I_{n-i}\end{pmatrix},\\
Q_\star^HM_iQ_\star
&=\begin{pmatrix}0&0\\0&I_{n-i}\end{pmatrix}
  \begin{pmatrix}T_{11}&T_{12}\\0&T_{22}\end{pmatrix}\\
&=\begin{pmatrix}0&0\\0&T_{22}\end{pmatrix}.
\end{aligned}
$$

This upper triangular matrix has the claimed diagonal entries. The first column of $T_{22}$ has $\lambda_{i+1}$ at the top and zeros below. Thus its first coordinate direction is an eigenvector; transforming back to the original coordinates gives $q_{\star,i+1}$.
:::

## Power Iteration on the Remaining Subspace

Under the strict magnitude ordering above, $\lambda_{i+1}$ is the unique dominant eigenvalue of $M_i$. Choose a random vector $g$ and initialize

$$
v_0=\frac{P_ig}{\|P_ig\|_2}.
$$

Redraw if the projection is zero. For fixed $i$, repeat

$$
\begin{aligned}
w_{\ell+1}
&=Av_\ell
-Q_{\star,i}\bigl(Q_{\star,i}^HAv_\ell\bigr),\\
v_{\ell+1}
&=\frac{w_{\ell+1}}{\|w_{\ell+1}\|_2}.
\end{aligned}
$$

Each step multiplies by $A$, subtracts the component in the known subspace, and normalizes. There is no need to form either $P_i$ or $M_i$. Here $i$ counts known directions and $\ell$ counts power iterations for the next direction.

As in the power method, $v_0$ must have a nonzero coefficient in the dominant direction when expanded in an eigenvector basis of $M_i$. Projecting a random vector drawn from a continuous distribution satisfies this condition, and gives a nonzero initial projection, with probability one in exact arithmetic.

Then the line spanned by $v_\ell$ converges to the line spanned by $q_{\star,i+1}$. For $i\le n-2$, the geometric convergence factor is governed by $|\lambda_{i+2}/\lambda_{i+1}|$. After $n-1$ directions have been found, their one-dimensional orthogonal complement supplies the last direction without further power iteration. Signs or complex phases do not affect these subspaces.

This construction assumes the preceding directions are exact. In computation, they are still approximate. [Orthogonal iteration](orthogonal_iteration.md) updates all directions together, using QR to perform the successive projections. There, $Q_k$ denotes the basis at iteration $k$, while $Q_\star$ remains the exact Schur target. Equal eigenvalue magnitudes may require following several directions as a subspace rather than one direction at a time.
