# Orthogonal Matrices

Orthogonal matrices describe linear transformations of real Euclidean space that preserve lengths and angles. Rotations and reflections are examples. These transformations fix the origin; translations are not represented by orthogonal matrices.

## Definition and Core Properties

A real square matrix $Q\in\mathbb{R}^{n\times n}$ is **orthogonal** if

$$Q^TQ=I_n.$$

Writing its columns as $q_1,\ldots,q_n$, we have $(Q^TQ)_{ij}=q_i^Tq_j$. Thus the definition says exactly that the columns are **orthonormal**: each column has length one, and distinct columns are perpendicular.

The columns are therefore linearly independent, so $Q$ is invertible. Multiplying $Q^TQ=I_n$ on the right by $Q^{-1}$ gives

$$Q^{-1}=Q^T,\qquad QQ^T=I_n.$$

The second identity says that the rows are orthonormal as well. Conversely, $QQ^T=I_n$ implies the same properties, so for square matrices either identity is sufficient.

### Rectangular Matrices with Orthonormal Columns

A matrix $Q\in\mathbb{R}^{m\times n}$ with $m>n$ can also satisfy $Q^TQ=I_n$. We will call this a **matrix with orthonormal columns**, reserving *orthogonal matrix* for the square case. Such a matrix preserves lengths when mapping $\mathbb{R}^n$ into $\mathbb{R}^m$, but it has no two-sided inverse and $QQ^T\neq I_m$.

For example,

$$
Q=\begin{pmatrix}1&0\\0&1\\0&0\end{pmatrix}
\quad\Longrightarrow\quad
Q^TQ=I_2,\qquad
QQ^T=\begin{pmatrix}1&0&0\\0&1&0\\0&0&0\end{pmatrix}.
$$

In general, $QQ^T$ is the orthogonal projection onto the column space of $Q$. It leaves vectors in that space unchanged and sends vectors perpendicular to it to zero. We will discuss this in the next section on [projections](projections.md).

### Unitary Matrices

The complex counterpart is a **unitary matrix**: a square matrix $Q\in\mathbb{C}^{n\times n}$ satisfying

$$Q^HQ=QQ^H=I_n,\qquad Q^{-1}=Q^H,$$

where $Q^H$ denotes the conjugate transpose. The conjugation is essential for preserving complex inner products and lengths.

## Key Results and Mathematical Formulas

### Length, Distance, and Angle Preservation

For real vectors $x,y$ and an orthogonal matrix $Q$,

$$
\begin{aligned}
(Qx)^T(Qy)
&=x^TQ^TQy \\
&=x^Ty.
\end{aligned}
$$

Thus inner products are preserved. Taking $y=x$ gives

$$\|Qx\|_2^2=x^Tx=\|x\|_2^2,$$

so $\|Qx\|_2=\|x\|_2$. Applying this to $x-y$ also gives

$$\|Qx-Qy\|_2=\|x-y\|_2.$$

Because both inner products and lengths are unchanged, the angle between any two nonzero vectors is unchanged. A distance-preserving map is called an **isometry**. These arguments also apply to rectangular matrices with orthonormal columns. For unitary matrices, replace each transpose by a conjugate transpose.

In particular, the operator 2-norm of an orthogonal or unitary matrix is $1$:

$$\|Q\|_2=\max_{x\neq 0}\frac{\|Qx\|_2}{\|x\|_2}=1.$$

### Products and Inverses

The product of two orthogonal matrices of the same size is orthogonal:

$$
\begin{aligned}
(Q_1Q_2)^T(Q_1Q_2)
&=Q_2^TQ_1^TQ_1Q_2 \\
&=Q_2^TQ_2=I_n.
\end{aligned}
$$

The inverse $Q^{-1}=Q^T$ is also orthogonal, since $QQ^T=I_n$. Together with the identity matrix, these properties make the orthogonal matrices a group, denoted $O(n)$. The corresponding statements hold for unitary matrices using conjugate transposes.

### Determinant and Orientation

For a real orthogonal matrix,

$$
\begin{aligned}
1=\det(Q^TQ)
&=\det(Q^T)\det(Q) \\
&=\det(Q)^2.
\end{aligned}
$$

Hence $\det(Q)=\pm1$. Both cases preserve volume, because $|\det(Q)|=1$.

- If $\det(Q)=1$, the transformation preserves orientation. These matrices form the **special orthogonal group**, denoted $SO(n)$. In two dimensions they are rotations about the origin; in three dimensions they are rotations about an axis through the origin.
- If $\det(Q)=-1$, the transformation reverses orientation. It can be a reflection, but need not be a single reflection across a hyperplane. For example, $-I_3$ reverses all three coordinate directions and has determinant $-1$; a reflection across a plane leaves that plane fixed.

For a complex unitary matrix, the corresponding conclusion is

$$
1=\det(Q^HQ)
=\overline{\det(Q)}\det(Q)
=|\det(Q)|^2.
$$

Thus $|\det(Q)|=1$, although the determinant need not be real.

## Why Orthogonal Matrices Are Useful in Numerical Linear Algebra

Orthogonal transformations preserve the size of an error already present in a vector. If $x$ is replaced by $x+e$, then, in exact arithmetic,

$$\|Q(x+e)-Qx\|_2=\|Qe\|_2=\|e\|_2.$$

Computing the transformation can introduce additional rounding errors, so this identity alone does not guarantee the accuracy of an algorithm. We will study those questions in the next chapter, [Solving Linear Systems](solving_linear_systems.md).

Orthogonal matrices also make changes of coordinates convenient. Applying the inverse means applying $Q^T$, without computing a general matrix inverse. For a square matrix $A$, changing to an orthonormal basis gives $Q^TAQ$, which is similar to $A$. Later chapters use orthogonal transformations in QR factorization, the singular value decomposition, and eigenvalue algorithms.

## Reflections Across Hyperplanes

Let $w\in\mathbb{R}^n$ be nonzero. The reflection across the hyperplane perpendicular to $w$ is

$$H=I_n-2\frac{ww^T}{w^Tw}.$$

To see its action, write $x=x_{\perp}+\alpha w$, where $w^Tx_{\perp}=0$. Then

$$Hx=x_{\perp}-\alpha w.$$

Thus $H$ leaves the hyperplane fixed and reverses the perpendicular component. Applying it twice returns the original vector, so $H^2=I_n$. The formula also gives $H^T=H$, and therefore

$$H^TH=H^2=I_n.$$

Hence $H$ is orthogonal and $H^{-1}=H$. In an orthonormal basis consisting of $w/\|w\|_2$ and vectors perpendicular to $w$, its matrix is $\operatorname{diag}(-1,1,\ldots,1)$, so $\det(H)=-1$.

These are **Householder reflections**. Their computational use is developed later in the section on [Householder transformations](householder_reflections.md).

## The Cartan–Dieudonné Theorem

The [Cartan–Dieudonné theorem](https://www.cis.upenn.edu/~jean/math-deep.pdf) shows that hyperplane reflections suffice to construct every real orthogonal transformation.

````{prf:theorem} Cartan–Dieudonné theorem
:label: thm:cartan_dieudonne
Every orthogonal transformation of $\mathbb{R}^n$ is a composition of at most $n$ reflections across hyperplanes through the origin. Equivalently,

$$Q=H_1H_2\cdots H_k,\qquad k\leq n,$$

where each $H_i$ is a hyperplane reflection. For $k=0$, the product is the identity matrix.
````

### Proof by Induction

````{prf:proof}
For $n=1$, an orthogonal matrix is either $[1]$ or $[-1]$. These require zero or one reflection, respectively.

Assume the theorem holds in dimension $n-1$, and let $Q$ be orthogonal on $\mathbb{R}^n$. Write $e_1=(1,0,\ldots,0)^T$ and $v=Qe_1$. Length preservation gives $\|v\|_2=1$.

**Case 1: $v=e_1$.** Then $Q$ fixes $e_1$. The subspace

$$S=\{x\in\mathbb{R}^n:e_1^Tx=0\}$$

is preserved by $Q$, because for $x\in S$,

$$e_1^TQx=(Qe_1)^TQx=e_1^Tx=0.$$

The restriction of $Q$ to $S$ is therefore an orthogonal transformation of an $(n-1)$-dimensional Euclidean space. By induction, it is a product of at most $n-1$ reflections on $S$.

Each such reflection has a normal vector $w\in S$. Using the same formula $I_n-2ww^T/(w^Tw)$ extends it to a reflection on $\mathbb{R}^n$ that fixes $e_1$. The product of these extended reflections agrees with $Q$ on both $S$ and $e_1$, so it equals $Q$ on all of $\mathbb{R}^n$.

**Case 2: $v\neq e_1$.** Set $w=v-e_1$ and define

$$H=I_n-2\frac{ww^T}{w^Tw}.$$

We verify that $Hv=e_1$. Since $v^Tv=1$, writing $v_1=e_1^Tv$ gives

$$
\begin{aligned}
w^Tv&=1-v_1, \\
w^Tw&=2(1-v_1)>0.
\end{aligned}
$$

Consequently,

$$Hv=v-2w\frac{1-v_1}{2(1-v_1)}=v-w=e_1.$$

The matrix $HQ$ is orthogonal and fixes $e_1$. By Case 1, it is a product of at most $n-1$ reflections. Since $H^{-1}=H$, multiplying by $H$ expresses $Q$ as a product of at most $n$ reflections.
````

Because each reflection has determinant $-1$, a product of $k$ reflections has determinant $(-1)^k$. Thus orientation-preserving orthogonal transformations require an even number of reflections, and orientation-reversing ones require an odd number.
