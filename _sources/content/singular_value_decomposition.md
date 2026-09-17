# Singular Value Decomposition

The singular value decomposition (SVD) describes a matrix through orthonormal input directions, orthonormal output directions, and the factors by which those directions are stretched. It exists for every matrix, including rectangular and rank-deficient matrices. This makes it useful for understanding linear transformations, solving least-squares problems, and approximating data by a matrix of lower rank.

## The Full SVD

````{prf:theorem} Singular Value Decomposition
:label: thm:svd
Every matrix $A\in\mathbb{R}^{m\times n}$ can be written as

$$
A=U\Sigma V^T,
$$

where:

- $U\in\mathbb{R}^{m\times m}$ is orthogonal. Its columns $u_1,\ldots,u_m$ are **left singular vectors**.
- $V\in\mathbb{R}^{n\times n}$ is orthogonal. Its columns $v_1,\ldots,v_n$ are **right singular vectors**.
- $\Sigma\in\mathbb{R}^{m\times n}$ is rectangular diagonal: its only possibly nonzero entries are $\Sigma_{ii}=\sigma_i$ for $1\leq i\leq p$, where $p=\min(m,n)$.

The **singular values** are ordered so that

$$
\sigma_1\geq\sigma_2\geq\cdots\geq\sigma_p\geq0.
$$
````

The factorization gives the paired relations

$$
Av_i=\sigma_i u_i,
\qquad
A^Tu_i=\sigma_i v_i,
\qquad 1\leq i\leq p.
$$

Thus, $A$ maps the input direction $v_i$ to the output direction $u_i$, with scale factor $\sigma_i$. An eigenvector describes a direction that stays on the same line under a square matrix. The SVD allows the input and output directions to differ, and even to belong to spaces of different dimensions.

For a complex matrix, the corresponding factorization is $A=U\Sigma V^H$, with $U$ and $V$ unitary and the same real, nonnegative singular values. We develop the real case below; the results extend to complex matrices by using conjugate transposes and complex inner products.

## Geometry: Directions and Stretching

Let $r$ be the number of positive singular values. Multiplying out the SVD gives the **outer-product expansion**

$$
A=\sum_{i=1}^r\sigma_i u_i v_i^T.
$$

Consequently, for any input $x$,

$$
Ax=\sum_{i=1}^r\sigma_i(v_i^Tx)u_i.
$$

The scalar $v_i^Tx$ is the component of $x$ along $v_i$. The matrix scales this component by $\sigma_i$ and sends it along $u_i$. Components along $v_{r+1},\ldots,v_n$ are sent to zero.

This gives a geometric interpretation of the three factors:

1. $V^T$ expresses the input in an orthonormal coordinate system.
2. $\Sigma$ scales the first $p$ coordinates. If $n>m$, it discards the remaining input coordinates; if $m>n$, it appends zeros.
3. $U$ expresses the result in the output coordinate system.

The image of the unit ball is an ellipsoid in $R(A)$ with **semiaxis lengths** $\sigma_1,\ldots,\sigma_r$ and directions $u_1,\ldots,u_r$. If $r<m$, this ellipsoid lies in a lower-dimensional subspace of the output space.

### A Rectangular Example

Consider

$$
A=\begin{pmatrix}
0&3&0&0\\
0&0&-2&0\\
0&0&0&0
\end{pmatrix},
\qquad
Ax=\begin{pmatrix}3x_2\\-2x_3\\0\end{pmatrix}.
$$

The second input coordinate is stretched by $3$, the third by $2$ with a reversal of direction, and the first and fourth are lost. Writing $e_i$ for standard basis vectors in the appropriate space, take

$$
U=[e_1,-e_2,e_3],
\qquad
V=[e_2,e_3,e_1,e_4],
$$

and

$$
\Sigma=\begin{pmatrix}
3&0&0&0\\
0&2&0&0\\
0&0&0&0
\end{pmatrix}.
$$

Then $A=U\Sigma V^T$. The singular values are $3,2,0$, and the rank is $2$. The unit ball in $\mathbb{R}^4$ maps to a filled ellipse in the first two output coordinates, with semiaxes $3$ and $2$.

## Why the SVD Always Exists

The proof uses the spectral theorem for real symmetric matrices, established in [Normal Matrices](normal_matrices.md).

````{prf:proof} Existence of the SVD.
The matrix $A^TA$ is symmetric and positive semidefinite because

$$
x^TA^TAx=\|Ax\|_2^2\geq0.
$$

It therefore has an orthonormal eigenbasis $v_1,\ldots,v_n$, with eigenvalues $\mu_1\geq\cdots\geq\mu_n\geq0$. Let $r$ be the number of positive eigenvalues. For $1\leq i\leq r$, define

$$
\sigma_i=\sqrt{\mu_i},
\qquad
u_i=\frac{Av_i}{\sigma_i}.
$$

These vectors are orthonormal, since

$$
u_i^Tu_j
=\frac{v_i^TA^TAv_j}{\sigma_i\sigma_j}
=\frac{\mu_j v_i^Tv_j}{\sigma_i\sigma_j}
=\delta_{ij}.
$$

In particular, $r\leq m$. Complete $u_1,\ldots,u_r$ to an orthonormal basis of $\mathbb{R}^m$ and form $U=[u_1,\ldots,u_m]$ and $V=[v_1,\ldots,v_n]$.

For $i>r$, we have

$$
\|Av_i\|_2^2=v_i^TA^TAv_i=\mu_i=0,
$$

so $Av_i=0$. Thus the columns of $AV$ are $\sigma_i u_i$ for $i\leq r$ and zero otherwise. Equivalently, $AV=U\Sigma$, where $\Sigma$ is rectangular diagonal with these $r$ positive entries and zeros elsewhere. Multiplying by $V^T$ gives $A=U\Sigma V^T$.

Since $U$ and $V$ are invertible, $\operatorname{rank}(A)=\operatorname{rank}(\Sigma)=r$. The argument also covers $A=0$: there are no positive singular values, and any orthonormal bases give an SVD with $\Sigma=0$.
````

This is an existence proof. Methods for computing the SVD will be developed later in the course.

## Singular Values and Eigenvalues

Multiplying the factors gives

$$
\begin{aligned}
A^TA&=V(\Sigma^T\Sigma)V^T,\\
AA^T&=U(\Sigma\Sigma^T)U^T.
\end{aligned}
$$

The matrices $\Sigma^T\Sigma$ and $\Sigma\Sigma^T$ are diagonal, but their sizes differ: they are $n\times n$ and $m\times m$, respectively. Both have positive eigenvalues $\sigma_1^2,\ldots,\sigma_r^2$. The remaining eigenvalues are zero: $n-r$ for $A^TA$ and $m-r$ for $AA^T$.

Thus the right singular vectors are eigenvectors of $A^TA$, and the left singular vectors are eigenvectors of $AA^T$. For a positive singular value, the two vectors are linked by $u_i=Av_i/\sigma_i$; they cannot be chosen independently.

The singular values are uniquely determined, but the singular vectors need not be. A positive pair $(u_i,v_i)$ may have both signs reversed. Within a repeated positive singular value, the left and right bases may undergo the same orthogonal change of basis. The bases for the two null spaces may be chosen independently.

### A Symmetric Block Matrix

There is also a direct connection with the symmetric matrix

$$
H=\begin{pmatrix}0&A\\A^T&0\end{pmatrix}.
$$

For every positive singular value, the paired relations give

$$
\begin{aligned}
H\frac{1}{\sqrt{2}}\begin{pmatrix}u_i\\v_i\end{pmatrix}
&=\sigma_i\frac{1}{\sqrt{2}}\begin{pmatrix}u_i\\v_i\end{pmatrix},\\
H\frac{1}{\sqrt{2}}\begin{pmatrix}u_i\\-v_i\end{pmatrix}
&=-\sigma_i\frac{1}{\sqrt{2}}\begin{pmatrix}u_i\\-v_i\end{pmatrix}.
\end{aligned}
$$

The remaining eigenvectors are

$$
\begin{pmatrix}u_i\\0\end{pmatrix}\quad(r<i\leq m),
\qquad
\begin{pmatrix}0\\v_j\end{pmatrix}\quad(r<j\leq n),
$$

all with eigenvalue zero. Together these form an orthonormal basis of $\mathbb{R}^{m+n}$. Hence the eigenvalues of $H$ are $\pm\sigma_1,\ldots,\pm\sigma_r$, together with $m+n-2r$ zeros. Including these zero directions is essential for rectangular or rank-deficient matrices.

## The Four Fundamental Subspaces

The SVD provides orthonormal bases for all [four fundamental subspaces](four_fundamental_subspaces.md):

| Subspace | Dimension | Orthonormal basis |
| :--- | :--- | :--- |
| Column space $R(A)$ | $r$ | $u_1,\ldots,u_r$ |
| Left null space $N(A^T)$ | $m-r$ | $u_{r+1},\ldots,u_m$ |
| Row space $R(A^T)$ | $r$ | $v_1,\ldots,v_r$ |
| Null space $N(A)$ | $n-r$ | $v_{r+1},\ldots,v_n$ |

An empty list is the basis of the zero subspace. To see the null-space statement directly, use the orthonormality of the $u_i$ to obtain

$$
\|Ax\|_2^2=\sum_{i=1}^r\sigma_i^2(v_i^Tx)^2.
$$

This is zero exactly when $x$ is orthogonal to $v_1,\ldots,v_r$. The column-space statement follows from the outer-product expansion and $Av_i=\sigma_i u_i$. Applying the same reasoning to $A^T$ gives the other two spaces.

Write $U_r=[u_1,\ldots,u_r]$ and $V_r=[v_1,\ldots,v_r]$. Then

$$
P_{R(A)}=U_rU_r^T,
\qquad
P_{R(A^T)}=V_rV_r^T
$$

are the **orthogonal projections** onto the column and row spaces. Their complementary projections, $I_m-U_rU_r^T$ and $I_n-V_rV_r^T$, project orthogonally onto the left null space and null space.

## Reduced, Compact, and Truncated Forms

The full SVD includes complete bases for both spaces. Often only some columns are needed. For $0\leq k\leq p$, let $U_k$ and $V_k$ contain the first $k$ columns, and let $\Sigma_k=\operatorname{diag}(\sigma_1,\ldots,\sigma_k)$.

- The **reduced SVD** (also called the economy SVD) keeps $p=\min(m,n)$ columns:

  $$A=U_p\Sigma_pV_p^T.$$

  The factor sizes are $m\times p$, $p\times p$, and $n\times p$ for $U_p$, $\Sigma_p$, and $V_p$. This factorization is exact, including when some singular values are zero.

- The **compact SVD** keeps only the $r$ positive singular values:

  $$A=U_r\Sigma_rV_r^T.$$

  It is also exact. When $r<p$, it omits additional zero singular directions. The term *thin SVD* is used for either of these forms in different references, so check the stated dimensions.

- A **truncated SVD** retains just the first $k<r$ terms:

  $$A_k=U_k\Sigma_kV_k^T=\sum_{i=1}^k\sigma_i u_i v_i^T.$$

  It has rank $k$ and approximates $A$. We set $A_0=0$; once $k\geq r$, the sum equals $A$ exactly.

Each retained term describes one input-output direction pair. Storing the factors requires $k(m+n+1)$ numbers instead of $mn$, so a small $k$ can reduce storage substantially.

## Matrix Norms from Singular Values

The largest singular value measures the maximum stretching of a unit vector:

$$
\|A\|_2=\max_{\|x\|_2=1}\|Ax\|_2=\sigma_1.
$$

Indeed, the expression for $\|Ax\|_2^2$ above is at most $\sigma_1^2\|x\|_2^2$, with equality at $x=v_1$. This also holds for $A=0$.

The Frobenius norm combines the stretching in all directions. Using the trace identity from [Trace](trace.md),

$$
\|A\|_F^2
=\operatorname{tr}(A^TA)
=\sum_{i=1}^p\sigma_i^2.
$$

In the rectangular example, $\|A\|_2=3$ whereas $\|A\|_F=\sqrt{13}$. These norms measure different things.

The [Schatten norms](operator_norms.md) apply a vector norm to the singular values. Using $q$ for the norm exponent,

$$
\|A\|_{S,q}=\left(\sum_{i=1}^p\sigma_i^q\right)^{1/q},
\qquad 1\leq q<\infty.
$$

The important special cases are

$$
\begin{aligned}
\|A\|_{S,1}&=\sum_{i=1}^p\sigma_i=\|A\|_* &&\text{(nuclear norm)},\\
\|A\|_{S,2}&=\|A\|_F &&\text{(Frobenius norm)},\\
\|A\|_{S,\infty}&=\sigma_1=\|A\|_2 &&\text{(spectral norm)}.
\end{aligned}
$$

The subscript $S$ matters: the Schatten $2$-norm is the Frobenius norm, while the induced matrix $2$-norm is the spectral norm. Likewise, the nuclear norm is generally different from the induced matrix $1$-norm.

## Best Low-Rank Approximation

Keeping the largest singular values gives the best approximation of a prescribed rank in two common norms.

````{prf:theorem} Eckart–Young–Mirsky Theorem
:label: thm:eckart-young-mirsky
Let $A_k=\sum_{i=1}^k\sigma_i u_i v_i^T$, with $0\leq k<p$. Among all matrices $B\in\mathbb{R}^{m\times n}$ of rank at most $k$, $A_k$ minimizes both the spectral-norm error and the Frobenius-norm error:

$$
\begin{aligned}
\min_{\operatorname{rank}(B)\leq k}\|A-B\|_2
&=\|A-A_k\|_2=\sigma_{k+1},\\
\min_{\operatorname{rank}(B)\leq k}\|A-B\|_F
&=\|A-A_k\|_F
=\left(\sum_{i=k+1}^p\sigma_i^2\right)^{1/2}.
\end{aligned}
$$

For $k=p$, $A_p=A$ and both errors are zero.
````

````{prf:proof} Spectral-Norm Error.
Let $B$ have rank at most $k$. The restriction of $B$ to the $(k+1)$-dimensional space $\operatorname{span}(v_1,\ldots,v_{k+1})$ has a nonzero null vector. Normalize it to obtain a unit vector $z=\sum_{i=1}^{k+1}c_iv_i$ with $Bz=0$. Then

$$
\begin{aligned}
\|A-B\|_2^2
&\geq\|(A-B)z\|_2^2=\|Az\|_2^2\\
&=\sum_{i=1}^{k+1}\sigma_i^2c_i^2
\geq\sigma_{k+1}^2.
\end{aligned}
$$

On the other hand, $A-A_k=\sum_{i=k+1}^p\sigma_i u_i v_i^T$ has largest singular value $\sigma_{k+1}$. Thus $A_k$ attains the bound.
````

````{prf:proof} Frobenius-Norm Error.
Let $B$ have rank at most $k$. Its null space has dimension at least $n-k$, so choose orthonormal vectors $z_1,\ldots,z_{n-k}$ in $N(B)$ and complete them to an orthonormal basis of $\mathbb{R}^n$. The Frobenius norm is unchanged by an orthogonal change of coordinates, so its square is the sum of squared output lengths over this basis. Hence

$$
\|A-B\|_F^2
\geq\sum_{j=1}^{n-k}\|(A-B)z_j\|_2^2
=\sum_{j=1}^{n-k}\|Az_j\|_2^2.
$$

Set $d_i=\sigma_i^2$ for $i\leq p$ and $d_i=0$ for $p<i\leq n$. Expanding each $z_j$ in the right singular basis gives

$$
\sum_{j=1}^{n-k}\|Az_j\|_2^2
=\sum_{i=1}^n d_i w_i,
\qquad
w_i=\sum_{j=1}^{n-k}(v_i^Tz_j)^2.
$$

Here $0\leq w_i\leq1$: $w_i$ is the squared length of the orthogonal projection of the unit vector $v_i$ onto $\operatorname{span}(z_1,\ldots,z_{n-k})$. Also $\sum_iw_i=n-k$, because each $z_j$ has unit length.

Since $d_1\geq\cdots\geq d_n\geq0$, these weights give a sum no smaller than the sum of the $n-k$ smallest $d_i$. Explicitly,

$$
\begin{aligned}
\sum_{i=1}^n d_iw_i-\sum_{i=k+1}^n d_i
&=\sum_{i=1}^k d_iw_i-\sum_{i=k+1}^n d_i(1-w_i)\\
&\geq d_{k+1}\left(\sum_{i=1}^k w_i
-\sum_{i=k+1}^n(1-w_i)\right)\\
&=0.
\end{aligned}
$$

Therefore $\|A-B\|_F^2\geq\sum_{i=k+1}^p\sigma_i^2$. The discarded terms of the SVD give exactly this squared error for $A_k$, so the bound is attained.
````

In the rectangular example, the best rank-one approximation retains the entry $3$ and replaces the entry $-2$ by zero. Both approximation errors equal $2$. A rank-two approximation is already exact.

These formulas let us choose how many directions to retain for a desired error. A small discarded singular value means a small contribution to this matrix approximation; it does not by itself mean that the corresponding information is noise.

## Determinant and Volume

````{prf:theorem} Determinant and Singular Values
:label: thm:det-svd
For a square real or complex matrix $A\in\mathbb{C}^{n\times n}$,

$$
|\det(A)|=\prod_{i=1}^n\sigma_i.
$$
````

````{prf:proof} Determinant Formula.
From $A=U\Sigma V^H$ and the multiplicative property of determinants,

$$
|\det(A)|=|\det(U)|\,|\det(\Sigma)|\,|\det(V^H)|.
$$

The unitary factors have determinants of magnitude one, and $\det(\Sigma)=\prod_i\sigma_i$. This gives the result.
````

For a real square matrix, $|\det(A)|$ is the volume scaling factor: it is the product of the stretches in the singular directions. The singular values do not record the sign of the determinant, which describes orientation. In particular, a square matrix is invertible exactly when all its singular values are positive.

## Connections to Later Topics

The SVD will recur throughout the course:

- **Least squares:** in singular coordinates, the equations separate into scalar equations $\sigma_i y_i=c_i$. Positive singular values can be inverted, while zero singular values identify directions the matrix cannot determine. We will develop the pseudoinverse and minimum-norm solutions in the least-squares chapter.
- **Principal component analysis:** if the rows of $Z\in\mathbb{R}^{N\times n}$ are centered observations, the covariance matrix is $C=Z^TZ/N$. Its principal directions are the right singular vectors of $Z$, with variances $\sigma_i^2/N$ (and additional zeros when needed). This connects the SVD to the PCA discussion in [Applications of Eigenvalues](applications_of_eigenvalues.md).
- **Data compression:** a truncated SVD represents a matrix using a few direction pairs. The approximation theorem gives the error incurred by discarding the remaining pairs.

The algorithms for computing these factorizations and using them in numerical problems will be studied later.
