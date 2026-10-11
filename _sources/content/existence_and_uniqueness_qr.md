# Existence and Uniqueness of QR Factorization

The preceding algorithms construct QR in different ways. To understand whether they can produce the same factors, we need to separate three questions: whether QR exists, whether the signs of the basis vectors are fixed, and whether unused orthogonal directions remain to be chosen.

All statements here concern exact arithmetic.

## Existence

:::{prf:theorem} Existence of QR Factorization
:label: thm:existence_qr
Every real matrix $A\in\mathbb{R}^{m\times n}$ has a full QR factorization

$$
A=QR,\qquad Q\in\mathbb{R}^{m\times m},\qquad Q^TQ=I_m,
$$

where $R\in\mathbb{R}^{m\times n}$ is upper trapezoidal: $r_{ij}=0$ whenever $i>j$. No rank assumption is required.
:::

:::{prf:proof} Construction by Householder Reflections.
At step $k$, apply a Householder reflection to rows $k,\ldots,m$ to eliminate the entries below the diagonal in column $k$. If those entries are already zero, use the identity transformation.

This step preserves earlier columns: they are already zero in rows $k,\ldots,m$, so combining those rows does not change their zeros. Continue for $s=\min(n,m-1)$ steps. Writing $P_k$ for the full $m\times m$ transformation at step $k$ gives

$$
R=P_s\cdots P_1A.
$$

Each $P_k$ is orthogonal. Therefore, defining

$$
Q=P_1^T\cdots P_s^T
$$

gives $Q^TQ=I_m$ and $A=QR$. If no elimination is needed, the empty product is the identity. The construction never requires division by a matrix pivot or by the norm of a zero column.
:::

When $m\ge n$, partition the factors as

$$
Q=\begin{pmatrix}Q_1&Q_2\end{pmatrix},
\qquad
R=\begin{pmatrix}R_1\\0\end{pmatrix}.
$$

Then

$$
A=Q_1R_1,\qquad
Q_1\in\mathbb{R}^{m\times n},\qquad Q_1^TQ_1=I_n,
$$

with $R_1\in\mathbb{R}^{n\times n}$ upper triangular. This is the **reduced (thin) QR factorization**, which also exists for rank-deficient matrices. Moreover, $\operatorname{rank}(A)=\operatorname{rank}(R_1)$: multiplication by $Q_1$ preserves rank because $Q_1$ has a left inverse $Q_1^T$.

Basic Gram-Schmidt, however, stops when a residual is zero. Continuing that construction requires choosing a new unit vector perpendicular to the previous basis vectors. Such a direction exists because at step $k\le n\le m$, only $k-1<m$ orthonormal vectors have been chosen. Thus existence of QR does not imply that basic Gram-Schmidt completes without modification.

## Fixing the Signs

Even for independent columns, QR is not unique unless we fix the signs. If $D$ is diagonal with entries $d_i\in\{-1,1\}$, then

$$
A=Q_1R_1=(Q_1D)(DR_1).
$$

The columns of $Q_1D$ remain orthonormal, and $DR_1$ remains upper triangular. Changing $d_i$ reverses both column $i$ of $Q_1$ and row $i$ of $R_1$.

For full column rank, $R_1$ is invertible, so none of its diagonal entries is zero. We can fix these sign choices by requiring $r_{ii}>0$.

## Uniqueness for Full Column Rank

:::{prf:theorem} Unique Reduced QR Factorization
:label: thm:unique_thin_qr
Let $A\in\mathbb{R}^{m\times n}$ with $m\ge n$ have full column rank. There is exactly one reduced QR factorization $A=Q_1R_1$ with $Q_1^TQ_1=I_n$ and $R_1$ upper triangular with positive diagonal.
:::

:::{prf:proof} Uniqueness Through Cholesky Factorization.
The sign adjustment above shows that a factorization with positive diagonal exists. To prove uniqueness, take any such factorization. Then

$$
A^TA=R_1^TQ_1^TQ_1R_1=R_1^TR_1.
$$

The matrix $A^TA$ is positive definite: for every nonzero $x\in\mathbb{R}^n$, full column rank gives $Ax\ne0$, so

$$
x^TA^TAx=\|Ax\|_2^2>0.
$$

Thus $R_1$ is the upper triangular Cholesky factor of $A^TA$. The [uniqueness of Cholesky factorization](cholesky.md) with positive diagonal determines $R_1$ uniquely.

Since $R_1$ is invertible,

$$
Q_1=AR_1^{-1},
$$

so $Q_1$ is uniquely determined as well.
:::

Without the positive-diagonal convention, a full-column-rank real matrix has exactly $2^n$ reduced QR factorizations, corresponding to the sign choices in $D$.

### The Remaining Freedom in Full QR

For a full-column-rank matrix, fixing the positive diagonal determines $Q_1$ and $R_1$. If $m>n$, the remaining columns $Q_2$ form an orthonormal basis for the orthogonal complement of the column space of $A$.

Changing this basis does not change $A$. For any orthogonal $W\in\mathbb{R}^{(m-n)\times(m-n)}$,

$$
A=
\begin{pmatrix}Q_1&Q_2W\end{pmatrix}
\begin{pmatrix}R_1\\0\end{pmatrix}.
$$

Every other orthonormal basis of this complement has the form $Q_2W$. Consequently, with the positive-diagonal convention, full QR is unique when $m=n$, has exactly two choices when $m-n=1$, and has infinitely many choices when $m-n\ge2$.

## Rank-Deficient Matrices

Suppose $m\ge n$ and $A$ has dependent columns. The reduced factor $R_1$ is singular, so at least one diagonal entry $r_{kk}$ is zero. Strict positivity of every diagonal entry is therefore impossible.

Even requiring a **nonnegative** diagonal does not ensure uniqueness. Choose $D$ with $d_k=-1$ at a zero diagonal entry and $d_i=1$ elsewhere. The factorization $(Q_1D)(DR_1)$ has the same nonnegative diagonal as $Q_1R_1$, but a different $Q_1$.

For example,

$$
\begin{pmatrix}1&0\\0&0\end{pmatrix}
=
\begin{pmatrix}1&0\\0&1\end{pmatrix}
\begin{pmatrix}1&0\\0&0\end{pmatrix}
=
\begin{pmatrix}1&0\\0&-1\end{pmatrix}
\begin{pmatrix}1&0\\0&0\end{pmatrix}.
$$

With a nonnegative diagonal, the first column forces $r_{11}=1$ and $q_1=e_1$. The zero second column then forces $r_{12}=r_{22}=0$, since the columns of $Q_1$ are orthonormal. Thus $R_1$ is unique in this example, but $Q_1$ is not.

For complex matrices, replace orthogonal by unitary and $Q^T$ by $Q^H$. Requiring positive real diagonal entries fixes the phase choices and gives the same reduced-QR uniqueness result for full column rank.
