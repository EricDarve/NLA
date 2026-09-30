# The LU Decomposition Algorithm

An **LU factorization** writes a square matrix as

$$
A=LU,
$$

where $L$ is lower triangular and $U$ is upper triangular. We use the convention $l_{ii}=1$, so $L$ is **unit lower triangular**. This is also called the Doolittle convention.

Once the factors are available, solving $Ax=b$ requires two triangular solves:

$$
Ly=b,\qquad Ux=y.
$$

We first explain these solves, then derive LU as a sequence of rank-one updates. Throughout this section, we factor without interchanging rows. This requires nonzero pivots at the elimination steps; an invertible matrix need not satisfy this requirement.

## Solving Triangular Systems

### Forward Substitution

For a lower triangular system $Ly=b$, the first equation determines $y_1$, the second determines $y_2$, and so on. For example,

$$
\begin{pmatrix}
l_{11}&0&0\\
l_{21}&l_{22}&0\\
l_{31}&l_{32}&l_{33}
\end{pmatrix}
\begin{pmatrix}y_1\\y_2\\y_3\end{pmatrix}
=\begin{pmatrix}b_1\\b_2\\b_3\end{pmatrix}
$$

gives

$$
\begin{aligned}
y_1&=b_1/l_{11},\\
y_2&=(b_2-l_{21}y_1)/l_{22},\\
y_3&=(b_3-l_{31}y_1-l_{32}y_2)/l_{33}.
\end{aligned}
$$

The general formula is

$$
y_i=\frac{b_i-\sum_{j=1}^{i-1}l_{ij}y_j}{l_{ii}},
\qquad i=1,\ldots,n.
$$

All entries on the right have already been computed. We need $l_{ii}\neq0$; for the unit lower triangular factor in LU, every denominator is one.

### Backward Substitution

For an upper triangular system $Ux=y$, start with the last equation and work upward:

$$
x_i=\frac{y_i-\sum_{j=i+1}^n u_{ij}x_j}{u_{ii}},
\qquad i=n,n-1,\ldots,1.
$$

This requires $u_{ii}\neq0$. In both formulas, an empty sum is zero. The order of computation proves correctness: each step solves one equation whose other unknowns have already been determined.

The following functions assume square triangular NumPy arrays and matching one-dimensional right-hand sides. The output uses floating-point or complex storage even when the input right-hand side contains integers.

```python
import numpy as np

def forward_substitution(L: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Solve Ly = b for a nonsingular lower triangular L."""
    n = L.shape[0]
    y = np.zeros(n, dtype=np.result_type(L.dtype, b.dtype, np.float64))
    for i in range(n):
        if L[i, i] == 0:
            raise np.linalg.LinAlgError("Zero diagonal in lower triangular solve")
        y[i] = (b[i] - L[i, :i] @ y[:i]) / L[i, i]
    return y


def backward_substitution(U: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Solve Ux = y for a nonsingular upper triangular U."""
    n = U.shape[0]
    x = np.zeros(n, dtype=np.result_type(U.dtype, y.dtype, np.float64))
    for i in range(n - 1, -1, -1):
        if U[i, i] == 0:
            raise np.linalg.LinAlgError("Zero diagonal in upper triangular solve")
        x[i] = (y[i] - U[i, i + 1:] @ x[i + 1:]) / U[i, i]
    return x
```

For complex matrices, these are ordinary matrix products; no conjugation is needed.

### Cost of the Two Solves

We count each real addition, subtraction, multiplication, or division as one **floating-point operation** (flop). A multiply followed by an addition counts as two flops, even if a processor performs them in one fused instruction.

Row $i$ of forward substitution uses $i-1$ multiplications, $i-1$ subtractions, and one division. Thus a general triangular solve costs

$$
\sum_{i=1}^n\bigl(2(i-1)+1\bigr)=n^2
$$

flops. For unit lower triangular $L$, the divisions can be omitted, reducing the count to $n(n-1)$. Together, forward and backward substitution cost approximately $2n^2$ flops. Complex arithmetic has different constants but the same $O(n^2)$ cost.

The factors depend only on $A$. If several right-hand sides must be solved with the same matrix, we reuse $L$ and $U$ and repeat only the triangular solves.

## Deriving LU through Outer Products

Recall the [outer-product view of matrix multiplication](matrix_matrix_multiplication.md):

$$
A=LU=\sum_{k=1}^n l_{:,k}u_{k,:}.
$$

Here $l_{:,k}$ is column $k$ of $L$ and $u_{k,:}$ is row $k$ of $U$. Each product has rank at most one. The algorithm determines one such column-row pair at a time, then factors the remainder.

### The First Column and Row

Partition $A$ as

$$
A=\begin{pmatrix}a_{11}&r\\c&B\end{pmatrix},
$$

where $r$ is a row vector and $c$ is a column vector. If $a_{11}\neq0$, block multiplication verifies

$$
\begin{pmatrix}a_{11}&r\\c&B\end{pmatrix}
=
\begin{pmatrix}1&0\\c/a_{11}&I\end{pmatrix}
\begin{pmatrix}a_{11}&r\\0&S\end{pmatrix},
\qquad
S=B-\frac{c\,r}{a_{11}}.
$$

The scalar $a_{11}$ is the first **pivot**, and the entries of $c/a_{11}$ are the **elimination multipliers**. The matrix $S$ is the [Schur complement](block_matrices.md) of the pivot.

In outer-product terms, we have determined

$$
l_{:,1}=\begin{pmatrix}1\\c/a_{11}\end{pmatrix},
\qquad
u_{1,:}=\begin{pmatrix}a_{11}&r\end{pmatrix}.
$$

Subtracting their outer product leaves

$$
A-l_{:,1}u_{1,:}
=\begin{pmatrix}0&0\\0&S\end{pmatrix}.
$$

The first row and column have been accounted for. Only the smaller matrix $S$ remains to be factored. If $S=L_SU_S$, then

$$
L=\begin{pmatrix}1&0\\c/a_{11}&L_S\end{pmatrix},
\qquad
U=\begin{pmatrix}a_{11}&r\\0&U_S\end{pmatrix}
$$

satisfy $A=LU$. Repeating this argument proves the recursive construction whenever the required pivots are nonzero.

### Connection to Elimination and Projections

In Gaussian elimination, subtracting $(c_i/a_{11})$ times the first row from each lower row produces $S$ in the trailing block. LU records those multipliers in $L$ and the resulting pivot rows in $U$.

The rank-one remainder also has an [oblique projection](projections.md) interpretation. We work in $\mathbb{R}^n$, with $n\geq2$, for the geometric discussion. Set

$$
m=\frac{c}{a_{11}},\qquad l=l_{:,1}=\begin{pmatrix}1\\m\end{pmatrix},
\qquad P=I-le_1^T=\begin{pmatrix}0&0\\-m&I\end{pmatrix},
$$

where $e_1$ is the first coordinate vector. To see what $R$ does, apply it to a vector with first coordinate $s$:

$$
x=\begin{pmatrix}s\\z\end{pmatrix},\qquad
Px=\begin{pmatrix}0\\z-ms\end{pmatrix},\qquad
x-Px=s\begin{pmatrix}1\\m\end{pmatrix}=sl.
$$

Thus $Px$ lies in the hyperplane $H=\{x:x_1=0\}$, and the displacement from $x$ to $Px$ is parallel to $l$. Geometrically, follow the line through $x$ parallel to $l$ until it meets $H$.

#### Why This Is a Projection

The formulas above give

$$
\operatorname{range}(P)=H,\qquad
\operatorname{null}(P)=\operatorname{span}\{l\},\qquad
P^2=P.
$$

````{prf:proof}
Every output of $P$ has first coordinate zero, so its range is contained in $H$. Conversely, $P(0,z)^T=(0,z)^T$: every vector in $H$ is fixed by $P$ and is therefore in its range. Also, $Px=0$ exactly when $z=ms$, or $x=sl$, proving the null-space statement. In particular, $Pl=0$.

Since $P$ fixes every vector in its range, applying it twice has the same effect as applying it once. Algebraically, $e_1^Tl=1$ gives

$$
P^2=(I-le_1^T)^2
=I-2le_1^T+l(e_1^Tl)e_1^T
=I-le_1^T=P.
$$
````

This projection is orthogonal exactly when $l$ is perpendicular to $H$. Since $H^\perp=\operatorname{span}\{e_1\}$ and $l_1=1$, this happens exactly when $m=0$. Otherwise it is **oblique**: the direction we project along is tilted relative to the normal to $H$.

To recover the LU remainder, observe that $e_1^TA=u_{1,:}$ is the first row of $A$. Hence

$$
PA=A-l(e_1^TA)=A-lu_{1,:}
=\begin{pmatrix}0&0\\c-ma_{11}&B-mr\end{pmatrix}
=\begin{pmatrix}0&0\\0&S\end{pmatrix}.
$$

The first column vanishes because it is $a_{11}l$, a vector in the null space of $P$. Every column is projected into $H$, so the first row vanishes as well. The remaining block is the Schur complement.

#### Small Pivots and Long Shadows

Think of $H$ as the ground and of lines parallel to $l$ as rays of sunlight. The projected point is the shadow of $x$. When the rays are nearly perpendicular to the ground, the shadow stays close to the point. When they are nearly parallel to the ground, a small change in height can move the shadow a long distance.

Let $\theta\in(0,\pi/2]$ be the angle between the line through $l$ and $H$. Because the component of $l$ perpendicular to $H$ is $e_1$, of length one,

$$
\sin\theta=\frac{1}{\|l\|_2},\qquad
\|l\|_2=\sqrt{1+\|m\|_2^2}
=\sqrt{1+\frac{\|c\|_2^2}{|a_{11}|^2}}.
$$

The normalization $l_1=1$ matters here: rescaling a direction would not change its angle, but our LU convention fixes its scale. Thus a large $\|l\|_2$ means a small angle with $H$.

```{figure} ../_static/lu_oblique_projection.svg
:alt: Two projections of the same unit vector x=e1 onto H, drawn horizontally with identical scales. For l=(1,0.5), the shadow is at (0,-0.5) and the projection norm is about 1.12. For l=(1,4), the shadow is at (0,-4) and the norm is about 4.12. The orange displacement from Px to x equals l and makes a smaller angle with H in the second panel.

A good case and a bad case for amplification in one projection. In this two-dimensional drawing, $x_1$ is vertical and $x_2$ is horizontal, so $H$ is the ground. Both panels use the same scale and the same input $x=e_1$. The orange arrow is $x-Px=l$; the blue arrow is the shadow $Px$. A larger $\|l\|_2$ tilts the projection direction toward $H$ and lengthens the shadow.
```

The largest possible amplification of a vector is the operator norm of $P$, and it has a particularly simple formula:

$$
\boxed{\|P\|_2=\|l\|_2=\sqrt{1+\|m\|_2^2}=\frac{1}{\sin\theta}.}
$$

````{prf:proof}
Using the block form of $P$,

$$
PP^T=\begin{pmatrix}0&0\\0&I+mm^T\end{pmatrix}.
$$

If $m\neq0$, then $(I+mm^T)m=(1+\|m\|_2^2)m$, while every vector perpendicular to $m$ is left unchanged. Thus the largest eigenvalue of $PP^T$ is $1+\|m\|_2^2$. The same formula holds when $m=0$, since the trailing block is then $I$. The operator norm is the largest singular value, or the square root of this largest eigenvalue. This proves the formula.
````

For a fixed $P$, an input perturbation $\delta x$ changes the projection by $P\delta x$, with

$$
\|P\delta x\|_2\leq\|l\|_2\,\|\delta x\|_2.
$$

This bound is attained for a suitable perturbation direction. Even the simple choice $\delta x=\varepsilon e_1$ gives $P\delta x=(0,-\varepsilon m)^T$: a height change of size $|\varepsilon|$ moves the shadow by $|\varepsilon|\,\|m\|_2$. The bad case is therefore a pivot small **relative to** $\|c\|_2$, which makes the multipliers large. This is one reason for [row pivoting](lu_pivoting.md).

## The General Elimination Step

Let $A^{(0)}=A$ and define the mathematical remainder after $k$ steps by

$$
A^{(k)}=A-\sum_{j=1}^k l_{:,j}u_{j,:}.
$$

Its first $k$ rows and columns are zero. At step $k$, use the trailing block of $A^{(k-1)}$ to compute

$$
\begin{aligned}
u_{kj}&=a_{kj}^{(k-1)}, &&j=k,\ldots,n,\\
l_{kk}&=1,\\
l_{ik}&=a_{ik}^{(k-1)}/u_{kk}, &&i=k+1,\ldots,n.
\end{aligned}
$$

The remaining entries are updated by

$$
a_{ij}^{(k)}=a_{ij}^{(k-1)}-l_{ik}u_{kj},
\qquad i,j=k+1,\ldots,n.
$$

Thus the work at each step is a division to obtain the multipliers, followed by a rank-one update of the trailing block. The pivot is taken from the **current remainder**, not from the original diagonal of $A$.

At $k=n$, only $u_{nn}$ remains to be read off; there are no multipliers or trailing entries to update.

## In-Place Implementation

We can store both factors in one array:

- The entries on and above the diagonal hold $U$.
- The entries below the diagonal hold the multipliers in $L$.
- The unit diagonal of $L$ is implicit.

This packed array is different from the mathematical remainder $A^{(k)}$. After storing the multipliers and pivot row, we update **only the trailing block**, so the factors already computed are preserved.

```python
def lu_inplace(A: np.ndarray) -> np.ndarray:
    """Overwrite a square floating-point or complex array with packed LU.

    No row pivoting is performed. Return the same array.
    """
    if A.ndim != 2 or A.shape[0] != A.shape[1]:
        raise ValueError("A must be square")
    if not np.issubdtype(A.dtype, np.inexact):
        raise TypeError("A must have a floating-point or complex dtype")

    n = A.shape[0]
    for k in range(n - 1):
        if A[k, k] == 0:
            raise np.linalg.LinAlgError("Zero pivot: row pivoting may be needed")
        A[k + 1:, k] /= A[k, k]
        A[k + 1:, k + 1:] -= np.outer(A[k + 1:, k], A[k, k + 1:])
    return A
```

The function modifies its argument; pass a copy to retain the original matrix. Its zero-pivot check detects division by zero, but does not guarantee accuracy when pivots are nonzero. We examine the need for pivoting later in this chapter.

## A Worked Factorization and Solve

Consider

$$
A=\begin{pmatrix}2&1&1\\4&3&3\\8&7&9\end{pmatrix}.
$$

The first pivot is $2$, with multipliers $l_{21}=2$ and $l_{31}=4$. The trailing block becomes

$$
\begin{pmatrix}3&3\\7&9\end{pmatrix}
-\begin{pmatrix}2\\4\end{pmatrix}\begin{pmatrix}1&1\end{pmatrix}
=\begin{pmatrix}1&1\\3&5\end{pmatrix}.
$$

The second pivot is $1$, so $l_{32}=3$. The last pivot is $5-3\cdot1=2$. Therefore,

$$
L=\begin{pmatrix}1&0&0\\2&1&0\\4&3&1\end{pmatrix},
\qquad
U=\begin{pmatrix}2&1&1\\0&1&1\\0&0&2\end{pmatrix}.
$$

For $b=(3,7,17)^T$, forward substitution gives $y=(3,1,2)^T$, and backward substitution gives $x=(1,0,1)^T$.

Here is the complete computation using the functions above:

```python
A = np.array([[2., 1., 1.],
              [4., 3., 3.],
              [8., 7., 9.]])
b = np.array([3, 7, 17])

packed = lu_inplace(A.copy())
L = np.tril(packed, k=-1) + np.eye(A.shape[0])
U = np.triu(packed)
y = forward_substitution(L, b)
x = backward_substitution(U, y)

assert np.allclose(L @ U, A)
assert np.allclose(A @ x, b)
```

## Cost of LU Factorization

For a dense real matrix, step $k$ uses $n-k$ divisions for the multipliers and $2(n-k)^2$ flops for the trailing update. Summing gives

$$
\begin{aligned}
\text{factorization cost}
&=\sum_{k=1}^{n-1}\bigl((n-k)+2(n-k)^2\bigr)\\
&=\frac{2}{3}n^3+O(n^2).
\end{aligned}
$$

For $s$ right-hand sides, factoring once and reusing the factors costs approximately

$$
\frac{2}{3}n^3+2sn^2
$$

flops. The distinction is important: factorization costs $O(n^3)$, while each additional solve costs $O(n^2)$.

## What If a Pivot Is Zero?

The no-pivoting algorithm divides by $u_{kk}$ for $k=1,\ldots,n-1$. If one of these pivots is zero, the stated algorithm stops. This can happen even when $A$ is invertible. For example,

$$
A=\begin{pmatrix}0&1\\1&0\end{pmatrix}
$$

has determinant $-1$, but its first pivot is zero. Swapping the two rows removes the difficulty. Row interchanges lead to a factorization $PA=LU$, where $P$ is a permutation matrix; the triangular solves then use $Ly=Pb$ and $Ux=y$.

A zero **last** pivot is different: no elimination step divides by it. LU may still be completed, but $U$ is singular and backward substitution cannot produce a unique solution for arbitrary $b$. More generally, stopping on a zero pivot does not rule out all LU factorizations of a singular matrix.

The [next section](existence_lu.md) establishes the conditions for existence and uniqueness. We then study floating-point arithmetic and [row pivoting](lu_pivoting.md), including why nonzero pivots alone do not ensure an accurate computation.
