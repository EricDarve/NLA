# Givens Rotations

A **Givens rotation** introduces a zero by rotating just two coordinates. In QR factorization, this means combining two rows to eliminate one entry below the diagonal. The other rows remain unchanged. This local action is useful when a matrix has only a few entries to eliminate, or when we need to update an existing factorization.

As in the [Householder construction](householder_reflections.md), we work with $A\in\mathbb{R}^{m\times n}$, where $m\ge n$, and seek $A=QR$ with $Q$ orthogonal and $R$ upper trapezoidal. We first construct one rotation, then explain how a sequence of rotations produces QR.

The algorithmic distinction is how many coordinates each step changes: two coordinates at a time, or a whole active column segment. A two-coordinate step can use either a rotation or a reflection, as we show below.

## Rotating Two Coordinates

Consider the vector $x=(a,b)^T$. To eliminate its second entry, we want

$$
G\begin{pmatrix}a\\b\end{pmatrix}
=\begin{pmatrix}r\\0\end{pmatrix},
\qquad
G=\begin{pmatrix}c&s\\-s&c\end{pmatrix},
\qquad c^2+s^2=1.
$$

If $c=\cos\theta$ and $s=\sin\theta$, this matrix rotates vectors through the angle $-\theta$. For positive $\theta$, the rotation is clockwise. We use this sign convention throughout; no angle needs to be computed in the algorithm.

:::{figure} ../_static/givens_rotation.svg
:alt: A vector x=(a,b) points up and right. A clockwise rotation through theta moves it to r e_1 on the horizontal axis. The vectors have the same length, and a curved arrow shows the direction of rotation.
:width: 600px

The rotation moves $x=(a,b)^T$ onto the horizontal axis, eliminating $b$ while preserving the length $r=\sqrt{a^2+b^2}$.
:::

:::{prf:proof} Orthogonality and the Inverse.
Multiplying the transpose by the original matrix gives

$$
\begin{aligned}
G^TG&=
\begin{pmatrix}c&-s\\s&c\end{pmatrix}
\begin{pmatrix}c&s\\-s&c\end{pmatrix}\\
&=
\begin{pmatrix}c^2+s^2&cs-sc\\sc-cs&s^2+c^2\end{pmatrix}
=I.
\end{aligned}
$$

Thus $G$ is orthogonal, $G^{-1}=G^T$, and

$$
\|Gz\|_2^2=z^TG^TGz=\|z\|_2^2
$$

for every vector $z$. Its determinant is $c^2+s^2=1$, so it is a rotation rather than a reflection. Transposing $G$ changes $s$ to $-s$ and reverses the rotation.
:::

### Choosing $c$ and $s$ to Eliminate $b$

Expanding the matrix-vector product gives

$$
G\begin{pmatrix}a\\b\end{pmatrix}
=\begin{pmatrix}ca+sb\\-sa+cb\end{pmatrix}.
$$

The second entry vanishes if $sa=cb$. For $(a,b)\ne(0,0)$, choose

$$
\boxed{r=\sqrt{a^2+b^2},
\qquad c=\frac{a}{r},
\qquad s=\frac{b}{r}.}
$$

These formulas satisfy both requirements:

$$
\begin{aligned}
c^2+s^2&=\frac{a^2+b^2}{r^2}=1,\\
-sa+cb&=-\frac{ba}{r}+\frac{ab}{r}=0,\\
ca+sb&=\frac{a^2+b^2}{r}=r.
\end{aligned}
$$

If $a=b=0$, use $c=1$, $s=0$, and $r=0$, so $G=I$. If $b=0$, there is already no entry to eliminate, and the QR algorithm may skip the rotation. In that case the first entry remains $a$, which need not be positive. The signs of the diagonal entries of $R$ do not affect the validity of QR.

## Two-Coordinate Versus Global Transformations

A two-coordinate reflection can perform the same elimination. For the same nonzero pair $(a,b)$ and coefficients $c=a/r$, $s=b/r$, define

$$
H=\begin{pmatrix}c&s\\s&-c\end{pmatrix}.
$$

It sends the vector to the same target:

$$
H\begin{pmatrix}a\\b\end{pmatrix}
=\begin{pmatrix}ca+sb\\sa-cb\end{pmatrix}
=\begin{pmatrix}r\\0\end{pmatrix}.
$$

This matrix is symmetric, and multiplying it by itself gives

$$
H^2=
\begin{pmatrix}c^2+s^2&cs-sc\\sc-cs&s^2+c^2\end{pmatrix}
=I.
$$

Thus $H^TH=I$ and $H^{-1}=H$. Its determinant is $-c^2-s^2=-1$, so it is a reflection. It is exactly a two-dimensional [Householder reflector](https://nhigham.com/2020/09/15/what-is-a-householder-matrix/).

The rotation and reflection differ only in the sign of their second output:

$$
H=\begin{pmatrix}1&0\\0&-1\end{pmatrix}G.
$$

Both eliminate $b$. When embedded in a larger matrix, both change only two rows, preserve the other rows, and require the same arithmetic per pair of entries. Symmetry simplifies the inverse, but does not by itself improve numerical accuracy.

:::{important}
For QR, the useful distinction is **two-coordinate versus global transformations**. A two-coordinate rotation or reflection eliminates one selected entry. A global Householder reflection acts on the whole active column segment and eliminates all entries below its first entry in one step.

A Householder reflection can also act on just two coordinates. The computational difference comes from the number of coordinates involved in each step, rather than from choosing a rotation or a reflection. Rotations acting on more than two coordinates are also possible. Householder reflections are the standard global construction used here.
:::

The QR algorithm below uses the rotation $G$. Replacing each rotation by the corresponding two-coordinate reflection gives another valid QR algorithm. The zero-preservation argument is unchanged, and each reflection is its own inverse.

## Computing the Rotation Without Overflow

The formulas above explain the mathematics, but evaluating $a^2+b^2$ directly can fail. For example, squaring a large finite number may overflow even when $\sqrt{a^2+b^2}$ is representable. Squaring tiny numbers may underflow to zero.

We avoid these unnecessary squares by scaling the inputs. For a nonzero pair, set

$$
h=\max(|a|,|b|),
\qquad
\alpha=\frac{a}{h},
\qquad
\beta=\frac{b}{h},
\qquad
d=\sqrt{\alpha^2+\beta^2}.
$$

Both scaled inputs have magnitude at most $1$, and at least one has magnitude $1$. Therefore $1\le d\le\sqrt{2}$. Compute

$$
c=\frac{\alpha}{d},
\qquad
s=\frac{\beta}{d},
\qquad
r=hd.
$$

These are the same formulas as before, since $r=h\sqrt{(a/h)^2+(b/h)^2}$. Computing $c$ and $s$ from the scaled inputs also avoids dividing by a very large or very small $r$.

The same scaled coefficients can be used in the two-coordinate reflection $H$.

Scaling cannot prevent overflow if the true $r$ exceeds the largest representable number. Production routines also handle extreme values and underflow carefully. LAPACK's [DLARTG](https://www.netlib.org/lapack/explore-html/da/dd3/group__lartg_ga86f8f877eaea0386cdc2c3c175d9ea88.html) uses safe scaling. Its sign convention allows a negative $r$ so that $c\ge0$; this differs from our choice $r\ge0$, but achieves the same elimination.

## Applying a Rotation to Two Rows

For indices $i<j$, the $m\times m$ rotation matrix is

$$
G_{ij}(c,s)=
\begin{pmatrix}
I_{i-1}&0&0&0&0\\
0&c&0&s&0\\
0&0&I_{j-i-1}&0&0\\
0&-s&0&c&0\\
0&0&0&0&I_{m-j}
\end{pmatrix}
\quad
\begin{array}{l}
\text{rows }1,\ldots,i-1\\
\text{row }i\\
\text{rows }i+1,\ldots,j-1\\
\text{row }j\\
\text{rows }j+1,\ldots,m
\end{array}.
$$

The columns follow the same grouping as the rows. The identity blocks leave all coordinates except $i$ and $j$ unchanged. Blocks of size zero are omitted.

Multiplying $A$ on the left by this matrix changes only rows $i$ and $j$. For each column $\ell$,

$$
\begin{aligned}
a_{i\ell}^{\mathrm{new}}&=ca_{i\ell}+sa_{j\ell},\\
a_{j\ell}^{\mathrm{new}}&=-sa_{i\ell}+ca_{j\ell}.
\end{aligned}
$$

Both right-hand sides use the **old** entries. An in-place implementation saves $a_{i\ell}$ and $a_{j\ell}$ before overwriting either. This is the update performed by the BLAS routine [DROT](https://www.netlib.org/lapack/explore-html/dc/d23/drot_8f_source.html).

Each column requires four multiplications and two additions or subtractions. Updating two row segments of length $t$ therefore costs $6t$ flops. We store only $i,j,c,s$, not the full matrix $G_{ij}$.

## Building the QR Factorization

The order of elimination matters: a rotation must preserve the zeros already created. One simple choice uses adjacent rows and works from the bottom of each column upward.

```{figure} ../_static/givens_dense_qr.svg
:alt: Four diagrams show a dense three-by-three matrix becoming upper triangular. First, rows 2 and 3 are rotated to clear the bottom entry of column 1. Next, rows 1 and 2 are rotated to clear the remaining entry below its diagonal, leaving the bottom zero untouched. Finally, rows 2 and 3 are rotated to clear the entry below the diagonal in column 2, preserving both zeros in column 1. Blue cells mark the active row segments, orange circles mark entries to eliminate, and green zeros mark completed eliminations.
:width: 780px
:name: fig-givens-dense-qr

Adjacent-row elimination for a dense matrix. In the first column, begin at the bottom and move upward. Then move to the next column. Blue cells are mixed by the next rotation, the orange circle marks its target, and green zeros show earlier eliminations. In the third panel, both active rows are already zero in the first column, so those zeros survive the final rotation. Dots denote entries that may be nonzero.
```

Starting with $A^{(0)}=A$, process columns $k=1,\ldots,n$. Within column $k$, process rows $j=m,m-1,\ldots,k+1$:

1. Take $a$ and $b$ to be the current entries in positions $(j-1,k)$ and $(j,k)$.
2. If $b=0$, continue to the next row. Otherwise, construct $c,s$ to send $(a,b)^T$ to $(r,0)^T$.
3. Apply the rotation to rows $j-1,j$, updating columns $k$ through $n$. Store $r$ and zero in the two positions of column $k$, and retain the rotation parameters for later use.

Why does this preserve the earlier work? In every earlier column, both affected rows lie below its diagonal, so their entries are zero. Rotating the pair $(0,0)^T$ leaves it zero. Within the current column, the entries below row $j$ have already been eliminated, and this rotation does not touch those rows. Induction therefore shows that each finished column is zero below its diagonal.

Number the rotations in the order they are applied as $G_1,\ldots,G_N$. Then

$$
A^{(\ell)}=G_\ell A^{(\ell-1)},
\qquad
R=G_N\cdots G_2G_1A.
$$

Define

$$
Q^T=G_N\cdots G_2G_1,
\qquad
Q=G_1^TG_2^T\cdots G_N^T.
$$

Each rotation is orthogonal, so their product is orthogonal. Since $R=Q^TA$, multiplying by $Q$ gives $A=QR$. Together with the zero-preservation argument, this proves the factorization. If no rotations are needed, $Q=I$.

**Unlike a Householder reflection, a Givens rotation is generally not symmetric.** The transposes in the formula for $Q$ are necessary: they reverse each rotation.

### Applying $Q$ and $Q^T$ Without Forming Them

To compute $Q^Tb$, apply $G_1$ first, followed by $G_2,\ldots,G_N$. Each rotation changes only two entries of the vector. To compute $Qb$, apply $G_N^T$ first and continue in reverse order through $G_1^T$, changing $s$ to $-s$ for each rotation.

For least squares, partition the transformed right-hand side as

$$
Q^Tb=\begin{pmatrix}c_1\\c_2\end{pmatrix},
\qquad
R=\begin{pmatrix}R_1\\0\end{pmatrix},
\qquad c_1\in\mathbb{R}^n.
$$

Here $R_1\in\mathbb{R}^{n\times n}$ is upper triangular. As explained in the [Householder least-squares construction](householder_reflections.md), solve $R_1x=c_1$ by back substitution. For full column rank $A$, this gives the unique minimizer, and the minimum residual norm is $\|c_2\|_2$. The algorithm can transform $b$ alongside $A$ if that is the only right-hand side needed.

## Computational Cost for Dense Matrices

In a dense matrix, column $k$ has $m-k$ entries to eliminate. For the algorithm above, the number of rotations is at most

$$
N_{\max}=\sum_{k=1}^n(m-k)
=mn-\frac{n(n+1)}{2}.
$$

The maximum is attained when none of the entries selected for elimination is already zero. A rotation in column $k$ updates two row segments of length $n-k+1$, costing approximately $6(n-k+1)$ flops. Hence the leading cost of dense Givens QR is

$$
\begin{aligned}
\operatorname{cost}(QR)
&=6\sum_{k=1}^n(m-k)(n-k+1)+O(mn)\\
&=3mn^2-n^3+O(mn).
\end{aligned}
$$

To evaluate the sum, set $t=n-k+1$. Then $m-k=m-n-1+t$, so

$$
\begin{aligned}
6\sum_{t=1}^n(m-n-1+t)t
&=3(m-n-1)n(n+1)\\
&\qquad{}+n(n+1)(2n+1)\\
&=n(n+1)(3m-n-2).
\end{aligned}
$$

Constructing the rotations and assigning the eliminated entries adds only lower-order work. This cost does not include forming $Q$ explicitly.

For square matrices, the leading cost is $2n^3$. Householder QR costs $\frac43n^3$, so Givens QR uses about 50% more arithmetic. For tall matrices, the leading costs are $3mn^2$ and $2mn^2$, respectively. Blocked Householder algorithms also group updates into matrix multiplications, which usually makes them the preferred choice for dense QR.

Applying the stored rotations to one vector costs $6N=O(mn)$ flops. Back substitution adds $O(n^2)$ work. Storing a general dense rotation sequence requires $O(mn)$ space.

## When the Local Action Helps

### Hessenberg and Tridiagonal Matrices

An $n\times n$ **upper Hessenberg matrix** satisfies $a_{ij}=0$ whenever $i>j+1$. Each column has at most one entry below the diagonal. Rotate rows $k,k+1$ to eliminate that entry in column $k$.

```{figure} ../_static/givens_hessenberg.svg
:alt: Four diagrams show a four-by-four upper Hessenberg matrix becoming upper triangular. The first rotation mixes rows 1 and 2 and eliminates the entry below the first pivot. The next mixes rows 2 and 3, preserving that zero while eliminating the next subdiagonal entry. The last mixes rows 3 and 4. Blue cells show the active row segments, orange circles mark entries to eliminate, and green zeros show completed eliminations.
:width: 780px
:name: fig-givens-hessenberg

Upper Hessenberg QR uses adjacent row pairs, moving down the subdiagonal. Each rotation eliminates one circled entry. Blue cells show the row segments mixed by the next rotation, and green zeros show earlier eliminations. Dots denote entries that may be nonzero.
```

This does not create entries to eliminate farther down: both affected rows were already zero in earlier columns, and all lower rows remain untouched. At most $n-1$ rotations are needed. Each updates at most $O(n)$ entries, so computing $R$ and storing the rotations costs $O(n^2)$.

A **tridiagonal matrix** has $a_{ij}=0$ whenever $|i-j|>1$. The same adjacent-row rotations eliminate its subdiagonal, but now the affected row segments are short.

```{figure} ../_static/givens_tridiagonal.svg
:alt: Four diagrams show a four-by-four tridiagonal matrix becoming upper triangular. The same row pairs are rotated as in the Hessenberg example, but each rotation touches at most three columns. New blue dots appear in the second superdiagonal in positions formerly occupied by zeros. Green zeros mark the eliminated subdiagonal entries. The final factor has at most two nonzero superdiagonals.
:width: 780px
:name: fig-givens-tridiagonal

Tridiagonal QR uses the same rotation order, with at most three active columns per step. The colors have the same meaning as above; blue dots additionally mark **fill-in**, where a zero becomes nonzero. The new entries lie on the second superdiagonal, so the final triangular factor remains narrow.
```

Here the construction costs only $O(n)$, because each rotation updates a fixed number of entries. Before step $k$, row $k$ can be nonzero only through column $k+1$, and row $k+1$ only through column $k+2$. Their combination can therefore create a new entry at $(k,k+2)$, but none farther to the right. The resulting $R$ has at most two nonzero superdiagonals.

These costs refer to $R$ and an implicit representation of $Q$. Explicitly forming the full $n\times n$ matrix $Q$ costs $O(n^2)$ even in the tridiagonal case.

The same costs are possible with two-coordinate reflections: they update the same row pairs and create the same pattern of zeros. Householder transformations can also exploit short column segments, so their cost on these structured matrices need not be cubic. What matters here is the local two-row action.

### General Sparse Matrices and Updates

For sparse rows, a rotation needs to update only columns where at least one of the two rows is nonzero. However, mixing the rows can create new nonzeros in either row. This is called **fill-in**. Local elimination therefore does not automatically preserve sparsity or make Givens QR the best choice for every sparse matrix.

Rotations are also useful when a new row is added to a least-squares problem. Apply the stored QR transformations to the old rows, leaving the added row unchanged. The old rows become a triangular block above a block of zeros. Those zero rows need no further work, so the figure shows only the triangular block and the added row.

```{figure} ../_static/givens_row_update.svg
:alt: Four matrix diagrams show a triangular three-column block with an added row. In the first panel, row 1 and the added row are highlighted, and the added row's first entry is circled. After one rotation that entry is zero, and row 2 and the remaining added-row entries are highlighted. After two rotations the first two added-row entries are zero, and row 3 is paired with the added row. The final panel has a triangular top block and an all-zero added row.
:width: 780px
:name: fig-givens-row-update

Updating QR after adding a row, illustrated with three columns. Blue cells are the two row segments mixed by the next rotation. The orange circle marks the entry to eliminate; green zeros mark entries already eliminated. A dot denotes an entry that may be nonzero. The same pattern continues across any number of columns.
```

Start in the first column. Rotate the first row of the triangular block together with the added row to make the added row's first entry zero. Apply this rotation to all the remaining entries in those two rows as well.

Move to the second column and pair the second row with the added row. Use the **current entries**, since the first rotation changed the added row. Both rows now have zero in the first column, so this rotation preserves that zero. Continue from left to right, pairing each successive row of the triangular block with the added row. If the entry to eliminate is already zero, skip the rotation.

At each stage, both affected rows are zero in the earlier columns, and all other rows stay unchanged. The top block therefore remains triangular. When the added row is entirely zero, the top block is the updated triangular factor. For $n$ columns, at most $n$ rotations are needed, each updating at most $n$ entries in each of two rows. This costs $O(n^2)$ rather than refactoring all the original rows.

For an updated least-squares solve with full column rank, apply each rotation to the corresponding right-hand-side entries too: the entry for the pivot row and the added observation. Then solve the updated triangular system. The original residual block is unchanged.

## Numerical Stability

Safe computation of $c$ and $s$ prevents avoidable overflow and underflow. For the complete factorization, use the two backward errors introduced in the {ref}`Householder discussion <qr-backward-error>`: the relative perturbation $\|\Delta A\|_F/\|A\|_F$ in the factorization equation and the orthogonality defect $\|E\|_2$.

If the full $\widehat Q$ is formed from the stored rotations, Givens QR satisfies the same form of bounds:

$$
\boxed{
\begin{aligned}
\frac{\|\Delta A\|_F}{\|A\|_F}
&\le C_1(m,n)u,\\
\|E\|_2
&\le C_2(m,n)u.
\end{aligned}
}
$$

These hold under the arithmetic assumptions stated for Householder QR, with dimension-dependent factors that do not involve $\kappa_2(A)$. Thus Givens QR also controls both defining equations, even though the computed $\widehat Q$ is not exactly orthogonal.

Orthogonal transformations help because they preserve vector lengths and do not amplify existing errors: $\|Gz\|_2=\|z\|_2$. The analysis must still account for new rounding errors in constructing and applying each rotation.
