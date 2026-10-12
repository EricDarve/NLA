# LU, QR, and the Determinant

Both LU and QR reduce a matrix to triangular form. The difference is how they introduce zeros: LU uses elimination and row interchanges, while QR uses orthogonal transformations. This affects their cost and numerical accuracy. Their triangular factors also provide an efficient way to compute determinants.

## Two Ways to Reach Triangular Form

For a square real matrix $A\in\mathbb{R}^{n\times n}$, the factorizations are

$$
\begin{aligned}
PA&=LU,\\
A&=QR,\qquad Q^TQ=I_n.
\end{aligned}
$$

Here $P$ records the row interchanges, $L$ is unit lower triangular, and $U$ and $R$ are upper triangular. LU records the elimination multipliers in $L$. QR records its orthogonal transformations in $Q$, usually stored implicitly as reflections or rotations.

Both factorizations exist for every square matrix, including singular matrices, when row pivoting is allowed for LU. Without row interchanges, LU need not exist. QR also exists for every rectangular matrix; for $m\ge n$, its reduced form has an $m\times n$ matrix with orthonormal columns and an $n\times n$ upper triangular factor. Full column rank is needed to make that triangular factor invertible, not to make QR exist.

For a nonsingular square system $Ax=b$, LU gives two triangular solves. QR instead applies $Q^T$ to $b$ and solves one triangular system. QR also preserves the residual norm in a least-squares problem, which makes it useful for rectangular systems. The [following sections](LS_using_QR.md) develop that application.

## Cost and Numerical Accuracy

For dense square matrices, the main comparison is between LU with partial pivoting and Householder QR:

| Property | LU with partial pivoting | Householder QR |
| :--- | :--- | :--- |
| Leading factorization cost | $\frac23n^3$ flops | $\frac43n^3$ flops |
| Control of intermediate values | Multipliers have magnitude at most $1$, but matrix entries can grow | Orthogonal transformations preserve the 2-norm of each column |
| Backward-error guarantee | Depends on element growth during elimination | Small factorization and orthogonality errors, with dimension-dependent factors times unit roundoff |

The costs omit lower-order terms. The QR cost includes computing $R$ and storing the reflections, but not forming $Q$ explicitly. For dense square systems, LU is therefore usually the less expensive choice. For symmetric positive definite matrices, [Cholesky factorization](cholesky.md) uses the additional structure and costs still less.

Partial pivoting is usually accurate in practice, but large element growth can cause large backward error, as shown in [the LU analysis](lu_pivoting.md). Householder QR avoids this dependence on element growth. Its guarantees assume the usual rounding model, no overflow or significant underflow, and dimension factors times the unit roundoff small compared with $1$. These are guarantees about the computed factorization, not guarantees of a small forward error in the solution. The [QR summary](qr_summary.md) compares the accuracy of Householder, Givens, and Gram-Schmidt methods.

Neither method removes the sensitivity of the original problem. In exact arithmetic, $A=QR$ gives

$$
\|Ax\|_2=\|QRx\|_2=\|Rx\|_2.
$$

Thus $A$ and $R$ have the same singular values, and $\kappa_2(R)=\kappa_2(A)$ when $A$ is nonsingular. QR can provide a stronger backward-error guarantee while the solution remains sensitive to small changes in the data.

## Determinants from the Triangular Factors

The determinant of a triangular matrix is the product of its diagonal entries. Combining this with $\det(BC)=\det(B)\det(C)$ gives a determinant formula for each factorization. The identities below concern exact factors and apply to singular matrices as well.

### From LU

Taking determinants in $PA=LU$ gives

$$
\det(P)\det(A)=\det(L)\det(U)
=\prod_{i=1}^n u_{ii},
$$

because the diagonal of $L$ consists of ones. If the algorithm makes $s$ row interchanges, then $\det(P)=(-1)^s$. Since this sign is its own inverse,

$$
\boxed{\det(A)=(-1)^s\prod_{i=1}^n u_{ii}.}
$$

The row interchanges affect the sign, so they must be counted. If no interchanges are needed, the determinant is simply the product of the pivots.

### From QR

:::{prf:theorem} Determinant via QR Factorization
:label: thm:det_via_qr
For a real square matrix $A=QR$, with $Q$ orthogonal and $R$ upper triangular,

$$
\boxed{
\det(A)=\det(Q)\prod_{i=1}^n r_{ii},
\qquad
|\det(A)|=\prod_{i=1}^n |r_{ii}|.
}
$$
:::

:::{prf:proof}
Orthogonality gives $Q^TQ=I_n$. Taking determinants yields

$$
1=\det(Q^TQ)=\det(Q)^2,
$$

so $\det(Q)=\pm1$. Also,

$$
\det(A)=\det(Q)\det(R)
=\det(Q)\prod_{i=1}^n r_{ii}.
$$

Taking absolute values gives the second formula.
:::

To obtain the sign, we do not need to form $Q$. Each nonidentity Householder reflection has determinant $-1$: it reverses one direction and fixes every perpendicular direction. Thus a product of $h$ such reflections has determinant $(-1)^h$. A product of Givens rotations has determinant $+1$, whereas each two-coordinate reflection contributes $-1$. Any later column-sign changes to $Q$ must also be included.

Once either factorization is available, computing the diagonal product costs only $O(n)$ operations. Starting from a dense matrix, the factorization costs $O(n^3)$, far less than evaluating the cofactor expansion for large $n$. If only the determinant is needed, LU with partial pivoting is usually the less expensive method.

## Computing the Magnitude Safely

A direct product of diagonal entries can overflow or underflow with **either LU or QR**, even if the factorization is accurate. To avoid forming that product, store the sign and the logarithm of the magnitude separately. For a nonsingular matrix,

$$
\log|\det(A)|
=\sum_{i=1}^n\log|u_{ii}|
=\sum_{i=1}^n\log|r_{ii}|.
$$

The signs are

$$
\begin{aligned}
\operatorname{sgn}(\det(A))
&=(-1)^s\prod_{i=1}^n\operatorname{sgn}(u_{ii}),\\
&=\det(Q)\prod_{i=1}^n\operatorname{sgn}(r_{ii}).
\end{aligned}
$$

For a singular matrix, use sign $0$ and log magnitude $-\infty$. NumPy's [`numpy.linalg.slogdet`](https://numpy.org/doc/stable/reference/generated/numpy.linalg.slogdet.html) returns this pair using LU factorization. Reconstructing the determinant by exponentiating the log magnitude can still overflow or underflow, so retain the pair when the magnitude is outside the representable range.

This representation addresses the range of floating-point numbers. It does not remove sensitivity to errors in $A$. Near a singular matrix, a small input perturbation can cause a large relative change in the determinant. A small determinant alone is also not evidence of ill-conditioning: $A=\alpha I_n$ has determinant $\alpha^n$, which can be tiny, but $\kappa_2(A)=1$ for every nonzero $\alpha$.
