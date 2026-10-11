# Summary of QR Methods

For $A\in\mathbb{R}^{m\times n}$ with $m\ge n$, the reduced QR factorization is

$$
A=QR,\qquad Q\in\mathbb{R}^{m\times n},\qquad
Q^TQ=I_n,
$$

where $R\in\mathbb{R}^{n\times n}$ is upper triangular. The columns of $Q$ provide orthonormal directions, and $R$ contains the coefficients expressing the columns of $A$ in those directions. If $A$ has full column rank, these directions form a basis for its column space and $R$ is invertible.

## How the Methods Construct QR

**[Householder reflections](householder_reflections.md)** act on an entire active column segment. One reflection eliminates all entries below the diagonal in that column. The algorithm stores $Q$ implicitly as a product of reflections. For a dense matrix available in advance, Householder QR is the usual choice. Blocked implementations group reflections into updates that use matrix multiplication efficiently.

**[Givens rotations](givens_rotations.md)** combine two rows at a time, eliminating one selected entry. This local action is useful for structured matrices, such as upper Hessenberg and tridiagonal matrices, and for updating a factorization when a new row is added. Rotations can also create new nonzeros, so sparsity alone does not guarantee low cost. The algorithm stores $Q$ implicitly as a sequence of rotations.

A two-coordinate transformation can also be a reflection. The useful computational distinction is therefore **two-coordinate versus whole-column transformations**, rather than rotation versus reflection.

**[Modified Gram-Schmidt](modified_gram_schmidt.md)** constructs the columns of $Q$ explicitly, one at a time. Each projection uses the current residual after the previous components have been removed. This is convenient when an iterative method generates new vectors from the existing basis, as in Arnoldi iteration. Householder and Givens methods also support incremental constructions, but their basis is naturally stored as transformations.

## Accuracy and Computational Cost

Rounding affects both defining equations of QR. For computed factors, write

$$
A+\Delta A=\widehat Q\widehat R,\qquad
\widehat Q^T\widehat Q=I_n+E.
$$

Backward accuracy in these equations requires both a small relative factorization error $\|\Delta A\|_F/\|A\|_F$ and a small orthogonality error $\|E\|_2$. The computed $\widehat Q$ is generally not exactly orthonormal.

For $A\ne0$, completed Householder, Givens, and MGS factorizations give bounds of the form

$$
\frac{\|\Delta A\|_F}{\|A\|_F}\le C(m,n)u,
$$

where $u$ is the unit roundoff and $C(m,n)$ is a dimension-dependent factor. Their orthogonality guarantees differ:

| Method | Leading cost for dense $A$ | Bound on $\lVert E\rVert_2$ |
| :--- | :--- | :--- |
| Householder | $2mn^2-\frac23n^3$ | $C(m,n)u$ |
| Givens | $3mn^2-n^3$ | $C(m,n)u$ |
| One-pass MGS | $2mn^2$ | $C(m,n)u\kappa_2(A)$ |
| MGS with two passes | $4mn^2$ | $C(m,n)u$ |
| One-pass classical Gram-Schmidt | $2mn^2$ | Can be large |

Here $\kappa_2(A)=\sigma_{\max}(A)/\sigma_{\min}(A)$ for full-column-rank $A$. The constants can differ between methods. The bounds assume the usual rounding model, no overflow or significant underflow, and small dimension factors times $u$. The MGS and two-pass guarantees additionally require a dimension-dependent multiple of $u\kappa_2(A)$ to be sufficiently small.

For Householder and Givens, the listed costs compute $R$ and store $Q$ implicitly. Forming $Q$ explicitly adds work. Gram-Schmidt produces the reduced $Q$ directly. All costs omit lower-order terms; structured matrices can require much less work. For dense square matrices, Givens uses about 50% more arithmetic than Householder. Flop counts alone do not determine execution time.

MGS with a second pass is called **reorthogonalization**. It removes components left by rounding in the first pass and can restore orthogonality close to working precision. One-pass MGS can lose substantial orthogonality when the columns are nearly dependent, even while reconstructing $A$ accurately. Classical Gram-Schmidt computes its coefficients from the original column and can lose much more orthogonality.

The preceding sections give the derivations and assumptions. [Higham's QR overview](https://nhigham.com/2020/11/10/what-is-a-qr-factorization/) also summarizes the accuracy guarantees.
