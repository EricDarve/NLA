# Symmetric and Nonsymmetric QR Iteration

For a real $n\times n$ matrix $A$, eigenvalue computation has two stages: reduce the matrix by orthogonal similarity, then apply shifted QR iteration. Symmetry makes both stages cheaper and ensures that all eigenvalues are real.

## Reduction: Hessenberg or Tridiagonal

As in [Hessenberg reduction](upper_hessenberg.md), first compute

$$
H=Q_{\mathrm{red}}^TAQ_{\mathrm{red}},
$$

where $Q_{\mathrm{red}}$ is orthogonal. For a general real matrix, $H$ is upper Hessenberg and may have $O(n^2)$ nonzero entries.

If $A=A^T$, then

$$
H^T=(Q_{\mathrm{red}}^TAQ_{\mathrm{red}})^T=Q_{\mathrm{red}}^TAQ_{\mathrm{red}}=H.
$$

Thus $H$ is both symmetric and upper Hessenberg, so it is **tridiagonal**. It has only a diagonal and matching first subdiagonal and superdiagonal, requiring $O(n)$ storage. Exploiting symmetry also reduces the work of the initial Householder reduction.

We keep the notation $H$ for the current matrix, which becomes smaller as eigenvalues are deflated.

## The Nonsymmetric Case

The [previous section](qr_iteration_with_shifts.md) describes shift selection and deflation for this case. The eigenvalues of the trailing $2\times2$ block supply shift estimates. Nonreal shifts are used in conjugate pairs, and the Francis double-shift step allows real arithmetic throughout. Its efficient implementation costs $O(n^2)$ per step; the implementation details are beyond this summary.

The target is **real Schur form**: an upper block triangular matrix with $1\times1$ diagonal blocks for real eigenvalues and $2\times2$ blocks for nonreal conjugate pairs. Once a $2\times2$ block separates from the rest of the matrix, its two eigenvalues can be computed directly. Its internal subdiagonal entry need not vanish.

Trailing-block shifts can stagnate, so practical methods also use exceptional shifts and other safeguards. Rapid local convergence does not guarantee convergence from every starting matrix.

## The Symmetric Case

Because every eigenvalue is real, a single real shift is sufficient. Choose a shift $\mu$ and take the same step as before:

$$
\boxed{
\begin{aligned}
H-\mu I&=UR,\\
H_{\mathrm{new}}&=RU+\mu I=U^THU.
\end{aligned}}
$$

Here $U$ is orthogonal, $R$ is upper triangular, and $H_{\mathrm{new}}$ is the updated matrix. Orthogonal similarity preserves symmetry, and the adjacent-rotation construction preserves Hessenberg form. Consequently, $H_{\mathrm{new}}$ is again symmetric tridiagonal.

### Why a QR Step Costs Only $O(n)$

The narrow band makes both the factorization and the reverse product inexpensive. To factor $H-\mu I$, use at most $n-1$ Givens rotations, eliminating the subdiagonal from top to bottom. Each rotation updates only a fixed number of potentially nonzero entries. The resulting $R$ has at most two superdiagonals.

Store these rotations and apply them from the right to compute $RU$. The intermediate matrices also have a band of fixed width, so each right rotation again requires only a fixed amount of work. Both parts of the step therefore cost $O(n)$. This cost assumes that $U$ is stored as rotations and that the reverse product is computed by applying them, rather than by dense matrix multiplication.

### Choosing the Wilkinson Shift

The bottom-right entry $h_{nn}$ is an inexpensive shift, but it can cycle even for a symmetric matrix. For example, $H=\begin{pmatrix}0&1\\1&0\end{pmatrix}$ has bottom-right entry zero. With this shift, the magnitude of its subdiagonal entry remains 1.

The **Wilkinson shift** uses the trailing symmetric block

$$
B=\begin{pmatrix}a&\varepsilon\\\varepsilon&d\end{pmatrix},
\qquad d=h_{nn}.
$$

Here $\varepsilon$ is the off-diagonal entry. Compute the two eigenvalues of $B$ and choose the one closer to $d$, using a fixed rule to break ties. Explicitly, the candidates are

$$
\mu_\pm=\frac{a+d}{2}
\pm\sqrt{\left(\frac{a-d}{2}\right)^2+\varepsilon^2}.
$$

They are always real. In the example above they are $\pm1$; either choice gives deflation in one step. If $\varepsilon=0$, the last eigenvalue is already separated and no QR step is needed. Numerical implementations use scaled formulas to avoid overflow and unnecessary cancellation.

**For real symmetric tridiagonal matrices, Wilkinson-shift QR converges to deflation from every starting matrix in exact arithmetic.** The bottom subdiagonal tends to zero, so it eventually meets any positive deflation tolerance. Repeating the procedure on the remaining blocks computes all the eigenvalues. This global convergence guarantee is a main difference from the usual nonsymmetric trailing-block shift rules; see [Leite, Saldanha, and Tomei](https://arxiv.org/abs/0912.3376).

### Why Convergence Is Often Cubic

To see the effect of symmetry, use the simpler bottom-right shift $\mu=d$ on the $2\times2$ matrix $H=B$ above. Assume $\varepsilon\ne0$, and let $\varepsilon_{\mathrm{new}}$ denote its subdiagonal entry after the step. Since both off-diagonal entries equal $\varepsilon$, the previous section's $2\times2$ calculation gives

$$
\begin{aligned}
|\varepsilon_{\mathrm{new}}|
&=\frac{|\varepsilon|^3}{(a-d)^2+\varepsilon^2}\\
&\approx\frac{|\varepsilon|^3}{(a-d)^2}
\qquad\text{when }|\varepsilon|\ll|a-d|.
\end{aligned}
$$

In this regime, the new off-diagonal entry is approximately proportional to the cube of the old one. This is **cubic local convergence**: once the off-diagonal entry is small and the diagonal entries are separated, it can decrease very rapidly.

The Wilkinson shift also usually gives cubic local convergence, but some exceptional cases converge only quadratically. The cubic rate is therefore not a universal guarantee; see [Leite, Saldanha, and Tomei's analysis of shift strategies](https://arxiv.org/abs/0912.3376).

When a subdiagonal entry is set to zero, symmetry makes the matching superdiagonal entry zero as well. The matrix splits into independent diagonal blocks. Repeated scalar deflation leads to a diagonal matrix $D$, so the final factorization is the eigendecomposition

$$
A=QDQ^T,
$$

where $Q$ is the accumulated orthogonal transformation. Its columns are eigenvectors of $A$.

## Cost Comparison

The following costs are for real arithmetic and **eigenvalues only**. Reduction costs cover $H$ and the stored Householder reflections. The total QR costs assume that only a few steps are needed per eigenvalue on average; more iterations increase these costs.

| Quantity | Symmetric | Nonsymmetric |
| :--- | :--- | :--- |
| Reduced form | Tridiagonal | Upper Hessenberg |
| Reduction cost | About $\frac43n^3$ flops | About $\frac{10}{3}n^3$ flops |
| QR step on a block of size $n$ | $O(n)$ | $O(n^2)$ |
| Typical total QR-phase cost | $O(n^2)$ | $O(n^3)$ |
| Total cost including reduction | $O(n^3)$ | $O(n^3)$ |
| Final real Schur form | Diagonal | Upper block triangular |

For a dense symmetric matrix, the initial reduction usually dominates the eigenvalue-only cost. In the nonsymmetric case, both reduction and QR iteration contribute at order $n^3$.

If eigenvectors are also required, the symmetric QR phase becomes more expensive. Each Givens rotation must update two full columns of the accumulated matrix $Q$, costing $O(n)$ per rotation. A full-size step then costs $O(n^2)$, and the usual total cost of accumulating all eigenvectors is $O(n^3)$. The $O(n)$ step cost applies to eigenvalues alone. [LAPACK's symmetric eigenproblem guide](https://www.netlib.org/lapack/lug/node48.html) distinguishes these outputs and also describes alternatives to QR iteration.
