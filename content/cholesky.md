---
jupytext:
  formats: md:myst
  text_representation:
    extension: .md
    format_name: myst
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Cholesky Factorization

For a symmetric positive definite matrix, elimination needs only one triangular factor. The upper triangular factor is the transpose of the lower one:

$$
\boxed{A=LL^T,\qquad l_{ii}>0.}
$$

Here $L$ is lower triangular. This is the **Cholesky factorization**. It requires about half the arithmetic of general LU factorization and can store the matrix and factor in a single triangle. Positive definiteness ensures that every pivot is positive in exact arithmetic, so no pivoting is needed.

The block derivation below explains why the factorization exists and how to compute it. The error analysis then explains its numerical advantage: under the stated rounding assumptions, a completed factorization has backward error of order $u$, up to a dimension-dependent factor. As with LU, solution accuracy also depends on conditioning.

## Symmetric Positive Definite Matrices

````{prf:definition} Symmetric Positive Definite Matrix
A real matrix $A\in\mathbb{R}^{n\times n}$ is **symmetric positive definite (SPD)** if

$$
A^T=A,
\qquad
x^TAx>0\quad\text{for every }x\ne0.
$$
````

For a symmetric matrix, positive definiteness is equivalent to having only positive eigenvalues. Indeed, write $A=Q\Lambda Q^T$ and set $z=Q^Tx$. Then

$$
x^TAx=z^T\Lambda z=\sum_{i=1}^n\lambda_i z_i^2.
$$

If every $\lambda_i>0$, this sum is positive for every nonzero $x$. Conversely, choosing $x$ to be an eigenvector gives $x^TAx=\lambda_i\|x\|_2^2>0$, so $\lambda_i>0$.

In particular, an SPD matrix is nonsingular, and each diagonal entry is positive because $a_{ii}=e_i^TAe_i>0$. Positive diagonal entries alone are not sufficient: for example, $\begin{pmatrix}1&2\\2&1\end{pmatrix}$ has eigenvalues $3$ and $-1$.

For complex matrices, the corresponding condition is **Hermitian positive definiteness**: $A^H=A$ and $x^HAx>0$ for every nonzero complex $x$. The factorization is then $A=LL^H$, again with positive real diagonal entries in $L$. The derivation and code below use real arithmetic; conjugate transposes give the analogous formulas in the complex case.

## Deriving the Factorization

As in LU, computing one column reduces the factorization to a smaller matrix. Separate the first row and column of the SPD matrix $A$ and its proposed factor $L$:

$$
A=\begin{pmatrix}a_{11}&c^T\\c&B\end{pmatrix},
\qquad
L=\begin{pmatrix}l_{11}&0\\\ell&L_S\end{pmatrix}.
$$

Multiplying the blocks gives

$$
LL^T=
\begin{pmatrix}
l_{11}^2&l_{11}\ell^T\\
l_{11}\ell&\ell\ell^T+L_SL_S^T
\end{pmatrix}.
$$

Matching these blocks with $A$ determines

$$
\begin{aligned}
l_{11}&=\sqrt{a_{11}},\\
\ell&=c/l_{11},\\
L_SL_S^T&=B-\ell\ell^T.
\end{aligned}
$$

Thus the remaining problem is to factor the **Schur complement**

$$
S=B-\frac{cc^T}{a_{11}}.
$$

To repeat this step, the remaining matrix must also be SPD. Symmetry is immediate from the formula for $S$; positive definiteness needs a proof. It ensures that the next pivot, and every later pivot, is positive.

### Why the Schur Complement Remains Positive Definite

For $x=(t,y^T)^T$, completing the square gives

$$
\begin{aligned}
x^TAx
&=a_{11}t^2+2t\,c^Ty+y^TBy\\
&=a_{11}\left(t+\frac{c^Ty}{a_{11}}\right)^2+y^TSy.
\end{aligned}
$$

Choose $t=-c^Ty/a_{11}$ to make the square vanish. If $y\ne0$, then $x\ne0$, so

$$
y^TSy=x^TAx>0.
$$

Also, $S$ is symmetric. Therefore $S$ is SPD: after eliminating one variable, we have a smaller problem with exactly the same structure.

### Existence and Uniqueness

````{prf:theorem} Existence and Uniqueness of Cholesky Factorization
:label: thm:cholesky_existence
A real symmetric matrix $A$ is positive definite if and only if it has a factorization $A=LL^T$ with $L$ lower triangular and $l_{ii}>0$. This factor is unique.
````

````{prf:proof} Existence and Uniqueness.
For $n=1$, the unique positive factor is $l_{11}=\sqrt{a_{11}}$.

Suppose the result holds in dimension $n-1$. For an SPD matrix $A$ of size $n$, the block equations above uniquely determine $l_{11}=\sqrt{a_{11}}$ and $\ell=c/l_{11}$. The Schur complement $S$ is SPD, so induction gives its unique factor $L_S$. Consequently,

$$
L=\begin{pmatrix}\sqrt{a_{11}}&0\\c/\sqrt{a_{11}}&L_S\end{pmatrix}
$$

is the unique lower triangular factor with positive diagonal satisfying $LL^T=A$.

Conversely, if such an $L$ exists, it is invertible. Hence, for $x\ne0$,

$$
x^TAx=x^TLL^Tx=\|L^Tx\|_2^2>0.
$$
````

This proof also explains why every pivot is positive in exact arithmetic. Positive **semidefinite** matrices can have zero pivots and require a modified algorithm; the divisions above are not valid at a zero pivot.

### Connection with LU

The Cholesky factor $L$ does not have unit diagonal. To recover the usual LU convention, let

$$
D=\operatorname{diag}(l_{11}^2,\ldots,l_{nn}^2),
\qquad
L_0=LD^{-1/2}.
$$

Then $L_0$ is unit lower triangular and

$$
A=L_0DL_0^T=L_0U,
\qquad U=DL_0^T.
$$

Thus Cholesky and unpivoted LU describe the same elimination, with different scalings of the triangular factors. The diagonal matrix $D$ contains the positive LU pivots. This also gives the **$LDL^T$ factorization**, here written with $L_0$ to distinguish its unit diagonal from that of the Cholesky factor.

## The In-Place Cholesky Algorithm

Repeating the block construction gives an algorithm that computes $L$ one column at a time. As in the [LU implementation](lu_pivoting.md), $a_{ij}$ denotes the **current value** in the array being overwritten. For $k=1,\ldots,n$:

1. Check that $a_{kk}>0$, and replace it by $\sqrt{a_{kk}}$.
2. Divide the entries below it by the new diagonal entry:

   $$
   a_{ik}\leftarrow a_{ik}/a_{kk},\qquad i>k.
   $$

3. Update only the lower triangle of the remaining block:

   $$
   a_{ij}\leftarrow a_{ij}-a_{ik}a_{jk},
   \qquad k<j\le i\le n.
   $$

At completion, the lower triangle contains $L$. The following real-arithmetic implementation reads and overwrites only that triangle; it leaves the upper triangle untouched. The lower triangle defines the symmetric input matrix, so the routine does not check agreement with the upper triangle.

```{code-cell} ipython3
import numpy as np

def cholesky_in_place(A: np.ndarray) -> np.ndarray:
    """Overwrite the lower triangle with L; return the same array.

    The lower triangle defines a real symmetric matrix A_original.
    On success, A_original = L @ L.T up to rounding error.
    """
    if A.ndim != 2 or A.shape[0] != A.shape[1]:
        raise ValueError("A must be square")
    if not np.issubdtype(A.dtype, np.floating):
        raise TypeError("A must have a real floating-point dtype")

    n = A.shape[0]
    for j in range(n):
        if not np.all(np.isfinite(A[j:, j])):
            raise ValueError("The lower triangle must contain finite entries")

    for k in range(n):
        if not np.isfinite(A[k, k]) or A[k, k] <= 0:
            raise np.linalg.LinAlgError(
                f"Nonpositive or nonfinite pivot at step {k + 1}"
            )
        A[k, k] = np.sqrt(A[k, k])
        A[k + 1:, k] /= A[k, k]
        for j in range(k + 1, n):
            A[j:, j] -= A[j:, k] * A[j, k]
    return A
```

A nonpositive computed pivot stops the routine. This can indicate that the input is not positive definite, or that rounding has obscured a very small positive pivot.

### Arithmetic and Storage

If the remaining block has size $m=n-k$, its lower triangle has $m(m+1)/2$ entries. Each update uses one multiplication and one subtraction, giving $m(m+1)$ flops. Thus the updates cost

$$
\sum_{m=1}^{n-1}m(m+1)=\frac{n^3-n}{3}.
$$

Including the divisions and $n$ square roots gives a leading cost of **$n^3/3$ flops**, compared with $2n^3/3$ for general LU. Updating both triangles would duplicate work and lose this saving. Half the leading arithmetic does not imply exactly half the elapsed time; that also depends on the implementation and hardware.

Only $n(n+1)/2$ entries are needed to store the input and factor. The NumPy code above still uses a full $n\times n$ array, however. To realize the storage saving, the triangle must be stored in a packed format. In either layout, overwriting the input avoids allocating a separate array for the factor.

## Solving a Linear System

Once $A=LL^T$ is available, solve $Ax=b$ with two triangular solves:

$$
Ly=b,\qquad L^Tx=y.
$$

Each new right-hand side costs about $2n^2$ flops for the pair of solves; the factorization can be reused. There is no need to form $A^{-1}$.

## Backward Stability and Its Limits

In the LU analysis, the main obstacle to a small backward-error bound was growth in the factors. Cholesky avoids that obstacle because positive definiteness controls the factor entries. This gives the route to a solution-error bound: first bound the error in the factorization, then include the two triangular solves, and finally apply the condition number.

These error bounds assume that the computation completes. Cancellation can still cause a nonpositive computed pivot; that limitation is discussed after the error analysis.

### Why the Factor Entries Stay Controlled

In exact arithmetic, the diagonal of $A=LL^T$ gives

$$
a_{ii}=\sum_{k=1}^i l_{ik}^2.
$$

Every term is nonnegative, so

$$
\boxed{|l_{ik}|\le\sqrt{a_{ii}}.}
$$

The LU error bound involved $|L|\,|U|$, which measures the sizes of products before cancellation. For Cholesky, the corresponding matrix is $|L|\,|L|^T$. Cauchy–Schwarz bounds each of its entries:

$$
\begin{aligned}
(|L|\,|L|^T)_{ij}
&=\sum_k|l_{ik}l_{jk}|\\
&\le\left(\sum_k l_{ik}^2\right)^{1/2}
     \left(\sum_k l_{jk}^2\right)^{1/2}
=\sqrt{a_{ii}a_{jj}}.
\end{aligned}
$$

Thus these products remain controlled by the original diagonal entries. There is no counterpart to the exponential growth in the LU example.

### Backward Error of the Factorization

The preceding estimates concern the exact factor $L$. Rounding produces a different factor $\widehat{L}$, so those estimates cannot simply be applied to it. The rounding-error theorem first expresses $\widehat{L}$ as the exact factor of a perturbed matrix. Its diagonal equations then let us control the size of $\widehat{L}$.

````{prf:theorem} Backward Error of Cholesky Factorization
:label: thm:backward_error_cholesky
Assume real arithmetic with rounding to nearest, no overflow or underflow, and $(n+1)u<1$. If Cholesky factorization of an SPD matrix $A$ completes with a positive diagonal, the computed factor satisfies

$$
A+E=\widehat{L}\widehat{L}^T,
\qquad
|E|\le\gamma_{n+1}|\widehat{L}|\,|\widehat{L}|^T,
$$

Here $u$ is the unit roundoff, $\gamma_m=mu/(1-mu)$, and $E$ is symmetric. Absolute values and inequalities are entrywise; $|\widehat{L}|\,|\widehat{L}|^T$ is an ordinary matrix product of nonnegative matrices.
````

This is the Cholesky counterpart of the [LU backward-error theorem](lu_pivoting.md); the componentwise bound is stated by [Rump and Jeannerod](https://doi.org/10.1137/130927231). The goal is a normwise bound that depends only on $A$, $n$, and $u$.

To make the bound useful, the size of the computed factor must be bounded in terms of the original input $A$. Write $\gamma=\gamma_{n+1}$ and assume $\gamma<1$. The diagonal equations imply

$$
\sum_k\widehat{l}_{ik}^2=a_{ii}+e_{ii}
\le a_{ii}+\gamma\sum_k\widehat{l}_{ik}^2.
$$

Moving the last term to the left and summing over $i$ bounds the sum of squares of all factor entries:

$$
\|\widehat{L}\|_F^2\le\frac{\operatorname{tr}(A)}{1-\gamma}.
$$

Since $\bigl\||\widehat{L}|\bigr\|_2\le\|\widehat{L}\|_F$ and $\operatorname{tr}(A)\le n\|A\|_2$, we obtain

$$
\begin{aligned}
\|E\|_2
&\le\gamma\bigl\||\widehat{L}|\,|\widehat{L}|^T\bigr\|_2\\
&\le\gamma\|\widehat{L}\|_F^2
\le\frac{n\gamma}{1-\gamma}\|A\|_2.
\end{aligned}
$$

```{admonition} Key result: Backward stability of Cholesky
:class: important

Under the theorem's assumptions, a completed factorization is the exact Cholesky factor of a nearby SPD matrix. If $\gamma_{n+1}<1$, then

$$
\frac{\|E\|_2}{\|A\|_2}
\le\frac{n\gamma_{n+1}}{1-\gamma_{n+1}}.
$$

For small $nu$, this bound is approximately $n(n+1)u$, so it is small when $n^2u\ll1$. The dimension factor is a worst-case bound; there is no additional factor from uncontrolled element growth.
```

### Backward Error of the Computed Solution

A small factorization error is only part of the accuracy argument: the two triangular solves introduce further rounding errors. To assess the complete solve, use the backward error $\eta(\widehat{x})$ defined in the [LU section](lu_pivoting.md). It is the smallest common bound on relative changes in $A$ and $b$ that make $\widehat{x}$ exact. For $b\ne0$, it can be computed from the residual:

$$
\eta(\widehat{x})=
\frac{\|b-A\widehat{x}\|_2}
{\|A\|_2\,\|\widehat{x}\|_2+\|b\|_2}.
$$

Assume that factorization and both solves complete, with rounding to nearest, no overflow or underflow, and $(3n+1)u<1$. Combining their errors as in the [LU solve analysis](lu_pivoting.md) gives

$$
(A+\Delta A)\widehat{x}=b,
\qquad
|\Delta A|\le\gamma_{3n+1}|\widehat{L}|\,|\widehat{L}|^T.
$$

The factorization contributes $\gamma_{n+1}$; the two solves add $2\gamma_n+\gamma_n^2$. Their sum is at most $\gamma_{3n+1}$. Unlike the factorization perturbation $E$, the combined perturbation $\Delta A$ need not be symmetric, because the two solves introduce different rounding errors.

The pair $(\Delta A,0)$ is one admissible input perturbation. Since $\eta$ is the minimum over all admissible pairs, this pair gives an upper bound on $\eta$:

$$
\begin{aligned}
\eta(\widehat{x})
&\le\frac{\|\Delta A\|_2}{\|A\|_2}\\
&\le\gamma_{3n+1}
\frac{\bigl\||\widehat{L}|\,|\widehat{L}|^T\bigr\|_2}{\|A\|_2}\\
&\le\gamma_{3n+1}\frac{\|\widehat{L}\|_F^2}{\|A\|_2}.
\end{aligned}
$$

Substituting the factor estimate $\|\widehat{L}\|_F^2\le n\|A\|_2/(1-\gamma_{n+1})$ bounds the backward error entirely in terms of $n$ and $u$.

```{admonition} Key result: Backward error of a Cholesky solve
:class: important

Under the assumptions above,

$$
\boxed{\eta(\widehat{x})\le\beta_n,
\qquad
\beta_n=\frac{n\gamma_{3n+1}}{1-\gamma_{n+1}}.}
$$

For small $nu$, $\beta_n\approx n(3n+1)u$. This bounds the error of the complete solve using only the dimension and arithmetic precision. The actual $\eta(\widehat{x})$, measured from the residual, can be much smaller than this worst-case bound.
```

### Forward Error of the Computed Solution

Let $x$ be the exact solution of $Ax=b$. The bound on $\eta(\widehat{x})$ controls the input perturbation needed to explain the computed answer $\widehat{x}$. The condition number determines how much that perturbation can affect the solution. Substituting $\eta(\widehat{x})\le\beta_n$ into the [forward-error bound for linear systems](lu_pivoting.md) gives, when $\kappa_2(A)\beta_n<1$,

$$
\frac{\|\widehat{x}-x\|_2}{\|x\|_2}
\le\frac{2\kappa_2(A)\beta_n}{1-\kappa_2(A)\beta_n}.
$$

Thus $\kappa_2(A)\beta_n\ll1$ guarantees a small relative solution error. For SPD matrices,

$$
\kappa_2(A)=\frac{\lambda_{\max}(A)}{\lambda_{\min}(A)}.
$$

An eigenvalue that is small relative to the largest one can therefore make the solution sensitive, despite the backward-error guarantee.

### Completion in Floating-Point Arithmetic

Here $A$ denotes the matrix already stored in working precision. If data are converted from a higher precision, positive definiteness must be assessed after that conversion.

For a stored SPD matrix, standard Cholesky is guaranteed to run to completion when $\kappa_2(A)$ is safely below $u^{-1}$. If $\kappa_2(A)$ is comparable to $u^{-1}$, the matrix is numerically singular at the working precision and rounding can produce a nonpositive computed pivot. Thus the backward-error bounds above assume successful completion; see [Higham's discussion of Cholesky factorization](https://nhigham.com/2020/08/11/what-is-a-cholesky-factorization/).

## Choosing a Direct Solver

For a general dense nonsingular system, use LU with partial pivoting. When the matrix is known to be symmetric positive definite, Cholesky exploits that structure to reduce arithmetic and storage. In both cases, factor once, solve triangular systems for each right-hand side, and interpret the residual together with the condition number when assessing solution accuracy.

## Optional Advanced Material: Numerical Breakdown

```{warning}
Reading this section is **not required**. It goes beyond the course and is included only for students who are curious about the limits of the numerical stability guarantees.

This is not covered on the exams or graded assignments.
```

How well conditioned can a stored SPD matrix be if Cholesky nevertheless produces a zero or negative pivot? This section gives a family with condition number of order $1/(n^{3/2}u)$. Together with a sufficient condition for success, it shows that the power $n^{3/2}$ cannot be improved uniformly in dimension and precision for the algorithm considered here. The exact best multiplicative constant is not determined.

### A Guarantee of Successful Completion

Wilkinson's bound guarantees successful completion when

$$
20n^{3/2}u\kappa_2(A)\le 1.
$$

See [Higham and Pranesh, p. A261](https://epubs.siam.org/doi/pdf/10.1137/19M1298263?download=true). A refinement gives the sufficient condition

$$
\kappa_2(A)\le\frac{1}{3.9n^{3/2}u},
\qquad n>10,\qquad 3nu<0.1;
$$

see [Terenin et al., Result 1, p. 4](https://jmlr.org/papers/volume25/22-1170/22-1170.pdf), citing Kiełbasiński. These results assume the usual floating-point error model, with no underflow or overflow.

To state the extremal question precisely, fix an idealized binary floating-point system with a $p$-bit significand, rounding to nearest, and an unbounded exponent range, so that $u=2^{-p}$ and underflow and overflow are excluded. For a fixed algorithm, define

$$
c^*_{n,u}
=\inf\left\{
u\kappa_2(A):
A\text{ is a stored SPD matrix on which the algorithm fails}
\right\}.
$$

Then $\kappa_2(A)<c^*_{n,u}/u$ guarantees success in this model. For a finite-exponent implementation such as binary64, the guarantee applies only to computations without underflow or overflow. The threshold depends on the arithmetic and the order of operations, not just on $n$.

The construction below uses the in-place algorithm given earlier: each update performs one rounded multiplication followed by one rounded subtraction. Changing the accumulation order can change the outcome. This family also breaks down with fused multiply-add updates, as explained below.

### An Exactly Representable Matrix Family

Assume binary rounding to nearest, with unit roundoff $u=2^{-p}$, and no underflow or overflow. Choose

$$
k=4^q\ge16,\qquad k^3u\le10^{-4},\qquad n=9k.
$$

Let $H_k$ be the Walsh–Hadamard matrix, defined recursively by

$$
H_1=[1],\qquad
H_{2m}=\begin{bmatrix}H_m&H_m\\H_m&-H_m\end{bmatrix},
\qquad Q=\frac{H_k}{\sqrt{k}}.
$$

Thus $Q^TQ=I_k$. Because $k$ is a power of four, every entry of $Q$ is a signed power of two and is exactly representable. Write $J_{a,b}$ for the $a\times b$ all-ones matrix, and $J_m=J_{m,m}$.

Set

$$
h=\frac{2u}{\sqrt{k}},\qquad
\eta=6kh=12\sqrt{k}\,u,\qquad
\delta=(k-2)\eta.
$$

Here $\eta$ is a construction parameter, unrelated to the backward error $\eta(\widehat{x})$ used earlier. Obtain $t$ by rounding $\sqrt h$ and then advancing twice to the next floating-point number toward $+\infty$. This choice ensures

$$
h<\operatorname{fl}(t^2)\le h(1+12u),\qquad
0<t^2-h\le12uh.
$$

Define

$$
T=\begin{bmatrix}
tJ_{6k,2k}\\[1mm]
\frac32[I_k\;\;Q]
\end{bmatrix},
\qquad
G=\frac94\begin{bmatrix}I_k&Q\\Q^T&I_k\end{bmatrix},
$$

and

$$
\boxed{
A_{9k}=\begin{bmatrix}
4I_{7k}&2T\\
2T^T&G+\eta J_{2k}+\delta I_{2k}
\end{bmatrix}.
}
$$

All entries of this matrix are exactly representable under the stated size restriction. In particular, $\eta$ and $\delta$ lie on the floating-point grids of the entries to which they are added. The small entries occur first in $T$; their order is essential to the example.

### Positive Definiteness of the Stored Matrix

The exact Schur complement of $4I_{7k}$ is particularly simple:

$$
\begin{aligned}
S
&=G+\eta J_{2k}+\delta I_{2k}-T^TT\\
&=\delta I_{2k}-\varepsilon J_{2k},
\qquad \varepsilon=6k(t^2-h).
\end{aligned}
$$

The all-ones matrix $J_{2k}$ has eigenvalues $2k$ and $0$, so

$$
\lambda_{\min}(S)
=\delta-2k\varepsilon
\ge\delta-24ku\eta>0.
$$

Therefore $A_{9k}$ is SPD in exact arithmetic. Positive definiteness holds for the floating-point entries themselves; it is not a property lost while forming or storing the input.

### How Rounding Produces a Nonpositive Pivot

The factorization of the leading block is exact: $4I_{7k}=(2I_{7k})(2I_{7k})^T$, and the divisions below it recover $T$. The first $6k$ rows therefore apply repeated updates of size $\operatorname{fl}(t^2)$ to the trailing block.

In its cross block, the entries start at $\pm9/(4\sqrt{k})+\eta$. Their floating-point spacing is $2h$. Each subtraction of $\operatorname{fl}(t^2)$, which is slightly larger than $h$, rounds to a decrease of $2h$. After $6k$ such updates, those entries are

$$
\pm\frac{9}{4\sqrt{k}}+\eta-6k(2h)
=\pm\frac{9}{4\sqrt{k}}-\eta.
$$

The true decrease would be only slightly larger than $\eta$. Rounding has almost doubled it. The diagonal entries, whose spacing is larger, remain unchanged during these tiny updates.

The next $k$ leading rows cancel $G$. The resulting stored trailing matrix has the form

$$
\widehat S=
\begin{bmatrix}
dI_k+e(J_k-I_k)&-\eta J_k\\
-\eta J_k&dI_k
\end{bmatrix},
\qquad d=(k-1)\eta,
$$

where $e\le0$ and $|e|\le16ku\eta$. The small residual $e$ does not affect the sign argument. For the normalized all-ones vector $x\in\mathbb R^{2k}$,

$$
x^T\widehat Sx
=-\eta+\frac{k-1}{2}e
\le-\eta.
$$

````{prf:proof} Details of the stored Schur complement and subsequent breakdown
Let $y=\operatorname{fl}(t^2)$. A within-group off-diagonal entry starts at $\eta$ and undergoes $6k$ rounded subtractions of $y$. Comparing this recurrence with the exactly representable sequence $\eta-jh$, and using $y>h$ and monotonicity of rounding, shows that its final value $e$ is nonpositive.

The summation error bound controls the distance from $e$ to the exact sum $\eta-6ky$:

$$
|e-(\eta-6ky)|\le\gamma_{6k}(\eta+6ky).
$$

We must also bound that exact sum. Since $\eta=6kh$ and $h<y\le h(1+12u)$, we have $|\eta-6ky|\le12u\eta$. Combining the two estimates gives

$$
|e|\le\left[12u+\gamma_{6k}(2+12u)\right]\eta
<16ku\eta,
$$

where the last inequality uses $k\ge16$ and $k^3u\le10^{-4}$.

During the remaining $k$ updates, the left within-group off-diagonal entries receive only zero products, so they retain $e$. In the right group, the first nonzero product has magnitude $9/(4k)$. The residual $e$ is too small to affect its rounding: the condition $k^3u\le10^{-4}$ gives

$$
|e|\le192k^{3/2}u^2<\frac{u}{k}.
$$

The subsequent signed products are dyadic numbers whose sums are exact and cancel by Hadamard orthogonality. The cross-block and diagonal cancellations are also exact. This gives the displayed formula for $\widehat S$.

An indefinite intermediate matrix alone does not prove that a floating-point routine must report failure. To finish the argument, suppose the remaining Cholesky factorization completed with positive pivots. Its rounding-error analysis would give

$$
\widehat L\widehat L^T=\widehat S+E,
\qquad
|E|\le\gamma_{2k+1}|\widehat L|\,|\widehat L|^T.
$$

This algebraic error bound holds whenever the computation completes, even if the input is indefinite. The diagonal equations and the trace argument used earlier imply

$$
\|E\|_2
\le\frac{\gamma_{2k+1}}{1-\gamma_{2k+1}}
\operatorname{tr}(\widehat S)
<10^{-3}\eta,
$$

because $\operatorname{tr}(\widehat S)=2k(k-1)\eta$ and $k^3u\le10^{-4}$. But then

$$
x^T\widehat L\widehat L^Tx
=x^T\widehat Sx+x^TEx<0,
$$

which is impossible. The computation must encounter a zero or negative pivot.

For fused multiply-add updates, each tiny update is $\operatorname{fl}(a-t^2)$. Because $t^2$ is slightly larger than $h$, the cross-block entries still decrease by exactly $2h$. Replacing $y$ by the exact $t^2$ in the residual estimates gives the same bounds on $e$; the remaining nonzero products are exact. The breakdown argument therefore also applies to fused updates.
````

### Condition Number and the Remaining Gap

To bound the condition number, set

$$
b=\frac92+2k\eta(1+12u),\qquad
s=\delta-24ku\eta.
$$

Then $\|T\|_2^2\le b$ and $\lambda_{\min}(S)\ge s$. The block factorization and its inverse give

$$
\|A_{9k}\|_2\le4+b+\delta,
\qquad
\|A_{9k}^{-1}\|_2
\le\frac14+\frac{1+b/4}{s}.
$$

Consequently, under the stated restrictions,

$$
\boxed{
\kappa_2(A_{9k})
\le(4+b+\delta)\left(\frac14+\frac{1+b/4}{s}\right)
<\frac{47}{n^{3/2}u}.
}
$$

To prove a matching lower bound, write $A=A_{9k}$ and

$$
M=[2I_{7k}\;\;T],\qquad
N=\begin{bmatrix}-T/2\\I_{2k}\end{bmatrix}.
$$

Here $X\succeq Y$ means that $X-Y$ is positive semidefinite. Since $A\succeq M^TM$ and $MM^T=4I_{7k}+TT^T$,

$$
\lambda_{\max}(A)\ge4+\|T\|_2^2\ge\frac{17}{2}.
$$

Also, $0\prec S\preceq\delta I_{2k}$ and the block inverse satisfies

$$
A^{-1}=\begin{bmatrix}I_{7k}/4&0\\0&0\end{bmatrix}
+NS^{-1}N^T
\succeq\frac{1}{\delta}NN^T.
$$

Therefore

$$
\|A^{-1}\|_2\ge\frac{1+\|T\|_2^2/4}{\delta}
\ge\frac{17}{8\delta},
\qquad
\kappa_2(A)\ge\frac{289}{16\delta}.
$$

For fixed $k$ as $u\to0$, the upper bound above has the same leading term: $b\to9/2$, $s/\delta\to1$, and $\delta\to0$. Thus the two bounds prove

$$
u\kappa_2(A_{9k})
\sim\frac{2601}{64}\frac{k}{k-2}\,n^{-3/2}.
$$

The coefficient approaches $2601/64=40.640625$ as $k$ grows and the precision increases sufficiently. Thus, for this algorithm and the admissible dimensions,

$$
\frac{1}{3.9n^{3/2}}\le c^*_{n,u}<\frac{47}{n^{3/2}}.
$$

These bounds establish the power of $n$, but leave a gap in the multiplicative constant. The claim concerns a joint regime of increasing dimension and decreasing roundoff, with $k^3u\le10^{-4}$; it is not an assertion about arbitrarily large $n$ at fixed precision. Other sufficiently large dimensions can be covered by padding a smaller member with $4I$, preserving its condition number and breakdown but increasing the constant in the bound expressed using the padded dimension.

### A Binary64 Example

The following values were obtained with $u=2^{-53}$ and the same order of operations as `cholesky_in_place`. Positive definiteness was checked using exact rational arithmetic. The condition numbers shown are rounded values from analytic bounds evaluated at 80-digit precision, rather than from a binary64 condition estimator.

| Size $n$ | Approximate $\kappa_2(A_n)$ | Failing pivot (counting from 1) | Computed pivot |
|---:|---:|---:|---:|
| 144 | $2.42102\times10^{14}$ | 143 | $-1.19904\times10^{-12}$ |
| 576 | $2.73341\times10^{13}$ | 575 | $-4.23022\times10^{-11}$ |
| 2304 | $3.33605\times10^{12}$ | 2303 | $-1.38609\times10^{-9}$ |

This code constructs the smallest example in the table and applies the earlier `cholesky_in_place` function. The exact-fraction check establishes that the stored matrix is SPD before the numerical factorization begins.

```python
import math
from fractions import Fraction
import numpy as np


def cholesky_breakdown_matrix(k=16):
    u = 2.0**-53
    root_k = math.isqrt(k)
    assert k >= 16 and root_k * root_k == k and k & (k - 1) == 0
    assert k**3 * u <= 1e-4

    H = np.ones((1, 1))
    while H.shape[0] < k:
        H = np.block([[H, H], [H, -H]])
    Q = H / root_k

    h = 2 * u / root_k
    eta = 6 * k * h
    delta = (k - 2) * eta
    t = math.nextafter(math.nextafter(math.sqrt(h), math.inf), math.inf)

    T = np.vstack((np.full((6 * k, 2 * k), t),
                   1.5 * np.hstack((np.eye(k), Q))))
    G = 2.25 * np.block([[np.eye(k), Q], [Q.T, np.eye(k)]])
    C = G + eta * np.ones((2 * k, 2 * k)) + delta * np.eye(2 * k)

    # Check every stored entry of C, grouping equal base values in G.
    for base in np.unique(G):
        expected = Fraction(float(base)) + Fraction(eta)
        if base == 2.25:  # Only diagonal entries have this base value.
            expected += Fraction(delta)
        for stored in np.unique(C[G == base]):
            assert Fraction(float(stored)) == expected

    epsilon = 6 * k * Fraction(t)**2 - Fraction(eta)
    assert epsilon > 0
    assert Fraction(delta) - 2 * k * epsilon > 0
    return np.block([[4 * np.eye(7 * k), 2 * T], [2 * T.T, C]])


A = cholesky_breakdown_matrix()
try:
    cholesky_in_place(A.copy())
except np.linalg.LinAlgError as error:
    print(error)
# Nonpositive or nonfinite pivot at step 143
```
