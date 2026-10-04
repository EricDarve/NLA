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

Here $L$ is lower triangular. This is the **Cholesky factorization**. It requires about half the arithmetic of general LU factorization. Only one triangle is needed to store the matrix and its factor. Positive definiteness ensures that every pivot is positive in exact arithmetic, so no pivoting is needed.

## Symmetric Positive Definite Matrices

````{prf:definition} Symmetric Positive Definite Matrix
A real matrix $A\in\mathbb{R}^{n\times n}$ is **symmetric positive definite (SPD)** if

$$
A^T=A,
\qquad
x^TAx>0\quad\text{for every }x\ne0.
$$
````

For a symmetric matrix, positive definiteness is equivalent to having only positive eigenvalues. To see this, write $A=Q\Lambda Q^T$ and set $z=Q^Tx$. Then

$$
x^TAx=z^T\Lambda z=\sum_{i=1}^n\lambda_i z_i^2.
$$

If every $\lambda_i>0$, this sum is positive for every nonzero $x$. Conversely, choosing $x$ to be an eigenvector gives $x^TAx=\lambda_i\|x\|_2^2>0$, so $\lambda_i>0$.

An SPD matrix is therefore nonsingular. Each diagonal entry is also positive because $a_{ii}=e_i^TAe_i>0$. Positive diagonal entries alone do not imply positive definiteness: $\begin{pmatrix}1&2\\2&1\end{pmatrix}$ has eigenvalues $3$ and $-1$.

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

To repeat the elimination step, we need $S$ to be SPD. Its symmetry follows from the formula. Proving positive definiteness will ensure that the next pivot, and every later pivot, is positive.

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

Together with symmetry, this proves that $S$ is SPD. We can therefore repeat the same elimination step on the smaller matrix.

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

At completion, the lower triangle contains $L$. The implementation below reads and overwrites only that triangle. It interprets these entries as a symmetric matrix, without checking or changing the upper triangle.

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

A nonpositive computed pivot stops the routine. The input may not be positive definite, or rounding may have changed the sign of a pivot that is positive in exact arithmetic.

### Arithmetic and Storage

If the remaining block has size $m=n-k$, its lower triangle has $m(m+1)/2$ entries. Each update uses one multiplication and one subtraction, giving $m(m+1)$ flops. Thus the updates cost

$$
\sum_{m=1}^{n-1}m(m+1)=\frac{n^3-n}{3}.
$$

Including the divisions and $n$ square roots gives a leading cost of **$n^3/3$ flops**, compared with $2n^3/3$ for general LU. Updating both triangles would duplicate work and lose this saving. The elapsed time also depends on the implementation and hardware, so it need not be exactly half that of LU.

Only $n(n+1)/2$ entries are needed to store the input and factor. The NumPy code above uses a full $n\times n$ array. Storing just the triangle in a packed format would save space. In either layout, overwriting the input avoids allocating a separate array for the factor.

## Solving a Linear System

Once $A=LL^T$ is available, solve $Ax=b$ with two triangular solves:

$$
Ly=b,\qquad L^Tx=y.
$$

Each new right-hand side costs about $2n^2$ flops for the pair of solves; the factorization can be reused. There is no need to form $A^{-1}$.

## Backward Stability and Its Limits

In the LU analysis, large entries in the factors could amplify rounding errors. For Cholesky, positive definiteness bounds the factor entries in terms of the input matrix. We use this bound to estimate the factorization error, then include the errors from the triangular solves. The condition number converts the resulting backward-error bound into a bound on solution accuracy.

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

The preceding estimates concern the exact factor $L$. To bound the computed factor $\widehat{L}$, we first express it as the exact factor of a perturbed matrix. Let $u$ be the unit roundoff and write $\gamma_m=mu/(1-mu)$.

````{prf:theorem} Backward Error of Cholesky Factorization
:label: thm:backward_error_cholesky
Assume real arithmetic with rounding to nearest, no overflow or underflow, and $(n+1)u<1$. If Cholesky factorization of an SPD matrix $A$ completes with a positive diagonal, the computed factor satisfies

$$
A+E=\widehat{L}\widehat{L}^T,
\qquad
|E|\le\gamma_{n+1}|\widehat{L}|\,|\widehat{L}|^T,
$$

The perturbation $E$ is symmetric. Absolute values and inequalities are entrywise; $|\widehat{L}|\,|\widehat{L}|^T$ is an ordinary matrix product of nonnegative matrices.
````

This is the Cholesky counterpart of the [LU backward-error theorem](lu_pivoting.md). The componentwise bound is stated by [Rump and Jeannerod](https://doi.org/10.1137/130927231). It still depends on $\widehat{L}$, so we need to bound that factor using $A$. Write $\gamma=\gamma_{n+1}$ and assume $\gamma<1$. The diagonal equations give

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

The two triangular solves introduce rounding errors in addition to those from the factorization. To assess their combined effect, use the backward error $\eta(\widehat{x})$ defined in the [LU section](lu_pivoting.md). It is the smallest common bound on relative changes in $A$ and $b$ that make $\widehat{x}$ exact. For $b\ne0$, the residual gives

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

The perturbation pair $(\Delta A,0)$ changes $A$ and leaves $b$ unchanged. It makes $\widehat{x}$ an exact solution, so its size bounds the minimum perturbation $\eta(\widehat{x})$:

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

For small $nu$, $\beta_n\approx n(3n+1)u$. This bound accounts for both the factorization and the triangular solves. The actual backward error, measured from the residual, can be much smaller.
```

### Forward Error of the Computed Solution

Let $x$ be the exact solution of $Ax=b$. The backward error tells us how much the input must change to make $\widehat{x}$ exact. The condition number tells us how much that change can affect the solution. When $\kappa_2(A)\beta_n<1$, substituting $\eta(\widehat{x})\le\beta_n$ into the [forward-error bound for linear systems](lu_pivoting.md) gives

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

For a stored SPD matrix, Wilkinson's sufficient condition for successful completion is

$$
20n^{3/2}u\kappa_2(A)\le1,
$$

assuming the usual floating-point error model and no underflow or overflow. The dimension factor matters: failure does not by itself imply that $\kappa_2(A)$ is comparable to $u^{-1}$. The backward-error bounds above assume successful completion. See the {ref}`optional advanced section <cholesky-numerical-breakdown>` for the source, a sharper sufficient condition, and examples of numerical breakdown.

## Choosing a Direct Solver

For a general dense nonsingular system, use LU with partial pivoting. For a known SPD matrix, use Cholesky to reduce arithmetic and storage. In either case, reuse the factorization for each new right-hand side. Assess solution accuracy using both the residual and the condition number.

(cholesky-numerical-breakdown)=
## Optional Advanced Material: Numerical Breakdown

```{warning}
Reading this section is **not required**. It goes beyond the course and is included only for students who are curious about the limits of the numerical stability guarantees.

*This is not covered on the exams or graded assignments.*
```

To assess how sharp the completion guarantee is, we need SPD matrices on which Cholesky fails. The family below has condition numbers bounded by a constant times $\max\{1,1/(n^{3/2}u)\}$. Combined with a sufficient condition for success, it determines the order of the smallest condition number permitting failure. The best multiplicative constant remains undetermined.

### A Guarantee of Successful Completion

Wilkinson's bound guarantees successful completion when

$$
20n^{3/2}u\kappa_2(A)\le 1.
$$

See [Higham and Pranesh, p. A261](https://epubs.siam.org/doi/pdf/10.1137/19M1298263?download=true). A sharper result comes from {download}`Kiełbasiński, A Note on Rounding-Error Analysis of Cholesky Factorization (1987) <../addl_material/1-s2.0-0024379587901212-main.pdf>`, Theorem 2, p. 491. Write $\rho_n$ for the paper's $p_n$, to distinguish it from the significand precision:

$$
\rho_1=2,\qquad \rho_2=3,\qquad
\rho_n=1+\sqrt{\rho_{n-1}^2+4}\quad(n>2).
$$

The theorem guarantees completion if

$$
3u\rho_n\|A^{-1}\|_2\|A\|_F\le1,
$$

and gives the Frobenius backward-error bound

$$
\widehat L\widehat L^T=A+E,\qquad
\|E\|_F<u(1+3u)^n\rho_n\|A\|_F.
$$

Since $\|A\|_F\le\sqrt n\|A\|_2$, a convenient sufficient condition is

$$
\kappa_2(A)\le\frac{1}{3.9n^{3/2}u},
\qquad n\ge13.
$$

To obtain the factor $3.9$, use the recurrence to get $\rho_{13}<16.8<1.3\cdot13$. The inequality $1+\sqrt{(1.3n)^2+4}\le1.3(n+1)$ then gives $\rho_n\le1.3n$ by induction for $n\ge13$. All dimensions used below satisfy this range. These results assume the usual floating-point error model, with no underflow or overflow.

The extra assumption $3nu\le0.1$ in the paper's Corollary 2 simplifies the backward-error bound; it is not an additional hypothesis of Theorem 2's completion guarantee. The paper's discussion of near-sharpness on p. 494 concerns the Frobenius backward error of a successful factorization. It does not construct matrices on which Cholesky fails, or establish failure with condition numbers close to $1$.

To compare a success guarantee with a failing example, we need a precise threshold. Fix binary arithmetic with a $p$-bit significand and rounding to nearest, so $u=2^{-p}$. Assume an unbounded exponent range to exclude underflow and overflow. For a fixed algorithm, define

$$
c^*_{n,u}
=\inf\left\{
u\kappa_2(A):
A\text{ is a stored SPD matrix on which the algorithm fails}
\right\}.
$$

Then $\kappa_2(A)<c^*_{n,u}/u$ guarantees success in this model. For a finite-exponent implementation such as binary64, the guarantee applies only to computations without underflow or overflow. The threshold depends on the arithmetic and the order of operations, not just on $n$.

The construction below uses the in-place algorithm given earlier: each update performs one rounded multiplication followed by one rounded subtraction. Changing the accumulation order can change the outcome. This family also breaks down with fused multiply-add updates, as explained below.

### What It Means to Keep the Precision Fixed

Fix $u$ throughout the construction. The relevant target across all dimensions is

$$
\kappa_2(A_n)\le C\max\left\{1,\frac{1}{n^{3/2}u}\right\},
$$

with a constant $C$ independent of both $n$ and $u$. At fixed $u$, the expression $1/(n^{3/2}u)$ eventually drops below $1$, while every condition number is at least $1$. The maximum accounts for this limit. At dimensions of order $u^{-2/3}$ and beyond, the target is a failing matrix with a bounded condition number.

The basic blocks below satisfy $k^{3/2}u\le1/256$. Choosing the largest admissible block gives a condition-number bound independent of precision. Adding a diagonal block then produces failing matrices of every larger size without changing the condition number.

### An Exactly Representable Matrix Family

We use a diagonal leading block so that its factorization is exact. Rounding in the updates to the remaining block will cause the failure.

Assume binary rounding to nearest, with fixed unit roundoff $u=2^{-p}$, and no underflow or overflow. Choose an integer $q$ such that

$$
q\ge3,\qquad p+q\text{ is even},\qquad 3q+8\le p,
\qquad k=4^q,\qquad n=9k.
$$

The restriction $3q+8\le p$ is equivalent to $k^{3/2}u\le1/256$. Requiring $p+q$ to be even makes the small entries below exactly representable powers of two times an integer.

Let $H_k$ be the Walsh–Hadamard matrix, defined recursively by

$$
H_1=[1],\qquad
H_{2m}=\begin{bmatrix}H_m&H_m\\H_m&-H_m\end{bmatrix},
\qquad Q=\frac{H_k}{\sqrt{k}}.
$$

Thus $Q^TQ=I_k$, and every entry of $Q$ is an exactly representable signed power of two. Write $J_{a,b}$ for the $a\times b$ all-ones matrix, and $J_m=J_{m,m}$.

Set

$$
h=\frac{2u}{\sqrt{k}},\qquad
\eta=6kh=12\sqrt{k}\,u,\qquad
\delta=k\eta,\qquad
t=3\,2^{-(p+q+2)/2}.
$$

Here $\eta$ is a construction parameter, unrelated to the backward error $\eta(\widehat{x})$ used earlier. The choice of $t$ gives the exact identity

$$
t^2=\frac98h.
$$

Both $t$ and $t^2$ are floating-point numbers; no rounded square root is needed to form them. Define

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

All entries are exactly representable. The parameters $\eta$ and $\delta$ are multiples of $4u$, the spacing near the diagonal value $9/4$. The shifts in the cross block are multiples of $4u/\sqrt{k}$. The size restriction keeps these entries in intervals with the stated spacing. The small rows occur first in $T$ so that their updates accumulate rounding errors before the larger rows are eliminated.

### Positive Definiteness of the Stored Matrix

The exact Schur complement of $4I_{7k}$ is

$$
\begin{aligned}
S
&=G+\eta J_{2k}+\delta I_{2k}-T^TT\\
&=\delta I_{2k}-\frac{\eta}{8}J_{2k}.
\end{aligned}
$$

Since $J_{2k}$ has eigenvalues $2k$ and $0$, and $\delta=k\eta$,

$$
\boxed{
\lambda_{\min}(S)=\frac34k\eta>0,
\qquad \lambda_{\max}(S)=k\eta.
}
$$

Both $4I_{7k}$ and its exact Schur complement are positive definite, so the stored matrix is SPD.

### How Rounding Changes the Schur Complement

The factorization of the leading block is exact: $4I_{7k}=(2I_{7k})(2I_{7k})^T$, and the divisions below it recover $T$. The first $6k$ rows apply repeated updates of size $t^2=9h/8$ to the trailing block.

In the cross block, entries start at $\pm9/(4\sqrt{k})+\eta$ and have floating-point spacing $2h$. To determine how each subtraction rounds, we must check that this spacing stays unchanged throughout the updates.

The size restriction gives $ku\le1/2048$, so $\eta\sqrt{k}=12ku\le3/512<1/4$. The path from $\pm9/(4\sqrt{k})+\eta$ to $\pm9/(4\sqrt{k})-\eta$ therefore stays strictly between $2/\sqrt{k}$ and $4/\sqrt{k}$ in magnitude. These endpoints are consecutive powers of two, so all entries stay in the same *binade*, where the spacing is $2h$. Each subtraction of $9h/8$ rounds to a decrease of $2h$, giving after $6k$ updates

$$
\pm\frac{9}{4\sqrt{k}}+\eta-6k(2h)
=\pm\frac{9}{4\sqrt{k}}-\eta.
$$

The diagonal entries have larger spacing and remain unchanged during these tiny updates.

The off-diagonal entries within each group start at $\eta$. After $j$ updates their exact value is $(48k-9j)h/8$, with $0\le j\le6k$. The integer coefficient has magnitude at most $48k<2^{2q+6}$. Each value therefore requires at most $2q+6\le p$ significand bits. Every subtraction is exact, leaving

$$
\eta-6k\frac98h=-\frac{\eta}{8}.
$$

The next $k$ leading rows cancel $G$ exactly and retain the residual $-\eta/8$ in both groups. Because $q\ge3$, this residual is a multiple of $4u$, as are the Hadamard products $9/(4k)$. All intermediate sums have magnitude less than $4$, so they are exactly representable.

Consequently the actual stored trailing matrix is

$$
\boxed{
\widehat S=
\begin{bmatrix}
aI_k-\frac{\eta}{8}J_k&-\eta J_k\\
-\eta J_k&aI_k-\frac{\eta}{8}J_k
\end{bmatrix},
\qquad a=\left(k+\frac98\right)\eta.
}
$$

Its extreme eigenvalues are

$$
\lambda_{\min}(\widehat S)=-\frac{k-9}{8}\eta,
\qquad
\lambda_{\max}(\widehat S)=\frac{15k+9}{8}\eta.
$$

To force breakdown, the negative eigenvalue must be large enough to survive the remaining rounding errors. Here its magnitude is a fixed fraction of the matrix norm: since $k\ge64$,

$$
\frac{-\lambda_{\min}(\widehat S)}{\|\widehat S\|_2}
=\frac{k-9}{15k+9}>\frac1{18}.
$$

### Why a Nonpositive Pivot Must Follow

Rounding could still allow the routine to complete on an indefinite intermediate matrix. We rule this out by using the negative off-diagonal entries to bound the remaining rounding error.

````{prf:proof} Breakdown from the sign pattern
Suppose the remaining Cholesky factorization completes with positive pivots, producing $\widehat L$. Every off-diagonal entry of $\widehat S$ is negative. Each off-diagonal update subtracts a nonnegative product, so every off-diagonal entry of $\widehat L$ is nonpositive.

Let $D$ be the diagonal part of $\widehat L$. Then

$$
|\widehat L|=2D-\widehat L,
\qquad
\bigl\||\widehat L|\bigr\|_2
\le2\|D\|_2+\|\widehat L\|_2
\le3\|\widehat L\|_2.
$$

The componentwise backward-error bound, with $\gamma=\gamma_{2k+1}$, gives

$$
\widehat L\widehat L^T=\widehat S+E,
\qquad
\|E\|_2\le\gamma\bigl\||\widehat L|\bigr\|_2^2
\le9\gamma\|\widehat S+E\|_2.
$$

Hence

$$
\|E\|_2
\le\frac{9\gamma}{1-9\gamma}\|\widehat S\|_2
<\frac1{100}\|\widehat S\|_2.
$$

For the last inequality, use $ku\le1/2048$, which follows from $k^{3/2}u\le1/256$ and $k\ge64$. Also $p\ge3q+8\ge17$, so $u\le2^{-17}$. Together these give $(2k+1)u\le129/131072$, and hence

$$
\frac{9\gamma}{1-9\gamma}
=\frac{9(2k+1)u}{1-10(2k+1)u}
\le\frac{1161}{129782}<\frac1{100}.
$$

This perturbation is too small to remove an eigenvalue below $-\|\widehat S\|_2/18$. It contradicts $\widehat L\widehat L^T\succeq0$, so the factorization must encounter a zero or negative pivot.

All update products during the leading elimination are exactly representable. A fused multiply-add rounds the product and subtraction only once, but here it gives the same stored $\widehat S$. The sign and error estimates also imply breakdown for that version of the algorithm.
````

### Condition Number and the Remaining Gap

Put $x=k^{3/2}u$, so that $\delta=12x$ and $\lambda_{\min}(S)=9x$. Since

$$
T^TT=G+\frac98\eta J_{2k},
$$

we have $\|T\|_2^2\le9/2+27x$. The block factorization and its inverse give

$$
\|A_{9k}\|_2\le\frac{17}{2}+39x,
\qquad
\|A_{9k}^{-1}\|_2\le1+\frac{17}{72x}.
$$

Therefore

$$
\boxed{
\kappa_2(A_{9k})
\le\left(\frac{17}{2}+39x\right)
\left(1+\frac{17}{72x}\right)
<\frac{57}{n^{3/2}u}.
}
$$

To check the constant $57$, use $n^{3/2}u=27x$. For $0<x\le1/256$, the scaled upper bound satisfies

$$
27\left(\frac{17}{2}+39x\right)
\left(x+\frac{17}{72}\right)
\le\frac{3674685}{65536}<57.
$$

We also need a lower bound to see how small this family's condition number can be. Write $A=A_{9k}$ and

$$
M=[2I_{7k}\;\;T],\qquad
N=\begin{bmatrix}-T/2\\I_{2k}\end{bmatrix}.
$$

Here $X\succeq Y$ means that $X-Y$ is positive semidefinite. Since $T^TT\succeq G$ and $\lambda_{\max}(G)=9/2$, we have $\|T\|_2^2\ge9/2$. Now $A\succeq M^TM$ and $MM^T=4I_{7k}+TT^T$, so $\lambda_{\max}(A)\ge4+9/2=17/2$. Also,

$$
A^{-1}=\begin{bmatrix}I_{7k}/4&0\\0&0\end{bmatrix}
+NS^{-1}N^T
\succeq\frac1{\delta}NN^T.
$$

It follows that $\|A^{-1}\|_2\ge17/(8\delta)$ and thus

$$
\frac{2601/64}{n^{3/2}u}
\le\kappa_2(A_{9k})
<\frac{57}{n^{3/2}u}.
$$

The lower and upper bounds have the same dependence on $n$ and $u$, but their constants differ. They do not establish an asymptotic equivalence. Combining the upper bound with the sufficient condition for success gives

$$
\frac{1}{3.9n^{3/2}}\le c^*_{n,u}<\frac{57}{n^{3/2}}
$$

for the admissible dimensions. The exact best multiplicative constant remains undetermined.

### Bounded Condition Numbers at Fixed Precision

To obtain a condition-number bound independent of precision, choose the largest admissible basic block. Assume $p\ge20$, let $q_*$ be the largest integer satisfying the construction's restrictions, and set $k_*=4^{q_*}$. The next candidate, $q_*+2$, violates $3q+8\le p$. Since $p+q_*$ is even, we get

$$
r=p-3q_*\in\{8,10,12\},\qquad
x_*=k_*^{3/2}u=2^{-r}\in\{2^{-8},2^{-10},2^{-12}\}.
$$

The smallest possible $x_*$ is therefore $2^{-12}$, regardless of the precision. Substituting these three possible values into the condition-number estimate gives

$$
\kappa_2(A_{9k_*})\le
\left(\frac{17}{2}+39x_*\right)
\left(1+\frac{17}{72x_*}\right)
\le\frac{303691615}{36864}<8239.
$$

Expanding the upper bound as $289/(144x_*)+425/24+39x_*$ shows that its largest value occurs at $x_*=2^{-12}$. Thus $8239$ bounds the condition number independently of both dimension and precision.

For every $N\ge9k_*$, define

$$
B_N=\operatorname{diag}(A_{9k_*},4I_{N-9k_*}).
$$

Cholesky fails within the original block, before reaching the added diagonal block. The new eigenvalues are all $4$, which lies between the smallest and largest eigenvalues of the original block. Thus $\kappa_2(B_N)=\kappa_2(A_{9k_*})<8239$ for every larger size, with $u$ fixed.

For binary64, $p=53$, $q_*=15$, and $x_*=2^{-8}$. Hence

$$
N\ge9\cdot2^{30}
\quad\Longrightarrow\quad
\kappa_2(B_N)\le\frac{1224895}{2304}<532.
$$

This binary64 block is far too large to store as a dense matrix in practice. For binary32, the largest admissible block has size $2304$. The numerical example below verifies its failure, and its condition-number bound is below $8239$.

### Sharpness for Decreasing and Bounded Condition Numbers

At fixed $u$, the target bound decreases like $1/(N^{3/2}u)$ until $N^{3/2}u=1$, then stays constant. To cover every dimension, choose the largest admissible block that fits and add a diagonal block as needed.

For any $p\ge20$ and $N\ge2304$, choose the largest admissible $q$ for which $m=9\cdot4^q\le N$, and form

$$
B_N=\operatorname{diag}(A_m,4I_{N-m}).
$$

If $q=q_*$, the preceding bound gives $\kappa_2(B_N)<8239$. Otherwise the next admissible size is $16m>N$, and therefore

$$
\kappa_2(B_N)<\frac{57}{m^{3/2}u}
<\frac{3648}{N^{3/2}u}.
$$

Combining the two cases proves the uniform bound

$$
\boxed{\kappa_2(B_N)<8239\max\left\{1,\frac{1}{N^{3/2}u}\right\}.}
$$

The sufficient condition for success and the inequality $\kappa_2(A)\ge1$ give lower bounds on $c^*_{N,u}$. The failing matrices $B_N$ give an upper bound. Together,

$$
\boxed{
\max\left\{u,\frac{1}{3.9N^{3/2}}\right\}
\le c^*_{N,u}
<8239\max\left\{u,\frac{1}{N^{3/2}}\right\},
\qquad p\ge20,\quad N\ge2304.
}
$$

Here $N\ge2304\ge13$, so the sufficient condition derived from Theorem 2 applies. When $1/(3.9N^{3/2})<u$, the lower bound $c^*_{N,u}\ge u$ suffices.

Dividing by $u$ shows that the smallest condition number permitting failure has order $\max\{1,1/(N^{3/2}u)\}$ in the specified model. The constants are independent of dimension and precision, but they are not claimed to be optimal.

### Operation Order and Condition Numbers Near One

The breakdown proved here depends on the order of operations. Our algorithm subtracts each product from the matrix entry as it is computed. An alternative is to accumulate the products first and then subtract their sum. These evaluations can round differently:

$$
\operatorname{fl}\!\left(\operatorname{fl}(a-b)-c\right)
\quad\text{can differ from}\quad
\operatorname{fl}\!\left(a-\operatorname{fl}(b+c)\right).
$$

In both numerical examples below, successive updates lead to failure, whereas accumulating each inner product sequentially from zero and then subtracting it succeeds. The counterexamples therefore do not establish failure for every implementation of Cholesky.

Our construction cannot produce condition numbers arbitrarily close to $1$. For all admissible parameters, the lower bound proved above gives

$$
\kappa_2(A_{9k})\ge\frac{289}{192x}
\ge\frac{1156}{3}>385,
\qquad x=k^{3/2}u\le\frac1{256}.
$$

Padding preserves the condition number, so increasing the dimension this way cannot close the gap. Whether a different construction can produce numerical breakdown with $\kappa_2(A)\to1$ is not resolved by this analysis.

### Numerical Examples

These examples use the same order of operations as `cholesky_in_place`. The stored entries and the exact Schur complement were checked using rational arithmetic; the computed trailing block was also checked entry by entry against the formula for $\widehat S$. Condition-number bounds in the table are rounded upward.

| Arithmetic | Size $n$ | Upper bound on $\kappa_2(A_n)$ | Failing pivot (counting from 1) | Computed pivot |
|---|---:|---:|---:|---:|
| Binary64, $u=2^{-53}$ | 576 | $3.531\times10^{13}$ | 565 | $-1.66255\times10^{-12}$ |
| Binary32, $u=2^{-24}$ | 2304 | $8239$ | 2252 | $-1.14055\times10^{-3}$ |

For the binary32 example, $k^{3/2}u=1/4096$, so the construction's size restriction is satisfied at the fixed precision $u=2^{-24}$.

The following code constructs either example and applies the earlier `cholesky_in_place` function.

```python
import math
from fractions import Fraction
import numpy as np


def cholesky_breakdown_matrix(q=3, dtype=np.float64):
    assert dtype in (np.float32, np.float64)
    p = np.finfo(dtype).nmant + 1
    u = 2.0**-p
    assert q >= 3 and (p + q) % 2 == 0
    assert 3 * q + 8 <= p  # Equivalent to k**1.5 * u <= 1/256.
    k = 4**q
    root_k = 2**q

    H = np.ones((1, 1), dtype=dtype)
    while H.shape[0] < k:
        H = np.block([[H, H], [H, -H]])
    Q = H / dtype(root_k)

    h = 2 * u / root_k
    eta = 6 * k * h
    delta = k * eta
    t = math.ldexp(3.0, -(p + q + 2) // 2)
    assert Fraction(t)**2 == Fraction(9, 8) * Fraction(h)

    I = np.eye(k, dtype=dtype)
    T = np.vstack((np.full((6 * k, 2 * k), t, dtype=dtype),
                   dtype(1.5) * np.hstack((I, Q))))
    G = dtype(2.25) * np.block([[I, Q], [Q.T, I]])
    C = G + dtype(eta) * np.ones((2 * k, 2 * k), dtype=dtype)
    C += dtype(delta) * np.eye(2 * k, dtype=dtype)

    # Check every stored entry of C, grouping equal base values in G.
    for base in np.unique(G):
        expected = Fraction(float(base)) + Fraction(eta)
        if base == 2.25:  # Only diagonal entries have this base value.
            expected += Fraction(delta)
        for stored in np.unique(C[G == base]):
            assert Fraction(float(stored)) == expected

    assert Fraction(delta) - Fraction(k, 4) * Fraction(eta) > 0
    return np.block([[dtype(4) * np.eye(7 * k, dtype=dtype), dtype(2) * T],
                     [dtype(2) * T.T, C]])


A = cholesky_breakdown_matrix()
# For the binary32 example: A = cholesky_breakdown_matrix(4, np.float32)
try:
    cholesky_in_place(A.copy())
except np.linalg.LinAlgError as error:
    print(error)
# Binary64: Nonpositive or nonfinite pivot at step 565
# Binary32: Nonpositive or nonfinite pivot at step 2252
```

### Failure Is Not Monotone in Precision

Changing $u$ in our construction changes the matrix. A smaller $u$ makes its positive Schur complement smaller as well. This does not describe what happens when we factor the same SPD matrix at higher precision. For a fixed matrix and dimension, sufficiently small $u$ guarantees success.

Success at one precision does not guarantee success at the next higher precision. For example, keep the following matrix fixed:

$$
A=\frac1{32}\begin{bmatrix}4&5\\5&7\end{bmatrix}.
$$

Its determinant is $3/1024>0$, and all entries are exactly representable in both of the following binary formats. With separate rounding of square root, division, multiplication, and subtraction, Cholesky gives

| Significand precision | $\widehat l_{11}$ | $\widehat l_{21}$ | Computed second pivot |
|---|---:|---:|---:|
| $p=3$, $u=1/8$ | $3/8$ | $7/16$ | $1/32>0$ |
| $p=4$, $u=1/16$ | $11/32$ | $15/32$ | $0$ |

Here the coarser calculation succeeds and the finer one fails. Individual rounding errors can cancel differently, so failure at one precision does not imply failure at every coarser precision.
