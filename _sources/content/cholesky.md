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

For a symmetric positive definite matrix, symmetry and positive definiteness allow elimination to use a single triangular factor:

$$
\boxed{A=LL^T,\qquad l_{ii}>0.}
$$

This is the **Cholesky factorization**. It requires about half the arithmetic of general LU factorization and can store the matrix and factor in a single triangle. No pivoting is needed in exact arithmetic. In floating-point arithmetic, a completed Cholesky factorization has a small backward error; as with LU, the accuracy of the solution also depends on conditioning.

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

The key to continuing the algorithm is that $S$ is also SPD.

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

For a concrete example, take

$$
A=\begin{pmatrix}4&2&2\\2&5&3\\2&3&6\end{pmatrix}.
$$

The first column of $L$ is $(2,1,1)^T$, and the first Schur complement is

$$
\begin{pmatrix}5&3\\3&6\end{pmatrix}
-\begin{pmatrix}1\\1\end{pmatrix}\begin{pmatrix}1&1\end{pmatrix}
=\begin{pmatrix}4&2\\2&5\end{pmatrix}.
$$

Continuing gives

$$
L=\begin{pmatrix}2&0&0\\1&2&0\\1&1&2\end{pmatrix}.
$$

The next cell factors $A$ and solves for $b=(6,5,8)^T$ using forward and backward substitution. Copies and an explicit $L$ are used here so that we can check the result against the original data.

```{code-cell} ipython3
A = np.array([[4., 2., 2.],
              [2., 5., 3.],
              [2., 3., 6.]])
b = np.array([6., 5., 8.])

packed = cholesky_in_place(A.copy())
L = np.tril(packed)
n = len(b)
y = np.zeros(n)
for i in range(n):
    y[i] = (b[i] - L[i, :i] @ y[:i]) / L[i, i]
xhat = np.zeros(n)
for i in range(n - 1, -1, -1):
    xhat[i] = (y[i] - L[i + 1:, i] @ xhat[i + 1:]) / L[i, i]

r = b - A @ xhat
eta = np.linalg.norm(r) / (
    np.linalg.norm(A, 2) * np.linalg.norm(xhat) + np.linalg.norm(b)
)
assert np.allclose(L @ L.T, A)
assert np.allclose(xhat, [1., 0., 1.])
print("L =")
print(L)
print("Computed solution:", xhat)
print(f"Relative backward error: {eta:.3e}")
```

The factorization and triangular solves in this example are exact in binary arithmetic, so the backward error is zero. For general data, rounding errors enter both stages. Their effect depends on the sizes of the terms formed during the computation.

## Backward Stability and Its Limits

### Why the Factor Entries Stay Controlled

In exact arithmetic, the diagonal of $A=LL^T$ gives

$$
a_{ii}=\sum_{k=1}^i l_{ik}^2.
$$

Every term is nonnegative, so

$$
\boxed{|l_{ik}|\le\sqrt{a_{ii}}.}
$$

More generally, Cauchy–Schwarz bounds the products that enter the backward-error analysis:

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

### Cancellation and Breakdown

Bounded factor entries do not rule out cancellation during the computation. A diagonal step computes

$$
l_{kk}=\sqrt{a_{kk}-\sum_{j<k}l_{kj}^2},
$$

where $a_{kk}$ in this formula is the original diagonal entry. The difference inside the square root can be small relative to its terms. For example,

$$
A=\begin{pmatrix}1&1\\1&1+\epsilon\end{pmatrix},
\qquad
L=\begin{pmatrix}1&0\\1&\sqrt{\epsilon}\end{pmatrix},
\qquad \epsilon>0.
$$

The second pivot is obtained by subtracting $1$ from $1+\epsilon$. If $\epsilon$ is small enough, input rounding alone can replace $1+\epsilon$ by $1$, making the stored matrix singular. Even for an SPD stored matrix, errors during factorization can produce a nonpositive pivot when the matrix is sufficiently close to singularity. Exact existence and successful floating-point completion are separate statements; see [Higham's discussion of Cholesky factorization](https://nhigham.com/2020/08/11/what-is-a-cholesky-factorization/).

### Backward Error of the Factorization

When the algorithm completes, the same control of the factor entries leads to a backward-error bound. The bound must account for rounding in the **computed** factor $\widehat{L}$.

````{prf:theorem} Backward Error of Cholesky Factorization
:label: thm:backward_error_cholesky
Assume real arithmetic with rounding to nearest, no overflow or underflow, and $(n+1)u<1$. If Cholesky factorization of an SPD matrix $A$ completes with a positive diagonal, the computed factor satisfies

$$
A+E=\widehat{L}\widehat{L}^T,
\qquad
|E|\le\gamma_{n+1}|\widehat{L}|\,|\widehat{L}|^T,
$$

where $E$ is symmetric, absolute values and inequalities are entrywise, and $\gamma_m=mu/(1-mu)$.
````

This is the Cholesky counterpart of the [LU backward-error theorem](lu_pivoting.md); the componentwise bound is stated by [Rump and Jeannerod](https://doi.org/10.1137/130927231). Its right-hand side involves the computed factor. Bounding that factor in terms of $A$ gives a normwise error bound without a growth-factor assumption.

Write $\gamma=\gamma_{n+1}$ and assume $\gamma<1$. The diagonal equations imply

$$
\sum_k\widehat{l}_{ik}^2=a_{ii}+e_{ii}
\le a_{ii}+\gamma\sum_k\widehat{l}_{ik}^2.
$$

Summing over $i$ yields

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

For small $nu$, this bound is approximately $n(n+1)u$. The dimension factor is a worst-case bound; there is no additional factor from uncontrolled element growth.
```

### Backward Error of the Computed Solution

The backward error of the computed solution $\widehat{x}$ accounts for rounding errors in the factorization and both triangular solves. For $b\ne0$, it is

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

The factorization contributes $\gamma_{n+1}$; the two solves add $2\gamma_n+\gamma_n^2$. Their sum is at most $\gamma_{3n+1}$. Here $\Delta A$ includes all three stages and need not be symmetric.

Since this perturbation makes $\widehat{x}$ exact without changing $b$, the definition of backward error gives

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

For small $nu$, $\beta_n\approx n(3n+1)u$. This gives an explicit backward-error guarantee for the complete solve, without a growth-factor assumption.
```

### Forward Error of the Computed Solution

The bound on $\eta(\widehat{x})$ controls the input perturbation needed to explain the computed answer. The condition number determines how much that perturbation can affect the solution. Substituting $\eta(\widehat{x})\le\beta_n$ into the [forward-error bound for linear systems](lu_pivoting.md) gives, when $\kappa_2(A)\beta_n<1$,

$$
\frac{\|\widehat{x}-x\|_2}{\|x\|_2}
\le\frac{2\kappa_2(A)\beta_n}{1-\kappa_2(A)\beta_n}.
$$

Thus $\kappa_2(A)\beta_n\ll1$ guarantees a small relative solution error. For SPD matrices,

$$
\kappa_2(A)=\frac{\lambda_{\max}(A)}{\lambda_{\min}(A)}.
$$

An eigenvalue that is small relative to the largest one can therefore make the solution sensitive, despite the backward-error guarantee.

## Choosing a Direct Solver

For a general dense nonsingular system, use LU with partial pivoting. When the matrix is known to be symmetric positive definite, Cholesky exploits that structure to reduce arithmetic and storage. In both cases, factor once, solve triangular systems for each right-hand side, and interpret the residual together with the condition number when assessing solution accuracy.
