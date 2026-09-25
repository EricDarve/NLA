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

# LU Factorization with Row Pivoting

LU factorization can fail at a zero pivot even when the matrix is invertible. A nonzero pivot can also cause trouble if it is much smaller than the entries below it: elimination then forms large intermediate quantities whose cancellation may expose rounding errors.

**Partial pivoting chooses the largest available entry in the current column as the pivot.** This guarantees that the elimination multipliers have magnitude at most 1 in exact arithmetic. It also guarantees that a row-permuted LU factorization exists for every square matrix.

We first describe the algorithm and its implementation. We then introduce backward error and sensitivity to assess the accuracy of the computed solution. Finally, we examine what partial pivoting controls and construct an example where it produces a large backward error.

## Why Interchange Rows?

Recall the example from [the discussion of LU existence](existence_lu.md):

$$
A_\epsilon=\begin{pmatrix}\epsilon&1\\1&\pi\end{pmatrix},
\qquad 0<\epsilon\ll1.
$$

Without pivoting, its exact factors are

$$
L=\begin{pmatrix}1&0\\1/\epsilon&1\end{pmatrix},
\qquad
U=\begin{pmatrix}\epsilon&1\\0&\pi-1/\epsilon\end{pmatrix}.
$$

The $(2,2)$ entry of $A_\epsilon$ is reconstructed by the cancellation

$$
\pi=\frac{1}{\epsilon}+\left(\pi-\frac{1}{\epsilon}\right).
$$

The two terms have magnitude about $1/\epsilon$, but their sum has magnitude about 1. This is the scale comparison from [floating-point arithmetic](floating_point.md): errors introduced at the large intermediate scale can overwhelm the small result. For example, at $\epsilon=10^{-18}$, double precision loses the contribution of $\pi$ when forming $\pi-1/\epsilon$.

Interchanging the rows gives

$$
P=\begin{pmatrix}0&1\\1&0\end{pmatrix},
\qquad
PA_\epsilon
=\begin{pmatrix}1&\pi\\\epsilon&1\end{pmatrix}.
$$

Now the factors are

$$
L=\begin{pmatrix}1&0\\\epsilon&1\end{pmatrix},
\qquad
U=\begin{pmatrix}1&\pi\\0&1-\epsilon\pi\end{pmatrix}.
$$

All entries remain moderate. The row swap avoids the large intermediate quantities and the subsequent cancellation in this example.

## Partial Pivoting and the Factorization $PA=LU$

We perform elimination **in place**, overwriting the entries of the array $A$ as we proceed. In the algorithm below, $a_{ij}$ always means the current value of that entry, including updates and row interchanges from earlier steps. For $k=1,\ldots,n-1$, perform the following steps:

1. Choose a row $p\in\{k,\ldots,n\}$ such that

   $$
   |a_{pk}|=\max_{i=k,\ldots,n}|a_{ik}|.
   $$

2. Interchange the entire rows $k$ and $p$, including multipliers stored in earlier columns.
3. If the pivot is nonzero, compute the multipliers, store them below the pivot, and update the trailing block:

   $$
   \begin{aligned}
   l_{ik}&=a_{ik}/a_{kk}, &&i>k,\\
   a_{ik}&\leftarrow l_{ik}, &&i>k,\\
   a_{ij}&\leftarrow a_{ij}-a_{ik}a_{kj}, &&i,j>k.
   \end{aligned}
   $$

The overwritten entries $a_{ik}$ store the multipliers in $L$; the pivot row supplies the corresponding row of $U$. Because the pivot was chosen to have largest magnitude in its active column,

$$
|l_{ik}|\le1.
$$

If the selected pivot is zero, every entry below it in that column is also zero. No elimination is needed there: set those multipliers to zero and continue. This case matters for singular matrices. For a nonsingular matrix, partial pivoting encounters no zero pivots in exact arithmetic.

The product of the row interchanges is a permutation matrix $P$, and the resulting factorization is

$$
\boxed{PA=LU.}
$$

In this identity, $A$ denotes the **original input matrix**. Its storage now contains $L$ below the diagonal and $U$ on and above it, with the unit diagonal of $L$ implicit. To solve $Ax=b$, apply the same permutation to the right-hand side and solve

$$
Ly=Pb,\qquad Ux=y.
$$

The original unknowns have not been reordered, so there is no permutation to undo in $x$.

### Why the Factorization Always Exists

````{prf:theorem} Existence of LU Factorization with Row Pivoting
:label: thm:existence_lu_pivoting
For every $A\in\mathbb{F}^{n\times n}$, where $\mathbb{F}=\mathbb{R}$ or $\mathbb{C}$, there exist a permutation matrix $P$, a unit lower triangular matrix $L$, and an upper triangular matrix $U$ such that $PA=LU$. No assumption of nonsingularity is needed.
````

````{prf:proof} Existence with Row Pivoting.
We use induction on $n$. For $n=1$, take $P=L=[1]$ and $U=A$.

For $n>1$, choose a largest-magnitude entry in the first column and move it to the first row with a permutation $P_1$. Write

$$
P_1A=\begin{pmatrix}\alpha&w\\v&B\end{pmatrix},
$$

where $w$ is a row vector and $v$ is a column vector. If $\alpha\ne0$, set $\ell=v/\alpha$. If $\alpha=0$, the pivot choice implies $v=0$, and we set $\ell=0$. In either case, $v=\alpha\ell$. Define the remaining block by $S=B-\ell w$.

By induction, $QS=L_SU_S$ for a permutation $Q$ and triangular factors $L_S,U_S$. Set

$$
\begin{aligned}
P&=\begin{pmatrix}1&0\\0&Q\end{pmatrix}P_1,\\
L&=\begin{pmatrix}1&0\\Q\ell&L_S\end{pmatrix},
\qquad
U=\begin{pmatrix}\alpha&w\\0&U_S\end{pmatrix}.
\end{aligned}
$$

Then

$$
LU
=\begin{pmatrix}\alpha&w\\Qv&Q\ell w+QS\end{pmatrix}
=\begin{pmatrix}\alpha&w\\Qv&QB\end{pmatrix}
=PA.
$$

The matrices have the required triangular structure, completing the induction.
````

The factor $Q\ell$ explains an implementation detail: later row interchanges must also reorder the multipliers already stored in $L$. For singular $A$, the factorization still exists, but $U$ has a zero diagonal entry and ordinary backward substitution cannot produce a unique solution for every right-hand side.

### In-Place Implementation

As in the [unpivoted algorithm](lu_decomposition.md), the factors share the original array's $n^2$ storage locations. Avoiding separate arrays for $A$, $L$, and $U$ matters for large matrices: those extra copies may exceed the available memory. The following function modifies a square NumPy array in place and returns $P$.

For clarity, this teaching implementation forms a full permutation matrix and a temporary array for the outer product. Implementations designed to limit memory use store row indices and update the trailing matrix in blocks.

```{code-cell} ipython3
import numpy as np

def lu_factorization_with_row_pivoting(A: np.ndarray) -> np.ndarray:
    """Overwrite A with packed LU; return P such that P @ A_original = L @ U."""
    if A.ndim != 2 or A.shape[0] != A.shape[1]:
        raise ValueError("A must be square")
    if not np.issubdtype(A.dtype, np.inexact):
        raise TypeError("A must have a floating-point or complex dtype")
    if not np.all(np.isfinite(A)):
        raise ValueError("A must contain finite entries")

    n = A.shape[0]
    P = np.eye(n)
    for k in range(n - 1):
        p = k + np.argmax(np.abs(A[k:, k]))
        if p != k:
            # Swap entire rows, including multipliers from earlier steps.
            A[[k, p], :] = A[[p, k], :]
            P[[k, p], :] = P[[p, k], :]

        if A[k, k] == 0:
            # The entire active column is zero; no elimination is needed.
            continue

        A[k + 1:, k] /= A[k, k]
        A[k + 1:, k + 1:] -= np.outer(A[k + 1:, k], A[k, k + 1:])
    return P
```

Pass a copy to preserve the original matrix, and convert integer data to floating-point storage before calling this in-place routine. The zero-pivot test permits a factorization of a singular matrix; it is not a test for numerical rank or for a well-conditioned solve.

The small example below keeps a copy of $A$ and extracts $L$ and $U$ to check the factorization, using the triangular-solve functions from [LU decomposition](lu_decomposition.md). These extra arrays are for demonstration; a solve can use the packed factors directly.

```{code-cell} ipython3
:tags: [hide-input]

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

```{code-cell} ipython3
A = np.array([[2., 1., 1.],
              [4., 3., 3.],
              [8., 7., 9.]])
b = np.array([3., 7., 17.])

packed = A.copy()
P = lu_factorization_with_row_pivoting(packed)
L = np.tril(packed, k=-1) + np.eye(A.shape[0])
U = np.triu(packed)
y = forward_substitution(L, P @ b)
x = backward_substitution(U, y)

assert np.allclose(P @ A, L @ U)
assert np.allclose(x, [1., 0., 1.])
print("Computed solution:", x)
```

This example swaps rows at both elimination steps, so the second swap must move the multipliers from the first step. Pivot searches and row swaps add $O(n^2)$ work; the leading arithmetic cost remains $\frac23n^3$ flops for a dense real factorization.

## Forward Error, Backward Error, and the Residual

We can now compute a solution, but how accurate is it? For the error analysis, assume $A$ is nonsingular and $b\ne0$. Let $x$ be the exact solution of $Ax=b$, and let $\widehat{x}$ be the computed solution. The hat marks a computed quantity that may contain rounding error. Unless a subscript specifies otherwise, we use the Euclidean vector norm and its induced matrix norm, $\|\cdot\|_2$.

### Forward Error

The **forward error** measures the difference between the computed and exact answers:

$$
\begin{aligned}
\text{absolute forward error}&=\|\widehat{x}-x\|,\\
\text{relative forward error}&=\frac{\|\widehat{x}-x\|}{\|x\|}.
\end{aligned}
$$

**Forward error is what we ultimately want to control:** it tells us how accurate our answer is. But we usually do not know $x$, and directly tracking how rounding errors affect the final answer is difficult. Instead, we first derive a backward-error bound, which is often easier to obtain. Sensitivity analysis then relates this change in the input to the forward error in the solution.

### Backward Error

The **backward error** measures how much we must change the input data to make the computed answer exact. For a linear system, the data are $A$ and $b$. We define the relative backward error $\eta(\widehat{x})$ as the smallest $\eta$ for which

$$
(A+E)\widehat{x}=b+e,
\qquad
\|E\|\le\eta\|A\|,
\qquad
\|e\|\le\eta\|b\|.
$$

Thus $\eta(\widehat{x})=10^{-8}$ means that relative changes of at most $10^{-8}$ in the data suffice to make $\widehat{x}$ exact.

We can compute this quantity without knowing the true solution $x$. First form the **residual**, which measures how far the computed answer is from satisfying the equations:

$$
r=b-A\widehat{x}.
$$

Then the [relative backward error](https://eprints.maths.manchester.ac.uk/2562/1/paper.pdf#page=6) is exactly

$$
\boxed{\eta(\widehat{x})=
\frac{\|r\|}{\|A\|\,\|\widehat{x}\|+\|b\|}.}
$$

The denominator accounts for changes in both inputs: a perturbation to $A$ acts on $\widehat{x}$, while a perturbation to $b$ changes the right-hand side directly. This scaled residual measures the size of the input changes needed to explain the error in the equations.

```{admonition} Key definition: Backward stability
:class: important

**A method is backward stable if it guarantees $\eta(\widehat{x})$ of order the unit roundoff $u$**, allowing for a modest factor depending on the dimension.
```

Computing the residual costs $O(n^2)$ operations for a dense system. A small backward error tells us that the computed answer solves a nearby problem. To decide whether it is close to the answer we wanted, we must ask how much the solution changes when the data change.

## Sensitivity: Relating Backward Error to Forward Error

For a general problem $x=f(d)$, let $d$ denote the input data and $f(d)$ the exact answer. **Sensitivity measures how much the answer changes when the input changes.** The **local absolute sensitivity** is

$$
S_f(d)=\lim_{\varepsilon\to0^+}
\sup_{0<\|\delta d\|\le\varepsilon}
\frac{\|f(d+\delta d)-f(d)\|}{\|\delta d\|}.
$$

The ratio measures the change in the answer per unit change in the input. The supremum selects the perturbation with the largest amplification, and the limit restricts attention to small changes near $d$. If $f$ is differentiable at $d$, then $S_f(d)=\|Df(d)\|$, the induced norm of its derivative. For a scalar function, this is simply $|f'(d)|$.

Suppose backward error analysis shows that $\widehat{x}=f(d+\delta d)$. For differentiable $f$, the definition gives the local bound

$$
\underbrace{\|\widehat{x}-x\|}_{\text{forward error}}
\le S_f(d)\underbrace{\|\delta d\|}_{\text{backward perturbation}}
+o(\|\delta d\|).
$$

Thus, to first order,

$$
\boxed{\text{forward error}\ \lesssim\ \text{sensitivity}\times\text{backward error}.}
$$

To compare relative changes, assuming $d\ne0$ and $f(d)\ne0$, define the **relative condition number**

$$
\kappa_f(d)=S_f(d)\frac{\|d\|}{\|f(d)\|}.
$$

It bounds the amplification of relative input errors into relative output errors for small perturbations. A problem with a large relative condition number is **ill-conditioned**; one with a modest condition number is **well-conditioned**.

Sensitivity is a property of the problem $f$ at the input $d$. Backward stability is a property of the algorithm. We need both to assess the accuracy of a computed answer.

## Conditioning of Linear Systems

We now apply this framework to $Ax=b$. The input data are $A$ and $b$, and the output is $x=A^{-1}b$. We first relate the residual to the solution error, then prove a bound for perturbations in both inputs.

### The Condition Number and Backward Error

Let $\delta x=\widehat{x}-x$. Since $r=b-A\widehat{x}$ and $Ax=b$,

$$
A\delta x=-r,
\qquad
\delta x=-A^{-1}r.
$$

Therefore $\|\delta x\|\le\|A^{-1}\|\,\|r\|$. Our backward-error formula gives

$$
\|r\|=\eta(\widehat{x})
\bigl(\|A\|\,\|\widehat{x}\|+\|b\|\bigr).
$$

Define the **matrix condition number** by

$$
\boxed{\kappa(A)=\|A\|\,\|A^{-1}\|.}
$$

Substituting the residual identity and using $\|b\|\le\|A\|\,\|x\|$ yields

$$
\frac{\|\delta x\|}{\|x\|}
\le \kappa(A)\eta(\widehat{x})
\left(1+\frac{\|\widehat{x}\|}{\|x\|}\right).
$$

```{admonition} Key result: Conditioning and backward error
:class: important

Here the product **$\kappa(A)\eta(\widehat{x})$** connects the problem's sensitivity to the algorithm's backward error. When $\widehat{x}$ is close to $x$, the factor in parentheses is close to 2, giving the first-order bound

$$
\frac{\|\widehat{x}-x\|}{\|x\|}
\lesssim 2\kappa(A)\eta(\widehat{x}).
$$
```

The **condition number** satisfies $\kappa(A)\ge1$, since $1=\|I\|\le\|A\|\,\|A^{-1}\|$. In the 2-norm, the SVD gives

$$
\kappa_2(A)=\frac{\sigma_{\max}(A)}{\sigma_{\min}(A)}.
$$

Geometrically, $\kappa_2(A)$ is the ratio of the largest to the smallest stretching of a unit vector by $A$. A large value means that small relative changes in the data can produce large relative changes in the solution.

For example, take $A=\operatorname{diag}(1,\epsilon)$ with $0<\epsilon\ll1$, and $b=(1,0)^T$. The solution is $x=(1,0)^T$. Changing $b$ by $(0,\epsilon)^T$ changes the solution by $(0,1)^T$: a relative input change of $\epsilon$ causes a relative solution change of 1. Here $\kappa_2(A)=1/\epsilon$.

Row permutations preserve singular values, so $\kappa_2(PA)=\kappa_2(A)$. Pivoting changes the elimination process; it does not remove the sensitivity of the original problem.

### A General Forward-Error Bound

We next consider finite perturbations $E$ and $e$ in $(A+E)\widehat{x}=b+e$. The first step is to show that $A+E$ remains invertible when $E$ is sufficiently small. The following lemma supplies both invertibility and a bound on the inverse.

````{prf:lemma} Banach Lemma
:label: lem:banach
If $\|X\|<1$ in an induced matrix norm, then $I+X$ is invertible and

$$
\|(I+X)^{-1}\|\le\frac{1}{1-\|X\|}.
$$
````

````{prf:proof} Banach Lemma.
The series $S=\sum_{j=0}^{\infty}(-X)^j$ converges because $\|X^j\|\le\|X\|^j$. Its partial sums satisfy

$$
(I+X)\sum_{j=0}^{m}(-X)^j=I-(-X)^{m+1}.
$$

Taking the limit gives $(I+X)S=I$, so $S=(I+X)^{-1}$. Finally,

$$
\|S\|\le\sum_{j=0}^{\infty}\|X\|^j
=\frac{1}{1-\|X\|}.
$$
````

````{prf:theorem} Perturbation Bound for Linear Systems
:label: thm:perturbation_bound
Let $A$ be nonsingular, $b\ne0$, and $Ax=b$. If $\|A^{-1}\|\,\|E\|<1$, then $A+E$ is nonsingular, and the solution of $(A+E)\widehat{x}=b+e$ satisfies

$$
\frac{\|\widehat{x}-x\|}{\|x\|}
\le
\frac{\kappa(A)}{1-\kappa(A)\|E\|/\|A\|}
\left(\frac{\|E\|}{\|A\|}+\frac{\|e\|}{\|b\|}\right).
$$
````

````{prf:proof} Perturbation Bound.
Write $A+E=A(I+A^{-1}E)$. Since $\|A^{-1}E\|\le\|A^{-1}\|\,\|E\|<1$, the Banach lemma proves invertibility and gives

$$
\|(A+E)^{-1}\|
\le\frac{\|A^{-1}\|}{1-\|A^{-1}\|\,\|E\|}.
$$

Let $\delta x=\widehat{x}-x$. Subtracting the two systems gives

$$
(A+E)\delta x=e-Ex.
$$

Consequently,

$$
\frac{\|\delta x\|}{\|x\|}
\le\frac{\|A^{-1}\|}{1-\|A^{-1}\|\,\|E\|}
\left(\frac{\|e\|}{\|x\|}+\|E\|\right).
$$

Use $\|b\|\le\|A\|\,\|x\|$ and $\kappa(A)=\|A\|\,\|A^{-1}\|$ to obtain the stated bound.
````

To express the result using our backward error, choose perturbations with
$\|E\|/\|A\|\le\eta(\widehat{x})$ and $\|e\|/\|b\|\le\eta(\widehat{x})$. Substitution gives the following bound.

```{admonition} Key result: From backward error to forward error
:class: important

**If $\kappa(A)\eta(\widehat{x})<1$, then**

$$
\frac{\|\widehat{x}-x\|}{\|x\|}
\le\frac{2\kappa(A)\eta(\widehat{x})}
{1-\kappa(A)\eta(\widehat{x})}.
$$

When $\kappa(A)\eta(\widehat{x})\ll1$, the denominator is close to 1, recovering the first-order estimate above. Thus a small product of condition number and backward error guarantees a small relative forward error.
```

The factor 2 comes from allowing relative changes of size $\eta$ in each of $A$ and $b$.

## Backward Error of LU and the Role of Cancellation

The preceding results apply to any computed solution. We now ask what backward error LU produces. Let $\widehat{L}$ and $\widehat{U}$ be the computed factors; their hats distinguish them from exact factors.

The following theorem gives separate bounds for the factorization and the complete solve. We state the bounds for real arithmetic; complex arithmetic has analogous bounds with different constants.

````{prf:theorem} Backward Error of LU Factorization and Solution
:label: thm:backward_error_lu
Assume real arithmetic with rounding to nearest, no overflow or underflow, and $3nu<1$. If LU with partial pivoting completes, its computed factors satisfy

$$
PA+F=\widehat{L}\widehat{U},
\qquad
|F|\le\gamma_n|\widehat{L}|\,|\widehat{U}|,
$$

where $\gamma_m=mu/(1-mu)$.

If the computed $\widehat{U}$ has no zero diagonal entries, the solution $\widehat{x}$ obtained by forward and backward substitution satisfies

$$
(A+E)\widehat{x}=b,
\qquad
|PE|\le\gamma_{3n}|\widehat{L}|\,|\widehat{U}|.
$$

Absolute values and inequalities are entrywise. The product $|\widehat{L}|\,|\widehat{U}|$ is an ordinary matrix product of nonnegative matrices. The factor $P$ in $|PE|$ puts the perturbation in the same row ordering as the factors.
````

The solve result fits our definition of $\eta(\widehat{x})$ with the perturbations $E$ and $e=0$. Since $\eta(\widehat{x})$ is the smallest admissible relative perturbation, this particular choice gives $\eta(\widehat{x})\le\|E\|/\|A\|$. Taking norms in the theorem and using $\|PE\|=\|E\|$ yields

$$
\boxed{
\eta(\widehat{x})
\le\frac{\|E\|}{\|A\|}
\le\gamma_{3n}
\frac{\bigl\||\widehat{L}|\,|\widehat{U}|\bigr\|}{\|A\|}.}
$$

For $3nu\ll1$, we have $\gamma_{3n}\approx3nu$. Thus, if $\bigl\||\widehat{L}|\,|\widehat{U}|\bigr\|$ is moderate relative to $\|A\|$, the theorem guarantees $\eta(\widehat{x})$ of order $u$, up to a modest dimension-dependent factor. This is the backward-stability criterion above.

To interpret the bound, compare

$$
\begin{aligned}
(\widehat{L}\widehat{U})_{ij}
&=\sum_k\widehat{l}_{ik}\widehat{u}_{kj},\\
(|\widehat{L}|\,|\widehat{U}|)_{ij}
&=\sum_k|\widehat{l}_{ik}\widehat{u}_{kj}|.
\end{aligned}
$$

The first sum allows cancellation; the second measures the sizes of all terms before cancellation. This is the same distinction as in the [summation error bound](floating_point.md), {prf:ref}`thm:summation_error`. Large terms that cancel to reproduce a much smaller matrix entry can carry rounding errors that are large relative to that entry.

For example, in the unpivoted $A_\epsilon$ factorization, the two terms reconstructing $a_{22}$ have total magnitude about $2/\epsilon$, although $a_{22}=\pi$. After swapping rows, the corresponding products stay moderate, giving a much smaller backward-error bound.

The full rounding-error derivation is given in [Higham, *Accuracy and Stability of Numerical Algorithms*, Chapter 9](https://epubs.siam.org/doi/10.1137/1.9780898718027.ch9); the solve bound is also stated in [Carson and Higham, equation (7.1)](https://eprints.maths.manchester.ac.uk/2562/1/paper.pdf#page=13).

## Element Growth: The Limitation of Partial Pivoting

The LU error bound depends on both factors. Partial pivoting bounds the multipliers in $L$, but entries in $U$ can still grow. To measure this growth, write $a_{ij}^{(k)}$ for the value of entry $(i,j)$ at the start of step $k$, with $a_{ij}^{(1)}$ the original input entries. The superscript records values at different steps; the algorithm still uses a single array. For $A\ne0$, define the **growth factor**

$$
\rho_n
=\frac{\displaystyle\max_{1\le k\le n}\max_{k\le i,j\le n}|a_{ij}^{(k)}|}
{\displaystyle\max_{1\le i,j\le n}|a_{ij}^{(1)}|}.
$$

This definition includes all intermediate trailing matrices, beginning with $A$. It excludes the stored multipliers in a packed implementation. In particular, it bounds the size of every entry eventually placed in $U$.

In exact arithmetic, an update satisfies

$$
|a_{ij}-l_{ik}a_{kj}|
\le |a_{ij}|+|a_{kj}|
$$

because $|l_{ik}|\le1$. Each elimination step can therefore at most double the largest active entry. There are $n-1$ updates, so

$$
\boxed{\rho_n\le2^{n-1}.}
$$

### A Matrix That Attains the Bound

For any $n\ge2$, consider the $n\times n$ matrix

$$
A=\begin{pmatrix}
1&0&\cdots&0&1\\
-1&1&\cdots&0&1\\
\vdots&\ddots&\ddots&\vdots&\vdots\\
-1&\cdots&-1&1&1\\
-1&\cdots&-1&-1&1
\end{pmatrix}.
$$

Its last column consists of ones. In each of the first $n-1$ columns, the diagonal entry is 1, entries below it are $-1$, and entries above it are zero.

At each elimination step the pivot candidates have equal magnitude. If ties are resolved by choosing the first candidate, as in the implementation above, no row swaps occur. Every multiplier is $-1$, so elimination adds the pivot row to each row below it. The last-column entries therefore double at each step. The factors are

$$
L=\begin{pmatrix}
1&0&\cdots&0\\
-1&1&\ddots&\vdots\\
\vdots&\ddots&\ddots&0\\
-1&\cdots&-1&1
\end{pmatrix},
\qquad
U=\begin{pmatrix}
1&0&\cdots&0&1\\
0&1&\cdots&0&2\\
\vdots&\ddots&\ddots&\vdots&\vdots\\
0&\cdots&0&1&2^{n-2}\\
0&\cdots&0&0&2^{n-1}
\end{pmatrix}.
$$

Thus $u_{in}=2^{i-1}$ and $\rho_n=2^{n-1}$: the growth is exponential in the dimension.

### A Right-Hand Side That Produces a Large Backward Error

We can make the **actual backward error** large for this same matrix. Assume binary arithmetic with rounding to nearest, unit roundoff $u$, and no overflow or underflow. Choose

$$
b=e_1+\delta_n e_n,
\qquad
\delta_n=u\,2^{n-3}.
$$

Here $e_i$ is the $i$th coordinate vector. Both nonzero entries of $b$ are represented exactly, as are the computed LU factors. The error will arise in forward substitution.

For the forward solve $Ly=b$, the first $n-1$ components are

$$
y_1=1,
\qquad
y_i=2^{i-2},\quad 2\le i\le n-1.
$$

The last component should be

$$
y_n=\delta_n+\sum_{j=1}^{n-1}y_j
=\delta_n+2^{n-2}.
$$

Accumulate the sum in increasing index order. Each partial sum is a power of two, so this sum is computed exactly. But the next representable floating-point number above $2^{n-2}$ is

$$
2^{n-2}+2u\,2^{n-2}=2^{n-2}+4\delta_n.
$$

There is no representable number strictly between these two values. The exact sum $2^{n-2}+\delta_n$ lies only one quarter of the way from the lower value to the upper one. Rounding to nearest therefore returns the lower value, losing the contribution of $\delta_n$:

$$
\widehat{y}_n=\operatorname{fl}(2^{n-2}+\delta_n)=2^{n-2}.
$$

Backward substitution now gives

$$
\widehat{x}_n=\frac{2^{n-2}}{2^{n-1}}=\frac12,
\qquad
\widehat{x}_1=\frac12,
\qquad
\widehat{x}_i=0\quad(2\le i<n).
$$

These operations are exact. Half the sum of the first and last columns of $A$ is $e_1$, so $A\widehat{x}=e_1$. The residual is therefore

$$
r=b-A\widehat{x}=\delta_n e_n,
$$

and, since $\|\widehat{x}\|_2=1/\sqrt{2}$ and $\|b\|_2=\sqrt{1+\delta_n^2}$, the backward error is

$$
\boxed{
\eta(\widehat{x})
=\frac{\delta_n}
{\|A\|_2/\sqrt{2}+\sqrt{1+\delta_n^2}},
\qquad \delta_n=u\,2^{n-3}.}
$$

This is the backward error of the computed solution, not an upper bound on it.

The matrix norm grows only linearly with $n$:

$$
\lfloor n/2\rfloor\le\|A\|_2\le\|A\|_F\le n.
$$

For the lower bound, use the square block of $-1$ entries in the bottom-left corner, of size $\lfloor n/2\rfloor$; its 2-norm is $\lfloor n/2\rfloor$.

Consequently, while $u2^n\ll n$, the actual backward error grows like $u2^n/n$, up to constant factors. Once $\delta_n$ is much larger than $n$, it approaches 1. It cannot grow indefinitely: the scaled residual defining $\eta$ is always at most 1.

In double precision, $u=2^{-53}$. The following cell computes the backward errors using the pivoted LU routine and triangular-solve functions above. The explicit forward-solve loop fixes the summation order used in the derivation. The output shows how the actual backward error grows with $n$.

```{code-cell} ipython3
u = 2.0**-53
for n in (10, 20, 30, 40, 50, 60, 70):
    A = np.eye(n) - np.tril(np.ones((n, n)), k=-1)
    A[:, -1] = 1.0
    b = np.zeros(n)
    b[0], b[-1] = 1.0, u * 2.0**(n - 3)

    packed = A.copy()
    P = lu_factorization_with_row_pivoting(packed)
    L = np.tril(packed, k=-1) + np.eye(n)
    U = np.triu(packed)
    pb = P @ b
    y = np.zeros(n)
    for i in range(n):
        s = 0.0
        for j in range(i):
            s += L[i, j] * y[j]
        y[i] = pb[i] - s  # L has unit diagonal.
    xhat = backward_substitution(U, y)

    r = b - A @ xhat
    eta = np.linalg.norm(r) / (
        np.linalg.norm(A, 2) * np.linalg.norm(xhat) + np.linalg.norm(b)
    )
    print(f"n = {n:2d}, eta = {eta:.3e}")
```

Here the residual is computed exactly: the entries of $\widehat{x}$ are zeros and halves, so $A\widehat{x}=e_1$ involves only exactly represented sums. The large backward errors in the output therefore come from the solve, not from inaccurate residual evaluation.

**The mechanism is the same as in the floating-point examples:** a contribution that matters to the original equations is lost when added to a much larger intermediate quantity. Partial pivoting keeps every multiplier bounded by 1, yet it does not prevent this failure.

Partial pivoting usually gives modest growth and is widely used for dense linear systems. Large growth can nevertheless occur, including in matrices from applications; it is not restricted to specially constructed examples. See [Higham's discussion of the growth factor](https://nhigham.com/2020/07/14/what-is-the-growth-factor-for-gaussian-elimination/).

## What Pivoting Does and Does Not Guarantee

- **Existence:** Row pivoting guarantees a factorization $PA=LU$, including for singular matrices. A unique solve additionally requires nonsingularity.
- **Rounding error:** Partial pivoting keeps the multipliers small. The remaining concern is growth in $U$ and cancellation among the products that reconstruct $PA$.
- **Solution accuracy:** Small backward error leads to small forward error when the problem is sufficiently well-conditioned. Pivoting does not change that conditioning.

In practice, assess a computed solution using a scaled residual and, when forward accuracy matters, information about the condition number. These checks address the two separate sources of difficulty: the algorithm's rounding errors and the problem's sensitivity.
