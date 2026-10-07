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

**Forward error is what we ultimately want to control:** it tells us how accurate our answer is.

### Backward Error

Directly estimating forward error is difficult: the exact solution $x$ is usually unknown, and tracking how rounding errors propagate to the final answer can become complicated. **Backward error is much easier to compute or bound.** It can be determined from the known inputs $A$ and $b$ and the computed answer $\widehat{x}$, by measuring how well that answer satisfies the equations.

This is the reason for taking a detour through backward error. First, bound the input changes needed to explain the computed answer. Then use the problem's sensitivity to bound the resulting change in the solution. Together, these two steps give a forward-error bound.

**1. Interpret the computed answer as an exact solution.** The computed vector $\widehat{x}$ generally does not satisfy $A\widehat{x}=b$ exactly. Backward error starts by keeping $\widehat{x}$ fixed and considering changes to the input data, $A$ and $b$, that would make it an exact solution.

Write the changed inputs as $A+E$ and $b+e$. The matrix $E$ changes the coefficients of the system, and the vector $e$ changes the right-hand side. They must satisfy

$$
(A+E)\widehat{x}=b+e.
$$

This equation is exact: $\widehat{x}$ solves the system with the changed inputs. It is a way to interpret the answer already computed; it does not require running the algorithm again.

Such changes always exist. For example, keeping $A$ unchanged and replacing $b$ by $A\widehat{x}$ makes the equation hold. What matters is **how much the inputs need to change**. Small changes mean that the computed answer solves a problem close to the one we intended to solve.

**2. Measure the size of a particular change.** For a chosen pair $(E,e)$, the relative changes in the two inputs are $\|E\|/\|A\|$ and $\|e\|/\|b\|$. Dividing by the input norms makes each change relative to the size of the data it modifies. To describe the pair by one number, take the larger relative change:

$$
\tau=\max\left\{\frac{\|E\|}{\|A\|},
\frac{\|e\|}{\|b\|}\right\}.
$$

For example, if the two relative changes are $10^{-8}$ and $3\times10^{-8}$, then $\tau=3\times10^{-8}$. Both inputs have changed by at most this fraction of their original size.

**3. Find the smallest change that works.** Many pairs $(E,e)$ can make $\widehat{x}$ exact, and they can have different sizes. The **relative backward error** $\eta(\widehat{x})$ is the smallest value of $\tau$ among all pairs that satisfy the perturbed equation:

$$
\boxed{
\eta(\widehat{x})=
\min_{\substack{E,e}} \;
\max\left\{\frac{\|E\|}{\|A\|},
\frac{\|e\|}{\|b\|}\right\}, \; \text{with} \; (A+E)\widehat{x}=b+e \,.
}
$$

The minimum is taken over the pair: changes to $A$ and $b$ can work together to make $\widehat{x}$ exact. Thus $\eta(\widehat{x})=10^{-8}$ means that relative changes of at most $10^{-8}$ in both inputs suffice, and no smaller common bound is possible. A different pair may have much larger changes; $\eta$ measures the best possible pair.

**4. Compute the minimum from the residual.** There is no need to search over possible perturbations or to know the true solution $x$. Form the **residual**

$$
r=b-A\widehat{x}.
$$

It measures how far $\widehat{x}$ is from satisfying the original equations. The [relative backward error](https://eprints.maths.manchester.ac.uk/2562/1/paper.pdf#page=6) is exactly

$$
\boxed{\eta(\widehat{x})=
\frac{\|r\|}{\|A\|\,\|\widehat{x}\|+\|b\|}.}
$$

The denominator accounts for both ways to correct the residual: changing $A$ affects the product $A\widehat{x}$, while changing $b$ affects the right-hand side directly. The formula gives the smallest common relative size of these changes.

The proof has two parts: every pair that makes $\widehat{x}$ exact must be at least this large, and a specific pair achieves this size.

````{prf:proof} Residual Formula for Backward Error.
Let $d=\|A\|\,\|\widehat{x}\|+\|b\|$, which is positive because $b\ne0$.

**Lower bound.** For any pair $(E,e)$ satisfying $(A+E)\widehat{x}=b+e$, set

$$
\tau=\max\left\{\frac{\|E\|}{\|A\|},
\frac{\|e\|}{\|b\|}\right\}.
$$

The perturbed equation gives $r=E\widehat{x}-e$. Hence

$$
\begin{aligned}
\|r\|&\le\|E\|\,\|\widehat{x}\|+\|e\|\\
&\le\tau\bigl(\|A\|\,\|\widehat{x}\|+\|b\|\bigr)
=\tau d.
\end{aligned}
$$

Every such pair therefore has $\tau\ge\|r\|/d$, so $\eta(\widehat{x})\ge\|r\|/d$.

**Attaining the bound.** If $\widehat{x}\ne0$, choose

$$
E=\frac{\|A\|}{d}\,\frac{r\widehat{x}^{H}}{\|\widehat{x}\|},
\qquad
e=-\frac{\|b\|}{d}\,r.
$$

Here $H$ denotes conjugate transpose, or ordinary transpose for real vectors. These choices satisfy

$$
E\widehat{x}-e
=\frac{\|A\|\,\|\widehat{x}\|+\|b\|}{d}\,r=r,
$$

so $(A+E)\widehat{x}=b+e$. In the Euclidean norm, the outer product satisfies $\|r\widehat{x}^H\|_2=\|r\|_2\,\|\widehat{x}\|_2$. Consequently,

$$
\frac{\|E\|}{\|A\|}
=\frac{\|e\|}{\|b\|}
=\frac{\|r\|}{d}.
$$

This pair attains the lower bound, proving the formula. If $\widehat{x}=0$, then $r=b$ and every admissible pair must have $e=-b$. Taking $E=0$ gives $\eta(0)=1=\|r\|/d$, so the formula also holds in this case.
````

The unit roundoff $u$ sets the scale of relative error in a single rounded arithmetic operation. It provides a natural benchmark for judging the backward error of an algorithm.

```{admonition} Key definition: Backward stability
:class: important

**A method is backward stable if it guarantees $\eta(\widehat{x})$ of order the unit roundoff $u$**, allowing for a modest factor depending on the dimension.
```

Backward error describes one computed answer; backward stability requires a guarantee for the method across its intended inputs. Computing the residual costs $O(n^2)$ operations for a dense system, so it provides a practical check on an individual answer. A small backward error tells us that the computed answer solves a nearby problem. To decide whether it is close to the answer we wanted, we must ask how much the solution changes when the data change.

## Sensitivity: Relating Backward Error to Forward Error

For a general problem $x=f(d)$, let $d$ denote the input data and $f(d)$ the exact answer. **Sensitivity measures how much the answer changes when the input changes.** The **local absolute sensitivity** is

$$
S_f(d)=\lim_{\varepsilon\to0^+}
\sup_{0<\|\delta d\|\le\varepsilon}
\frac{\|f(d+\delta d)-f(d)\|}{\|\delta d\|}.
$$

The ratio measures the change in the answer per unit change in the input. The supremum selects the perturbation with the largest amplification, and the limit restricts attention to small changes near $d$.

If $f$ is differentiable at $d$, its derivative, denoted $Df(d)$, is the linear map that predicts the change in the output caused by a small input change $\delta d$:

$$
f(d+\delta d)-f(d)\approx Df(d)\,\delta d.
$$

For vector inputs and outputs, $Df(d)$ is the **Jacobian matrix**: its $(i,j)$ entry is the partial derivative $\partial f_i/\partial d_j$, evaluated at $d$. Its induced norm measures the largest amplification of an input change by this linear map. Consequently,

$$
S_f(d)=\|Df(d)\|.
$$

For a scalar function of one variable, $Df(d)$ is the ordinary derivative $f'(d)$, so $S_f(d)=|f'(d)|$.

Suppose backward error analysis shows that $\widehat{x}=f(d+\delta d)$. For differentiable $f$, the definition gives the local bound

$$
\underbrace{\|\widehat{x}-x\|}_{\text{forward error}}
\le S_f(d) 
\hspace{-1.5em} \underbrace{\|\delta d\|}_{\text{backward perturbation}} \hspace{-1.5em} 
+ o(\|\delta d\|).
$$

The remainder $o(\|\delta d\|)$ becomes negligible compared with $\|\delta d\|$ as the perturbation tends to zero. Thus, to first order,

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

For $Ax=b$, the input data are $A$ and $b$, and the output is $x=A^{-1}b$. The matrix condition number measures the sensitivity relevant to this problem. Combining it with backward error gives a simple first-order estimate of solution error, followed by a rigorous bound that removes the approximation.

### The Matrix Condition Number

With $A$ fixed, changing $b$ by $\delta b$ changes the solution by $A^{-1}\delta b$. Thus $\|A^{-1}\|$ bounds the amplification of absolute changes. Converting to relative changes uses $\|b\|/\|x\|\le\|A\|$. This leads to the **matrix condition number**

$$
\boxed{\kappa(A)=\|A\|\,\|A^{-1}\|.}
$$

The **condition number** satisfies $\kappa(A)\ge1$, since $1=\|I\|\le\|A\|\,\|A^{-1}\|$. In the 2-norm, the SVD gives

$$
\kappa_2(A)=\frac{\sigma_{\max}(A)}{\sigma_{\min}(A)}.
$$

Geometrically, $\kappa_2(A)$ is the ratio of the largest to the smallest stretching of a unit vector by $A$. A large value means that small relative changes in the data can produce large relative changes in the solution.

For example, take $A=\operatorname{diag}(1,\epsilon)$ with $0<\epsilon\ll1$, and $b=(1,0)^T$. The solution is $x=(1,0)^T$. Changing $b$ by $(0,\epsilon)^T$ changes the solution by $(0,1)^T$: a relative input change of $\epsilon$ causes a relative solution change of 1. Here $\kappa_2(A)=1/\epsilon$.

Row permutations preserve singular values, so $\kappa_2(PA)=\kappa_2(A)$. Pivoting changes the elimination process; it does not remove the sensitivity of the original problem.

### A First-Order Forward-Error Estimate

The goal is to bound the relative forward error using the backward error $\eta(\widehat{x})$ and the condition number $\kappa(A)$. The connection is the residual: it measures how much the computed answer fails to satisfy the equations, and $A^{-1}$ converts that discrepancy into an error in the solution.

**1. Relate the solution error to the residual.** Let $x$ be the exact solution, let $\widehat{x}$ be the computed solution, and write their difference as $\delta x=\widehat{x}-x$. Recall that $Ax=b$ and $r=b-A\widehat{x}$. Multiplying the solution error by $A$ gives

$$
\begin{aligned}
A\delta x
&=A(\widehat{x}-x)\\
&=A\widehat{x}-Ax\\
&=A\widehat{x}-b=-r.
\end{aligned}
$$

Since $A$ is nonsingular, multiplying by $A^{-1}$ yields

$$
\delta x=-A^{-1}r.
$$

Thus the residual and the solution error are related, but they are generally different vectors with different magnitudes.

**2. Bound the relative solution error.** By the defining property of an induced matrix norm, $\|A^{-1}r\|\le\|A^{-1}\|\,\|r\|$. Taking norms in the identity above therefore gives

$$
\|\delta x\|=\|A^{-1}r\|
\le\|A^{-1}\|\,\|r\|.
$$

To measure the error relative to the exact answer, divide by $\|x\|$. This is valid because $b\ne0$ implies $x\ne0$:

$$
\frac{\|\delta x\|}{\|x\|}
\le\frac{\|A^{-1}\|\,\|r\|}{\|x\|}.
$$

**3. Express the residual in terms of backward error.** The formula proved earlier is

$$
\eta(\widehat{x})=
\frac{\|r\|}{\|A\|\,\|\widehat{x}\|+\|b\|}.
$$

Multiplying by its denominator gives

$$
\|r\|=\eta(\widehat{x})
\bigl(\|A\|\,\|\widehat{x}\|+\|b\|\bigr).
$$

Substitute this expression for $\|r\|$ into the relative-error bound:

$$
\frac{\|\delta x\|}{\|x\|}
\le\frac{\|A^{-1}\|\,\eta(\widehat{x})}{\|x\|}
\bigl(\|A\|\,\|\widehat{x}\|+\|b\|\bigr).
$$

**4. Bring in the condition number.** To combine the two terms in parentheses, use the exact equation $b=Ax$ to bound

$$
\|b\|=\|Ax\|\le\|A\|\,\|x\|.
$$

Replacing $\|b\|$ by this upper bound and factoring out $\|A\|$ gives

$$
\begin{aligned}
\frac{\|\delta x\|}{\|x\|}
&\le\frac{\|A^{-1}\|\,\eta(\widehat{x})}{\|x\|}
\bigl(\|A\|\,\|\widehat{x}\|+\|A\|\,\|x\|\bigr)\\
&=\|A^{-1}\|\,\|A\|\,\eta(\widehat{x})
\frac{\|\widehat{x}\|+\|x\|}{\|x\|}\\
&=\kappa(A)\eta(\widehat{x})
\left(1+\frac{\|\widehat{x}\|}{\|x\|}\right).
\end{aligned}
$$

The last line uses $\kappa(A)=\|A\|\,\|A^{-1}\|$. Every step so far is rigorous. The remaining obstacle is the unknown $\|x\|$ on the right-hand side.

**5. Make the first-order approximation.** If the computed answer is close to the exact answer, their norms are close as well. Assuming $\|\widehat{x}\|\approx\|x\|$, we have

$$
\frac{\|\widehat{x}\|}{\|x\|}\approx1,
\qquad
1+\frac{\|\widehat{x}\|}{\|x\|}\approx2.
$$

This is the only approximation in the derivation. It removes the unknown solution norm from the right-hand side.

```{admonition} Key result: First-order forward-error estimate
:class: important

Under the assumption $\|\widehat{x}\|\approx\|x\|$,

$$
\frac{\|\widehat{x}-x\|}{\|x\|}
\lesssim 2\kappa(A)\eta(\widehat{x}).
$$

This is an approximate estimate. It identifies **$\kappa(A)\eta(\widehat{x})$** as the quantity that controls forward error: the problem's sensitivity multiplies the algorithm's backward error. The next result removes the assumption about the solution norms and gives a rigorous guarantee.
```

### A Rigorous Forward-Error Bound

The first-order estimate assumes $\|\widehat{x}\|\approx\|x\|$. A rigorous bound removes this assumption and gives a condition that can be checked using $\kappa(A)$ and $\eta(\widehat{x})$:

```{admonition} Key result: From backward error to forward error
:class: important

**If $\kappa(A)\eta(\widehat{x})<1$, then**

$$
\frac{\|\widehat{x}-x\|}{\|x\|}
\le\frac{2\kappa(A)\eta(\widehat{x})}
{1-\kappa(A)\eta(\widehat{x})}.
$$

This is a rigorous bound under the stated condition. Its numerator is the same as in the first-order estimate; the denominator $1-\kappa(A)\eta(\widehat{x})$ accounts for the effect of finite perturbations. When $\kappa(A)\eta(\widehat{x})\ll1$, that denominator is close to 1, and the bound reduces to the earlier estimate to first order. This justifies using the first-order estimate when that product is small.
```

The factor 2 comes from allowing relative changes of size $\eta$ in each of $A$ and $b$.

To prove the bound, consider the perturbed system $(A+E)\widehat{x}=b+e$. Changes in $A$ also change its inverse, so we first need to control $\|(A+E)^{-1}\|$. The following lemma provides this control and explains where the denominator comes from.

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

Applying the lemma with $X=A^{-1}E$ gives a bound in terms of the perturbations $E$ and $e$. Choosing a pair that minimizes backward error will then give the result in terms of $\eta(\widehat{x})$.

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

For a perturbation pair attaining the minimum in the definition of backward error, both
$\|E\|/\|A\|\le\eta(\widehat{x})$ and $\|e\|/\|b\|\le\eta(\widehat{x})$. When $\kappa(A)\eta(\widehat{x})<1$, the theorem applies to this pair. Substituting these inequalities gives the stated bound:

$$
\frac{\|\widehat{x}-x\|}{\|x\|}
\le\frac{2\kappa(A)\eta(\widehat{x})}
{1-\kappa(A)\eta(\widehat{x})}.
$$

## Backward Error of LU and the Role of Cancellation

The residual formula measures the backward error of an answer already computed. To understand the reliability of LU, we also need a bound on the backward error its rounding errors can produce. Combined with the preceding forward-error bound, this will connect the factorization to solution accuracy.

Let $\widehat{L}$ and $\widehat{U}$ be the computed factors; their hats distinguish them from exact factors.

The following theorem gives separate bounds for the factorization and the complete solve. We state the bounds for real arithmetic; complex arithmetic has analogous bounds with different constants.

````{prf:theorem} Backward Error of LU Factorization and Solution
:label: thm:backward_error_lu
Assume real arithmetic with rounding to nearest, no overflow or underflow, and $3nu<1$.

**Factorization.** If LU with partial pivoting completes, its computed factors satisfy

$$
PA+F=\widehat{L}\widehat{U},
\qquad
|F|\le\gamma_n|\widehat{L}|\,|\widehat{U}|,
$$

where $\gamma_m=mu/(1-mu)$.

**Complete solve.** If the computed $\widehat{U}$ has no zero diagonal entries, the solution $\widehat{x}$ obtained by forward and backward substitution satisfies

$$
(A+E)\widehat{x}=b,
\qquad
|PE|\le\gamma_{3n}|\widehat{L}|\,|\widehat{U}|.
$$

Here $F$ accounts for errors in the factorization, while $E$ accounts for errors in both the factorization and the two triangular solves.

Absolute values and inequalities are entrywise. The product $|\widehat{L}|\,|\widehat{U}|$ is an ordinary matrix product of nonnegative matrices. The factor $P$ in $|PE|$ puts the perturbation in the same row ordering as the factors.
````

The theorem supplies one admissible pair $(E,0)$. For this pair, the maximum of the two relative changes is $\|E\|/\|A\|$. Since $\eta(\widehat{x})$ is the minimum over all admissible pairs, $\eta(\widehat{x})\le\|E\|/\|A\|$. The theorem's pair need not attain the minimum. Taking norms in the theorem and using $\|PE\|=\|E\|$ yields

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

For an overview of growth factors and how they depend on the pivoting strategy, see Nick Higham's note [*What Is the Growth Factor for Gaussian Elimination?*](https://nhigham.com/2020/07/14/what-is-the-growth-factor-for-gaussian-elimination/).

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

```{warning}
Reading this section is **optional. This is not covered on the exams or graded assignments.**
```

Large element growth alone does not establish that a computed solution has a large error. To show that **actual backward error** can be large, choose a right-hand side and follow the rounding in the solve. Assume binary arithmetic with rounding to nearest, unit roundoff $u$, and no overflow or underflow. Choose

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
    print(f"n = {n:2d}, eta = {eta:.1e}")
```

Here the residual is computed exactly: the entries of $\widehat{x}$ are zeros and halves, so $A\widehat{x}=e_1$ involves only exactly represented sums. The large backward errors in the output therefore come from the solve, not from inaccurate residual evaluation.

**The mechanism is the same as in the floating-point examples:** a contribution that matters to the original equations is lost when added to a much larger intermediate quantity. Partial pivoting keeps every multiplier bounded by 1, yet it does not prevent this failure.

Partial pivoting usually gives modest growth and is widely used for dense linear systems. Large growth can nevertheless occur, including in matrices from applications; it is not restricted to specially constructed examples. See [Higham's discussion of the growth factor](https://nhigham.com/2020/07/14/what-is-the-growth-factor-for-gaussian-elimination/).

## What Pivoting Does and Does Not Guarantee

- **Existence:** Row pivoting guarantees a factorization $PA=LU$, including for singular matrices. A unique solve additionally requires nonsingularity.
- **Rounding error:** Partial pivoting keeps the multipliers small. The remaining concern is growth in $U$ and cancellation among the products that reconstruct $PA$.
- **Solution accuracy:** Small backward error leads to small forward error when the problem is sufficiently well-conditioned. Pivoting does not change that conditioning.

In practice, assess a computed solution using a scaled residual and, when forward accuracy matters, information about the condition number. These checks address the two separate sources of difficulty: the algorithm's rounding errors and the problem's sensitivity.
