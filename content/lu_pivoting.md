# LU Factorization with Row Pivoting

LU factorization can fail at a zero pivot even when the matrix is invertible. A nonzero pivot can also cause trouble if it is much smaller than the entries below it: elimination then forms large intermediate quantities whose cancellation may expose rounding errors.

**Partial pivoting chooses the largest available entry in the current column as the pivot.** This guarantees that the elimination multipliers have magnitude at most 1 in exact arithmetic. It also guarantees that a row-permuted LU factorization exists for every square matrix.

We first describe the algorithm, then explain how to assess its accuracy. Two questions must be kept separate: how much error the algorithm introduces, and how sensitive the original linear system is to changes in its data.

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

All entries remain moderate. The row swap avoids the large intermediate quantities and the subsequent cancellation in this example. A small pivot or a large multiplier alone does not prove that a computation is inaccurate; we must examine how the resulting errors enter the answer.

## Partial Pivoting and the Factorization $PA=LU$

At step $k$, let $b_{ij}$ denote the entries of the current matrix after the preceding eliminations and row interchanges. The pivot is chosen from the **current column**, not from the original matrix $A$.

1. Choose a row $p\in\{k,\ldots,n\}$ such that

   $$
   |b_{pk}|=\max_{i=k,\ldots,n}|b_{ik}|.
   $$

2. Interchange rows $k$ and $p$.
3. If the pivot is nonzero, compute the multipliers and update the trailing block:

   $$
   \begin{aligned}
   l_{ik}&=b_{ik}/b_{kk}, &&i>k,\\
   b_{ij}&\leftarrow b_{ij}-l_{ik}b_{kj}, &&i,j>k.
   \end{aligned}
   $$

Because the pivot has largest magnitude in its active column,

$$
|l_{ik}|\le1.
$$

If the selected pivot is zero, every entry below it in that column is also zero. No elimination is needed there: set those multipliers to zero and continue. This case matters for singular matrices. For a nonsingular matrix, partial pivoting encounters no zero pivots in exact arithmetic.

The product of the row interchanges is a permutation matrix $P$, and the resulting factorization is

$$
\boxed{PA=LU.}
$$

Here $L$ is unit lower triangular and $U$ is upper triangular. To solve $Ax=b$, apply the same permutation to the right-hand side and solve

$$
Ly=Pb,\qquad Ux=y.
$$

The original unknowns have not been reordered, so there is no permutation to undo in $x$. Some software instead writes $A=PLU$; its permutation is the transpose of the one used here.

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

As in the [unpivoted algorithm](lu_decomposition.md), store $L$ below the diagonal and $U$ on and above it, with the unit diagonal of $L$ implicit. The following function modifies a square NumPy array in place and returns $P$.

```python
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

For example, using the triangular-solve functions from [LU decomposition](lu_decomposition.md):

```python
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
```

This example swaps rows at both elimination steps, so the second swap must move the multipliers from the first step. Pivot searches and row swaps add $O(n^2)$ work; the leading arithmetic cost remains $\frac23n^3$ flops for a dense real factorization. Implementations generally store the permutations as row indices rather than a full matrix $P$.

## Forward Error, Backward Error, and the Residual

We now turn to accuracy. Assume $A$ is nonsingular and $b\ne0$. Let $x$ be the exact solution of $Ax=b$, and let $\widehat{x}$ be the computed solution. The hat marks a computed quantity that may contain rounding error. Unless a subscript specifies otherwise, we use the Euclidean vector norm and its induced matrix norm, $\|\cdot\|_2$.

### Forward Error

The **forward error** measures the difference between the computed and exact answers:

$$
\text{absolute forward error}=\|\widehat{x}-x\|,
\qquad
\text{relative forward error}=\frac{\|\widehat{x}-x\|}{\|x\|}.
$$

We usually do not know $x$, so we cannot measure this error directly. We can still bound it without knowing $x$, as the perturbation analysis below will show.

### Backward Error

The **backward error** asks how much the input data must change to make $\widehat{x}$ an exact solution. We seek small perturbations $E$ and $e$ such that

$$
(A+E)\widehat{x}=b+e.
$$

The changes are measured relative to the data, through $\|E\|/\|A\|$ and $\|e\|/\|b\|$. A method is **backward stable** if it guarantees perturbations of order $u$, allowing for a modest factor depending on the dimension. Here $u$ is the unit roundoff. This is a guarantee about the algorithm over its stated class of inputs; a small error on one example does not establish it.

The **residual** provides a way to assess a particular computed solution:

$$
r=b-A\widehat{x}.
$$

Since $A\widehat{x}=b-r$, the computed vector exactly solves a problem with the same matrix and right-hand side $b-r$. Thus $\|r\|/\|b\|$ measures a relative backward error when only $b$ is allowed to change.

When both $A$ and $b$ may change, the smallest common relative bound on their perturbations is

$$
\eta(\widehat{x})
=\frac{\|b-A\widehat{x}\|}
{\|A\|\,\|\widehat{x}\|+\|b\|}.
$$

More precisely, this is the smallest $\eta$ for which some $E,e$ satisfy the perturbed system with $\|E\|\le\eta\|A\|$ and $\|e\|\le\eta\|b\|$. This [normwise backward-error formula](https://eprints.maths.manchester.ac.uk/2562/1/paper.pdf#page=6) gives a scale-aware residual check without requiring the exact solution.

A dense residual costs $O(n^2)$ operations. It is itself computed in floating point, so a very small residual may need higher precision to be measured reliably. A small residual indicates a nearby solved problem; whether its solution is close to $x$ depends on the sensitivity of the original system.

## Conditioning: How Changes in the Data Affect the Solution

**Conditioning is a property of the problem.** It describes how much the exact solution can change when the input data change. Stability describes the errors introduced by an algorithm used to solve that problem.

### Perturbing the Right-Hand Side

If $A(x+\delta x)=b+\delta b$, subtracting $Ax=b$ gives

$$
A\delta x=\delta b,
\qquad
\delta x=A^{-1}\delta b.
$$

Therefore $\|\delta x\|\le\|A^{-1}\|\,\|\delta b\|$. Using $\|b\|\le\|A\|\,\|x\|$ yields

$$
\frac{\|\delta x\|}{\|x\|}
\le \kappa(A)\frac{\|\delta b\|}{\|b\|},
\qquad
\boxed{\kappa(A)=\|A\|\,\|A^{-1}\|.}
$$

The **condition number** satisfies $\kappa(A)\ge1$, since $1=\|I\|\le\|A\|\,\|A^{-1}\|$. In the 2-norm, the SVD gives

$$
\kappa_2(A)=\frac{\sigma_{\max}(A)}{\sigma_{\min}(A)}.
$$

This is the ratio of the largest to the smallest stretching of a unit vector by $A$. A large condition number means that some small relative changes in the data can produce large relative changes in the solution. It is a worst-case bound; not every perturbation is amplified by this factor.

For example, take $A=\operatorname{diag}(1,\epsilon)$ with $0<\epsilon\ll1$, and $b=(1,0)^T$. The solution is $x=(1,0)^T$. Changing $b$ by $(0,\epsilon)^T$ changes the solution by $(0,1)^T$: a relative input change of $\epsilon$ causes a relative solution change of 1. Here $\kappa_2(A)=1/\epsilon$.

Row permutations preserve singular values, so $\kappa_2(PA)=\kappa_2(A)$. Pivoting changes the elimination process; it does not remove the sensitivity of the original problem.

### Perturbing Both $A$ and $b$

To treat changes in $A$ rigorously, we first establish when a small perturbation preserves invertibility.

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

When $\kappa(A)\|E\|/\|A\|\ll1$, the denominator is close to 1. The practical interpretation is

$$
\begin{aligned}
\text{relative forward error}
&\lesssim\kappa(A)\times{}\\
&\qquad\text{sum of relative data perturbations}.
\end{aligned}
$$

Even if the backward error is of order $u$, this simplification needs $\kappa(A)u$ to be small, with the dimension-dependent factors included. Backward stability alone does not guarantee small forward error for an ill-conditioned system.

As a useful special case, take $E=0$ and $e=-r$. Then

$$
\frac{\|\widehat{x}-x\|}{\|x\|}
\le\kappa(A)\frac{\|r\|}{\|b\|}.
$$

This bounds the forward error using the residual and the condition number, without knowing the exact solution.

## Backward Error of LU and the Role of Cancellation

Let $\widehat{L}$ and $\widehat{U}$ be the computed factors; their hats distinguish them from exact factors. The following standard result separates errors in the factorization from those in the complete solve. We state the scalar-operation bounds for real arithmetic; complex arithmetic has analogous bounds with different constants.

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

Absolute values and inequalities are entrywise. The product $|\widehat{L}|\,|\widehat{U}|$ is an ordinary matrix product of nonnegative matrices.
````

The full rounding-error derivation is given in [Higham, *Accuracy and Stability of Numerical Algorithms*, Chapter 9](https://epubs.siam.org/doi/10.1137/1.9780898718027.ch9); the solve bound is also stated in [Carson and Higham, equation (7.1)](https://eprints.maths.manchester.ac.uk/2562/1/paper.pdf#page=13). The factor $P$ in $|PE|$ puts the perturbation in the same row ordering as the factors.

To interpret the bound, compare

$$
(\widehat{L}\widehat{U})_{ij}
=\sum_k\widehat{l}_{ik}\widehat{u}_{kj}
\quad\text{with}\quad
(|\widehat{L}|\,|\widehat{U}|)_{ij}
=\sum_k|\widehat{l}_{ik}\widehat{u}_{kj}|.
$$

The first sum allows cancellation; the second measures the sizes of all terms before cancellation. This is the same distinction as in the [summation error bound](floating_point.md), {prf:ref}`thm:summation_error`. Large terms that cancel to reproduce a much smaller matrix entry can carry rounding errors that are large relative to that entry.

For example, in the unpivoted $A_\epsilon$ factorization, the two terms reconstructing $a_{22}$ have total magnitude about $2/\epsilon$, although $a_{22}=\pi$. After swapping rows, the corresponding products stay moderate. The bounds explain the advantage of pivoting without asserting that every large multiplier must cause a large error.

Since row permutations preserve the 2-norm, the solve bound implies

$$
\frac{\|E\|}{\|A\|}
\le\gamma_{3n}
\frac{\bigl\||\widehat{L}|\,|\widehat{U}|\bigr\|}{\|A\|}.
$$

Thus moderate size of $|\widehat{L}|\,|\widehat{U}|$ relative to $A$ gives a small backward error. This is a normwise statement; it does not require every perturbation entry to be small relative to the corresponding entry of $A$, which may be zero. Conversely, a large upper bound signals possible loss of accuracy, not proof that a large error occurred.

## Element Growth: The Limitation of Partial Pivoting

Partial pivoting bounds the multipliers in $L$, but entries in $U$ can still grow. To measure this, let $B^{(k)}$ be the active trailing matrix at the start of step $k$, before elimination in that column. For $A\ne0$, define the **growth factor**

$$
\rho_n
=\frac{\displaystyle\max_{k}\max_{i,j}|b_{ij}^{(k)}|}
{\displaystyle\max_{i,j}|a_{ij}|}.
$$

This definition includes all intermediate trailing matrices, beginning with $A$. It excludes the stored multipliers in a packed implementation. In particular, it bounds the size of every entry eventually placed in $U$.

In exact arithmetic, an update satisfies

$$
|b_{ij}-l_{ik}b_{kj}|
\le |b_{ij}|+|b_{kj}|
$$

because $|l_{ik}|\le1$. Each elimination step can therefore at most double the largest active entry. There are $n-1$ updates, so

$$
\boxed{\rho_n\le2^{n-1}.}
$$

### A Matrix That Attains the Bound

Consider

$$
A=\begin{pmatrix}
1&0&0&1\\
-1&1&0&1\\
-1&-1&1&1\\
-1&-1&-1&1
\end{pmatrix}.
$$

At each step the pivot candidates have equal magnitude. If ties are resolved by choosing the first candidate, as in the implementation above, no row swaps occur. The factors are

$$
L=\begin{pmatrix}
1&0&0&0\\
-1&1&0&0\\
-1&-1&1&0\\
-1&-1&-1&1
\end{pmatrix},
\qquad
U=\begin{pmatrix}
1&0&0&1\\
0&1&0&2\\
0&0&1&4\\
0&0&0&8
\end{pmatrix}.
$$

The last-column entries double at each step, giving $\rho_4=8$. The same construction in dimension $n$ gives $\rho_n=2^{n-1}$.

This example demonstrates growth, not necessarily rounding error: its displayed factorization is computed exactly in binary arithmetic. It shows that bounding the multipliers does not, by itself, ensure a small backward-error bound. Nearby data or other right-hand sides can expose the risks of such growth.

Partial pivoting usually gives modest growth and is widely used for dense linear systems. Large growth can nevertheless occur, including in matrices from applications; it is not restricted to specially constructed examples. See [Higham's discussion of the growth factor](https://nhigham.com/2020/07/14/what-is-the-growth-factor-for-gaussian-elimination/).

## What Pivoting Does and Does Not Guarantee

- **Existence:** Row pivoting guarantees a factorization $PA=LU$, including for singular matrices. A unique solve additionally requires nonsingularity.
- **Rounding error:** Partial pivoting keeps the multipliers small. The remaining concern is growth in $U$ and cancellation among the products that reconstruct $PA$.
- **Solution accuracy:** Small backward error leads to small forward error when the problem is sufficiently well-conditioned. Pivoting does not change that conditioning.

In practice, assess a computed solution using a scaled residual and, when forward accuracy matters, information about the condition number. These checks address the two separate sources of difficulty: the algorithm's rounding errors and the problem's sensitivity.
