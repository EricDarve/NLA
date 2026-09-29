# Existence and Uniqueness of LU Factorization

An invertible matrix need not have an LU factorization without row interchanges. A singular matrix may have one, and that factorization may or may not be unique.

Unless stated otherwise, **LU factorization** means $A=LU$, where $L$ is **unit lower triangular** and $U$ is upper triangular. Fixing the diagonal of $L$ is essential when discussing uniqueness. We work over $\mathbb R$ or $\mathbb C$.

## An Invertible Matrix with a Zero Pivot

Consider the matrix

$$
A=\begin{pmatrix}
1&6&1&0\\
0&1&9&3\\
1&6&1&1\\
0&0&1&9
\end{pmatrix}.
$$

The first pivot is $1$. Subtracting the first row from the third gives

$$
\begin{pmatrix}
1&6&1&0\\
0&1&9&3\\
0&0&0&1\\
0&0&1&9
\end{pmatrix}.
$$

The second pivot is also $1$, and the entries below it are already zero. At the third step, however, elimination would require dividing the entry $1$ in position $(4,3)$ by the zero pivot in position $(3,3)$.

The matrix itself is invertible: its determinant is $-1$. Swapping the third and fourth rows at this step allows elimination to continue. The obstruction concerns the chosen row order, not invertibility of the whole matrix.

## Leading Principal Minors and Pivots

Let

$$
A_k=A[1:k,1:k],\qquad \Delta_k=\det(A_k),
\qquad \Delta_0=1,
$$

where the index ranges include both endpoints. The matrix $A_k$ is the **leading principal submatrix** of order $k$, and $\Delta_k$ is its **leading principal minor**.

If $A=LU$, triangularity gives $A_k=L_kU_k$, where $L_k$ and $U_k$ are the corresponding leading blocks. Consequently,

$$
\Delta_k=\det(L_k)\det(U_k)=\prod_{j=1}^k u_{jj}.
$$

Whenever $\Delta_{k-1}\neq0$, we can divide two consecutive identities to obtain

$$
\boxed{u_{kk}=\frac{\Delta_k}{\Delta_{k-1}}.}
$$

This also describes the pivot reached after the first $k-1$ elimination steps have succeeded. It does not define a pivot when the denominator is zero.

For the opening example,

$$
\Delta_1=1,\qquad\Delta_2=1,\qquad\Delta_3=0,
\qquad\Delta_4=-1.
$$

Thus the third pivot is zero. In fact, no LU factorization in this row order can exist: $\Delta_3=0$ would force one of the first three diagonal entries of $U$ to vanish, contradicting $\det(A)=\prod_{j=1}^4u_{jj}\neq0$.

The determinant formula is useful after factorization. We do not compute determinants to choose pivots in the algorithm; the elimination updates already give them.

## When Is the Factorization Unique?

````{prf:theorem} Existence and Uniqueness of LU Factorization
:label: thm:existence_lu
An $n\times n$ matrix $A$ has a unique factorization $A=LU$, with $L$ unit lower triangular and $U$ upper triangular, if and only if

$$
\Delta_k\neq0,\qquad k=1,\ldots,n-1.
$$

No condition on $\Delta_n=\det(A)$ is needed for this statement.
````

````{prf:proof} Existence and Uniqueness.
**Sufficiency.** We use induction on $n$. For $n=1$, the unique factorization is $[a]=[1][a]$, including when $a=0$.

For $n>1$, partition

$$
A=\begin{pmatrix}B&b\\c&d\end{pmatrix},
$$

where $B=A_{n-1}$ and $c$ is a row vector. By induction, $B$ has a unique factorization $B=L_BU_B$. Also $\det(B)\neq0$, so both factors are invertible. Seek factors of the form

$$
L=\begin{pmatrix}L_B&0\\\ell&1\end{pmatrix},
\qquad
U=\begin{pmatrix}U_B&u\\0&\alpha\end{pmatrix}.
$$

Block multiplication requires

$$
L_Bu=b,\qquad \ell U_B=c,\qquad \alpha=d-\ell u.
$$

These equations uniquely determine $u$, $\ell$, and $\alpha$, giving a unique LU factorization of $A$. The last pivot $\alpha$ is allowed to be zero.

**Necessity.** Suppose an LU factorization exists but some $\Delta_k=0$ with $k<n$. The product formula for $\Delta_k$ implies that the leading $(n-1)\times(n-1)$ block of $U$ is singular. There is therefore a nonzero row vector $w$ such that

$$
wU_{n-1}=0.
$$

For any scalar $t$, define the unit lower triangular matrix

$$
M_t=\begin{pmatrix}I_{n-1}&0\\tw&1\end{pmatrix}.
$$

Then

$$
A=(LM_t)(M_t^{-1}U).
$$

The first factor remains unit lower triangular. The second is upper triangular: it changes only the last row of $U$, and its first $n-1$ entries remain zero because $wU_{n-1}=0$. Different values of $t$ give different first factors, since $L$ is invertible and $w\neq0$. Hence the factorization is not unique.
````

For a **nonsingular** matrix, any LU factorization is automatically unique: $\det(A)=\prod_j u_{jj}\neq0$ forces every pivot, and hence every leading principal minor, to be nonzero. The minor condition is unchanged; in this case it characterizes existence alone as well as existence and uniqueness.

The [no-pivoting algorithm](lu_decomposition.md) needs nonzero pivots only through step $n-1$. A nonzero final pivot is additionally required to solve arbitrary systems uniquely by backward substitution.

## What Changes for Singular Matrices?

A singular matrix can still have a unique LU factorization. For example,

$$
\begin{pmatrix}1&2\\3&6\end{pmatrix}
=
\begin{pmatrix}1&0\\3&1\end{pmatrix}
\begin{pmatrix}1&2\\0&0\end{pmatrix}.
$$

Here $\Delta_1=1$, so the theorem guarantees uniqueness even though the final pivot is zero.

If an earlier leading principal minor vanishes, LU may still exist, but it cannot be unique. To illustrate the construction in the proof, start with

$$
A=\begin{pmatrix}0&1\\0&1\end{pmatrix},
\qquad L=I,\qquad U=A.
$$

Here $U_{n-1}=[0]$, so choose $w=[1]$. The matrices from the proof are

$$
M_t=\begin{pmatrix}1&0\\t&1\end{pmatrix},
\qquad
M_t^{-1}=\begin{pmatrix}1&0\\-t&1\end{pmatrix}.
$$

Thus, for every scalar $t$,

$$
A=(LM_t)(M_t^{-1}U)
=
\begin{pmatrix}1&0\\t&1\end{pmatrix}
\begin{pmatrix}0&1\\0&1-t\end{pmatrix}.
$$

The standard algorithm would attempt a division by zero at its first step, yet these are valid factorizations. In contrast, the singular matrix

$$
\begin{pmatrix}0&0\\1&0\end{pmatrix}
$$

has no factorization with unit lower triangular $L$: $u_{11}=a_{11}=0$ would imply $a_{21}=l_{21}u_{11}=0$, a contradiction.

### Two Existence Criteria

For singular matrices, leading principal minors alone cannot decide whether an LU factorization exists. The criterion depends on whether $L$ is required to have unit diagonal. Both results below concern $A=LU$ in the original row and column ordering.

Write $\operatorname{null}(B)=\dim\{x:Bx=0\}$ for the dimension of the null space of a matrix $B$, and define the leading column and row blocks

$$
C_k=A[1:n,1:k],\qquad R_k=A[1:k,1:n].
$$

Thus $C_k$ and $R_k^T$ each have $k$ columns. Their nullities count linear dependencies among the first $k$ columns and the first $k$ rows of $A$, respectively.

````{prf:theorem} Existence of General LU
:label: thm:existence_lu_general
An $n\times n$ matrix $A$ has a factorization $A=LU$, with $L$ lower triangular and $U$ upper triangular, allowing zero diagonal entries in either factor, if and only if

$$
\operatorname{null}(A_k)
\le\operatorname{null}(C_k)+\operatorname{null}(R_k^T),
\qquad k=1,\ldots,n-1.
$$
````

Requiring $L$ to be unit lower triangular gives a stronger condition:

````{prf:theorem} Existence of LU with Unit Lower Triangular L
:label: thm:existence_lu_existence
A factorization $A=LU$ with $L$ unit lower triangular and $U$ upper triangular, allowing zero diagonal entries in $U$, exists if and only if

$$
\operatorname{null}(A_k)=\operatorname{null}(C_k),
\qquad k=1,\ldots,n-1.
$$
````

Since $C_kx=0$ implies $A_kx=0$, the second condition means that the two null spaces coincide: every linear dependency among the columns of $A_k$ remains valid when those columns are extended to all $n$ rows. In both criteria, $k=n$ adds only a trivial identity or inequality.

The distinction matters for the singular matrix considered above:

$$
\begin{pmatrix}0&0\\1&0\end{pmatrix}
=\underbrace{\begin{pmatrix}0&0\\1&0\end{pmatrix}}_{L}
\underbrace{\begin{pmatrix}1&0\\0&1\end{pmatrix}}_{U}.
$$

It has a general LU factorization, but no factorization with unit lower triangular $L$.

We state these results without proof; see [Property 2.1 and Theorem 2.5 for general LU](https://arxiv.org/pdf/2601.07791#page=5) and [Theorem 3.1 for unit lower triangular $L$](https://arxiv.org/pdf/2601.07791#page=10).

These existence criteria do not justify continuing the standard elimination algorithm through a zero pivot. A zero pivot cannot eliminate a nonzero entry below it; even when those entries are already zero, the choice of multipliers can affect later steps. The cited construction for unit lower triangular $L$ uses internal column permutations and proves that the resulting factors retain the required triangular structure in the original ordering.

## Existence Does Not Guarantee Numerical Accuracy

The results above concern exact arithmetic. Even when LU exists uniquely, the product of the computed factors need not approximate $A$ accurately. Consider

$$
A_\epsilon=\begin{pmatrix}\epsilon&1\\1&\pi\end{pmatrix},
\qquad \epsilon>0.
$$

Its exact factors are

$$
L=\begin{pmatrix}1&0\\1/\epsilon&1\end{pmatrix},
\qquad
U=\begin{pmatrix}\epsilon&1\\0&\pi-1/\epsilon\end{pmatrix}.
$$

When $1/\epsilon$ is much larger than $\pi$, rounding the subtraction can lose the contribution of $\pi$. For example, with $\epsilon=10^{-18}$, double-precision arithmetic rounds $\pi-1/\epsilon$ to the same value as $-1/\epsilon$. The $(2,2)$ entry reconstructed from these rounded factor entries is then zero instead of $\pi$.

The issue is the loss of a contribution that matters to reconstructing $A_\epsilon$. A small pivot or a large multiplier alone does not prove that a result is inaccurate, and the effect on a computed solution also depends on the right-hand side.

Swapping the two rows puts $1$ in the first pivot position and avoids the large multiplier in this example. The next sections explain [floating-point arithmetic](floating_point.md) and [LU with row pivoting](lu_pivoting.md). These separate the question of whether a factorization exists from whether we can compute it accurately.
