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

# Modified Gram-Schmidt

[Householder reflections](householder_reflections.md) and [Givens rotations](givens_rotations.md) compute QR by introducing zeros in a matrix. Gram-Schmidt takes a different approach: it builds an orthonormal basis for the column space, one vector at a time. The coefficients used to express the original columns in this basis form the triangular factor $R$.

This approach is useful when vectors arrive sequentially. For example, an iterative method may generate its next vector using the basis vectors already computed. Gram-Schmidt can extend that basis without changing the earlier vectors. We will use this construction later in [Arnoldi iteration](arnoldi_iteration.md). Householder QR can also be extended, but it stores the basis implicitly as a product of transformations.

Throughout this page, $A\in\mathbb{R}^{m\times n}$ with $m\ge n$ has linearly independent columns. We seek the **reduced QR factorization**

$$
A=QR,\qquad Q\in\mathbb{R}^{m\times n},\qquad
Q^TQ=I_n,
$$

where $R\in\mathbb{R}^{n\times n}$ is upper triangular with positive diagonal entries. Since $R$ is invertible, $Q=AR^{-1}$. This viewpoint is sometimes called **triangular orthogonalization**, in contrast to applying orthogonal transformations to triangularize $A$. The algorithm does not form $R^{-1}$.

## Removing One Direction at a Time

Let $a_j$ denote the original $j$-th column of $A$. The first column of $A=QR$ satisfies $a_1=r_{11}q_1$, so

$$
r_{11}=\|a_1\|_2,\qquad q_1=\frac{a_1}{r_{11}}.
$$

For every remaining column, remove its component along $q_1$:

$$
r_{1j}=q_1^Ta_j,\qquad
a_j^{(1)}=a_j-q_1r_{1j},\qquad j>1.
$$

Here $a_j^{(1)}$ is the remaining part of $a_j$ after the first projection. It is perpendicular to $q_1$, because

$$
q_1^Ta_j^{(1)}
=q_1^Ta_j-(q_1^Tq_1)r_{1j}=0.
$$

Normalize $a_2^{(1)}$ to obtain $q_2$, then remove the $q_2$ component from each remaining column. Continuing in this way gives **modified Gram-Schmidt (MGS)**. Each projection coefficient is computed from the current residual column, after the preceding components have been removed.

At step $k$, write $a_j^{(k-1)}$ for this residual, with $a_j^{(0)}=a_j$. The formulas are

$$
\begin{aligned}
r_{kk}&=\|a_k^{(k-1)}\|_2,\\
q_k&=\frac{a_k^{(k-1)}}{r_{kk}},\\
r_{kj}&=q_k^Ta_j^{(k-1)},\qquad j>k,\\
a_j^{(k)}&=a_j^{(k-1)}-q_kr_{kj},\qquad j>k.
\end{aligned}
$$

The update is an orthogonal projection onto the space perpendicular to $q_k$:

$$
a_j^{(k)}=(I-q_kq_k^T)a_j^{(k-1)}.
$$

There is also an outer-product interpretation. If $r_k^T$ denotes row $k$ of $R$, then

$$
A=\sum_{k=1}^n q_kr_k^T.
$$

Each step extracts one rank-one term. Since $r_{kj}=0$ for $j<k$, only the remaining columns need to be updated.

## The In-Place Algorithm

The superscripts above explain the derivation. In an implementation, one working array is enough. Its column $k$ is normalized and stored as $q_k$; its later columns hold the current residuals. Initialize $R$ to zero.

~~~text
for k = 1, ..., n:
    r[k,k] = norm(A[:,k], 2)
    if r[k,k] == 0:
        stop: the current column adds no independent direction
    A[:,k] = A[:,k] / r[k,k]
    for j = k+1, ..., n:
        r[k,j] = A[:,k]^T A[:,j]
        A[:,j] = A[:,j] - A[:,k] r[k,j]
~~~

On completion, the working array $A$ contains $Q$, and $R$ is stored separately. The input matrix is overwritten; retain a copy if it is needed later.

The zero test describes breakdown in exact arithmetic. In floating-point arithmetic, a tiny residual may contain mostly rounding error even when its norm is nonzero. Detecting numerical dependence requires a tolerance relative to the input scale and desired accuracy. A small $r_{kk}$ alone does not establish that the original columns are dependent.

## Why the Algorithm Works

The proof needs to establish three facts: the vectors are orthonormal, normalization never divides by zero under our rank assumption, and the stored coefficients reconstruct $A$.

:::{prf:proof} Correctness of Modified Gram-Schmidt in Exact Arithmetic.
Suppose that $q_1,\ldots,q_{k-1}$ are orthonormal and span the first $k-1$ original columns. At the start of step $k$, each remaining column satisfies

$$
a_j=\sum_{i=1}^{k-1}q_ir_{ij}+a_j^{(k-1)},
\qquad q_i^Ta_j^{(k-1)}=0\quad(i<k).
$$

For $k=1$, the sum is empty and these statements hold with $a_j^{(0)}=a_j$.

The residual $a_k^{(k-1)}$ cannot be zero. Otherwise, the displayed identity would place $a_k$ in the span of its preceding columns, contradicting linear independence. Thus $r_{kk}>0$, and normalizing gives a unit vector $q_k$ perpendicular to the previous vectors. Also,

$$
a_k=\sum_{i=1}^{k-1}q_ir_{ik}+q_kr_{kk},
$$

so the first $k$ basis vectors span the first $k$ original columns.

For $j>k$, the update removes the component along $q_k$:

$$
q_k^Ta_j^{(k)}
=q_k^Ta_j^{(k-1)}-(q_k^Tq_k)r_{kj}=0.
$$

It also preserves perpendicularity to every earlier $q_i$, since

$$
q_i^Ta_j^{(k)}
=q_i^Ta_j^{(k-1)}-(q_i^Tq_k)r_{kj}=0,
\qquad i<k.
$$

Finally, substituting $a_j^{(k-1)}=q_kr_{kj}+a_j^{(k)}$ into the reconstruction identity extends its sum through $i=k$. This proves the induction step.

At completion, every original column satisfies

$$
a_j=\sum_{i=1}^j q_ir_{ij}.
$$

These column identities give $A=QR$, with $Q^TQ=I_n$ and $R$ upper triangular with positive diagonal.
:::

If the columns are dependent, the basic algorithm stops at the first zero residual. QR still exists, but completing an orthonormal basis requires an additional choice. This is one reason to distinguish existence of QR from successful completion of this particular algorithm.

## Comparison with Classical Gram-Schmidt

**Classical Gram-Schmidt (CGS)** computes all projection coefficients for a column from the original column:

$$
r_{ij}=q_i^Ta_j\quad(i<j),\qquad
v=a_j-\sum_{i=1}^{j-1}q_ir_{ij}.
$$

It then sets $r_{jj}=\|v\|_2$ and $q_j=v/r_{jj}$. MGS instead uses the successive residuals. For one column, its computation can be written as

$$
\begin{aligned}
v&=a_j,\\
r_{ij}&=q_i^Tv,\qquad v\leftarrow v-q_ir_{ij},
\quad i=1,\ldots,j-1,\\
r_{jj}&=\|v\|_2,\qquad q_j=v/r_{jj}.
\end{aligned}
$$

This column-by-column ordering computes the same MGS projections as the in-place algorithm above. It is convenient when new vectors are generated one at a time.

In exact arithmetic, CGS and MGS agree. Removing a component along $q_\ell$ does not change the inner product with another orthogonal vector $q_i$:

$$
q_i^T(v-q_\ell r_{\ell j})
=q_i^Tv-(q_i^Tq_\ell)r_{\ell j}=q_i^Tv,
\qquad i\ne\ell.
$$

In floating-point arithmetic, the computed basis vectors are not exactly orthogonal, so this equality no longer holds exactly. The choice of which vector enters each dot product then matters.

## Numerical Accuracy: The Condition-Number Factor

Use the two backward-error measures introduced for {ref}`Householder QR <qr-backward-error>`. For MGS, $\widehat Q$ is the computed reduced factor, so its orthogonality equation uses $I_n$.

For a full-column-rank matrix, let

$$
\kappa_2(A)=\frac{\sigma_{\max}(A)}{\sigma_{\min}(A)},
$$

where $\sigma_{\max}$ and $\sigma_{\min}$ are the largest and smallest singular values. Under the same rounding assumptions, MGS gives the bounds

$$
\boxed{
\begin{aligned}
\frac{\|A-\widehat Q\widehat R\|_F}{\|A\|_F}
&\le C_1(m,n)u,\\
\|\widehat Q^T\widehat Q-I_n\|_2
&\le C_2(m,n)u\kappa_2(A).
\end{aligned}
}
$$

The $C_i(m,n)$ are dimension-dependent factors independent of $A$. The orthogonality bound additionally requires $C_0(m,n)u\kappa_2(A)$ to be sufficiently small, so rounding errors remain below the scale of the smallest singular value.

**The new factor is $\kappa_2(A)$ in the orthogonality bound.** MGS still reconstructs $A$ accurately, but nearly dependent columns can cause substantial loss of orthogonality. Householder and Givens have no condition-number factor in either bound.

The reason is cancellation in the projection steps. A nearly dependent column leaves a residual much smaller than the quantities subtracted from it. Errors in those quantities can be large relative to the residual, and normalization turns them into errors in the new basis direction. MGS usually reduces this effect compared with one-pass CGS, but does not eliminate it. This motivates a second orthogonalization pass.

See [Higham's QR overview](https://nhigham.com/2020/11/10/what-is-a-qr-factorization/) for the comparison, and Chapter 19 of his [Accuracy and Stability of Numerical Algorithms](https://epubs.siam.org/doi/10.1137/1.9780898718027.ch19) for the analysis.

## Reorthogonalization

If accurate orthogonality is needed, project each residual against the previous basis vectors a second time before normalizing. This is called **reorthogonalization**. The second pass removes components left by rounding errors in the first pass.

For column $j$, first perform the column-by-column MGS pass above. Then, using its remaining vector $v$, compute

$$
\begin{aligned}
\delta r_{ij}&=q_i^Tv,\\
v&\leftarrow v-q_i\delta r_{ij},\\
r_{ij}&\leftarrow r_{ij}+\delta r_{ij},
\qquad i=1,\ldots,j-1.
\end{aligned}
$$

Finally set $r_{jj}=\|v\|_2$ and $q_j=v/r_{jj}$. Adding the correction to $r_{ij}$ is essential: both passes remove components that must be included in the identity $a_j=\sum_{i<j}q_ir_{ij}+q_jr_{jj}$.

For numerically independent columns, two passes can restore orthogonality to a dimension-dependent multiple of $u$. This guarantee requires rounding errors to be small relative to the independent part of each new column; a second pass cannot recover a direction that has already been lost in rounding. See [Giraud, Langou, and Rozloznik's analysis of reorthogonalization](https://www.cerfacs.fr/algor/reports/2002/TR_PA_02_33.pdf) for the precise assumptions and bounds.

### A Numerical Example

Consider three almost parallel columns:

$$
A_\varepsilon=
\begin{pmatrix}
1&1&1\\
\varepsilon&0&0\\
0&\varepsilon&0\\
0&0&\varepsilon
\end{pmatrix},\qquad \varepsilon>0.
$$

We have $A_\varepsilon^TA_\varepsilon=\boldsymbol{1}\boldsymbol{1}^T+\varepsilon^2I_3$, where $\boldsymbol{1}=(1,1,1)^T$. The rank-one matrix $\boldsymbol{1}\boldsymbol{1}^T$ has eigenvalue $3$ along $\boldsymbol{1}$ and eigenvalue $0$ on its perpendicular complement. Adding $\varepsilon^2I_3$ and taking square roots gives the singular values $\sqrt{3+\varepsilon^2},\varepsilon,\varepsilon$. Thus the columns are independent, but $\kappa_2(A_\varepsilon)=\sqrt{3+\varepsilon^2}/\varepsilon$ becomes large as $\varepsilon$ decreases.

The following experiment compares one-pass CGS, one-pass MGS, MGS with a second pass, and Householder QR. It uses the column-by-column ordering so the second pass is easy to identify.

~~~{code-cell} ipython3
import numpy as np

def gram_schmidt(A, classical=False, passes=1):
    m, n = A.shape
    Q = np.zeros((m, n))
    R = np.zeros((n, n))
    for j in range(n):
        v = A[:, j].copy()
        for _ in range(passes):
            if classical:
                # All coefficients use the same vector.
                coefficients = Q[:, :j].T @ v
                R[:j, j] += coefficients
                v -= Q[:, :j] @ coefficients
            else:
                # Each coefficient uses the updated residual.
                for i in range(j):
                    coefficient = Q[:, i] @ v
                    R[i, j] += coefficient
                    v -= coefficient * Q[:, i]
        R[j, j] = np.linalg.norm(v)
        if R[j, j] == 0:
            raise ValueError("Zero residual: cannot normalize this column.")
        Q[:, j] = v / R[j, j]
    return Q, R

epsilon = 1e-8
A = np.vstack((np.ones((1, 3)), epsilon * np.eye(3)))
factorizations = [
    ("CGS", gram_schmidt(A, classical=True)),
    ("MGS", gram_schmidt(A)),
    ("MGS, two passes", gram_schmidt(A, passes=2)),
    ("Householder", np.linalg.qr(A, mode="reduced")),
]
print(f"{'Method':<18} {'Factorization error':>20} {'Orthogonality error':>20}")
for name, (Q, R) in factorizations:
    residual = np.linalg.norm(A - Q @ R, "fro") / np.linalg.norm(A, "fro")
    orthogonality = np.linalg.norm(np.eye(3) - Q.T @ Q, 2)
    print(f"{name:<18} {residual:20.1e} {orthogonality:20.1e}")
~~~

All four methods reconstruct this input accurately. Their orthogonality errors differ substantially: MGS improves on CGS, while the second pass brings MGS close to Householder accuracy. The exact printed values can vary slightly with the numerical library.

## Computational Cost and Choice of Method

For each pair $k<j$, the dot product costs about $2m$ flops, and the residual update costs another $2m$. There are $n(n-1)/2$ such pairs. Computing a norm and normalizing each column costs about $3m$ flops, so the total is

$$
\begin{aligned}
\text{cost}
&=4m\frac{n(n-1)}{2}+O(mn)\\
&=\boxed{2mn^2+O(mn).}
\end{aligned}
$$

CGS has the same leading arithmetic cost. A second MGS pass doubles the projection work, giving $4mn^2+O(mn)$ flops when used for every column.

Householder QR costs $2mn^2-\tfrac23n^3+O(mn)$ flops to compute $R$ and store $Q$ implicitly; explicitly forming the reduced $Q$ requires additional work. MGS produces the reduced $Q$ directly. Thus being tall and thin does not by itself make MGS the preferred method. Householder QR is the usual choice for a dense matrix available in advance; MGS with reorthogonalization is useful when an explicit basis must be built sequentially.

For complex matrices, replace transposes in the inner products by conjugate transposes. The projection formulas and exact-arithmetic proof then apply with $Q^HQ=I_n$.
