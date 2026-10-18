# The Orthogonal Iteration Algorithm

[Deflation](deflation.md) finds a new Schur direction by projecting out the directions already known. Its explanation assumed those directions were exact. In computation, we would have to wait for each power iteration to converge before starting the next one.

**Orthogonal iteration** updates all the directions together. At each step, it multiplies the current basis vectors by $A$, then uses QR to project out earlier directions and normalize. The projections use the current approximations rather than exact Schur vectors.

## The Exact Target and the Current Basis

For $A\in\mathbb C^{n\times n}$, write an ordered Schur decomposition as $A=Q_\star TQ_\star^H$, with $|\lambda_1|\ge\cdots\ge|\lambda_n|$, as in deflation. Choose how many directions to follow: $1\le p\le n$.

| Notation | Meaning |
| :--- | :--- |
| $Q_\star$ | The fixed, exact $n\times n$ Schur basis. |
| $Q_k$ | The $n\times p$ orthonormal basis at iteration $k$. |
| $\boldsymbol q_{\star,j}$, $\boldsymbol q_{k,j}$ | Column $j$ of the exact and current bases. |
| $Q_{\star,i}$, $Q_{k,i}$ | The first $i$ columns of those bases. |

Thus $k$ counts iterations and $i$ or $j$ counts columns. The target $p$-dimensional invariant subspace is

$$
\mathcal S_\star=\operatorname{range}(Q_{\star,p}).
$$

The exact basis describes what we want to find.

## Multiply, Then Orthonormalize

Start with $Q_0\in\mathbb C^{n\times p}$ whose columns are orthonormal. For $k=0,1,\ldots$, compute

$$
\boxed{Z_{k+1}=AQ_k,\qquad Z_{k+1}=Q_{k+1}R_{k+1}.}
$$

The second step is an **unpivoted thin QR factorization**. The columns of $Q_{k+1}$ are orthonormal, and $R_{k+1}$ is $p\times p$ upper triangular. Column order matters: each new direction is made perpendicular to the directions before it. Normalizing the columns separately would allow all of them to approach the same dominant direction.

For the analysis, assume that $AQ_k$ has rank $p$ at every step, so each $R_{k+1}$ is invertible. The convergence theorem below gives conditions that guarantee this.

## How QR Performs Approximate Deflation

Write the QR step as Gram–Schmidt in exact arithmetic, choosing positive diagonal entries in $R_{k+1}$. Other sign or phase choices give the same subspaces.

The first column is simply power iteration:

$$
\boldsymbol q_{k+1,1}
=\frac{A\boldsymbol q_{k,1}}{\|A\boldsymbol q_{k,1}\|_2}.
$$

For the second column, multiply $\boldsymbol q_{k,2}$ by $A$, subtract its component along the **new first column** $\boldsymbol q_{k+1,1}$, and normalize. Continue in the same way. For column $i+1$, the update is

$$
\begin{aligned}
\boldsymbol w_{k+1,i+1}
&=\bigl(I-Q_{k+1,i}Q_{k+1,i}^H\bigr)
  A\boldsymbol q_{k,i+1},\\
\boldsymbol q_{k+1,i+1}
&=\frac{\boldsymbol w_{k+1,i+1}}{\|\boldsymbol w_{k+1,i+1}\|_2}.
\end{aligned}
$$

The projector uses the newly computed first $i$ columns because QR processes the columns of $AQ_k$ in order. If they spanned the exact leading Schur subspace, then

$$
I-Q_{k+1,i}Q_{k+1,i}^H
=I-Q_{\star,i}Q_{\star,i}^H=P_i.
$$

The update would then be the projected power step from deflation. Its target satisfies

$$
P_iA\boldsymbol q_{\star,i+1}
=\lambda_{i+1}\boldsymbol q_{\star,i+1}.
$$

In the actual iteration, the first $i$ directions are still improving while we update direction $i+1$. This explains the mechanism, but it is not yet a convergence proof: the projector changes at every step. We will account for those changes by following the column spaces.

### The Same Two-Dimensional Example

Return to the matrix from deflation and start with a rotated orthonormal basis:

$$
A=\begin{pmatrix}4&1\\0&2\end{pmatrix},
\qquad
Q_0=\frac{1}{\sqrt2}\begin{pmatrix}1&-1\\1&1\end{pmatrix}.
$$

One multiplication and QR step give

$$
AQ_0=\frac{1}{\sqrt2}\begin{pmatrix}5&-3\\2&2\end{pmatrix},
\qquad
Q_1=\frac{1}{\sqrt{29}}\begin{pmatrix}5&-2\\2&5\end{pmatrix}.
$$

The first column moves closer to the horizontal dominant direction. The second stays perpendicular to the new first column, so it moves closer to the vertical second Schur direction. Neither direction has to converge before the other is updated.

:::{figure} ../_static/orthogonal_iteration.svg
:alt: Four diagrams show both orthonormal columns after zero, one, three, and six orthogonal iterations for the matrix with rows (4,1) and (0,2). The blue first column approaches the horizontal axis. The purple second column approaches the vertical axis while remaining perpendicular to the first column at every step. Dashed arrows show the two exact Schur directions.
:width: 850px
:name: fig-orthogonal-iteration

Both directions evolve together. The blue column follows power iteration; QR keeps the purple column perpendicular to it. The two columns span the whole plane from the start. What improves is the first column's line and, consequently, the perpendicular second direction.
:::

(sec:oi-qr-ak)=
## The Connection to Powers of $A$

The changing projections make a column-by-column convergence proof awkward. A simpler observation removes them from the analysis: **QR changes the basis, but preserves the column space**.

Define the accumulated triangular factor

$$
\widehat R_0=I_p,
\qquad
\widehat R_k=R_kR_{k-1}\cdots R_1.
$$

:::{prf:theorem} QR of $A^kQ_0$
:label: thm:AkQR

For every $k\ge0$, orthogonal iteration satisfies

$$
A^kQ_0=Q_k\widehat R_k.
$$
:::

:::{prf:proof}
For $k=0$, the identity is $Q_0=Q_0I_p$. If it holds at step $k$, then

$$
\begin{aligned}
A^{k+1}Q_0
&=A(Q_k\widehat R_k)\\
&=(AQ_k)\widehat R_k\\
&=Q_{k+1}R_{k+1}\widehat R_k\\
&=Q_{k+1}\widehat R_{k+1}.
\end{aligned}
$$

This proves the identity by induction.
:::

Because $\widehat R_k$ is upper triangular and invertible, this is a thin QR factorization of $A^kQ_0$. With positive diagonal entries in each $R_j$, it is the unique thin QR factorization with positive diagonal.

For $p=1$, $\widehat R_k$ is the product of the normalization factors in the power method. For full iteration with $Q_0=I_n$, the identity becomes

$$
A^k=Q_k\widehat R_k.
$$

The algorithm obtains these factors without forming $A^k$.

More importantly, invertibility gives

$$
\boxed{\operatorname{range}(Q_k)
=\operatorname{range}(A^kQ_0).}
$$

Upper triangularity gives the same result for every leading set of columns. Let $\widehat R_{k,i}$ be the leading $i\times i$ block of $\widehat R_k$. The first $i$ columns of the factorization satisfy

$$
A^kQ_{0,i}=Q_{k,i}\widehat R_{k,i},
\qquad
\operatorname{range}(Q_{k,i})
=\operatorname{range}(A^kQ_{0,i}).
$$

Thus each leading column space undergoes a power iteration of its own. QR keeps its basis orthonormal without changing that space.

## Convergence of a Leading Subspace

Orthogonal iteration aims to approximate the leading $p$-dimensional Schur subspace. The theorem below shows how quickly that approximation improves. This is the main convergence result for the method: it extends the power method's convergence theory from one direction to a space of $p$ directions.

:::{admonition} Key result: convergence of orthogonal iteration
:class: important

Under the assumptions below, the space spanned by the columns of $Q_k$ approaches the space spanned by the first $p$ exact Schur vectors. Its error is bounded by a constant times $\rho_p^k$, where

$$
\boxed{\rho_p=\frac{|\lambda_{p+1}|}{|\lambda_p|}<1,
\qquad p<n.}
$$

For $p=1$, the target is the dominant eigenvector's line and $\rho_1=|\lambda_2/\lambda_1|$: the familiar power-method convergence factor. For $p>1$, the result controls the accuracy of an entire approximate Schur subspace.
:::

To understand the ratio, recall that $\operatorname{range}(Q_k)=\operatorname{range}(A^kQ_0)$. Assuming $A$ is diagonalizable, multiplication by $A^k$ scales eigenvector component $i$ by $\lambda_i^k$. We want to retain the first $p$ eigenvector directions and suppress the rest. Their span is the same as that of the first $p$ Schur vectors.

Among the wanted components, $|\lambda_p|^k$ is the smallest scaling factor. Among the unwanted components, $|\lambda_{p+1}|^k$ is the largest. Comparing the largest unwanted scaling with the smallest wanted scaling gives

$$
\frac{|\lambda_{p+1}|^k}{|\lambda_p|^k}=\rho_p^k.
$$

This is a relative comparison, so it applies whether the components grow or decay. A smaller $\rho_p$ gives faster separation of the wanted space from the unwanted directions. Only the gap between these two groups must be strict; eigenvalue magnitudes within either group may be equal. The subspace can therefore converge even while its individual basis vectors continue to change within it.

For example, suppose the three largest eigenvalue magnitudes are $10$, $9.9$, and $1$, and the starting-space conditions hold. For the first direction, the error bound decays like $0.99^k$. For the first two columns together, it decays like $(1/9.9)^k$, with a convergence factor of about $0.10$. Thus a two-dimensional Schur subspace can become accurate much sooner than its individual Schur directions are resolved.

There is also a starting-space requirement. The power method needs a nonzero coefficient in the wanted eigenvector direction. For $p$ directions, the starting vectors must contain $p$ independent combinations of the wanted eigenvector components. Otherwise, repeated multiplication cannot recover all of them. The theorem makes this requirement precise.

### Measuring the Subspace Error

To measure the accuracy of the whole subspace, look at how far its unit vectors extend outside the target space. Changes of basis within the target should not count as error. Let $\Pi_\star=Q_{\star,p}Q_{\star,p}^H$ project onto that target. Define

$$
\begin{aligned}
d_k&=\|(I-\Pi_\star)Q_k\|_2\\
&=\max_{\substack{\boldsymbol z\in\operatorname{range}(Q_k)\\
                  \|\boldsymbol z\|_2=1}}
  \|(I-\Pi_\star)\boldsymbol z\|_2.
\end{aligned}
$$

The second expression follows because the columns of $Q_k$ are orthonormal. It measures the largest component outside the target among unit vectors in the current subspace. It is zero exactly when the two subspaces coincide and does not depend on which orthonormal bases represent them. For $p=1$, it measures the misalignment of the iterate with the dominant eigenvector's line.

### The Convergence Theorem

:::{prf:theorem} Convergence of Orthogonal Iteration
:label: thm:orth-it-nh-rigorous

Suppose $A=X\Lambda X^{-1}$, and for some $1\le p<n$,

$$
|\lambda_1|\ge\cdots\ge|\lambda_p|
>|\lambda_{p+1}|\ge\cdots\ge|\lambda_n|.
$$

Partition the eigenvector basis and the starting coordinates as

$$
X=[X_1,X_2],
\qquad
X^{-1}Q_0=\begin{pmatrix}C_1\\C_2\end{pmatrix},
$$

where $X_1$ has $p$ columns and $C_1$ is $p\times p$. If $C_1$ is invertible, orthogonal iteration is well-defined in exact arithmetic and its column space converges to

$$
\mathcal S_\star=\operatorname{range}(X_1)
=\operatorname{range}(Q_{\star,p}).
$$

For $k\ge1$,

$$
d_k\le K\rho_p^k,
\qquad
K=\|X_2\|_2\,\|X^{-1}\|_2\,\|C_2C_1^{-1}\|_2.
$$
:::

Invertibility of $C_1$ means that the starting vectors contain enough independent components to recover all $p$ wanted directions. For $p=1$, it reduces to $c_1\ne0$ in the power method. A random orthonormal start satisfies it with probability one in exact arithmetic.

These are coordinates in the eigenvector basis. For a normal matrix, we can take $X=Q_\star$, so $C_1=Q_{\star,p}^HQ_0$. For a nonnormal matrix, ordinary orthogonal overlap with the target does not replace the condition on $X^{-1}Q_0$.

:::{prf:proof} The Eigenvalue-Ratio Bound.
Partition $\Lambda=\operatorname{diag}(\Lambda_1,\Lambda_2)$ consistently with $X$. Repeated multiplication gives

$$
A^kQ_0=X_1\Lambda_1^kC_1+X_2\Lambda_2^kC_2.
$$

The gap implies that $\Lambda_1$ is invertible. The first $p$ rows of $X^{-1}A^kQ_0$ form the invertible matrix $\Lambda_1^kC_1$. Hence $A^kQ_0$ has rank $p$ at every step. Inductively, $AQ_k$ has the same rank as $A^{k+1}Q_0$, so every QR factor $R_k$ is invertible.

**Separate the wanted coordinates from the unwanted ones.** Right multiplication by $C_1^{-1}\Lambda_1^{-k}$ changes the basis but not the column space:

$$
\begin{aligned}
A^kQ_0C_1^{-1}\Lambda_1^{-k}
&=X_1+X_2F_k,\\
F_k&=\Lambda_2^kC_2C_1^{-1}\Lambda_1^{-k}.
\end{aligned}
$$

Thus $\operatorname{range}(Q_k)=\operatorname{range}(X_1+X_2F_k)$. The matrix $F_k$ contains the unwanted coordinates after the wanted coordinates have been set to the identity. Since both eigenvalue blocks are diagonal,

$$
\begin{aligned}
\|F_k\|_2
&\le\|\Lambda_2^k\|_2\,\|C_2C_1^{-1}\|_2\,
       \|\Lambda_1^{-k}\|_2\\
&=\|C_2C_1^{-1}\|_2
  \left(\frac{|\lambda_{p+1}|}{|\lambda_p|}\right)^k.
\end{aligned}
$$

**Translate this into geometric error.** Any unit vector in the current subspace has the form

$$
\boldsymbol z=(X_1+X_2F_k)\boldsymbol c,
\qquad
X^{-1}\boldsymbol z
=\begin{pmatrix}\boldsymbol c\\F_k\boldsymbol c\end{pmatrix}.
$$

Therefore,

$$
\|\boldsymbol c\|_2
\le\|X^{-1}\boldsymbol z\|_2
\le\|X^{-1}\|_2.
$$

Projection removes the wanted part because $(I-\Pi_\star)X_1=0$. It follows that

$$
\begin{aligned}
\|(I-\Pi_\star)\boldsymbol z\|_2
&=\|(I-\Pi_\star)X_2F_k\boldsymbol c\|_2\\
&\le\|X_2\|_2\,\|F_k\|_2\,\|X^{-1}\|_2.
\end{aligned}
$$

Maximizing over unit vectors $\boldsymbol z$ and substituting the bound on $F_k$ gives $d_k\le K\rho_p^k$.
:::

A ratio near one gives slow separation; a smaller ratio gives faster separation. The constant $K$ depends on the initial space and the eigenvector basis and can be large for nonnormal matrices. The theorem gives an upper bound; particular starting spaces can converge faster.

This proof assumes diagonalizability. Orthogonal iteration can also be applied to matrices that are not diagonalizable, but their convergence bounds may include polynomial factors in $k$. [Arbenz's notes on subspace iteration](https://people.inf.ethz.ch/arbenz/ewp/Lnotes/chapter8.pdf) give further analysis.

## From Nested Subspaces to Schur Directions

For full iteration, $p=n$, the columns of $Q_k$ always span $\mathbb C^n$. Convergence of that full space says nothing about the individual directions. The two-dimensional example illustrates what matters: the **leading column spaces** converge.

Assume now that

$$
|\lambda_1|>\cdots>|\lambda_n|>0,
$$

and that the starting-space condition holds for each leading set of columns: the first $i$ rows of $X^{-1}Q_{0,i}$ form an invertible $i\times i$ matrix. The nonzero eigenvalues ensure that full iteration has invertible QR factors. Applying the theorem separately for each $i<n$ gives

$$
\operatorname{range}(Q_{k,i})
\longrightarrow\operatorname{range}(Q_{\star,i}),
$$

with error bounded by a constant times $|\lambda_{i+1}/\lambda_i|^k$.

Column $i$ is the unit direction in the first $i$ columns that is perpendicular to the first $i-1$. To express this precisely, let

$$
\Pi_{k,i}=Q_{k,i}Q_{k,i}^H,
\qquad
\Pi_{\star,i}=Q_{\star,i}Q_{\star,i}^H.
$$

For equal-dimensional subspaces, their orthogonal projectors satisfy

$$
\|\Pi_{k,i}-\Pi_{\star,i}\|_2
=\|(I-\Pi_{\star,i})Q_{k,i}\|_2.
$$

Thus their projectors converge too. Subtracting two consecutive projectors isolates column $i$:

$$
\begin{aligned}
\boldsymbol q_{k,i}\boldsymbol q_{k,i}^H
&=\Pi_{k,i}-\Pi_{k,i-1}\\
&\longrightarrow\Pi_{\star,i}-\Pi_{\star,i-1}\\
&=\boldsymbol q_{\star,i}\boldsymbol q_{\star,i}^H.
\end{aligned}
$$

Here $\Pi_{k,0}=\Pi_{\star,0}=0$, and $\Pi_{k,n}=\Pi_{\star,n}=I$. This proves convergence of every Schur direction, including the last one.

As in the power method, signs or complex phases can keep changing. Nevertheless, the diagonal entries of

$$
T_k=Q_k^HAQ_k
$$

approach the ordered eigenvalues, and the entries below the diagonal approach zero. Entries above the diagonal can continue to change with those phases.

For partial iteration, a boundary gap $|\lambda_p|>|\lambda_{p+1}|$ alone guarantees convergence of the $p$-dimensional space, not its individual basis directions. Additional internal gaps and starting-space conditions give convergence of the corresponding leading spaces. A nonreal conjugate pair must be treated together in real arithmetic because its two eigenvalues have equal magnitude.

## A Basic Implementation

The following code performs a fixed number of steps under the full-column-rank assumption above. NumPy may choose different QR signs or phases; these do not affect the computed subspaces.

~~~python
import numpy as np

def orthogonal_iteration(A, p, num_iterations):
    n = A.shape[0]
    rng = np.random.default_rng()
    G = rng.standard_normal((n, p))
    if np.iscomplexobj(A):
        G = G + 1j * rng.standard_normal((n, p))
    Q, _ = np.linalg.qr(G, mode="reduced")

    for _ in range(num_iterations):
        Z = A @ Q
        Q, R = np.linalg.qr(Z, mode="reduced")

    return Q
~~~

For dense $A$, each step costs $O(n^2p)$ for multiplication and $O(np^2)$ for QR. Using fewer columns is useful when only a few eigenvalues are wanted; sparse matrix products can also reduce the multiplication cost.

For $p<n$, form

$$
T_k=Q_k^HAQ_k,
\qquad
E_k=AQ_k-Q_kT_k.
$$

A small $\|E_k\|_2$ means that applying $A$ nearly keeps vectors within the current subspace. It can be used as a stopping criterion relative to the scale of $A$. As the subspace converges, the eigenvalues of $T_k$, called **Ritz values**, approach the wanted eigenvalues. A small residual alone does not identify which invariant subspace has been found.

For full iteration, $E_k=0$ at every step because $Q_k$ is square and unitary. Under the strict-gap assumptions, the entries below the diagonal of $T_k$ measure progress toward Schur form. The [next section](orthogonal_iteration_and_eigenvalues.md) develops the eigenvalue interpretation and the role of shifts.
