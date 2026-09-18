# The Orthogonal Iteration Algorithm

[Deflation](deflation.md) finds successive Schur directions by projecting out those already known. **Orthogonal iteration** performs these power and projection steps together: multiply a set of vectors by $A$, then orthonormalize them. The first column follows the power method; each later column removes the directions represented by the preceding columns.

We first explain this connection to deflation. Then we prove that the leading subspaces converge, with a rate controlled by the ratio of the first unwanted eigenvalue to the last wanted eigenvalue.

## Keep the Target Separate from the Iterates

We use the notation introduced in the deflation section:

| Notation | Meaning |
| :--- | :--- |
| $Q_\star$ | A fixed, exact Schur basis: $A=Q_\star TQ_\star^H$. |
| $Q_k$ | The computed $n\times p$ orthonormal basis at iteration $k$. |
| $\boldsymbol q_{\star,j}$, $\boldsymbol q_{k,j}$ | Column $j$ of the exact and computed bases. |
| $Q_{\star,i}$, $Q_{k,i}$ | The first $i$ columns of those bases. |

The iteration index is always $k$; $i$ or $j$ counts directions. We order the Schur form so that $|\lambda_1|\geq\cdots\geq|\lambda_n|$. The desired $p$-dimensional invariant subspace is

$$
\mathcal S_\star=\operatorname{range}(Q_{\star,p}).
$$

The Schur basis is a target for the analysis, not an input to the algorithm.

## The Algorithm

Choose $Q_0\in\mathbb C^{n\times p}$ with orthonormal columns, where $1\leq p\leq n$. For $k=0,1,2,\ldots$, compute

$$
\boxed{Z_{k+1}=AQ_k,\qquad Z_{k+1}=Q_{k+1}R_{k+1}.}
$$

The second step is an **unpivoted thin QR factorization**: $Q_{k+1}$ has $p$ orthonormal columns and $R_{k+1}$ is $p\times p$ upper triangular. Keeping the column order is important for the connection to successive deflation.

For now, assume $AQ_k$ has rank $p$, so $R_{k+1}$ is invertible. The convergence theorem below gives conditions that guarantee this. Full iteration ($p=n$) requires $A$ to be nonsingular for this argument.

A basic implementation is:

```python
import numpy as np

# A is n x n; choose 1 <= p <= n and a number of steps.
n = A.shape[0]
rng = np.random.default_rng()
G = rng.standard_normal((n, p))
if np.iscomplexobj(A):
    G = G + 1j * rng.standard_normal((n, p))
Q, _ = np.linalg.qr(G, mode="reduced")

for _ in range(num_iterations):
    Z = A @ Q
    Q, R = np.linalg.qr(Z, mode="reduced")

T = Q.conj().T @ A @ Q  # A represented in the computed subspace
```

QR preserves the span of the columns while keeping them orthonormal. Normalizing each column separately would not remove the dominant directions from the later columns.

## QR Performs the Deflation Steps

To see what QR does, write it column by column as Gram–Schmidt in exact arithmetic. Take the diagonal entries of $R_{k+1}$ to be positive; other phase choices give the same subspaces.

The first column satisfies

$$
\boldsymbol q_{k+1,1}
=\frac{A\boldsymbol q_{k,1}}{\|A\boldsymbol q_{k,1}\|_2}.
$$

This is the power method. Under its usual gap and starting-vector assumptions, its direction converges to $\boldsymbol q_{\star,1}$.

For column $i+1$, QR subtracts the components along the **newly computed** first $i$ columns:

$$
\begin{aligned}
\boldsymbol w_{k+1,i+1}
&=\bigl(I-Q_{k+1,i}Q_{k+1,i}^H\bigr)
  A\boldsymbol q_{k,i+1},\\
\boldsymbol q_{k+1,i+1}
&=\frac{\boldsymbol w_{k+1,i+1}}{\|\boldsymbol w_{k+1,i+1}\|_2}.
\end{aligned}
$$

Suppose for a moment that these first $i$ columns span the exact leading Schur subspace. Their complementary projector is then

$$
I-Q_{k+1,i}Q_{k+1,i}^H
=I-Q_{\star,i}Q_{\star,i}^H=P_i.
$$

The update for column $i+1$ becomes power iteration on $P_iA$. From deflation,

$$
P_iA\boldsymbol q_{\star,i+1}
=\lambda_{i+1}\boldsymbol q_{\star,i+1}.
$$

If $\lambda_{i+1}$ is the unique dominant eigenvalue of $P_iA$ and the starting direction contains the required eigenvector component, this power iteration finds the next Schur direction. This explains the inductive structure: the first column finds the first direction, the second projects it out to find the next, and so on.

In the actual algorithm, the earlier columns are still changing. The argument above explains the mechanism; the proof below establishes convergence while accounting for all columns evolving together.

## The Key Identity: QR Changes the Basis, Not the Subspace

Repeatedly substituting $AQ_k=Q_{k+1}R_{k+1}$ gives

$$
A^kQ_0=Q_k\bigl(R_kR_{k-1}\cdots R_1\bigr).
$$

The factor on the right is invertible, so

$$
\boxed{\operatorname{range}(Q_k)
=\operatorname{range}(A^kQ_0).}
$$

The product is also upper triangular. Therefore, for every $i\leq p$, its first $i$ columns involve only the first $i$ columns of $Q_k$, giving

$$
\operatorname{range}(Q_{k,i})
=\operatorname{range}(A^kQ_{0,i}).
$$

Each leading subspace undergoes its own power iteration. We can analyze these subspaces without tracking the individual projections inside QR. The [next section](orthogonal_and_power_iteration.md) develops this QR identity further.

## Convergence of a Leading Subspace

As in the power-method analysis, assume $A$ is diagonalizable. The theorem allows equal eigenvalue magnitudes within the wanted group or within the unwanted group; only the gap between the two groups must be strict.

### Measuring the Error

Let $\Pi_\star=Q_{\star,p}Q_{\star,p}^H$ be the orthogonal projection onto $\mathcal S_\star$. Define

$$
\begin{aligned}
d_k&=\|(I-\Pi_\star)Q_k\|_2\\
&=\max_{\substack{\boldsymbol z\in\operatorname{range}(Q_k)\\
                  \|\boldsymbol z\|_2=1}}
  \|(I-\Pi_\star)\boldsymbol z\|_2.
\end{aligned}
$$

Thus $d_k$ is the largest component outside the target subspace among unit vectors in the current subspace. It is zero when the two subspaces coincide. Geometrically, it is the sine of their largest principal angle. This measures the subspaces themselves, independent of the orthonormal bases used to represent them.

### The Convergence Theorem

````{prf:theorem} Convergence of a Dominant Subspace
:label: thm:orth-it-nh-rigorous
Suppose $A=X\Lambda X^{-1}$, and for some $1\leq p<n$,

$$
|\lambda_1|\geq\cdots\geq|\lambda_p|
>|\lambda_{p+1}|\geq\cdots\geq|\lambda_n|.
$$

Partition the eigenvector basis and the initial coordinates as

$$
X=[X_1,X_2],
\qquad
X^{-1}Q_0=\begin{pmatrix}C_1\\C_2\end{pmatrix},
$$

where $X_1$ has $p$ columns and $C_1$ is $p\times p$. Assume $C_1$ is invertible. Then orthogonal iteration is well-defined in exact arithmetic, and its subspace converges to

$$
\operatorname{range}(X_1)=\mathcal S_\star.
$$

For $k\geq1$, its error satisfies

$$
d_k\leq K\left(\frac{|\lambda_{p+1}|}{|\lambda_p|}\right)^k,
$$

where one possible constant is

$$
K=\|X_2\|_2\,\|X^{-1}\|_2\,\|C_2C_1^{-1}\|_2.
$$
````

The condition on $C_1$ says that the starting vectors contain enough independent components to capture all $p$ wanted eigenvector directions. For $p=1$, it is exactly the power method's requirement $c_1\neq0$. A random orthonormal start satisfies it with probability one in exact arithmetic.

For a normal matrix, we may take $X=Q_\star$, so the condition becomes invertibility of $Q_{\star,p}^HQ_0$. For a nonnormal matrix, the eigenvector coordinates $X^{-1}Q_0$ matter; ordinary orthogonal overlap with the target is not the same condition.

````{prf:proof} The Eigenvalue-Ratio Bound.
Partition $\Lambda=\operatorname{diag}(\Lambda_1,\Lambda_2)$ consistently with $X$. The gap implies that every diagonal entry of $\Lambda_1$ is nonzero. Then

$$
A^kQ_0=X_1\Lambda_1^kC_1+X_2\Lambda_2^kC_2.
$$

Its leading eigenvector-coordinate block, $\Lambda_1^kC_1$, is invertible. Thus $A^kQ_0$ has rank $p$ for every $k$, which ensures that all QR steps retain $p$ independent columns.

**Separate the desired directions from the error.** Right multiplication by $C_1^{-1}\Lambda_1^{-k}$ changes the basis but not the range. Using the QR identity,

$$
\operatorname{range}(Q_k)
=\operatorname{range}(X_1+X_2F_k),
$$

where

$$
F_k=\Lambda_2^kC_2C_1^{-1}\Lambda_1^{-k}.
$$

The matrix $F_k$ gives the unwanted eigenvector coordinates when the wanted coordinates are set to the identity. Since $\Lambda_1$ and $\Lambda_2$ are diagonal,

$$
\begin{aligned}
\|F_k\|_2
&\leq\|\Lambda_2^k\|_2\,
       \|C_2C_1^{-1}\|_2\,
       \|\Lambda_1^{-k}\|_2\\
&=\|C_2C_1^{-1}\|_2
  \left(\frac{|\lambda_{p+1}|}{|\lambda_p|}\right)^k.
\end{aligned}
$$

**Translate this into subspace error.** Any unit vector in the computed subspace has the form

$$
\boldsymbol z=(X_1+X_2F_k)\boldsymbol c,
\qquad
X^{-1}\boldsymbol z
=\begin{pmatrix}\boldsymbol c\\F_k\boldsymbol c\end{pmatrix}.
$$

Hence $\|\boldsymbol c\|_2\leq\|X^{-1}\|_2$. Also, $(I-\Pi_\star)X_1=0$, so

$$
\begin{aligned}
\|(I-\Pi_\star)\boldsymbol z\|_2
&=\|(I-\Pi_\star)X_2F_k\boldsymbol c\|_2\\
&\leq\|X_2\|_2\,\|F_k\|_2\,\|X^{-1}\|_2.
\end{aligned}
$$

Maximizing over unit vectors $\boldsymbol z$ and applying the bound on $F_k$ proves the result.
````

The ratio $|\lambda_{p+1}/\lambda_p|$ compares the fastest-growing unwanted component with the slowest-growing wanted component. A ratio near one gives slow separation; a smaller ratio gives faster separation. The theorem is an upper bound, and particular starting spaces may converge faster. For nonnormal matrices, the constant $K$ can be large.

For example, suppose the three largest eigenvalue magnitudes are $10$, $9.9$, and $1$, and the starting-space conditions hold. The error bound for the first vector decays like $0.99^k$, but the bound for the first two columns together decays like $(1/9.9)^k$. We can capture their two-dimensional invariant subspace long before the individual directions are resolved.

This proof assumes diagonalizability. The algorithm also applies to matrices that are not diagonalizable, but their convergence analysis can include additional polynomial factors in $k$.

## From Subspaces to the Schur Vectors

To recover all Schur directions, consider full iteration ($p=n$) and suppose now that

$$
|\lambda_1|>\cdots>|\lambda_n|>0,
$$

and that the starting-space condition holds for each leading set of columns being analyzed: the first $i$ rows of $X^{-1}Q_{0,i}$ form an invertible $i\times i$ matrix. Apply the theorem separately to each $i<n$. It gives

$$
\operatorname{range}(Q_{k,i})
\longrightarrow\operatorname{range}(Q_{\star,i}),
$$

with an error bounded by a constant times $|\lambda_{i+1}/\lambda_i|^k$.

Why does this identify column $i$? It is the unit direction in the first $i$ columns that is orthogonal to the first $i-1$. More explicitly, let

$$
\Pi_{k,i}=Q_{k,i}Q_{k,i}^H,
\qquad
\Pi_{\star,i}=Q_{\star,i}Q_{\star,i}^H.
$$

For equal-dimensional subspaces, the error used above also equals the distance between their orthogonal projections:

$$
\|\Pi_{k,i}-\Pi_{\star,i}\|_2
=\|(I-\Pi_{\star,i})Q_{k,i}\|_2.
$$

Thus these projections converge. Subtract the projections for two consecutive subspaces:

$$
\begin{aligned}
\boldsymbol q_{k,i}\boldsymbol q_{k,i}^H
&=\Pi_{k,i}-\Pi_{k,i-1}\\
&\longrightarrow\Pi_{\star,i}-\Pi_{\star,i-1}\\
&=\boldsymbol q_{\star,i}\boldsymbol q_{\star,i}^H.
\end{aligned}
$$

Here $\Pi_{k,0}=\Pi_{\star,0}=0$. For full iteration, $\Pi_{k,n}=\Pi_{\star,n}=I$, so the last column follows too. This proves convergence of all Schur directions. Signs or complex phases may continue to change, as in the power method; the lines converge.

For $p=n$, the full column space is always $\mathbb C^n$. It is the convergence of these **nested leading subspaces**, not the full space, that produces Schur form. Under the assumptions above, the diagonal entries of $Q_k^HAQ_k$ approach the ordered eigenvalues and the entries below the diagonal approach zero.

For partial iteration ($p<n$), if only $|\lambda_p|>|\lambda_{p+1}|$ holds, the $p$-dimensional subspace can converge without individual columns settling into Schur directions. This is useful when wanted eigenvalues have equal or nearly equal magnitudes. For real matrices, a nonreal conjugate pair has equal magnitudes and must be treated together when working in real arithmetic.

## Using the Computed Subspace

For $p<n$, form the projected matrix and its residual:

$$
T_k=Q_k^HAQ_k,
\qquad
E_k=AQ_k-Q_kT_k.
$$

A small residual means that applying $A$ nearly keeps vectors within the computed subspace. The eigenvalues of $T_k$, called **Ritz values**, then provide approximations to the wanted eigenvalues as the subspace converges. A residual alone does not identify which invariant subspace has been found.

For full iteration, this residual is always zero because $Q_k$ is square and unitary. Under the strict-gap assumptions, we instead monitor the entries below the diagonal of $T_k$. The following sections develop the connection to eigenvalue computation and the QR algorithm.
