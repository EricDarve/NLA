# The Method of Deflation

The [power method](power_method.md) finds the dominant eigenvector. To find the next direction, we can remove the component along the one already found and apply the power method again. This is **deflation**.

For a general matrix, this process constructs **Schur vectors**, which need not be eigenvectors of the original matrix. The key fact is that the next Schur vector becomes an eigenvector of the deflated matrix.

## Notation and the Target

Write an ordered Schur decomposition as

$$
A=Q_\star TQ_\star^H,
\qquad
Q_\star=[\boldsymbol q_{\star,1},\ldots,\boldsymbol q_{\star,n}],
$$

where $Q_\star$ is unitary and $T$ is upper triangular with diagonal entries $\lambda_1,\ldots,\lambda_n$. The star in $Q_\star$ marks an **exact Schur basis**. In the next section, $Q_k$ will denote the **computed basis at iteration $k$**.

For the sequential power-method argument, assume

$$
|\lambda_1|>|\lambda_2|>\cdots>|\lambda_n|>0.
$$

These strict gaps let us find one direction at a time. The projection identities below do not require these gaps.

The Schur relation $AQ_\star=Q_\star T$ says

$$
A\boldsymbol q_{\star,j}
=\sum_{\ell=1}^{j-1}t_{\ell j}\boldsymbol q_{\star,\ell}
+\lambda_j\boldsymbol q_{\star,j}.
$$

Thus the first $j$ Schur vectors span an invariant subspace: applying $A$ keeps us within that subspace. Only the first Schur vector is necessarily an eigenvector of $A$.

## Removing the First Direction

Suppose the power method has found the first direction exactly. The matrix

$$
P_1=I-\boldsymbol q_{\star,1}\boldsymbol q_{\star,1}^H
$$

is the **orthogonal projection onto the complement** of that direction. For the second Schur vector,

$$
A\boldsymbol q_{\star,2}
=t_{12}\boldsymbol q_{\star,1}
+\lambda_2\boldsymbol q_{\star,2}.
$$

Projection removes the first term and leaves the second unchanged:

$$
P_1A\boldsymbol q_{\star,2}
=\lambda_2\boldsymbol q_{\star,2}.
$$

This is the central idea. Although $\boldsymbol q_{\star,2}$ need not be an eigenvector of $A$, it **is** an eigenvector of $P_1A$. Also, $P_1A\boldsymbol q_{\star,1}=0$. The dominant direction has been removed.

## Removing Several Directions

Suppose the first $i$ Schur vectors are known. Collect them in

$$
Q_{\star,i}=[\boldsymbol q_{\star,1},\ldots,\boldsymbol q_{\star,i}],
\qquad
P_i=I-Q_{\star,i}Q_{\star,i}^H.
$$

Here the second subscript $i$ counts columns. The deflated matrix is $M_i=P_iA$.

````{prf:theorem} Exact Schur Deflation
:label: thm:exact-schur-deflation
For $1\leq i<n$, the eigenvalues of $M_i$ are $i$ zeros followed by $\lambda_{i+1},\ldots,\lambda_n$, counted with multiplicity. Moreover,

$$
M_i\boldsymbol q_{\star,i+1}
=\lambda_{i+1}\boldsymbol q_{\star,i+1}.
$$
````

````{prf:proof} Deflation in Schur Coordinates.
Partition the Schur form after its first $i$ rows and columns:

$$
T=\begin{pmatrix}T_{11}&T_{12}\\0&T_{22}\end{pmatrix}.
$$

In this basis, $P_i$ removes the first $i$ coordinates. Therefore,

$$
\begin{aligned}
Q_\star^HP_iQ_\star
&=\begin{pmatrix}0&0\\0&I_{n-i}\end{pmatrix},\\
Q_\star^HM_iQ_\star
&=\begin{pmatrix}0&0\\0&I_{n-i}\end{pmatrix}
  \begin{pmatrix}T_{11}&T_{12}\\0&T_{22}\end{pmatrix}\\
&=\begin{pmatrix}0&0\\0&T_{22}\end{pmatrix}.
\end{aligned}
$$

This upper triangular matrix has the claimed eigenvalues. Its first remaining coordinate is an eigenvector with eigenvalue $\lambda_{i+1}$. Equivalently, applying $P_i$ to the Schur relation for $\boldsymbol q_{\star,i+1}$ removes all earlier vectors and leaves $\lambda_{i+1}\boldsymbol q_{\star,i+1}$.
````

The leading Schur subspace is invariant, so $P_iAQ_{\star,i}=0$. Consequently, $P_iA=P_iAP_i$: the deflated matrix discards the known subspace and acts within its orthogonal complement.

## Power Iteration on the Remaining Subspace

Under the strict-gap assumption, $\lambda_{i+1}$ is the unique dominant eigenvalue of $M_i$. Start with a unit vector $\boldsymbol v_0$ orthogonal to the known Schur vectors, and repeat

$$
\begin{aligned}
\boldsymbol w_{\ell+1}
&=A\boldsymbol v_\ell
-Q_{\star,i}\bigl(Q_{\star,i}^HA\boldsymbol v_\ell\bigr),\\
\boldsymbol v_{\ell+1}
&=\frac{\boldsymbol w_{\ell+1}}{\|\boldsymbol w_{\ell+1}\|_2}.
\end{aligned}
$$

There is no need to form $P_i$ or $M_i$. We multiply by $A$, subtract the projection onto the known subspace, and normalize.

As in the power method, the starting vector must have a nonzero coefficient in the dominant direction when expanded in an eigenvector basis of $M_i$. A random start in the complement satisfies this with probability one in exact arithmetic. Then the line spanned by $\boldsymbol v_\ell$ converges to the line spanned by $\boldsymbol q_{\star,i+1}$. For $i\leq n-2$, the geometric convergence factor is governed by $|\lambda_{i+2}/\lambda_{i+1}|$. Once $n-1$ directions are known, their one-dimensional orthogonal complement supplies the last one.

This gives an inductive construction: find the first Schur direction by power iteration, project it out to find the second, and continue. Each new vector is orthogonal to the previous ones, and together they span a larger invariant subspace. The corresponding eigenvalue is $\lambda_j=\boldsymbol q_{\star,j}^HA\boldsymbol q_{\star,j}$.

In computation, the preceding vectors are approximate. [Orthogonal iteration](orthogonal_iteration.md) updates all of them at every step, using a QR factorization to perform these projections together. If eigenvalues have equal magnitudes, we may need to follow a subspace containing several directions instead of seeking one direction at a time.
