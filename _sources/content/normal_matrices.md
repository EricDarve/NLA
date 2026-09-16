# Normal Matrices

An eigendecomposition $A=X\Lambda X^{-1}$ uses a basis of eigenvectors, but that basis need not be orthogonal. **Normal matrices** are exactly the matrices for which an orthonormal eigenvector basis exists. This makes both their geometry and their time evolution easier to describe.

## Definition and Examples

A matrix $A\in\mathbb{C}^{n\times n}$ is **normal** if it commutes with its conjugate transpose:

$$A^HA=AA^H.$$

For a real matrix, the condition is $A^TA=AA^T$.

Several familiar classes satisfy this condition:

- **Hermitian matrices**, $A^H=A$, are normal because both products equal $A^2$. Real Hermitian matrices are precisely real symmetric matrices.
- **Skew-Hermitian matrices**, $A^H=-A$, are normal because both products equal $-A^2$. Their real counterparts are skew-symmetric matrices.
- **Unitary matrices**, $A^HA=AA^H=I$, are normal. Their real counterparts are orthogonal matrices.
- **Diagonal matrices** are normal, even with complex diagonal entries, because $A^HA=AA^H=\operatorname{diag}(|a_{11}|^2,\ldots,|a_{nn}|^2)$.

Normality does not require symmetry or preservation of lengths. For example,

$$A=\begin{pmatrix}1&-2\\2&1\end{pmatrix}$$

satisfies $A^TA=AA^T=5I$. It is normal, but neither symmetric nor orthogonal. Its eigenvalues are $1\pm2i$.

Conversely, a diagonalizable matrix need not be normal. The example from the [previous section](eigendecomposition.md),

$$B=\begin{pmatrix}2&1\\0&1\end{pmatrix},$$

has two distinct eigenvalues, but

$$
B^TB=\begin{pmatrix}4&2\\2&2\end{pmatrix},\qquad
BB^T=\begin{pmatrix}5&1\\1&1\end{pmatrix}.
$$

Its eigenvectors form a basis, but no orthonormal eigenvector basis exists.

## The Spectral Theorem for Normal Matrices

````{prf:theorem} Spectral Theorem
:label: thm:spectral_theorem
A complex square matrix $A$ is normal if and only if it is unitarily diagonalizable:

$$A=Q\Lambda Q^H,$$

where $Q$ is unitary and $\Lambda$ is diagonal. The columns of $Q$ form an orthonormal basis of eigenvectors, and the diagonal entries of $\Lambda$ are their corresponding eigenvalues.
````

````{prf:proof} Spectral Theorem.
First suppose $A=Q\Lambda Q^H$. Then $A^H=Q\Lambda^HQ^H$, so

$$
\begin{aligned}
A^HA&=Q\Lambda^H\Lambda Q^H,\\
AA^H&=Q\Lambda\Lambda^HQ^H.
\end{aligned}
$$

The diagonal matrices $\Lambda$ and $\Lambda^H$ commute. Hence $A$ is normal.

Conversely, suppose $A$ is normal. Its complex Schur decomposition is $A=QTQ^H$, with $T$ upper triangular. Normality is preserved by this unitary change of basis:

$$T^HT-TT^H=Q^H(A^HA-AA^H)Q=0.$$

We show that an upper triangular normal matrix must be diagonal. Comparing the $(1,1)$ entries gives

$$
\begin{aligned}
(T^HT)_{11}&=|t_{11}|^2,\\
(TT^H)_{11}&=|t_{11}|^2+\sum_{j=2}^n|t_{1j}|^2.
\end{aligned}
$$

Equality forces $t_{12}=\cdots=t_{1n}=0$. Therefore $T$ is block diagonal:

$$T=\begin{pmatrix}t_{11}&0\\0&T_2\end{pmatrix}.$$

The equation $T^HT=TT^H$ now implies $T_2^HT_2=T_2T_2^H$. Thus $T_2$ is also upper triangular and normal. Repeating the argument, or using induction with the scalar case as the base, proves that $T$ is diagonal. Set $\Lambda=T$.
````

This is the precise link between Schur form and eigendecomposition: **normality forces the triangular Schur factor to be diagonal**.

### Orthogonality of Eigenspaces

For a normal matrix, eigenvectors belonging to distinct eigenvalues are orthogonal. This is stronger than the linear independence guaranteed for a general matrix.

````{prf:proof} Orthogonal Eigenvectors.
The spectral theorem gives $A=Q\Lambda Q^H$ and $A^H=Q\Lambda^HQ^H$. If $Ax=\lambda x$, the coordinates $Q^Hx$ can be nonzero only at indices whose diagonal entry in $\Lambda$ equals $\lambda$. Consequently,

$$A^Hx=\overline{\lambda}x.$$

Now let $Ay=\mu y$, with $\lambda\ne\mu$. Then

$$
\begin{aligned}
\lambda x^Hy
&=(A^Hx)^Hy\\
&=x^HAy\\
&=\mu x^Hy.
\end{aligned}
$$

Thus $(\lambda-\mu)x^Hy=0$, so $x^Hy=0$.
````

Vectors in the same eigenspace need not be orthogonal. We can choose an orthonormal basis within each eigenspace, including when an eigenvalue is repeated. Normal matrices are always diagonalizable and therefore never defective.

## Hermitian, Skew-Hermitian, and Unitary Matrices

Normality allows arbitrary complex eigenvalues. The additional structure of these special classes restricts where their eigenvalues lie:

| Class | Defining condition | Eigenvalues |
|---|---|---|
| Hermitian | $A^H=A$ | Real |
| Skew-Hermitian | $A^H=-A$ | Purely imaginary, including zero |
| Unitary | $A^HA=I$ | $\lvert\lambda\rvert=1$ |

For a **normal** matrix, each eigenvalue condition also implies the corresponding matrix condition. Indeed, using $A=Q\Lambda Q^H$,

$$
\begin{aligned}
A^H=A&\quad\Longleftrightarrow\quad\Lambda^H=\Lambda,\\
A^H=-A&\quad\Longleftrightarrow\quad\Lambda^H=-\Lambda,\\
A^HA=I&\quad\Longleftrightarrow\quad\Lambda^H\Lambda=I.
\end{aligned}
$$

The normality assumption matters in these converses: the nonnormal matrix $B$ above has real eigenvalues but is not Hermitian.

The reality of Hermitian eigenvalues also has a short direct proof.

````{prf:proof} Real Eigenvalues.
Let $A^H=A$ and $Ax=\lambda x$, with $x\ne0$. Then

$$\lambda=\frac{x^HAx}{x^Hx}.$$

The denominator is positive and real. The numerator is real because

$$\overline{x^HAx}=x^HA^Hx=x^HAx.$$

Therefore $\lambda\in\mathbb{R}$.
````

For complex matrices, **symmetric** means $A^T=A$, whereas **Hermitian** means $A^H=A$. These are different conditions. For example, $\operatorname{diag}(i,1)$ is complex symmetric, but one of its eigenvalues is $i$. The real-eigenvalue conclusion applies to Hermitian matrices, including real symmetric matrices.

## The Spectral Theorem for Real Symmetric Matrices

````{prf:theorem} Spectral Theorem for Real Symmetric Matrices
:label: thm:spectral_theorem_real
A real square matrix $A$ is symmetric if and only if it is orthogonally diagonalizable:

$$A=Q\Lambda Q^T,$$

where $Q$ is real orthogonal and $\Lambda$ is real diagonal.
````

````{prf:proof} Real Spectral Theorem.
Suppose $A$ is real symmetric. It is Hermitian, so all its eigenvalues are real. Its real Schur decomposition therefore has only $1\times1$ diagonal blocks:

$$A=QTQ^T,$$

where $Q$ is real orthogonal and $T$ is upper triangular. Also,

$$T^T=(Q^TAQ)^T=Q^TA^TQ=T.$$

A matrix that is both symmetric and upper triangular must be diagonal. Taking $\Lambda=T$ proves the decomposition.

Conversely, if $A=Q\Lambda Q^T$ with real diagonal $\Lambda$, then

$$A^T=Q\Lambda^TQ^T=Q\Lambda Q^T=A.$$
````

A real normal matrix need not have real eigenvalues, so it need not be diagonalizable by a real orthogonal matrix. The normal matrix with eigenvalues $1\pm2i$ at the start of this section is one example. It has an orthonormal eigenvector basis over $\mathbb{C}$.

### Geometry and a Worked Example

For a real symmetric matrix, the factors in $Q\Lambda Q^T$ have a concrete meaning:

1. $Q^Tv$ gives the coordinates of $v$ in the orthonormal eigenvector basis.
2. $\Lambda$ scales each coordinate by its eigenvalue.
3. $Q$ converts the result back to the original coordinates.

Orthogonal coordinate changes can involve rotations or reflections. The scaling step stretches or compresses by $|\lambda_j|$, reverses the corresponding direction if $\lambda_j<0$, and removes it if $\lambda_j=0$.

For example,

$$
A=\begin{pmatrix}2&1\\1&2\end{pmatrix}
$$

has the orthonormal eigenvectors

$$
q_1=\frac{1}{\sqrt2}\begin{pmatrix}1\\1\end{pmatrix},\qquad
q_2=\frac{1}{\sqrt2}\begin{pmatrix}1\\-1\end{pmatrix},
$$

with eigenvalues $3$ and $1$. Thus $Q=[q_1,q_2]$ and $\Lambda=\operatorname{diag}(3,1)$. Along the line $v_1=v_2$, the matrix triples lengths; along the perpendicular line $v_1=-v_2$, it leaves vectors unchanged.

## Spectral Decomposition into Orthogonal Projections

Write $Q=[q_1,\ldots,q_n]$. Expanding $A=Q\Lambda Q^H$ gives

$$A=\sum_{j=1}^n\lambda_jq_jq_j^H.$$

Each $q_jq_j^H$ is the orthogonal projection onto the line spanned by $q_j$. Thus applying $A$ means separating a vector into orthogonal eigenvector components and scaling each component by its eigenvalue.

For a repeated eigenvalue, group all its eigenvectors together. Let $Q_\lambda$ contain an orthonormal basis of its eigenspace. Then $P_\lambda=Q_\lambda Q_\lambda^H$ is the **orthogonal projection onto the eigenspace associated with $\lambda$**. These projections satisfy

$$
\begin{aligned}
A&=\sum_{\lambda}\lambda P_\lambda,\\
I&=\sum_{\lambda}P_\lambda,\\
P_\lambda^2&=P_\lambda=P_\lambda^H,\\
P_\lambda P_\mu&=0\quad(\lambda\ne\mu).
\end{aligned}
$$

The sums run over distinct eigenvalues. These identities follow from orthonormality of the columns of $Q$. Although the basis within an eigenspace is not unique, its orthogonal projection $P_\lambda$ is unique.

The null space is spanned by the eigenvectors with eigenvalue zero, and the range by those with nonzero eigenvalues. Since $A^H$ has the same eigenvectors with conjugate eigenvalues,

$$
\begin{aligned}
N(A)&=N(A^H),\\
R(A)&=R(A^H)=N(A)^\perp.
\end{aligned}
$$

## Norms and Time Evolution

For a normal matrix, write $v=Qc$. Since $Q$ preserves the Euclidean norm,

$$
\begin{aligned}
\|v\|_2^2&=\sum_{j=1}^n|c_j|^2,\\
\|Av\|_2^2&=\sum_{j=1}^n|\lambda_j|^2|c_j|^2.
\end{aligned}
$$

It follows that

$$\|Av\|_2\le\left(\max_j|\lambda_j|\right)\|v\|_2.$$

Equality is attained by a unit eigenvector corresponding to an eigenvalue of largest magnitude. Hence the operator 2-norm is

$$\boxed{\|A\|_2=\max_j|\lambda_j|.}$$

This equality need not hold for a nonnormal matrix. For $B$ above, the unit vector $v=(1,1)^T/\sqrt2$ satisfies $\|Bv\|_2=\sqrt5>2$, even though the largest eigenvalue magnitude is $2$.

The same reasoning applies to repeated evolution $v_{k+1}=Av_k$. If $v_0=\sum_jc_jq_j$, then for $k\ge1$,

$$\|v_k\|_2^2=\sum_{j=1}^n|c_j|^2|\lambda_j|^{2k}.$$

For normal matrices, the orthogonality of the modes lets us add their squared magnitudes directly. If all $|\lambda_j|\le1$, the state norm never increases; if all $|\lambda_j|<1$, every state tends to zero. If $A$ is unitary, every $|\lambda_j|=1$ and the state norm is preserved.

## Positive Definite Matrices

A Hermitian matrix is **positive definite** if $v^HAv>0$ for every nonzero $v$, and **positive semidefinite** if $v^HAv\ge0$ for every $v$. These quantities are real because $A$ is Hermitian.

In orthonormal eigenvector coordinates, $v=Qc$ gives

$$v^HAv=c^H\Lambda c=\sum_{j=1}^n\lambda_j|c_j|^2.$$

Therefore, a Hermitian matrix is positive definite exactly when every eigenvalue is positive, and positive semidefinite exactly when every eigenvalue is nonnegative. The forward implications follow by choosing $v=q_j$; the reverse implications follow from the sum above.

For real symmetric matrices, replace $H$ by $T$. We will use these characterizations when studying [Cholesky factorization](cholesky.md) and [conjugate gradients](conjugate_gradient.md).
