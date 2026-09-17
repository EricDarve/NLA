# Eigenvalues and Singular Values

Eigenvalues and singular values answer different questions about a matrix. An eigenvalue describes how the matrix acts along an eigenvector: $Ax=\lambda x$. A singular value describes how it maps one orthonormal direction to another: $Av_i=\sigma_i u_i$.

Eigenvalues are defined for square matrices and may be complex. Singular values exist for rectangular matrices as well and are always real and nonnegative.

For an evolving system $x_{k+1}=Ax_k$, eigenvalues describe the growth, decay, and oscillation of individual modes. Singular values describe size and scale in space: they measure stretching in the Euclidean norm. The largest stretching factor is

$$
\sigma_1(A)=\max_{\|x\|_2=1}\|Ax\|_2.
$$

The distinction matters when eigenvectors are not orthogonal: a combination of modes can grow much more than the eigenvalue magnitudes suggest.

## An Example with Nearly Parallel Eigenvectors

Let $0<\epsilon<1$ and define

$$
u=\begin{pmatrix}1\\\epsilon\end{pmatrix},
\qquad
v=\begin{pmatrix}1\\-\epsilon\end{pmatrix}.
$$

These vectors are linearly independent, but become nearly parallel as $\epsilon$ decreases. We construct a matrix with eigenvector $u$ for eigenvalue $1$ and eigenvector $v$ for eigenvalue $-1$. Set

$$
X=\begin{pmatrix}1&1\\\epsilon&-\epsilon\end{pmatrix},
\qquad
D=\begin{pmatrix}1&0\\0&-1\end{pmatrix}.
$$

Then

$$
X^{-1}=\frac12
\begin{pmatrix}1&1/\epsilon\\1&-1/\epsilon\end{pmatrix},
$$

and multiplication gives

$$
A=XDX^{-1}
=\begin{pmatrix}0&1/\epsilon\\\epsilon&0\end{pmatrix}.
$$

The eigenvalues remain $1$ and $-1$ for every $\epsilon$ in this range. Yet some vectors can be stretched by a very large factor.

### What the Singular Values Show

We have

$$
A^TA=\begin{pmatrix}\epsilon^2&0\\0&1/\epsilon^2\end{pmatrix}.
$$

Since $0<\epsilon<1$, the singular values in decreasing order are

$$
\sigma_1=\frac1\epsilon,
\qquad
\sigma_2=\epsilon.
$$

The corresponding directions are especially simple:

$$
Ae_2=\frac1\epsilon e_1,
\qquad
Ae_1=\epsilon e_2.
$$

Thus $e_2$ is stretched by $1/\epsilon$ and sent along $e_1$, while $e_1$ is contracted by $\epsilon$ and sent along $e_2$. In particular, $\|A\|_2=1/\epsilon$ even though both eigenvalues have magnitude one.

### Why the Eigenvector Expansion Allows This

The unit vector $e_2$ can be written as

$$
e_2=\frac{u-v}{2\epsilon}.
$$

The two terms on the right are large when $\epsilon$ is small, but their first components cancel. Applying $A$ changes the sign of the second eigenvector contribution:

$$
Ae_2=\frac{Au-Av}{2\epsilon}
=\frac{u+v}{2\epsilon}
=\frac1\epsilon e_1.
$$

The cancellation has disappeared. Neither eigenmode has increased in magnitude, but their sum has become much larger. Orthogonal eigenvectors cannot produce this effect, because the squared norm of their sum is the sum of the squared component magnitudes.

### What Happens over Several Steps

Direct multiplication gives $A^2=I$. Therefore,

$$
A^{2j}=I,
\qquad
A^{2j+1}=A,
\qquad j\geq0.
$$

Starting from $x_0=e_2$, the sequence alternates between $e_2$ and $e_1/\epsilon$. It is bounded for each fixed $\epsilon$, but its size can be much larger than that of the initial state. In norm,

$$
\|A^{2j}\|_2=1,
\qquad
\|A^{2j+1}\|_2=\frac1\epsilon.
$$

This also shows why singular values cannot generally be raised to powers to find the singular values of $A^k$: here $\sigma_1(A^2)=1$, whereas $\sigma_1(A)^2=1/\epsilon^2$.

## Large Amplification before Eventual Decay

A small modification gives an example in which every initial state eventually tends to zero. Let

$$
B=\alpha A,
\qquad 0<\alpha<1.
$$

Its eigenvalues are $\alpha$ and $-\alpha$, and

$$
\begin{aligned}
B^{2j}&=\alpha^{2j}I,\\
B^{2j+1}&=\alpha^{2j+1}A.
\end{aligned}
$$

Both expressions tend to zero as $j\to\infty$. Nevertheless, the first step amplifies $e_2$ by $\alpha/\epsilon$, which can be arbitrarily large. For example, with $\alpha=1/2$ and $\epsilon=1/100$,

$$
e_2\ \longmapsto\ 50e_1\ \longmapsto\ \frac14e_2
\ \longmapsto\ \frac{25}{2}e_1\ \longmapsto\ \frac1{16}e_2\ \longmapsto\cdots.
$$

This is **transient growth**: a state becomes larger before eventually decaying. The eigenvalues describe the eventual decay of the modes; the singular values of $B^k$ give the largest possible amplification at each step.

As discussed in [Eigendecomposition](eigendecomposition.md), for a diagonalizable matrix, eigenvalues of magnitude less than one imply decay of every state. They do not imply that the norm decreases at every step. Also, eigenvalues of magnitude exactly one do not by themselves guarantee bounded powers when the matrix is not diagonalizable.

## When Do Eigenvalue Magnitudes Equal Singular Values?

Recall that a complex square matrix is **normal** if $A^HA=AA^H$. For a real matrix this becomes $A^TA=AA^T$.

````{prf:theorem} Singular Values of a Normal Matrix
If $A$ is normal, its singular values are the magnitudes of its eigenvalues, counted with multiplicity and placed in decreasing order.
````

````{prf:proof} The Normal Case.
By the spectral theorem, $A=Q\Lambda Q^H$ with $Q$ unitary. Consequently,

$$
A^HA=Q\Lambda^H\Lambda Q^H.
$$

The diagonal entries of $\Lambda^H\Lambda$ are $|\lambda_i|^2$. Taking their nonnegative square roots gives the singular values.
````

For a Hermitian positive semidefinite matrix, the eigenvalues are already nonnegative, so the eigenvalues and singular values coincide.

A normal matrix admits an orthonormal eigenvector basis over $\mathbb{C}$. This is an existence statement: arbitrary eigenvectors chosen within a repeated eigenspace need not be orthogonal. For a real normal matrix, an eigenvector basis may require complex vectors.

For normal $A$, the same unitary basis diagonalizes every power. Hence, for integers $k\geq1$,

$$
\|A^k\|_2=\max_i|\lambda_i|^k.
$$

If all eigenvalues have magnitude at most one, no initial state can increase in Euclidean norm. In particular, the growth seen above cannot occur in this case.

Our example is not normal: for $0<\epsilon<1$, the following matrices differ.

$$
\begin{aligned}
A^TA&=\begin{pmatrix}\epsilon^2&0\\0&1/\epsilon^2\end{pmatrix},\\
AA^T&=\begin{pmatrix}1/\epsilon^2&0\\0&\epsilon^2\end{pmatrix},
\end{aligned}
$$

## Bounds That Hold for Every Square Matrix

Even without normality, every eigenvalue satisfies

$$
\sigma_n(A)\leq|\lambda|\leq\sigma_1(A).
$$

````{prf:proof} Bounds on Eigenvalue Magnitudes.
For a square matrix, the SVD gives

$$
\sigma_n(A)\|x\|_2\leq\|Ax\|_2
\leq\sigma_1(A)\|x\|_2.
$$

Apply this to a unit eigenvector. Since $\|Ax\|_2=|\lambda|$, the result follows. For a real matrix with a complex eigenvalue, use its eigenvector in $\mathbb{C}^n$; the same norm inequalities hold.
````

The **spectral radius** is $\rho(A)=\max_i|\lambda_i|$. The upper bound says $\rho(A)\leq\|A\|_2$, with equality for normal matrices. In our example the ratio $\|A\|_2/\rho(A)=1/\epsilon$ can be arbitrarily large.

The distinction also appears under a change of coordinates. Similarity preserves eigenvalues, but a general change of basis can change Euclidean lengths and therefore singular values: $D$ has singular values $1,1$, while $A=XDX^{-1}$ has singular values $1/\epsilon,\epsilon$. An orthogonal or unitary change of basis preserves lengths and preserves both sets of values.

Eigenvalues describe time evolution along individual modes; singular values describe size and scale in space. The geometry of the eigenvectors explains why eigenvalue magnitudes alone may not capture how much a matrix can stretch a vector.
