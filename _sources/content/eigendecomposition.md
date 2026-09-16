# Eigendecomposition

An eigendecomposition expresses a matrix in a basis of eigenvectors. In that basis, applying the matrix amounts to multiplying each coordinate by a scalar. The central question is whether there are enough independent eigenvectors to form a basis.

Throughout this section, $A$ is an $n\times n$ matrix with $n\ge1$. We allow complex eigenvalues and eigenvectors, even when $A$ has real entries.

## Eigenvalues, Eigenvectors, and Eigenspaces

A scalar $\lambda$ is an **eigenvalue** of $A$ if there is a nonzero vector $x$ such that

$$Ax=\lambda x.$$

The vector $x$ is a corresponding **eigenvector**. The requirement $x\ne0$ is essential: the zero vector satisfies this equation for every $\lambda$.

For real $x$ and $\lambda$, the line through $x$ is mapped into itself. A positive eigenvalue preserves the direction along this line, a negative one reverses it, and $\lambda=0$ sends the line to zero. Lengths along the line are multiplied by $|\lambda|$. Complex eigenvalues require a complex vector space; the picture of stretching a real line does not apply to them.

If $x$ is an eigenvector, so is every nonzero scalar multiple of $x$. More generally, the **eigenspace** associated with $\lambda$ is

$$E_\lambda=N(A-\lambda I).$$

It contains zero and all the eigenvectors for that eigenvalue. In particular, when zero is an eigenvalue, its eigenspace is $N(A)$.

For example,

$$
A=\begin{pmatrix}2&1\\0&1\end{pmatrix}
$$

has eigenvectors

$$
\begin{aligned}
x_1&=\begin{pmatrix}1\\0\end{pmatrix}, & Ax_1&=2x_1,\\
x_2&=\begin{pmatrix}-1\\1\end{pmatrix}, & Ax_2&=x_2.
\end{aligned}
$$

These two eigenvectors are independent, but they are not orthogonal. We will use them to diagonalize $A$ below.

## Existence of Eigenvalues

The eigenvalue equation is equivalent to $(A-\lambda I)x=0$ with $x\ne0$. Thus

$$\lambda\text{ is an eigenvalue}\quad\Longleftrightarrow\quad
\det(A-\lambda I)=0.$$

The **characteristic polynomial** is

$$p_A(t)=\det(A-tI).$$

It has degree $n$ and leading coefficient $(-1)^n$. Some texts use $\det(tI-A)$ instead; the roots are the same. For a $2\times2$ matrix,

$$
\begin{aligned}
p_A(t)
&=\det\begin{pmatrix}a-t&b\\c&d-t\end{pmatrix}\\
&=(a-t)(d-t)-bc\\
&=t^2-(a+d)t+(ad-bc).
\end{aligned}
$$

````{prf:theorem} Existence of Eigenvalues
:label: thm:eigenvalue_existence
Every $n\times n$ complex matrix, with $n\ge1$, has $n$ eigenvalues in $\mathbb{C}$ when roots are counted with their multiplicities.
````

````{prf:proof} Using the Characteristic Polynomial.
By the Fundamental Theorem of Algebra, the degree-$n$ polynomial $p_A$ factors into $n$ linear factors over $\mathbb{C}$. Each root $\lambda$ makes $A-\lambda I$ singular, so $N(A-\lambda I)$ contains a nonzero vector. That vector is an eigenvector for $\lambda$.
````

Counting roots with multiplicity does **not** guarantee $n$ independent eigenvectors.

The existence of at least one eigenvalue can also be proved without determinants.

````{prf:proof} Using Linear Dependence.
Choose any nonzero $x\in\mathbb{C}^n$. The $n+1$ vectors

$$x,Ax,A^2x,\ldots,A^nx$$

are linearly dependent. Hence there are scalars $c_0,\ldots,c_n$, not all zero, such that

$$c_0x+c_1Ax+\cdots+c_nA^nx=0.$$

Define $P(t)=c_0+c_1t+\cdots+c_nt^n$. Then $P(A)x=0$, where the constant term in $P(A)$ is $c_0I$. This nonzero polynomial cannot be constant, since $c_0x=0$ with $x\ne0$ would force $c_0=0$.

Let $k\ge1$ be its degree, so $c_k\ne0$. By the Fundamental Theorem of Algebra,

$$P(t)=c_k\prod_{j=1}^k(t-\lambda_j).$$

Substituting $A$ and dividing by $c_k$ gives

$$(A-\lambda_1I)\cdots(A-\lambda_kI)x=0.$$

Apply these factors from right to left. More explicitly, set $v_k=x$ and define

$$v_{j-1}=(A-\lambda_jI)v_j,
\qquad j=k,k-1,\ldots,1.$$

We start with $v_k\ne0$ and end with $v_0=0$. At the first step that produces zero, the input $v_j$ is nonzero and satisfies

$$(A-\lambda_jI)v_j=0.$$

Thus $Av_j=\lambda_jv_j$: we have found an eigenvalue and a corresponding eigenvector.
````

### Real Matrices Can Have Complex Eigenvalues

The real rotation matrix

$$R=\begin{pmatrix}0&-1\\1&0\end{pmatrix}$$

has characteristic polynomial $t^2+1$. Its eigenvalues are $i$ and $-i$, with eigenvectors $(1,-i)^T$ and $(1,i)^T$, respectively. It has no real eigenvectors: a rotation by $90^\circ$ preserves no real line.

For any real matrix, nonreal eigenvalues occur in conjugate pairs. Indeed, taking complex conjugates of $Ax=\lambda x$ gives

$$A\overline{x}=\overline{\lambda}\,\overline{x}.$$

The corresponding roots have the same multiplicity because the characteristic polynomial has real coefficients.

## Diagonalization: A Change to an Eigenvector Basis

Suppose $A$ has $n$ independent eigenvectors $x_1,\ldots,x_n$, with corresponding eigenvalues $\lambda_1,\ldots,\lambda_n$. Form

$$
X=[x_1,\ldots,x_n],\qquad
\Lambda=\operatorname{diag}(\lambda_1,\ldots,\lambda_n).
$$

The eigenvalue equations combine into

$$AX=X\Lambda.$$

Since the columns of $X$ are independent, $X$ is invertible. Therefore,

$$\boxed{A=X\Lambda X^{-1}.}$$

This is an **eigendecomposition**, or **diagonalization**, of $A$.

Conversely, if $A=X\Lambda X^{-1}$ with $X$ invertible and $\Lambda$ diagonal, then $AX=X\Lambda$ shows that each column of $X$ is an eigenvector. Hence:

> A matrix is diagonalizable exactly when it has a basis of eigenvectors.

Over $\mathbb{R}$, this requires a basis of real eigenvectors and real eigenvalues. A real matrix can be diagonalizable over $\mathbb{C}$ without being diagonalizable over $\mathbb{R}$, as the rotation example illustrates.

### What the Three Factors Do

For any vector $v$, write

$$v=Xc=\sum_{j=1}^n c_jx_j.$$

The factors in $X\Lambda X^{-1}$ act from right to left:

1. $X^{-1}$ finds the coordinates $c$ in the eigenvector basis.
2. $\Lambda$ multiplies coordinate $c_j$ by $\lambda_j$.
3. $X$ converts back to the original coordinates.

Thus

$$Av=\sum_{j=1}^n\lambda_jc_jx_j.$$

The eigenvectors need not be orthogonal, so in general $X^{-1}\ne X^H$. Normalizing the columns of $X$ does not make them orthogonal.

### A Worked Diagonalization

For the earlier example, take

$$
A=\begin{pmatrix}2&1\\0&1\end{pmatrix},\qquad
X=\begin{pmatrix}1&-1\\0&1\end{pmatrix}.
$$

Then

$$
\Lambda=\begin{pmatrix}2&0\\0&1\end{pmatrix},\qquad
X^{-1}=\begin{pmatrix}1&1\\0&1\end{pmatrix},
$$

and direct multiplication verifies $A=X\Lambda X^{-1}$.

For $v=(0,1)^T$, the eigenvector coordinates are $c=X^{-1}v=(1,1)^T$. Hence $v=x_1+x_2$ and

$$Av=2x_1+x_2=\begin{pmatrix}1\\1\end{pmatrix}.$$

The decomposition is not unique. Eigenvectors can be rescaled, the eigenpairs can be reordered together, and any basis can be chosen within a repeated eigenvalue's eigenspace.

## When Are There Enough Eigenvectors?

### Distinct Eigenvalues Give Independent Eigenvectors

More generally, eigenspaces for distinct eigenvalues form a direct sum.

````{prf:proof}
Let $\lambda_1,\ldots,\lambda_s$ be distinct, and let $v_j\in E_{\lambda_j}$ satisfy

$$v_1+\cdots+v_s=0.$$

We prove by induction on $s$ that every $v_j$ is zero. The case $s=1$ is immediate. Apply $A-\lambda_sI$ to the relation:

$$\sum_{j=1}^{s-1}(\lambda_j-\lambda_s)v_j=0.$$

By the induction hypothesis, every term is zero. Since $\lambda_j\ne\lambda_s$, this gives $v_1=\cdots=v_{s-1}=0$, and the original relation then gives $v_s=0$.
````

In particular, a matrix with $n$ distinct eigenvalues is diagonalizable. Distinct eigenvalues are sufficient, but not necessary: the identity matrix has just one eigenvalue, yet every nonzero vector is an eigenvector.

### Repeated Eigenvalues: Two Different Multiplicities

For an eigenvalue $\lambda$:

- Its **algebraic multiplicity** $a_\lambda$ is its multiplicity as a root of $p_A$.
- Its **geometric multiplicity** $g_\lambda=\dim E_\lambda$ counts the independent eigenvectors available for that eigenvalue.

These satisfy

$$1\le g_\lambda\le a_\lambda.$$

````{prf:proof}
Let $g=g_\lambda$. If $g=n$, then $A=\lambda I$ and the result follows directly. Otherwise, extend a basis of $E_\lambda$ to a basis of the whole space, and put these basis vectors in the columns of an invertible matrix $P$. In this basis,

$$P^{-1}AP=\begin{pmatrix}\lambda I_g&B\\0&C\end{pmatrix}.$$

A change of basis preserves the characteristic polynomial, because

$$
\begin{aligned}
\det(P^{-1}AP-tI)
&=\det\bigl(P^{-1}(A-tI)P\bigr)\\
&=\det(A-tI).
\end{aligned}
$$

The block triangular determinant formula therefore gives

$$p_A(t)=(\lambda-t)^g\det(C-tI).$$

Thus $\lambda$ occurs at least $g$ times as a root.
````

Because eigenspaces for distinct eigenvalues form a direct sum, their bases together provide $\sum_\lambda g_\lambda$ independent eigenvectors. The algebraic multiplicities sum to $n$. Therefore, over $\mathbb{C}$,

$$A\text{ is diagonalizable}
\quad\Longleftrightarrow\quad
 g_\lambda=a_\lambda\text{ for every eigenvalue }\lambda.$$

### A Matrix That Cannot Be Diagonalized

Consider

$$J=\begin{pmatrix}1&1\\0&1\end{pmatrix}.$$

Its characteristic polynomial is $(1-t)^2$, so the eigenvalue $1$ has algebraic multiplicity two. But

$$
(J-I)x=0
\quad\Longleftrightarrow\quad x_2=0.
$$

Thus $E_1=\operatorname{span}\{(1,0)^T\}$ has dimension one. There is no eigenvector basis, so $J$ is not diagonalizable. Such a matrix is called **defective**.

Compare $J$ with $I_2$: both have eigenvalues $1,1$, but only $I_2$ is diagonalizable. Also, $J$ is invertible. Invertibility and diagonalizability are different properties.

## Matrix Powers and Repeated Application

Eigenvalues are useful in science and engineering because they describe how a system evolves in time. A vector $v_k$ can record its state at time step $k$: temperatures, displacements, or concentrations, for example. A linear model advances that state according to

$$v_{k+1}=Av_k.$$

Applying the same rule repeatedly gives $v_k=A^kv_0$. Thus understanding powers of $A$ tells us how the initial state changes over time.

Linear models also help us understand nonlinear systems. Near an equilibrium, the evolution of small deviations can often be approximated by a linear equation. In that setting, $v_k$ represents the deviation from the equilibrium, and the approximation describes local behavior.

### Evolution in an Eigenvector Basis

If $A=X\Lambda X^{-1}$, adjacent factors $X^{-1}X$ cancel when we multiply copies of $A$. For example,

$$
\begin{aligned}
A^2&=X\Lambda(X^{-1}X)\Lambda X^{-1}\\
&=X\Lambda^2X^{-1}.
\end{aligned}
$$

Continuing in the same way gives, for integers $k\ge1$,

$$A^k=X\Lambda^kX^{-1},\qquad
\Lambda^k=\operatorname{diag}(\lambda_1^k,\ldots,\lambda_n^k).$$

Write the initial state as $v_0=\sum_jc_jx_j$. Then

$$v_k=A^kv_0=\sum_{j=1}^n c_j\lambda_j^kx_j.$$

Each eigenvector describes a pattern in the state, often called a **mode**. Its eigenvalue tells us how that mode changes during one time step: its coefficient is multiplied by $\lambda_j$. The initial state determines which modes are present through the coefficients $c_j$.

### Example: Two Bodies Exchanging Heat

Consider two bodies with equal heat capacities that exchange heat with each other but not with their surroundings. In a simple discrete-time approximation, each temperature moves $10\%$ of the way toward the other body's previous temperature during each step:

$$
\begin{pmatrix}T_{1,k+1}\\T_{2,k+1}\end{pmatrix}
=\underbrace{\begin{pmatrix}0.9&0.1\\0.1&0.9\end{pmatrix}}_{A}
\begin{pmatrix}T_{1,k}\\T_{2,k}\end{pmatrix}.
$$

The two eigenvectors and eigenvalues are

$$
\begin{aligned}
x_1&=\begin{pmatrix}1\\1\end{pmatrix}, & \lambda_1&=1,\\
x_2&=\begin{pmatrix}1\\-1\end{pmatrix}, & \lambda_2&=0.8.
\end{aligned}
$$

The first mode gives both bodies the same temperature. Its eigenvalue is one, so this common temperature remains unchanged. The second mode raises one temperature and lowers the other by the same amount. Its eigenvalue is $0.8$, so the temperature difference decreases by $20\%$ per step.

Starting from $80$ and $20$ degrees,

$$
v_0=\begin{pmatrix}80\\20\end{pmatrix}=50x_1+30x_2.
$$

After $k$ steps,

$$
v_k=50x_1+30(0.8)^kx_2
=\begin{pmatrix}50+30(0.8)^k\\50-30(0.8)^k\end{pmatrix}.
$$

Both temperatures approach $50$ degrees. The eigenvalues explain both the preserved average and the rate at which the temperatures equalize.

### Growth, Decay, and Oscillation

For a diagonalizable matrix, the expansion above gives:

- If every $|\lambda_j|<1$, then $v_k\to0$ for every initial vector.
- If some $|\lambda_j|>1$, initial vectors with a nonzero component in a corresponding eigendirection grow without bound. Other initial vectors may still decay.
- If every $|\lambda_j|\le1$, the sequence is bounded, but components with $|\lambda_j|=1$ need not decay or converge.

A negative real eigenvalue makes its component alternate sign. For a complex eigenvalue $\lambda=re^{i\theta}$, the factor $\lambda^k=r^ke^{ik\theta}$ combines growth or decay with oscillation.

These conclusions follow here from the assumption that $A$ is diagonalizable. Without it, repeated eigenvalues can produce additional behavior. For the defective matrix above,

$$J^k=\begin{pmatrix}1&k\\0&1\end{pmatrix}.$$

Both eigenvalues equal one, yet $J^k(0,1)^T=(k,1)^T$ grows. Eigenvalues alone do not describe every aspect of repeated application.

### Continuous-Time Systems

The same idea applies to a linear differential equation

$$\frac{dv}{dt}=Av.$$

The matrix $A$ in this equation maps the state to its rate of change.

If $A=X\Lambda X^{-1}$ and $v(t)=Xc(t)$, then $c'(t)=\Lambda c(t)$. Each coordinate therefore satisfies $c_j'(t)=\lambda_jc_j(t)$, giving

$$v(t)=\sum_{j=1}^n c_j(0)e^{\lambda_jt}x_j.$$

In continuous time, the real part of $\lambda_j$ determines exponential growth or decay, and its imaginary part determines oscillation. In discrete time, the magnitude $|\lambda_j|$ determines growth or decay per step. In both cases, an eigenvector basis separates the evolution into individual modes.

## The Complex Schur Decomposition

Diagonalization requires an eigenvector basis. Every complex square matrix does, however, admit a triangular representation in an orthonormal basis:

$$\boxed{A=QTQ^H,}$$

where $Q$ is unitary and $T$ is upper triangular. This is the **complex Schur decomposition**.

````{prf:proof} Existence of the Complex Schur Decomposition.
We use induction on $n$. For $n=1$, take $Q=[1]$ and $T=A$.

Choose a unit eigenvector $q_1$ with eigenvalue $\lambda_1$ and extend it to an orthonormal basis. Let $U$ contain these basis vectors as columns. Since $Aq_1=\lambda_1q_1$,

$$
U^HAU=\begin{pmatrix}\lambda_1&w^H\\0&A_2\end{pmatrix}.
$$

By induction, $A_2=VT_2V^H$ for some unitary $V$ and upper triangular $T_2$. Set

$$Q=U\begin{pmatrix}1&0\\0&V\end{pmatrix}.$$

Then $Q$ is unitary, and block multiplication gives

$$Q^HAQ=\begin{pmatrix}\lambda_1&w^HV\\0&T_2\end{pmatrix}.$$

This matrix is upper triangular, completing the induction.
````

### What Schur Form Tells Us

Similarity preserves the characteristic polynomial, and the determinant of a triangular matrix is the product of its diagonal entries. Hence

$$p_A(t)=\det(T-tI)=\prod_{j=1}^n(t_{jj}-t).$$

The diagonal entries of $T$ are therefore the eigenvalues of $A$, counted with algebraic multiplicity. This also proves, even when $A$ is defective, that

$$\operatorname{tr}(A)=\sum_{j=1}^n\lambda_j,
\qquad \det(A)=\prod_{j=1}^n\lambda_j.$$

In particular, $A$ is invertible exactly when zero is not an eigenvalue.

The columns of $Q$ are called **Schur vectors**. They need not be eigenvectors. From $AQ=QT$,

$$Aq_j=\sum_{i=1}^j t_{ij}q_i.$$

Thus the first column is an eigenvector, and the span of the first $j$ columns is an **invariant subspace**: applying $A$ to a vector in that subspace keeps it there.

An eigendecomposition permits a general invertible basis matrix $X$. Schur form requires an orthonormal basis and permits nonzero entries above the diagonal of $T$. Even a diagonalizable matrix can have a nondiagonal Schur form: our example $\begin{pmatrix}2&1\\0&1\end{pmatrix}$ is already in Schur form with $Q=I$.

The stronger property of having an orthonormal eigenvector basis is called **unitary diagonalization**. The next section on [normal matrices](normal_matrices.md) characterizes exactly when this is possible.

## The Real Schur Decomposition

Every real square matrix has a decomposition

$$A=QSQ^T,$$

where $Q$ is real orthogonal and $S$ is real block upper triangular, with diagonal blocks of size $1\times1$ or $2\times2$. The $1\times1$ blocks are real eigenvalues; each $2\times2$ block has a nonreal conjugate pair as its eigenvalues. This is the **real Schur decomposition**.

For example, a $3\times3$ real Schur form with one real eigenvalue and one nonreal pair has the structure

$$
S=\begin{pmatrix}
\lambda&*&*\\
0&b_{11}&b_{12}\\
0&b_{21}&b_{22}
\end{pmatrix},
$$

where the stars denote unrestricted entries. The nonreal eigenvalues belong to the entire bottom-right block; they cannot be read from its diagonal entries alone.

````{prf:proof} Existence of the Real Schur Decomposition.
We again use induction, with the scalar case immediate.

If $A$ has a real eigenvalue, choose a real unit eigenvector. The complex Schur argument can then be carried out with real orthogonal matrices, leaving a smaller real matrix to which induction applies.

Otherwise, choose a nonreal eigenvalue $\lambda=a+ib$, with $b\ne0$, and an eigenvector $z=u+iv$, where $u,v$ are real. Separating the real and imaginary parts of $Az=\lambda z$ gives

$$
\begin{aligned}
Au&=au-bv,\\
Av&=bu+av.
\end{aligned}
$$

The vectors $u,v$ are independent. Otherwise, $z$ would be a nonzero complex multiple of a real vector, and the eigenvalue equation would force $\lambda$ to be real.

Therefore $W=\operatorname{span}\{u,v\}$ is a two-dimensional real invariant subspace. Choose an orthonormal basis for $W$ and extend it to an orthonormal basis of $\mathbb{R}^n$. With these vectors in the columns of $U$,

$$U^TAU=\begin{pmatrix}B&C\\0&D\end{pmatrix}.$$

The matrix $B$ represents the restriction of $A$ to $W$. In the basis $u,v$, that restriction has matrix

$$\begin{pmatrix}a&b\\-b&a\end{pmatrix},$$

so $B$ has eigenvalues $a\pm ib$. The change to an orthonormal basis can change the entries of this $2\times2$ matrix.

If $n=2$, the proof is complete. Otherwise, apply induction to $D=VS_2V^T$ and set $Q=U\operatorname{diag}(I_2,V)$. Then

$$Q^TAQ=\begin{pmatrix}B&CV\\0&S_2\end{pmatrix},$$

which has the required block upper triangular form.
````

A real $2\times2$ Schur block need not have off-diagonal entries of equal magnitude and opposite sign. For example,

$$B=\begin{pmatrix}0&-2\\1&0\end{pmatrix}$$

is a valid block with eigenvalues $\pm i\sqrt{2}$.

The distinction to retain is that **Schur form always exists, whereas an eigenvector basis may not**. Numerical methods for computing these decompositions will be developed in the later chapter on [eigenvalue computations](eigenvalues.md).
