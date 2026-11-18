# The Conjugate Gradient (CG) Method

The conjugate gradient method chooses an approximation from a growing Krylov subspace. Its main idea is to change the inner product used to define the best approximation. This makes the coefficients computable without knowing the exact solution. Symmetry and the Krylov structure then reduce the construction to a short recurrence.

Throughout this section, $A\in\mathbb R^{n\times n}$ is symmetric positive definite (SPD), and $x=A^{-1}b$ is the exact solution. As in the chapter introduction, $x_k$ is the approximation after $k$ steps. We write $e_k=x-x_k$ for its error and $r_k=b-Ax_k=Ae_k$ for its residual. We first take $x_0=0$, so $r_0=b$. If $b=0$, the initial guess is already the solution. All orthogonality and optimality statements below concern exact arithmetic.

## Nested Krylov Spaces and Search Directions

Starting from $b\ne0$, define

$$
\mathcal K_0=\{0\},
\qquad
\mathcal K_k=\operatorname{span}\{b,Ab,\ldots,A^{k-1}b\}.
$$

These spaces are nested:

$$
\mathcal K_0\subseteq\mathcal K_1\subseteq\mathcal K_2\subseteq\cdots.
$$

Each multiplication can add at most one new independent direction. While the spaces are growing, $\dim\mathcal K_k=k$. We seek directions $p_1,p_2,\ldots$ that follow this growth:

$$
\operatorname{span}\{p_1,\ldots,p_k\}=\mathcal K_k.
$$

With such a basis, an approximation in $\mathcal K_k$ has the form

$$
x_k=\sum_{j=1}^k\mu_jp_j.
$$

If the optimal coefficients of the old directions remain unchanged when we enlarge the space, only one new term is needed:

$$
x_k=x_{k-1}+\mu_kp_k.
$$

This is what an orthogonal basis can provide. We need to choose both the geometry that defines orthogonality and a normalization for each direction.

## A First Attempt: Minimize the Euclidean Error

A natural goal is to find the vector in $\mathcal K_k$ that minimizes $\|x-x_k\|_2$. Choose mutually orthogonal directions adapted to the nested spaces. The new direction must satisfy

$$
p_k\in\mathcal K_k,
\qquad p_k\perp\mathcal K_{k-1}.
$$

When $\dim\mathcal K_k=k$, these conditions specify a unique direction up to a nonzero scalar: the orthogonal complement of $\mathcal K_{k-1}$ inside $\mathcal K_k$ is one-dimensional. If we give each direction unit Euclidean norm, these are precisely the [Lanczos vectors](lanczos_iteration.md), up to their signs:

$$
p_k=q_k,
\qquad Q_k=[q_1,\ldots,q_k],
\qquad Q_k^TQ_k=I.
$$

Suppose $x$ first belongs to $\mathcal K_m$, where $m\le n$. The orthonormal basis of this space gives the expansion

$$
x=\sum_{j=1}^m\mu_jp_j,
\qquad \mu_j=p_j^Tx.
$$

The Euclidean best approximation in $\mathcal K_k$ is the partial sum through $k$. Orthogonality ensures that adding a new direction does not change the previous coefficients. The finite-convergence theorem below shows that some $m\le n$ exists.

The figure illustrates this process in two dimensions. The first approximation is the closest point to $x$ on the line $\mathcal K_1$. Adding the perpendicular direction $p_2$ captures the remaining component without changing the coefficient along $p_1$.

:::{figure} ../_static/cg_euclidean_projection.svg
:alt: Left, the exact solution projects perpendicularly onto the first Krylov line, giving the Euclidean best approximation x1. Right, adding the second orthogonal direction reaches the solution without changing the first component.
:width: 100%
:name: fig-cg-euclidean-projection

Euclidean projection for $A=\begin{pmatrix}2&-1\\-1&2\end{pmatrix}$ and $b=(3,0)^T$. Here $p_1=(1,0)^T$, $p_2=(0,1)^T$, and $x=(2,1)^T$. In the left panel, the circle centered at $x$ touches $\mathcal K_1$ at the nearest point $x_1=2p_1$, where the error is perpendicular to the line. In the right panel, $\mathcal K_2$ is the whole plane, so adding $p_2$ gives $x_2=2p_1+p_2=x$. Drawing these Euclidean projections requires the exact solution, which is unavailable in the actual solve.
:::

The difficulty is the coefficient formula: $p_j^Tx$ involves the unknown solution. Lanczos gives us the directions, but does not directly give these projection coefficients. We need an inner product that lets us use the known equation $Ax=b$.

## The Key Change: Use the A-Inner Product

Compare the two quantities

$$
p_j^Tx
\qquad\text{and}\qquad
p_j^TAx=p_j^Tb.
$$

The second can be evaluated from known data. This suggests replacing Euclidean orthogonality by orthogonality in the geometry defined by $A$.

:::{prf:definition} The A-Inner Product
:label: a_inner_prod

For an SPD matrix $A$, define

$$
\langle y,z\rangle_A=y^TAz,
\qquad
\|z\|_A=\sqrt{z^TAz}.
$$

Two vectors are **A-orthogonal**, or **conjugate**, if $y^TAz=0$.
:::

The expression $y^TAz$ is linear in each vector. Symmetry of $A$ makes it symmetric in $y$ and $z$, and positive definiteness gives $z^TAz>0$ for every nonzero $z$. Together, these properties make it a valid inner product. A general nonsymmetric or indefinite matrix does not provide this geometry.

This is the [energy norm](gauss_seidel_method.md) introduced for Gauss-Seidel. To see its geometry, write $A=V\Lambda V^T$, with $V$ orthogonal and $\Lambda$ positive diagonal. Then

$$
\begin{aligned}
\langle y,z\rangle_A
&=(\Lambda^{1/2}V^Ty)^T(\Lambda^{1/2}V^Tz)\\
&=(A^{1/2}y)^T(A^{1/2}z),\\
\|z\|_A&=\|A^{1/2}z\|_2.
\end{aligned}
$$

Thus the new inner product is an ordinary dot product after rotating into the eigenvector coordinates and scaling by $\sqrt{\lambda_i}$. The symmetric square root $A^{1/2}=V\Lambda^{1/2}V^T$ includes a final orthogonal change of coordinates, which preserves dot products. No square root is needed in the algorithm.

We can interpret this scaling as a **change of units** along the eigenvector directions. In those coordinates, the A-inner product weights the product of the $i$th components by $\lambda_i$. Equivalently, multiply each vector's $i$th component by $\sqrt{\lambda_i}$, then take the ordinary dot product. This is an ordinary dot product after changing the coordinate scales.

:::{figure} ../_static/cg_inner_product.svg
:alt: Two vectors before and after doubling their first coordinates. They are not perpendicular in the original Euclidean coordinates, but become perpendicular after rescaling, illustrating their A-orthogonality for A diagonal with entries 4 and 1.
:width: 100%
:name: fig-cg-inner-product

The same vectors measured in two coordinate systems. For $A=\operatorname{diag}(4,1)$, the axes already follow the eigenvector directions. Rescaling the first coordinate by $\sqrt4=2$ turns squares of the gray grid into rectangles. The transformed vectors have ordinary dot product zero, so the original vectors are A-orthogonal. A change of units preserves the weighted inner product through this transformation, but can change Euclidean angles: A-orthogonality need not look like a right angle in the original coordinates.
:::

### Computable Expansion Coefficients

From this point on, $p_j$ denotes an **A-orthogonal search direction**, replacing the Euclidean directions of the first attempt. We still require the first $k$ directions to span $\mathcal K_k$:

$$
p_i^TAp_j=0\quad(i\ne j),
\qquad
\operatorname{span}\{p_1,\ldots,p_k\}=\mathcal K_k.
$$

As before, the new direction is unique up to scaling: it lies in $\mathcal K_k$ and is A-orthogonal to $\mathcal K_{k-1}$. Taking the A-inner product of the full expansion with $p_i$ gives

$$
\begin{aligned}
p_i^Tb=p_i^TAx
&=\sum_{j=1}^m\mu_jp_i^TAp_j\\
&=\mu_i p_i^TAp_i.
\end{aligned}
$$

Consequently,

$$
\boxed{\mu_i=\frac{p_i^Tb}{p_i^TAp_i}.}
$$

The denominator is positive. Both numerator and denominator are computable once $p_i$ is known. Changing the inner product has removed the unknown solution from the coefficient calculation.

In matrix notation, define $P_k=[p_1,\ldots,p_k]$ and $d_j=p_j^TAp_j$. Then

$$
P_k^TAP_k=D_k=\operatorname{diag}(d_1,\ldots,d_k),
\qquad
D_k\begin{pmatrix}\mu_1\\\vdots\\\mu_k\end{pmatrix}=P_k^Tb.
$$

We leave the scale of each $p_j$ free. When we construct the directions from the residuals, we will choose a normalization that gives a short recurrence.

## Why the Partial Sum Is Optimal

The computable coefficients have a geometric meaning: they give the A-orthogonal projection of $x$ onto $\mathcal K_k$.

:::{prf:theorem} CG Optimality in the Energy Norm
:label: cg_energy_optimality

Let $p_1,\ldots,p_k$ be an A-orthogonal basis of $\mathcal K_k$, and set

$$
x_k=\sum_{j=1}^k\frac{p_j^Tb}{p_j^TAp_j}\,p_j.
$$

Then $x_k$ is the unique minimizer of $\|x-z\|_A$ over $z\in\mathcal K_k$. Equivalently,

$$
p_j^TA(x-x_k)=0\qquad(j\le k).
$$
:::

:::{prf:proof}
For each $j\le k$, A-orthogonality gives

$$
p_j^TA(x-x_k)
=p_j^Tb-\mu_jp_j^TAp_j=0.
$$

Any other candidate in $\mathcal K_k$ is $x_k+v$, with $v\in\mathcal K_k$. Expanding its squared error yields

$$
\begin{aligned}
\|x-(x_k+v)\|_A^2
&=\|x-x_k\|_A^2-2v^TA(x-x_k)+\|v\|_A^2\\
&=\|x-x_k\|_A^2+\|v\|_A^2.
\end{aligned}
$$

The last term is positive unless $v=0$, proving the unique minimum.
:::

Enlarging the space adds one A-orthogonal direction, so the old coefficients stay fixed. The update $x_k=x_{k-1}+\mu_kp_k$ is therefore also a minimization over the whole current space, not just along the latest line.

The same projection identity shows how the error decreases. Since
$x-x_{k-1}=(x-x_k)+\mu_kp_k$ and $p_k^TA(x-x_k)=0$,

$$
\|x-x_{k-1}\|_A^2
=\|x-x_k\|_A^2+\mu_k^2p_k^TAp_k.
$$

Each nonzero step removes an A-orthogonal component of the error.

As in Gauss-Seidel, define the quadratic energy

$$
\Phi(z)=\frac12z^TAz-b^Tz.
$$

Since $Ax=b$,

$$
\Phi(z)-\Phi(x)=\frac12\|z-x\|_A^2.
$$

Thus CG also minimizes this computable quadratic over $\mathcal K_k$. It differs from Gauss-Seidel by choosing its directions from the Krylov spaces rather than the coordinate axes.

## Use the Residual to Find the Next Direction

The coefficient formula is useful only if we can construct the A-orthogonal directions from known quantities. The residual provides the new Krylov direction and leads directly to an efficient recurrence.

One possibility is to start with $p_1=b$, generate a new vector $Ap_k$, and subtract its A-projections onto the previous directions. For $i\le k-2$, symmetry gives

$$
p_i^TA(Ap_k)=(Ap_i)^TAp_k=0,
$$

because $Ap_i\in\mathcal K_{i+1}\subseteq\mathcal K_{k-1}$ and $p_k$ is A-orthogonal to that space. Thus only the projections onto $p_{k-1}$ and $p_k$ remain. This gives a three-term recurrence based on $Ap_k$.

Using the **residual** is more efficient because only the latest search direction is needed, and the coefficients reduce to ratios of residual norms. Recall that

$$
r_k=b-Ax_k.
$$

Because $x_k\in\mathcal K_k$, multiplication by $A$ gives $Ax_k\in\mathcal K_{k+1}$. Also $b\in\mathcal K_{k+1}$, so

$$
\boxed{r_k\in\mathcal K_{k+1}.}
$$

The optimality condition tells us more.

:::{prf:theorem} Orthogonality of Residuals
:label: cg_residual_orthogonality

The residual of the energy-minimizing approximation satisfies

$$
r_k\perp\mathcal K_k
$$

in the ordinary Euclidean inner product.
:::

:::{prf:proof}
For any $v\in\mathcal K_k$, the projection property gives

$$
v^Tr_k=v^T(b-Ax_k)=v^TA(x-x_k)=0.
$$
:::

Notice the two different orthogonalities: the error is **A-orthogonal** to $\mathcal K_k$, while the residual is **Euclidean-orthogonal** to it.

If $r_k\ne0$, it belongs to $\mathcal K_{k+1}$ but cannot belong to $\mathcal K_k$, since it is orthogonal to that space. It therefore supplies precisely the new independent direction we need. If $r_k=0$, then $Ax_k=b$ and we stop.

### Constructing an A-Orthogonal Direction from the Residual

Suppose $p_1,\ldots,p_k$ are already known. There is one new A-orthogonal direction in $\mathcal K_{k+1}$. Call it $p_{k+1}$ and expand the nonzero residual in the enlarged basis:

$$
r_k=\sum_{\ell=1}^{k+1}u_{\ell,k+1}p_\ell.
$$

For each known direction $p_i$, taking its A-inner product with this equation isolates one coefficient:

$$
p_i^TAr_k=u_{i,k+1}p_i^TAp_i,
\qquad
u_{i,k+1}=\frac{p_i^TAr_k}{d_i}
\quad(i\le k).
$$

The coefficient of $p_{k+1}$ is nonzero because $r_k\notin\mathcal K_k$. We choose the scale of $p_{k+1}$ so that this coefficient is $1$. Rearranging then gives the construction

$$
\boxed{
p_{k+1}=r_k-\sum_{i=1}^k\frac{p_i^TAr_k}{d_i}\,p_i.
}
$$

This is A-orthogonalization of the residual. It produces a nonzero vector in $\mathcal K_{k+1}$ that is A-orthogonal to every previous direction. Starting with $p_1=r_0=b$, it constructs the entire sequence. The sum appears to require all previous directions, but symmetry and the nesting of the Krylov spaces reduce it to one term.

## Only the Latest Direction Is Needed

Most coefficients in the preceding sum vanish. To see why, use symmetry:

$$
p_i^TAr_k=(Ap_i)^Tr_k.
$$

For $i<k$, the nesting of the Krylov spaces gives

$$
Ap_i\in\mathcal K_{i+1}\subseteq\mathcal K_k.
$$

Since $r_k\perp\mathcal K_k$, it follows that $(Ap_i)^Tr_k=0$ for all $i<k$. Only the coefficient of $p_k$ remains:

$$
r_k=u_{k,k+1}p_k+p_{k+1},
\qquad
u_{k,k+1}=\frac{p_k^TAr_k}{d_k}.
$$

Thus

$$
\boxed{p_{k+1}=r_k-u_{k,k+1}p_k.}
$$

Although we required A-orthogonality to *all* previous directions, we need to subtract only the projection onto the latest one. The Krylov structure and symmetry make all the earlier projections zero automatically. This fixes the normalization: each new direction contains the current residual with coefficient $1$. Its A-norm is generally not $1$.

## The Residual Basis and Its Orthogonality

The same subspace argument also identifies the residual sequence.

:::{prf:theorem} Orthogonality of Residuals to Each Other
:label: cg_residuals_orthogonal

Before termination, the nonzero residuals are mutually Euclidean-orthogonal:

$$
r_i^Tr_j=0\qquad(i\ne j).
$$
:::

:::{prf:proof}
For $i<j$, we have $r_i\in\mathcal K_{i+1}\subseteq\mathcal K_j$, while $r_j\perp\mathcal K_j$. Hence $r_i^Tr_j=0$.
:::

:::{prf:theorem} Residuals as a Basis
:label: cg_residual_basis

If $r_0,\ldots,r_{k-1}$ are nonzero, then

$$
\operatorname{span}\{r_0,\ldots,r_{k-1}\}
=\mathcal K_k
=\operatorname{span}\{p_1,\ldots,p_k\}.
$$
:::

:::{prf:proof}
Each $r_i$ belongs to $\mathcal K_{i+1}\subseteq\mathcal K_k$. These $k$ residuals are nonzero and mutually orthogonal, so they are linearly independent. Since $\mathcal K_k$ is generated by $k$ vectors, its dimension is at most $k$. The $k$ independent residuals therefore form a basis of $\mathcal K_k$.
:::

The residuals form a Euclidean-orthogonal basis, while the search directions form an A-orthogonal basis. Both describe the same nested Krylov spaces. This distinction explains how the residual can generate a new search direction without itself being A-orthogonal to the previous directions.

## Deriving the Efficient Coefficients

The direction update now involves only $r_k$ and $p_k$. Its two scalar coefficients can also be written using quantities already computed during the iteration.

(cg-step-length-derivation)=
### The Step Length

From the expansion formula, $\mu_k=p_k^Tb/d_k$, where $d_k=p_k^TAp_k$. Since $b=r_{k-1}+Ax_{k-1}$ and $p_k$ is A-orthogonal to $\mathcal K_{k-1}$,

$$
\begin{aligned}
p_k^Tb
&=p_k^Tr_{k-1}+p_k^TAx_{k-1}\\
&=p_k^Tr_{k-1}.
\end{aligned}
$$

For $k>1$, the short direction formula gives $p_k=r_{k-1}-u_{k-1,k}p_{k-1}$. The residual is orthogonal to $p_{k-1}\in\mathcal K_{k-1}$, so

$$
p_k^Tr_{k-1}
=\|r_{k-1}\|_2^2-u_{k-1,k}p_{k-1}^Tr_{k-1}
=\|r_{k-1}\|_2^2.
$$

For $k=1$, the same identity follows from $p_1=r_0$. Therefore,

$$
\boxed{\mu_k=\frac{\|r_{k-1}\|_2^2}{p_k^TAp_k}.}
$$

The residual also updates cheaply. Substituting $x_k=x_{k-1}+\mu_kp_k$ into its definition gives

$$
\begin{aligned}
r_k&=b-A(x_{k-1}+\mu_kp_k)\\
&=r_{k-1}-\mu_kAp_k.
\end{aligned}
$$

The product $Ap_k$ is already needed for the denominator of $\mu_k$, so this update needs no additional matrix-vector product.

(cg-direction-coefficient-derivation)=
### The Direction Coefficient

The direction update still contains $u_{k,k+1}=p_k^TAr_k/d_k$. The residual update gives

$$
Ap_k=\frac{r_{k-1}-r_k}{\mu_k}.
$$

Multiplying the step-length formula by $d_k$ gives $\mu_kd_k=\|r_{k-1}\|_2^2$. Using this identity, symmetry, and residual orthogonality, we obtain

$$
\begin{aligned}
u_{k,k+1}
&=\frac{(Ap_k)^Tr_k}{d_k}\\
&=\frac{r_{k-1}^Tr_k-r_k^Tr_k}{\mu_kd_k}\\
&=-\frac{\|r_k\|_2^2}{\|r_{k-1}\|_2^2}.
\end{aligned}
$$

Define its negative by

$$
\tau_k=-u_{k,k+1}
=\frac{\|r_k\|_2^2}{\|r_{k-1}\|_2^2}.
$$

The direction update becomes

$$
\boxed{p_{k+1}=r_k+\tau_kp_k.}
$$

## The Conjugate Gradient Algorithm

Start with $x_0=0$ and $r_0=b$. If $r_0=0$, stop. Otherwise, set $p_1=r_0$. For $k=1,2,\ldots$, perform

$$
\begin{aligned}
w_k&=Ap_k, &&\text{one matrix-vector product},\\
\mu_k&=\frac{r_{k-1}^Tr_{k-1}}{p_k^Tw_k}, &&\text{step length},\\
x_k&=x_{k-1}+\mu_kp_k, &&\text{solution update},\\
r_k&=r_{k-1}-\mu_kw_k, &&\text{residual update}.
\end{aligned}
$$

If $r_k=0$, stop before forming another direction. Otherwise, compute

$$
\begin{aligned}
\tau_k&=\frac{r_k^Tr_k}{r_{k-1}^Tr_{k-1}},\\
p_{k+1}&=r_k+\tau_kp_k.
\end{aligned}
$$

Positive definiteness ensures $p_k^TAp_k>0$ for every nonzero search direction. Each iteration uses one matrix-vector product, two inner products when previously computed residual norms are reused, and a few vector updates. Only a fixed number of vectors is stored. The full basis matrices are not formed.

For a general initial guess, start with $r_0=b-Ax_0$ and use exactly the same recurrence. It searches in $x_0+\mathcal K_k(A,r_0)$ and minimizes the energy error there. This is equivalent to applying the zero-start derivation to the correction equation $A(x-x_0)=r_0$.

In floating-point arithmetic, the recursively updated residual can drift from $b-Ax_k$. Before accepting the answer, recompute the residual explicitly and check a scaled quantity such as

$$
\eta_k=\frac{\|b-Ax_k\|_2}{\|A\|_2\|x_k\|_2+\|b\|_2}.
$$

A small $\eta_k$ means that small relative perturbations to $A$ and $b$ make $x_k$ an exact solution. This is a small normwise backward error. Forward accuracy also depends on the conditioning of $A$.

## Connection to Lanczos and Matrix Relations

CG and Lanczos build the same Krylov spaces with different bases. The Euclidean directions from our first attempt reappear as normalized CG residuals. Both $q_k$ from Lanczos and $r_{k-1}$ satisfy

$$
q_k,\ r_{k-1}\in\mathcal K_k,
\qquad
q_k,\ r_{k-1}\perp\mathcal K_{k-1}.
$$

Before termination, this new orthogonal direction is one-dimensional. Thus, for signs $\sigma_k\in\{-1,1\}$ determined by the Lanczos convention,

$$
\boxed{q_k=\sigma_k\frac{r_{k-1}}{\|r_{k-1}\|_2}.}
$$

The residuals, rather than the A-orthogonal search directions, are therefore parallel to the Lanczos vectors.

### Orthogonality Relations in Matrix Notation

Let $Q_k=[q_1,\ldots,q_k]$, $P_k=[p_1,\ldots,p_k]$, and $\mathcal R_k=[r_0,\ldots,r_{k-1}]$. The matrix $\mathcal R_k$ collects residual vectors. Their properties are

$$
\begin{aligned}
Q_k^TQ_k&=I,\\
\mathcal R_k^T\mathcal R_k&=\operatorname{diag}(\|r_0\|_2^2,\ldots,\|r_{k-1}\|_2^2),\\
P_k^TAP_k&=D_k.
\end{aligned}
$$

The matrices $Q_k^TAQ_k=T_k$ and $\mathcal R_k^TA\mathcal R_k$ are symmetric tridiagonal. The first is the Lanczos projection. Since the normalized residuals are the Lanczos vectors up to signs, $\mathcal R_k=Q_kS_k$ for a nonsingular diagonal matrix $S_k$ of signed residual norms. Therefore,

$$
\mathcal R_k^TA\mathcal R_k=S_kT_kS_k.
$$

### The Small Tridiagonal System

Write the CG approximation in the Lanczos basis as $x_k=Q_ky_k$. Residual orthogonality then gives

$$
\begin{aligned}
0=Q_k^Tr_k&=Q_k^Tb-Q_k^TAQ_ky_k,\\
T_ky_k&=\|b\|_2e_1.
\end{aligned}
$$

The right-hand side follows from $q_1=b/\|b\|_2$ and the orthogonality of the remaining columns to $q_1$. Here $e_1$ is the first coordinate vector. Since $A$ is SPD, so is $T_k=Q_k^TAQ_k$, and this projected system has a unique solution.

CG computes $x_k=Q_ky_k$ implicitly through its short recurrence, without storing $Q_k$ or explicitly forming and solving each $T_k$. This explains both its connection to Lanczos and its small storage requirement.

## Why the Name Conjugate Gradient?

The search directions are **conjugate** because $p_i^TAp_j=0$ for $i\ne j$. This makes the coefficient along each direction final: adding a new conjugate direction does not undo the minimizations already performed.

The **gradient** comes from the quadratic energy. Its gradient is

$$
\nabla\Phi(z)=Az-b,
\qquad
\nabla\Phi(x_k)=-r_k.
$$

Equivalently, the squared-error objective $L(z)=\tfrac12\|x-z\|_A^2$ differs from $\Phi(z)$ by a constant and has the same gradient. Thus $r_k$ is the **negative gradient**, the direction of steepest decrease in the Euclidean geometry. Steepest descent would use $r_k$ alone. CG adds a multiple of $p_k$ to make the new direction conjugate to all previous ones.

## Convergence of the CG Method

The construction explains which approximation CG selects. Its usefulness depends on how quickly the Krylov spaces contain accurate approximations.

### Finite Convergence in Exact Arithmetic

:::{prf:theorem} Finite Convergence of CG
:label: cg_finite_convergence

For an SPD $n\times n$ matrix, CG finds the exact solution in at most $n$ iterations in exact arithmetic. It may terminate earlier.
:::

:::{prf:proof}
If $r_0,\ldots,r_{n-1}$ are all nonzero, they are $n$ mutually orthogonal vectors and hence form a basis of $\mathbb R^n$. By {prf:ref}`cg_residual_basis`, $\mathcal K_n=\mathbb R^n$. The next residual $r_n$ is orthogonal to this whole space, so it must be zero. If any earlier residual vanishes, the solution has already been found.
:::

The Krylov space can stop growing before it fills $\mathbb R^n$. If $\mathcal K_{m+1}=\mathcal K_m$, then $r_m\in\mathcal K_{m+1}=\mathcal K_m$ and $r_m\perp\mathcal K_m$, forcing $r_m=0$. Thus termination occurs when growth stops, and the full solution has an expansion in the directions already constructed. We never need to assume that $\mathcal K_n$ has dimension $n$ for every right-hand side.

For a large sparse problem, the practical goal is an accurate approximation in far fewer than $n$ steps. Finite termination is a theoretical property. In floating-point arithmetic, conjugacy and residual orthogonality can be lost, so the $n$-step guarantee does not generally hold.

### The Error Polynomial and the Condition Number

The [Chebyshev analysis](chebyshev_acceleration.md) described error reduction through polynomials of $A$. The same viewpoint applies to CG. Since $x_k\in\mathcal K_k(A,b)$ and $b=Ax$, its error can be written

$$
x-x_k=p(A)x,
\qquad \deg p\le k,
\qquad p(0)=1.
$$

Conversely, every polynomial satisfying these conditions gives a candidate in $\mathcal K_k$: because $p(0)=1$, we can write $p(\lambda)=1-\lambda s(\lambda)$ with $\deg s\le k-1$. Then $x-p(A)x=s(A)b\in\mathcal K_k$. Thus minimizing over these polynomials is exactly the same as minimizing over the Krylov space. CG chooses the candidate with the smallest energy error.

To bound this error, write $x=\sum_j c_jv_j$ in an orthonormal eigenvector basis of $A$. Then

$$
\begin{aligned}
\|p(A)x\|_A^2
&=\sum_j\lambda_j c_j^2|p(\lambda_j)|^2\\
&\le\left(\max_j|p(\lambda_j)|\right)^2\|x\|_A^2.
\end{aligned}
$$

Let $\lambda_{\min}$ and $\lambda_{\max}$ be the smallest and largest eigenvalues. Bounding all values on this enclosing interval gives

$$
\|x-x_k\|_A
\le
\min_{\substack{\deg p\le k\\p(0)=1}}
\max_{\lambda\in[\lambda_{\min},\lambda_{\max}]}|p(\lambda)|\,\|x-x_0\|_A.
$$

The same bound holds for a general $x_0$ by applying the argument to the correction equation. In particular, CG has no larger energy error than the Chebyshev approximation of the same degree and initial guess.

For SPD $A$, define the 2-norm condition number and the associated convergence factor by

$$
\kappa_2(A)=\frac{\lambda_{\max}(A)}{\lambda_{\min}(A)},
\qquad
\rho=\frac{\sqrt{\kappa_2(A)}-1}{\sqrt{\kappa_2(A)}+1}.
$$

Applying the Chebyshev bound from {prf:ref}`thm:chebyshev-convergence` gives, for $k\ge1$,

$$
\boxed{
\|x-x_k\|_A
\le2\rho^k\|x-x_0\|_A.
}
$$

If $\kappa_2(A)=1$, then $A$ is a positive multiple of the identity and CG solves the system in one step. When the condition number is close to one, the factor $\rho$ is small and the bound predicts rapid convergence. For large $\kappa_2(A)$, the factor approaches one and the guarantee weakens. The factor bounds the overall error decay. It does not say that the actual error decreases by exactly $\rho$ at every step.

### Eigenvalue Clustering

The condition-number bound uses only the extreme eigenvalues. The actual rate depends on the whole eigenvalue distribution and on the initial error.

If $A$ has $s$ distinct eigenvalues $\theta_1,\ldots,\theta_s$, consider

$$
p(\lambda)=\prod_{j=1}^s\left(1-\frac{\lambda}{\theta_j}\right).
$$

It equals $1$ at zero and vanishes at every eigenvalue, so $p(A)=0$. This candidate has zero error, and CG optimality proves exact convergence in at most $s$ steps. Only eigenvalues present in the initial error are needed in this argument, so convergence can occur even sooner.

Tight eigenvalue clusters can also help because a polynomial can be small on the clusters without being small throughout the entire interval between the extreme eigenvalues. This can make convergence much faster than the condition-number bound predicts. The widths and separation of the clusters matter. Having $s$ clusters does not generally imply termination in $s$ steps.

The following experiment compares two $120\times120$ SPD systems with the same condition number, $\kappa_2(A)=100$. One has evenly spread eigenvalues, while the other has three tight clusters. Both start from zero with the same initial-error coefficients in their eigenvector bases. The plotted errors come from actual CG iterations, using the known solution to measure the energy error.

:::{figure} ../_static/cg_convergence.svg
:alt: Left, 120 eigenvalues spread between 1 and 100 compared with 120 eigenvalues grouped near 1, 20, and 100. Right, computed relative CG energy errors: the clustered system reaches 10 to the minus 12 in 23 iterations and the spread system in 73, both below their shared condition-number bound.
:width: 100%
:name: fig-cg-convergence

Same condition number, different convergence. The matrices are diagonal, so their coordinate axes are eigenvector directions. The spread spectrum is evenly spaced in $[1,100]$. The clustered spectrum has 40 eigenvalues in each of $[1,1.1]$, $[19.95,20.05]$, and $[99.9,100]$. In both cases, $x=(1,\ldots,1)^T/\sqrt{120}$, $b=Ax$, and $x_0=0$. The right panel shows the computed relative energy error and their common bound $2(9/11)^k$. Each run stops once the measured error is below $10^{-12}$. Clustering gives faster convergence here, even though all 120 eigenvalues remain distinct.
:::

This motivates [preconditioning](preconditioning.md): transform the system into an equivalent problem for which an inexpensive iterative solve converges faster. For CG, the transformation must retain the SPD structure. We develop this in [preconditioned CG](pcg.md).

### Superlinear Acceleration

CG sometimes converges faster as the iterations proceed. In this setting, **superlinear convergence** refers to an improving effective error-reduction factor, rather than one fixed geometric rate throughout the calculation. Such acceleration is possible, but not guaranteed for every problem or at every step.

The Lanczos connection helps explain it. The **Ritz values**, which are the eigenvalues of the small matrix $T_k=Q_k^TAQ_k$, approximate some eigenvalues of $A$. Once well-separated extreme eigenvalues are approximated accurately, the CG error polynomial can strongly suppress their error components. The remaining components may occupy a much narrower spectral interval, so CG behaves as though the condition number were smaller. This acceleration does not require any component to become exactly zero.

## Key Results to Remember

Assume $A$ is SPD and work in exact arithmetic. As in the derivation, take $x_0=0$, write $x=A^{-1}b$, and define

$$
e_k=x-x_k,
\qquad r_k=b-Ax_k=Ae_k,
\qquad \mathcal K_k=\mathcal K_k(A,b).
$$

The geometry is defined by $\langle u,v\rangle_A=u^TAv$ and $\|v\|_A=\sqrt{v^TAv}$.

The formulas below use $x_0=0$. The basis statements apply while the listed residuals and search directions are nonzero. For a general initial guess, apply them to the correction equation $A(x-x_0)=r_0$: replace $b$ by $r_0$ and the approximations $x_k$ by $x_k-x_0$. The errors and residuals are unchanged.

### One Krylov Space, Several Bases

The search directions, residuals, and Lanczos vectors span the same nested spaces:

$$
\begin{aligned}
\mathcal K_k
&=\operatorname{span}\{b,Ab,\ldots,A^{k-1}b\}\\
&=\operatorname{span}\{p_1,\ldots,p_k\}\\
&=\operatorname{span}\{r_0,\ldots,r_{k-1}\}\\
&=\operatorname{span}\{q_1,\ldots,q_k\}.
\end{aligned}
$$

The approximation and the next residual satisfy

$$
x_k\in\mathcal K_k,
\qquad
r_k\in\mathcal K_{k+1}.
$$

### The Orthogonality Relations

The search directions are **A-orthogonal**, while the residuals are **Euclidean-orthogonal**:

$$
\boxed{
\begin{aligned}
p_i^TAp_j&=0 &&(i\ne j),\\
r_i^Tr_j&=0 &&(i\ne j).
\end{aligned}
}
$$

The error and residual are orthogonal to the current search space in different inner products:

$$
\boxed{
\begin{aligned}
v^TAe_k&=0 &&\text{for every }v\in\mathcal K_k,\\
v^Tr_k&=0 &&\text{for every }v\in\mathcal K_k.
\end{aligned}
}
$$

These two statements are equivalent because $r_k=Ae_k$. In particular, $p_i^Tr_k=0$ for $i\le k$. A nonzero $r_k$ is therefore a new direction outside $\mathcal K_k$, even though it is generally not A-orthogonal to all the previous search directions.

### Optimality over the Whole Search Space

CG produces the unique best approximation in the energy norm:

$$
\boxed{x_k=\underset{z\in\mathcal K_k}{\operatorname{argmin}}\,\|x-z\|_A.}
$$

Equivalently, it minimizes the computable quadratic energy:

$$
x_k=\underset{z\in\mathcal K_k}{\operatorname{argmin}}\,\Phi(z),
\qquad
\Phi(z)=\frac12z^TAz-b^Tz.
$$

The projection identity makes the optimality explicit. For every $z\in\mathcal K_k$,

$$
\|x-z\|_A^2=\|x-x_k\|_A^2+\|x_k-z\|_A^2.
$$

Thus the energy error decreases as the search space grows. The Euclidean residual norm, however, need not decrease at every step.

The same optimality statement can be written in terms of the residual. Since

$$
\|r_k\|_{A^{-1}}^2
:=r_k^TA^{-1}r_k
=e_k^TAe_k,
$$

CG minimizes $\|b-Az\|_{A^{-1}}$ over $z\in\mathcal K_k$. It does not generally minimize the Euclidean residual norm $\|b-Az\|_2$.

### Computable Coefficients and Short Updates

With the normalization $p_1=r_0$ and $p_{k+1}=r_k+\tau_kp_k$, the coefficients are

$$
\begin{aligned}
\mu_k
&=\frac{p_k^Tr_0}{p_k^TAp_k}
=\frac{\|r_{k-1}\|_2^2}{p_k^TAp_k},\\
\tau_k&=\frac{\|r_k\|_2^2}{\|r_{k-1}\|_2^2}.
\end{aligned}
$$

The iterates, residuals, and directions update by

$$
\boxed{
\begin{aligned}
x_k&=x_{k-1}+\mu_kp_k,\\
r_k&=r_{k-1}-\mu_kAp_k,\\
p_{k+1}&=r_k+\tau_kp_k.
\end{aligned}
}
$$

The step along $p_k$ minimizes the energy along that line. Conjugacy ensures that this step also gives the minimum over the entire enlarged space without changing the old expansion coefficients. Stop if $r_k=0$, before generating another direction. Each iteration uses one matrix-vector product and stores only a fixed number of vectors.

### The Lanczos Connection

The normalized residuals are the Lanczos vectors up to signs:

$$
q_k=\sigma_k\frac{r_{k-1}}{\|r_{k-1}\|_2},
\qquad \sigma_k\in\{-1,1\}.
$$

For $Q_k=[q_1,\ldots,q_k]$, $P_k=[p_1,\ldots,p_k]$, and $\mathcal R_k=[r_0,\ldots,r_{k-1}]$,

$$
\begin{aligned}
Q_k^TQ_k&=I,\\
P_k^TAP_k&=D_k\quad\text{(positive diagonal)},\\
\mathcal R_k^T\mathcal R_k&\text{ is positive diagonal},\\
Q_k^TAQ_k=T_k&\text{ and }\mathcal R_k^TA\mathcal R_k\text{ are symmetric tridiagonal}.
\end{aligned}
$$

With $q_1=b/\|b\|_2$, the same approximation can be written

$$
x_k=Q_ky_k,
\qquad T_ky_k=\|b\|_2e_1.
$$

CG solves this projected problem implicitly through the short updates.

### Convergence and Its Limits

CG finds the exact solution in at most $n$ steps, or at most $s$ steps if $A$ has $s$ distinct eigenvalues. Earlier termination is possible when fewer eigenvalues are present in the initial error.

Its optimality can also be expressed through error polynomials:

$$
\|e_k\|_A
=\min_{\substack{\deg p\le k\\p(0)=1}}\|p(A)e_0\|_A.
$$

With $\kappa_2(A)=\lambda_{\max}(A)/\lambda_{\min}(A)$, the condition-number bound for $k\ge1$ is

$$
\boxed{
\|e_k\|_A
\le2\left(
\frac{\sqrt{\kappa_2(A)}-1}{\sqrt{\kappa_2(A)}+1}
\right)^k\|e_0\|_A.
}
$$

Eigenvalue clustering can give faster convergence than this bound predicts, and convergence can accelerate as some error components are suppressed. These improvements depend on the problem. They are not guaranteed at every step.

In floating-point arithmetic, the orthogonality relations and finite-termination guarantee need not hold exactly. Check the explicitly computed residual before accepting a solution, and use preconditioning when the unpreconditioned iteration converges too slowly.

## Optional Reading: The Geometry of the CG Coefficients

The algebraic derivation above gives the two CG coefficients. The geometry below gives a second derivation of the same formulas and can be skipped without affecting the rest of the chapter.

Consider one step starting at $x_k$, with residual $r_k=b-Ax_k$. The coefficient $\tau_k$ sets the next search direction,

$$
p_{k+1}=r_k+\tau_kp_k,
$$

and $\mu_{k+1}$ sets the distance traveled along it. Assume exact arithmetic and a nonzero residual throughout this discussion.

### The Step Length on the Energy Curve

The {ref}`earlier derivation of the step length <cg-step-length-derivation>` found $\mu_{k+1}$ as the coefficient of $p_{k+1}$ in the A-orthogonal expansion of the solution. We can recover that coefficient by asking where the search line reaches its lowest energy. Substituting $x_k+tp_{k+1}$ into the energy $\Phi(z)=\tfrac12 z^TAz-b^Tz$ gives a quadratic in $t$,

$$
\begin{aligned}
f(t)&=\Phi(x_k+tp_{k+1})\\
&=\Phi(x_k)+t\,p_{k+1}^T(Ax_k-b)
  +\tfrac12 t^2\,p_{k+1}^TAp_{k+1}\\
&=\Phi(x_k)-t\,\|r_k\|_2^2
  +\tfrac12 t^2\,p_{k+1}^TAp_{k+1}.
\end{aligned}
$$

In the expansion, symmetry of $A$ lets us combine the two cross terms. The last line uses $Ax_k-b=-r_k$ and $p_{k+1}^Tr_k=\|r_k\|_2^2$. The latter identity holds because $p_{k+1}=r_k+\tau_kp_k$ with $r_k\perp p_k$ for $k\ge1$, and $p_1=r_0$ at the first step.

The parabola opens upward because $p_{k+1}^TAp_{k+1}>0$. Its minimum occurs where the slope vanishes. Differentiating with respect to $t$ gives

$$
f'(t)=-\|r_k\|_2^2+t\,p_{k+1}^TAp_{k+1}.
$$

Setting $f'(\mu_{k+1})=0$ recovers the step length and the corresponding iterate,

$$
\mu_{k+1}=\frac{\|r_k\|_2^2}{p_{k+1}^TAp_{k+1}},
\qquad x_{k+1}=x_k+\mu_{k+1}p_{k+1}.
$$

For example, take $A=\left(\begin{smallmatrix}2&-1\\-1&2\end{smallmatrix}\right)$, $b=(3,0)^T$, and $x_0=0$. Along $p_1=b$ the energy is $f(t)=-9t+9t^2$. Its slope $f'(t)=-9+18t$ vanishes at $t=1/2$, as shown below.

:::{figure} ../_static/cg_step_parabola.svg
:alt: The energy reaches its minimum at one half, where the downhill slope falls from nine to zero.
:align: center

The right panel plots the downhill slope $-f'(t)$. The shaded triangle has area $\tfrac12(1/2)(9)=2.25$, the energy lost over the first step in the left panel.
:::

The contour plot shows two consecutive CG steps in the plane through $x_{k-1}$, $x_k$, and $x_{k+1}$. The horizontal step follows $p_k$. At its endpoint $x_k$, the line search gives $p_k^Tr_k=0$. Since $r_k=-\nabla\Phi(x_k)$ is normal to the energy contour there, the horizontal search line is tangent to that contour.

:::{figure} ../_static/cg_step_in_plane.svg
:width: 80%
:alt: The step from x_{k-1} to x_k is horizontal and tangent to the contour at x_k. The projected previous residual is horizontal, the current residual is vertical, and the step to x_{k+1} is diagonal.
:align: center

The ellipses are energy contours in the plane of the two steps. The black arrows show the downhill direction in that plane at each iterate. The blue triangle shows how the search direction is built from the residual and a component along the previous step. At $x_{k-1}$ this is the projection of $r_{k-1}$. At $x_k$ the full residual $r_k=p_{k+1}-\tau_kp_k$ already lies in the plane.
:::

The residuals are perpendicular in the full space. Because $r_k$ lies in the drawing plane, projecting $r_{k-1}$ into that plane preserves their right angle. From $x_k$, CG follows $p_{k+1}=r_k+\tau_kp_k$ to $x_{k+1}$, the minimum of the energy in this slice.

### Why the Direction Coefficient Is a Squared Length Ratio

Now draw the last two residuals, $r_{k-1}$ and $r_k$, as arrows from a common origin. They are perpendicular, so we can use them as the horizontal and vertical directions in the drawing. The arrow joining their tips is parallel to $Ap_k$, by the residual update

$$
r_{k-1}-r_k=\mu_kAp_k.
$$

We want $p_{k+1}$ to be A-orthogonal to $p_k$. In this picture that means ordinary perpendicularity to the arrow joining the tips, since $p_{k+1}^T(r_{k-1}-r_k)=\mu_kp_{k+1}^TAp_k=0$.

Only the part of each search direction in the drawing is visible. Write $E_k=\operatorname{span}\{r_{k-1},r_k\}$ for this residual plane and $\operatorname{proj}_{E_k}$ for Euclidean orthogonal projection onto it. The part of $p_k$ inherited from earlier steps belongs to $\mathcal K_{k-1}$, perpendicular to both residuals. Consequently,

$$
\operatorname{proj}_{E_k}p_k=r_{k-1},
\qquad
\operatorname{proj}_{E_k}p_{k+1}=r_k+\tau_kr_{k-1}.
$$

:::{figure} ../_static/cg_direction_triangle.svg
:width: 70%
:alt: Two right triangles share a vertical leg. The line joining the residual tips is perpendicular to the projected next search direction.
:align: center

The red segment joins the residual tips. The blue diagonal shows the part of $p_{k+1}$ in their plane. Projection preserves its right angle with the red segment because the red segment lies entirely in the plane.
:::

Let $a=\|r_{k-1}\|_2$ and $b_k=\|r_k\|_2$ be the horizontal and vertical lengths. The red segment has slope $-b_k/a$. The blue diagonal rises by $b_k$ and runs to the right by $\tau_k a$, so its slope is $b_k/(\tau_k a)$. The {ref}`earlier derivation of the direction coefficient <cg-direction-coefficient-derivation>` used the residual update to evaluate the dot product $p_k^TAr_k$. Here the same conjugacy condition appears as a right angle. Setting the blue slope equal to the negative reciprocal of the red slope gives the same formula for $\tau_k$ directly from the drawing.

$$
\frac{b_k}{\tau_k a}=\frac{a}{b_k}
\qquad\Longrightarrow\qquad
\tau_k=\frac{b_k^2}{a^2}
=\frac{\|r_k\|_2^2}{\|r_{k-1}\|_2^2}.
$$

### A Direction Perpendicular to All the Residual Edges

The triangle uses only $p_{k+1}^TAp_k=0$. Earlier we proved that $p_{k+1}$ is A-orthogonal to *every* preceding search direction. Each of those conditions gives a right angle of its own. Draw $r_0,\ldots,r_k$ from a common origin and join consecutive tips. For each edge, the residual update and conjugacy give

$$
p_{k+1}^T(r_{i-1}-r_i)=\mu_i p_{k+1}^TAp_i=0,
\qquad 1\le i\le k.
$$

For two tips, $p_{k+1}$ is perpendicular to the line joining them. For three, it is perpendicular to the plane through them. In general, the **affine span** of the tips is the set reached by starting at one tip and adding linear combinations of these edges. The search direction is perpendicular to every direction in that span.

:::{figure} ../_static/cg_direction_corner.svg
:width: 90%
:alt: Three mutually orthogonal residuals end at the vertices of a triangle. The next search direction is perpendicular to its face.
:align: center

The drawing uses three consecutive residuals, ending at $r_{k+1}$. The red edges are multiples of $Ap_k$ and $Ap_{k+1}$, so the next direction $p_{k+2}$ is perpendicular to both. The drawing shows its projection onto the space of the three residuals.
:::

Perpendicularity to each edge means that $p_{k+1}$ has the same dot product with every residual from $r_0$ through $r_k$. We can evaluate the last one using $p_k^Tr_k=0$, which gives $(r_k+\tau_kp_k)^Tr_k=\|r_k\|_2^2$. Thus

$$
p_{k+1}^Tr_0
=p_{k+1}^Tr_1
=\cdots
=p_{k+1}^Tr_k
=\|r_k\|_2^2.
$$

This gives an energy interpretation. At any earlier iterate $x_j$, the initial rate of decrease along a direction $d$ is $-\nabla\Phi(x_j)^Td=r_j^Td$. Thus $p_{k+1}$ would give the same initial downhill rate at *every* iterate from $x_0$ through $x_k$. To see how the direction combines information from these iterates, expand it in the orthogonal residual basis. The coefficient of $r_j$ is $p_{k+1}^Tr_j/\|r_j\|_2^2$, so

$$
\frac{p_{k+1}}{\|r_k\|_2^2}
=\sum_{j=0}^k\frac{r_j}{\|r_j\|_2^2}
=\frac{p_k}{\|r_{k-1}\|_2^2}+\frac{r_k}{\|r_k\|_2^2}.
$$

Each term in the sum gives downhill rate one at its own iterate and zero at the others. The last equality separates the newest term. In the original derivation, orthogonality made all but the most recent projection coefficient vanish. In the residual picture, $p_k$ already contains every earlier residual with the required relative weights. Multiplying through by $\|r_k\|_2^2$ gives the same short update, $p_{k+1}=r_k+\tau_kp_k$.
