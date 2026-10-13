# The Method of Normal Equations

The **method of normal equations** transforms a linear least-squares problem into a square symmetric system. It connects the geometry of orthogonal projection to the Cholesky solver from the preceding chapter.

The linear least-squares problem is defined as finding the vector $x$ that minimizes the 2-norm of the residual:

$$
x^* \in \arg\min_{x \in \mathbb{R}^n} \| Ax - b \|_2,
$$

where $A\in\mathbb{R}^{m\times n}$ and $b\in\mathbb{R}^m$. The main case here is $m\ge n$ with full column rank. The geometric derivation also applies when the columns are dependent, so $x^*$ denotes any minimizer until uniqueness has been established.

## Geometric Derivation

Every fitted vector $Ax$ lies in the column space $\operatorname{range}(A)$. Thus least squares asks for the point in this subspace closest to $b$. From the [geometry of orthogonal projections](projections.md), this point is the orthogonal projection $p$ of $b$ onto the column space.

Because $p$ belongs to the column space, some coefficient vector $x^*$ satisfies $Ax^*=p$. The corresponding residual is

$$
r^*=b-p=b-Ax^*.
$$

The projection separates $b$ into a fitted part $p$ in the column space and a residual $r^*$ perpendicular to it.

```{figure} ../_static/normal_equations_projection.svg
:alt: A blue vector b points above the column space, drawn as a horizontal line. Its purple projection p equals A x star. The orange residual r star points from p to b and is perpendicular to the line. Another fitted vector A x lies farther to the right. The vector from A x to b is the hypotenuse of a right triangle with sides p minus A x and r star.
:width: 720px
:name: fig-normal-equations-projection

The fitted vector $p=Ax^*$ is closest to $b$. For any other fitted vector $Ax$, the discrepancy $b-Ax$ splits into the perpendicular residual $r^*$ and the displacement $p-Ax$ within the column space. The column space is drawn as a line; the same argument holds for subspaces of any dimension.
```

To see why this projection minimizes the residual, compare $p$ with any other fitted vector $Ax$. Write

$$
b-Ax=(b-p)+(p-Ax)=r^*+(p-Ax).
$$

The first term is perpendicular to the column space, while the second lies in it. The Pythagorean theorem therefore gives

$$
\boxed{\|b-Ax\|_2^2=\|r^*\|_2^2+\|p-Ax\|_2^2.}
$$

The first term is fixed. The second is nonnegative and is zero precisely when $Ax=p$. Moving away from $p$ within the column space can only increase the squared distance to $b$. This proves that $x^*$ minimizes the residual and that every minimizer has the same fitted vector $p$.

We can express the perpendicularity condition using the columns $a_1,\ldots,a_n$ of $A$. A vector is perpendicular to their span if and only if it is perpendicular to each column:

$$
a_j^Tr^*=0\quad(j=1,\ldots,n)
\quad\Longleftrightarrow\quad
A^Tr^*=0.
$$

Substituting $r^*=b-Ax^*$ gives $A^T(b-Ax^*)=0$, which yields a square system for the coefficients.

## The Normal Equations

The resulting system is known as the **normal equations**:

$$
\boxed{(A^T A)x^*=A^T b.}
$$

The word *normal* refers to the residual being perpendicular to the column space. These $n$ equations characterize all least-squares minimizers: any solution has a perpendicular residual, and the Pythagorean argument then proves it minimizes the distance. They remain valid when $A$ is rank deficient, although the coefficient vector is then not unique.

### Derivation by Minimizing the Quadratic Form

Minimizing the residual norm is equivalent to minimizing its square. The same equations follow by differentiating this quadratic function.

:::{prf:proof}
Expanding the squared norm gives

$$
\begin{aligned}
f(x)&=\|Ax-b\|_2^2\\
&=x^TA^TAx-2x^TA^Tb+b^Tb.
\end{aligned}
$$

The two cross terms combine because $b^TAx=x^TA^Tb$ is a scalar. Since $A^TA$ is symmetric, the gradient (the vector of partial derivatives) is

$$
\nabla f(x)=2A^TAx-2A^Tb.
$$

At a minimizer this gradient is zero, giving $A^TAx^*=A^Tb$. Conversely, a solution of these equations has the perpendicular residual established above, so it is a global minimizer.
:::

### Properties of $A^T A$

To choose a solver for the normal equations, we need the structure of $A^TA$. It is symmetric, since $(A^TA)^T=A^TA$, and positive semidefinite because

$$
z^TA^TAz=\|Az\|_2^2\ge0.
$$

If $A$ has full column rank, then $Az\ne0$ for every nonzero $z$, making this expression strictly positive. Thus $A^TA$ is **symmetric positive definite (SPD)** and the normal equations have a unique solution $x^*$. If the columns are dependent, $A^TA$ is singular. The fitted vector is still unique, but adding any vector in the null space of $A$ to $x^*$ gives another minimizer.

For full column rank, the exact solution can be written as

$$
x^*=(A^TA)^{-1}A^Tb.
$$

In computation, we solve using [Cholesky factorization](cholesky.md) rather than form the inverse. Its stability guarantees apply to the matrix supplied to the solver; forming $A^TA$ in floating-point arithmetic introduces a separate source of error.

## Algorithm Summary

For full column rank $A$, the method has three steps:

1.  **Form the Gram matrix:** Compute the $n \times n$ matrix $C = A^T A$.
2.  **Form the right-hand side:** Compute the $n$-vector $d = A^T b$.
3.  **Solve using Cholesky:** Factor $C=LL^T$, then solve $Ly=d$ and $L^Tx^*=y$.

## Analysis and Drawbacks

### Computational Cost

Counting multiplication and addition as separate flops, the leading costs are:

-   Computing one triangle of $A^TA$ and copying it to the other takes approximately $mn^2$ flops. Computing every entry independently would take about $2mn^2$.
-   Computing $A^T b$ takes approximately $2mn$ flops.
-   Cholesky factorization takes approximately $\frac13n^3$ flops; the two triangular solves add $O(n^2)$.

The total leading cost is $mn^2+\frac13n^3$, compared with $2mn^2-\frac23n^3$ for Householder QR. For $m\gg n$, the normal-equations method uses about half the factorization arithmetic. This saving is the reason to consider it, but it comes with an accuracy risk.

:::{admonition} Accuracy Risk: Squaring the Condition Number
:class: warning

For full column rank $A$, the eigenvalues of $A^TA$ are $\sigma_i(A)^2$. Therefore the condition number of the system we form is

$$
\kappa_2(A^T A) = (\kappa_2(A))^2
$$

Rounding errors in forming $A^TA$ can be large relative to its smallest eigenvalues. These errors can obscure the directions associated with small singular values of $A$, so a stable Cholesky solve does not guarantee an accurate solution of the original least-squares problem. If the columns are nearly dependent, the rounded Gram matrix may even lose positive definiteness and Cholesky can fail.
:::

For example, $\kappa_2(A)=10^6$ gives $\kappa_2(A^TA)=10^{12}$. This signals a greater risk of losing precision, but does not by itself determine the error in a particular computed solution. The error also depends on $b$ and the rounding errors in forming and solving the system.

The [QR method](LS_using_QR.md) avoids forming $A^TA$ and works with a triangular factor having the same singular values as $A$. It preserves the conditioning of that factorization step, while the original least-squares problem can still be sensitive to perturbations.
