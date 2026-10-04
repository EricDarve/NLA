# Least Squares Problems

In data fitting, a model often has fewer parameters than there are measurements. The vector $x$ contains the parameters, $Ax$ contains the model's predictions, and $b$ contains the observations. Noise or limitations of the model can make it impossible to match every observation exactly. **Linear least squares chooses the parameters that minimize the sum of squared discrepancies.**

## From Linear Equations to Least Squares

Let $A\in\mathbb{R}^{m\times n}$ and $b\in\mathbb{R}^m$. When $m>n$, the system $Ax=b$ has more equations than unknowns and is called **overdetermined**. An exact solution exists if and only if $b$ lies in the column space of $A$. Having more equations than unknowns does not by itself rule out a solution, but it is common for these equations to be inconsistent.

The **residual** $r=b-Ax$ measures the discrepancy between the observations and the predictions. The linear least-squares problem is to find a vector $x^*$ satisfying

$$
\boxed{\|Ax^*-b\|_2=\min_{x\in\mathbb{R}^n}\|Ax-b\|_2.}
$$

Minimizing this norm is equivalent to minimizing its square,

$$
\|Ax-b\|_2^2=\sum_{i=1}^m\bigl((Ax)_i-b_i\bigr)^2,
$$

which explains the name *least squares*. If the equations are consistent, the minimum is zero and the minimizers are exactly the solutions of $Ax=b$.

The definition applies to matrices of any shape or rank. The main case in this chapter is $m\ge n$ with linearly independent columns, also called **full column rank**. Later, the SVD will allow us to handle dependent columns and underdetermined systems, where $m<n$.

## Geometry, Existence, and Uniqueness

Every prediction $Ax$ lies in the column space of $A$. Thus least squares asks for the point in this subspace closest to $b$. From the [geometry of orthogonal projections](projections.md), that point is the orthogonal projection $p$ of $b$ onto the column space.

To see why this minimizes the residual, write

$$
b-Ax=(b-p)+(p-Ax).
$$

The first term is perpendicular to the column space, and the second lies in it. By the Pythagorean theorem,

$$
\|b-Ax\|_2^2=\|b-p\|_2^2+\|p-Ax\|_2^2.
$$

The first term is fixed, and the second is minimized at zero. Since $p$ belongs to the column space, some $x^*$ satisfies $Ax^*=p$. **A least-squares minimizer therefore always exists, and the fitted vector $Ax^*$ is unique.**

The coefficient vector $x^*$ is unique precisely when $A$ has full column rank. Indeed, two minimizers produce the same fitted vector, so their difference lies in the null space of $A$. If this null space contains a nonzero vector $z$, then $x^*+tz$ is another minimizer for every scalar $t$.

The projection also characterizes the minimizer through its residual:

$$
r^*=b-Ax^*\perp\operatorname{range}(A),
\qquad A^Tr^*=0.
$$

Equivalently, every minimizer satisfies the **normal equations**

$$
A^TAx^*=A^Tb.
$$

Conversely, these equations make the residual perpendicular to the column space, so they characterize all least-squares minimizers, including when the minimizer is not unique. This geometric condition connects the three solution methods below.

## Three Approaches to Computing a Solution

### QR Factorization

For full column rank $A$, the **reduced QR factorization** is

$$
A=Q_1R_1,
\qquad
Q_1\in\mathbb{R}^{m\times n},
\quad R_1\in\mathbb{R}^{n\times n}.
$$

The columns of $Q_1$ are orthonormal, and $R_1$ is upper triangular and nonsingular. Here $Q_1^TQ_1=I_n$; $Q_1$ is rectangular when $m>n$.

Since $p=Q_1Q_1^Tb$, the condition $Ax^*=p$ reduces to the triangular system

$$
R_1x^*=Q_1^Tb,
$$

which can be solved by backward substitution. The factor $R_1$ has the same singular values as $A$, so it can still be ill-conditioned. Orthogonal transformations preserve lengths and provide a way to compute the factorization without forming $A^TA$.

The chapter begins with three ways to construct QR:

- **[Householder reflections](householder_reflections.md)** give a backward-stable method for dense QR factorization.
- **[Givens rotations](givens_rotations.md)** eliminate selected entries by acting on pairs of rows.
- **[Modified Gram-Schmidt](modified_gram_schmidt.md)** builds an orthonormal basis one vector at a time. This makes it useful when vectors arrive sequentially, as in iterative methods. Rounding can cause a loss of orthogonality, which requires separate attention.

Their properties and error analysis lead to the [QR method for least squares](LS_using_QR.md).

### Normal Equations

The normal equations give a direct connection to the solvers from the preceding chapter. When $A$ has full column rank, $A^TA$ is symmetric positive definite because

$$
z^TA^TAz=\|Az\|_2^2>0\qquad\text{for }z\ne0.
$$

Thus we can form $A^TA$ and $A^Tb$ and solve using Cholesky factorization. The [normal-equations method](normal_equations.md) usually requires less arithmetic than Householder QR, but forming $A^TA$ can lose accuracy. For full column rank $A$, define

$$
\kappa_2(A)=\frac{\sigma_{\max}(A)}{\sigma_{\min}(A)}.
$$

Then

$$
\kappa_2(A^TA)=\kappa_2(A)^2.
$$

The singular values are squared when forming $A^TA$, so rounding errors in this matrix can overwhelm information associated with small singular values. This is why a stable Cholesky solve alone does not ensure an accurate least-squares answer. The chapter compares this method with QR in more detail.

### SVD and the Minimum-Norm Solution

If $A$ has dependent columns, there are infinitely many least-squares minimizers. To select one, we can ask for the minimizer with the smallest 2-norm. This **minimum-norm least-squares solution** is unique: it is the minimizer perpendicular to the null space of $A$.

The [SVD method](LS_using_SVD.md) computes this solution for any matrix, regardless of its shape or rank. It leads to the **Moore-Penrose pseudoinverse** $A^\dagger$, with

$$
x_{\min}=A^\dagger b.
$$

The two minimizations have a specific order: first minimize the residual, then minimize $\|x\|_2$ among all vectors attaining that residual. The SVD makes both steps explicit by separating the directions associated with positive and zero singular values.
