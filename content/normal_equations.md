# The Method of Normal Equations

The first approach we will study for solving the linear least-squares problem is the **method of normal equations**. It's an intuitive method that transforms the problem into a square linear system that we already know how to solve.

The linear least-squares problem is defined as finding the vector $x$ that minimizes the 2-norm of the residual:

$$
x^* = \arg\min_{x \in \mathbb{R}^n} \| Ax - b \|_2
$$

where $A$ is an $m \times n$ matrix and $b$ is an $m$-vector. We typically assume $m \ge n$.

## Geometric Derivation

Minimizing the Euclidean distance $\| Ax - b \|_2$ is equivalent to finding the vector $Ax$ in the column space of $A$, denoted $\text{span}(A)$, that is closest to $b$. This occurs when the residual vector, $r = b - Ax$, is **orthogonal** to the column space of $A$.

This orthogonality condition can be stated mathematically as:

$$
A^T (b - Ax) = 0
$$

Expanding this expression gives:

$$
A^T b - A^T A x = 0
$$

Rearranging the terms yields a linear system of equations for $x$.

## The Normal Equations

The resulting system is known as the **normal equations**:

$$
(A^T A) x = A^T b
$$

This is a fundamental equation in linear algebra and optimization. It transforms the original, often rectangular, least-squares problem into a square $n \times n$ linear system.

### Properties of $A^T A$

The matrix $A^T A$ has several important properties:
1.  It is always **symmetric**, since $(A^T A)^T = A^T (A^T)^T = A^T A$.
2.  If the columns of $A$ are linearly independent (i.e., $A$ has **full column rank**), then $A^T A$ is **symmetric positive definite (SPD)**. This guarantees that a unique solution $x$ exists.

Because $A^T A$ is SPD (under the common assumption of full rank), the system can be solved efficiently and stably using **Cholesky factorization**. The explicit solution can be written as:
$$
x = (A^T A)^{-1} A^T b
$$

## Algorithm Summary

The method can be summarized in three steps:
1.  **Form the Gram matrix:** Compute the $n \times n$ matrix $C = A^T A$.
2.  **Form the right-hand side:** Compute the $n$-vector $d = A^T b$.
3.  **Solve the system:** Solve the $n \times n$ symmetric system $Cx = d$ for $x$, typically using a Cholesky factorization of $C$.

## Analysis and Drawbacks

### Computational Cost
The computational cost is dominated by the formation of $A^T A$.
-   Computing $A^T A$ takes approximately $m n^2$ floating-point operations (flops).
-   Computing $A^T b$ takes approximately $2mn$ flops.
-   Solving the system with Cholesky factorization takes $\approx \frac{1}{3}n^3$ flops.

For a "tall and skinny" matrix where $m \gg n$, the total cost is approximately $O(mn^2)$.

:::{admonition} Potential for Severe Numerical Instability
:class: warning

The primary drawback of the normal equations method is its impact on the **condition number** of the problem. The condition number of the Gram matrix $A^T A$ is the square of the condition number of the original matrix $A$:

$$
\kappa_2(A^T A) = (\kappa_2(A))^2
$$

**Implication:** If $A$ is even moderately ill-conditioned (e.g., $\kappa_2(A) \approx 10^6$), the matrix $A^T A$ will be severely ill-conditioned ($\kappa_2(A^T A) \approx 10^{12}$). Solving a system with such a high condition number can lead to a significant loss of precision, potentially rendering the computed solution useless. This squaring of the condition number is the main reason why the normal equations are often avoided in practice in favor of more numerically robust methods like QR factorization.
:::

## Derivation by Minimizing the Quadratic Form

````{prf:proof}
There is another proof based on optimality conditions. It is less elegant than that one given above, but it is more general and can be extended to other problems. This approach treats the least-squares problem as an optimization problem.

The objective is to find the vector $x$ that minimizes the squared 2-norm of the residual:

$$
f(x) = \| Ax - b \|_2^2
$$

### Step 1: Expand the Objective Function

First, we rewrite the squared norm as a dot product (or transpose product):

$$
f(x) = (Ax - b)^T (Ax - b)
$$

Expanding this product gives:

$$
\begin{aligned}
f(x) &= (x^T A^T - b^T)(Ax - b) \\
&= x^T A^T Ax - x^T A^T b - b^T Ax + b^T b
\end{aligned}
$$

Since $b^T A x$ is a scalar quantity, it is equal to its own transpose: $b^T A x = (b^T A x)^T = x^T A^T b$. This allows us to combine the two middle terms:

$$
f(x) = x^T (A^T A) x - 2 x^T (A^T b) + \|b\|_2^2
$$

This expression is a **quadratic form** in $x$. It's a multidimensional parabola (a paraboloid). Our goal is to find the value of $x$ at the bottom of this bowl, which is its minimum. 

### Step 2: Find the Minimum using the Gradient

The minimum of a convex function occurs at a stationary point, where its gradient with respect to $x$ is the zero vector. Let's compute the gradient of $f(x)$, which we denote as $\nabla_x f(x)$.

Using the standard rules for matrix calculus:

* The gradient of $x^T C x$ is $2Cx$ (for a symmetric matrix $C$ like $A^TA$).
* The gradient of $c^T x$ (or $x^T c$) is $c$.

Applying these rules to our objective function:

$$
\begin{gathered}
\nabla_x f(x) = \nabla_x \left( x^T (A^T A) x - 2 x^T (A^T b) + \|b\|_2^2 \right) \\
= 2(A^T A)x - 2(A^T b) + 0
\end{gathered}
$$

To find the minimum, we set the gradient to zero:

$$
2(A^T A)x - 2(A^T b) = 0
$$

### Step 3: Solve for $x$

Dividing by 2 and rearranging the terms, we arrive at the previous **normal equations**:

$$
(A^T A)x = A^T b
$$

The solution $x$ to this system is the vector that minimizes the least-squares objective function $f(x)$.

````

:::{admonition} Confirming a Minimum with the Hessian
:class: note

To be certain that this stationary point is a minimum (and not a maximum or saddle point), we can examine the second derivative, or the **Hessian matrix** ($\nabla_x^2 f(x)$). The Hessian of our objective function is:

$$
\nabla_x^2 f(x) = 2(A^T A)
$$

If the matrix $A$ has full column rank, then $A^T A$ is **positive definite**. A positive definite Hessian confirms that the function is strictly convex, and therefore the stationary point we found is a unique global minimum. 👍
:::