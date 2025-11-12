# The Chebyshev Iteration Method

The Chebyshev iteration method, also known as a Chebyshev semi-iterative method, is a powerful technique for accelerating the convergence of iterative solvers for the linear system $A\mathbf{x} = \mathbf{b}$. It is a **nonstationary iterative method**, meaning it improves upon basic stationary schemes (like the Richardson iteration) by using acceleration parameters $\tau_k$ that change at each step $k$. This strategy, which leverages the properties of Chebyshev polynomials, can result in significantly faster convergence.

## From Stationary to Nonstationary Iteration

Classical iterative methods can be expressed using the preconditioned "residual update" form. Let's start with the familiar splitting-style approach:

$$
\begin{aligned}
A &= M - N \\
\mathbf{x}^{(k+1)} &= M^{-1} (\mathbf{b} + N \mathbf{x}^{(k)}) \\
\mathbf{x}^{(k+1)} &= M^{-1} (\mathbf{b} + (M - A) \mathbf{x}^{(k)})
\end{aligned}
$$

This implies that

$$
\mathbf{x}^{(k+1)} = \mathbf{x}^{(k)} + M^{-1} (\mathbf{b} - A \mathbf{x}^{(k)})
$$

This is the residual iteration form. If we define $\tilde{A} = M^{-1}A$ and $\tilde{\mathbf{b}} = M^{-1}\mathbf{b}$, the iteration becomes:

$$\mathbf{x}^{(k+1)} = \mathbf{x}^{(k)} + (\tilde{\mathbf{b}} - \tilde{A}\mathbf{x}^{(k)})$$

This is the basic **Richardson iteration**. A stationary method would apply a single, fixed acceleration parameter $\tau$. The key idea of a nonstationary method is to allow this parameter to vary at each step.

For the remainder of this analysis, we will assume such preconditioning $M^{-1}A$ has already been applied and (for notational simplicity) analyze the iteration for the system $A\mathbf{x} = \mathbf{b}$. The general nonstationary form is:

$$\mathbf{x}^{(k+1)} = \mathbf{x}^{(k)} + \tau_k (\mathbf{b} - A\mathbf{x}^{(k)}), \quad k=0, 1, \dots$$

Here, $\{\tau_k\}$ is the sequence of step-varying parameters.

## Error Propagation and Matrix Polynomials

We can analyze the method's convergence by tracking the error $\mathbf{e}^{(k)} = \mathbf{x} - \mathbf{x}^{(k)}$. We can derive the recurrence for the error as follows:

$$
\begin{aligned}
\mathbf{e}^{(k+1)} &= \mathbf{x} - \mathbf{x}^{(k+1)} \\
&= \mathbf{x} - \left[ \mathbf{x}^{(k)} + \tau_k (\mathbf{b} - A\mathbf{x}^{(k)}) \right] \\
&= (\mathbf{x} - \mathbf{x}^{(k)}) - \tau_k (A\mathbf{x} - A\mathbf{x}^{(k)}) \quad (\text{since } \mathbf{b} = A\mathbf{x}) \\
&= \mathbf{e}^{(k)} - \tau_k A \mathbf{e}^{(k)} \\
\mathbf{e}^{(k+1)} &= (I - \tau_k A) \mathbf{e}^{(k)}
\end{aligned}
$$

By unrolling this recurrence, we can express the error after $k$ steps as the result of a matrix polynomial applied to the initial error $\mathbf{e}^{(0)}$:

$$\mathbf{e}^{(k)} = (I - \tau_{k-1} A) \cdots (I - \tau_1 A) (I - \tau_0 A) \mathbf{e}^{(0)}$$

We define this product as the **matrix polynomial** $P_k(A)$:

$$\mathbf{e}^{(k)} = P_k(A) \mathbf{e}^{(0)} \quad \text{where} \quad P_k(A) = \prod_{l=0}^{k-1} (I - \tau_l A)$$

## The Polynomial Optimization Problem

The effectiveness of the iteration after $k$ steps is determined entirely by the properties of this $k$-th degree polynomial. This construction imposes two fundamental constraints on the corresponding scalar polynomial $P_k(\lambda) = \prod_{l=0}^{k-1} (1 - \tau_l \lambda)$:

1.  **Degree:** $P_k(\lambda)$ is a polynomial of degree $k$.
2.  **Normalization:** The polynomial must satisfy $P_k(0) = 1$. This is self-evident from the product form (setting $\lambda=0$) and is a necessary constraint for any method derived from this iterative form.

From this factorization, we also gain a key insight: the iteration parameters $\{\tau_l\}$ are the **inverses of the roots** of the scalar polynomial $P_k(\lambda)$, since $P_k(1/\tau_l) = 0$.

The central goal of the Chebyshev method is to choose the sequence of parameters $\{\tau_l\}$ to minimize the norm of the error $\mathbf{e}^{(k)}$ in the fastest possible way. Since $\mathbf{e}^{(k)} = P_k(A) \mathbf{e}^{(0)}$, this is equivalent to finding the polynomial $P_k(A)$ that is "as small as possible" in some sense.

If $A$ is diagonalizable, its eigenvectors $\mathbf{v}_i$ corresponding to eigenvalues $\lambda_i$ are mapped as:

$$P_k(A)\mathbf{v}_i = P_k(\lambda_i)\mathbf{v}_i$$

This shows that the error components in the direction of $\mathbf{v}_i$ are damped by a factor of $P_k(\lambda_i)$. The optimization problem is thus reduced to finding a polynomial $P_k(\lambda)$ of degree $k$ that satisfies:

1.  $P_k(0) = 1$
2.  $P_k(\lambda)$ is as small as possible across the spectrum of $A$, $\sigma(A)$.

This strategy—finding the optimal polynomial that is small on the spectrum while constrained to be 1 at the origin—is the foundation of the Chebyshev iteration method.

The convergence analysis of the Chebyshev iteration method is rooted in optimization theory, specifically solving a **minimax polynomial approximation problem**. The goal is to construct a polynomial that is maximally suppressed over the range of eigenvalues while satisfying a critical normalization condition.

## Error Representation and the Minimax Problem

Since the eigenvalues of $P_k(A)$ are $P_k(\lambda_i)$ for eigenvalues $\lambda_i$ of $A$, the convergence rate (determined by the spectral radius $\rho(P_k(A))$) is minimized by solving the following **minimax problem**:

$$\min_{P_k \in \Pi_k, P_k(0)=1} \max_{\lambda \in [a, b]} |P_k(\lambda)|$$

where $\Pi_k$ is the set of polynomials of degree at most $k$, and $[a, b]$ is the interval containing the eigenvalues of $A$, with $0 < a \le b$.

### Deriving the Optimal Polynomial $P_k(\lambda)$

The function that solves this specific minimax problem is a suitably scaled and shifted Chebyshev polynomial of the first kind, $T_k(z)$.

**Transformation to the Standard Interval**

The standard Chebyshev polynomial $T_k(z)$ minimizes its maximum magnitude over the interval $z \in [-1, 1]$. To map the eigenvalue interval $\lambda \in [a, b]$ to this standard interval, a linear transformation is required:

$$z(\lambda) = \frac{a + b - 2\lambda}{b - a}$$

This linear transformation maps $\lambda \in [a, b]$ to $z \in [-1, 1]$, reversing the order of the interval.

This transformation ensures that:
*   $\lambda = a \implies z = 1$
*   $\lambda = b \implies z = -1$

**Normalization**

The optimal polynomial $P_k(\lambda)$ is then constructed using $T_k(z(\lambda))$ but requires normalization to enforce the constraint $P_k(0)=1$. This is achieved by dividing by the polynomial's value at the origin:

$$P_k(\lambda) = \frac{T_k(z(\lambda))}{T_k(z(0))}$$

### Calculating the Exact Bound

The convergence bound arises from evaluating the magnitude of the optimal polynomial at its maximal points over the eigenvalue interval $[a, b]$.

**Evaluating the Maximum Modulus of the Numerator**

For $\lambda \in [a, b]$, the transformed variable $z(\lambda)$ lies in $[-1, 1]$. In this range, the Chebyshev polynomial's defining property ensures that its maximum modulus is unity:

$$\max_{\lambda \in [a, b]} |T_k(z(\lambda))| = \max_{z \in [-1, 1]} |T_k(z)| = 1$$

**Evaluating the Normalization Constant**

$$z(0) = \frac{b+a - 2(0)}{b-a} = \frac{b+a}{b-a}\,.$$

Since $a>0$, $b>a$, the argument $z(0) = \frac{b+a}{b-a}$ is greater than $1$. Substituting these results into the definition of the optimal polynomial gives the final exact bound for the maximum error reduction:

$$\min_{P_k \in \Pi_k, P_k(0)=1} \max_{a \le \lambda \le b} |P_k(\lambda)| = \frac{\max_{\lambda \in [a, b]} |T_k(z(\lambda))|}{|T_k(z(0))|} = \frac{1}{\Big|T_k(\frac{b+a}{b-a})\Big|}$$

### Asymptotic Rate

**Error Norm Application**

Recall that:

$$\mathbf{e}^{(k)} = P_k(A) \mathbf{e}^{(0)} = X P_k(\Lambda) X^{-1} \mathbf{e}^{(0)}$$

Assume that $A$ is symmetric positive definite (SPD), so it has an orthonormal eigenvector basis. Then:

$$\|\mathbf{e}^{(k)}\|_2 \le \|P_k(A)\|_2 \|\mathbf{e}^{(0)}\|_2 = \max_{a \le \lambda \le b} |P_k(\lambda)| \; \|\mathbf{e}^{(0)}\|_2$$

**Asymptotic Approximation and $\sigma$**

For large $k$, the Chebyshev polynomial $T_k(z)$ for $z > 1$ is approximated by its exponential form (hyperbolic cosine):

$$T_k(z) \approx \frac{1}{2} (z + \sqrt{z^2 - 1})^k$$

When $z_0 = \frac{b+a}{b-a}$, the term $z_0 + \sqrt{z_0^2 - 1}$ simplifies, defining the asymptotic reduction rate $\sigma$:

$$\|\mathbf{e}^{(k)}\|_2 \lessapprox 2 \cdot \left(\frac{1 - \sqrt{\frac{a}{b}}}{1 + \sqrt{\frac{a}{b}}}\right)^k \|\mathbf{e}^{(0)}\|_2 = 2 \sigma^k \|\mathbf{e}^{(0)}\|_2$$

where the reduction rate is defined as:

$$\sigma = \frac{1 - \sqrt{a/b}}{1 + \sqrt{a/b}}$$

This reduction rate $\sigma$ shows that the number of iterations required is proportional to $\ln(\epsilon) / \ln(\sigma)$, where $\epsilon$ is the desired error tolerance. Since $\sigma$ is approximately $1 - 2/\sqrt{\kappa}$ when the condition number $\kappa = b/a$ is large, the convergence rate is proportional to $\sqrt{\kappa}$, which explains the superior efficiency of Chebyshev iteration compared to linear convergence methods proportional to $\kappa$.

The convergence analysis hinges on finding the polynomial roots (the iteration parameters $\tau_l$) that distribute the suppression across the eigenvalue interval optimally, minimizing the largest residual component left after $k$ steps. The Chebyshev iteration achieves this by having the roots of $P_k(\lambda)$ correspond to the specific **zeros of $T_k(z)$** transformed back to the $\lambda$ axis.

## Limitations of the 1st Order Chebyshev Iterative Method

The first-order Chebyshev iterative method is defined by the one-step recurrence relation: 

$$x^{(k+1)} = x^{(k)} - \tau_k(Ax^{(k)} - b).$$

The limitations associated with this method are:

*   **A Priori Selection of Steps:** The method requires the user to choose the total number of iteration steps $k$ beforehand. This must be done using estimation formulas based on the desired accuracy and the eigenvalue bounds.
*   **Numerical Instability:** The major disadvantage is that the method can be **numerically unstable**. This instability arises because for certain indices $k$, the iteration matrices $I - \tau_k A$ can have a spectral radius much larger than $1$.
*   **Parameter Ordering Requirement:** To mitigate the instability, the iteration parameters $\tau_k$ must be selected and used in a specific order (so that large values of $\tau_k$ are balanced by small values).

### How to Build the 2nd Order Chebyshev Iterative Method

The 2nd order Chebyshev iterative method is developed to eliminate the disadvantages of the 1st order method—namely, the requirement to choose the total number of iteration steps $k$ beforehand and the potential numerical instability. The two-step version is **unconditionally stable**.

This method utilizes a two-step recurrence relation, often implemented with a preconditioning matrix $C$ (assuming $C$ and $A$ are symmetric positive definite, s.p.d.):

$$x^{(k+1)} = \alpha_k x^{(k)} + (1 - \alpha_k) x^{(k-1)} - \beta_k r^{(k)}, \quad k= 1, 2, \dots$$

where the iteration step starts with $x^{(1)} = x^{(0)} - \beta_0 r^{(0)}$, and $r^{(k)}$ is the preconditioned residual:

$$r^{(k)} = Ax^{(k)} - b$$

The convergence analysis relies on tracking the error vector $e^{(k)} = x - x^{(k)}$, which is related to the initial error $e^{(0)}$ by the error polynomial $P_k(A)$ such that $e^{(k)} = P_k(A)e^{(0)}$. The construction relates the error polynomials $P_k(A)$ generated by this recurrence to the standard three-term recurrence relation of the Chebyshev polynomials of the first kind:

$$T_{k+1}(z) = 2zT_k(z) - T_{k-1}(z), \quad k = 1, 2, \dots$$

To ensure the recursion for the error vector $e^{(k)}$ follows the structure derived from the normalized Chebyshev polynomial recurrence, the parameters $\alpha_k$ and $\beta_k$ must be chosen such that the error polynomial $P_k(\lambda)$ matches $T_k(z(\lambda))/T_k(z(0))$.

The choice of parameters relies on estimating the bounds of the eigenvalues of the preconditioned matrix $A$. Let $a$ and $b$ be the lower and upper bounds, respectively, of the eigenvalues $\lambda_j$ of $A$, such that $0 < a < \lambda_j < b$.

The optimal parameters are defined as:

$$
\alpha_k = \frac{a+b}{2}\,\beta_k, \qquad
\beta_k^{-1} = \frac{a+b}{2} - \left(\frac{b-a}{4}\right)^{\!2}\,\beta_{k-1}, \quad k \ge 2,
$$

with the initialization
$$
\beta_0 = \frac{4}{a+b}.
$$

<!-- This matches Algorithm 1 (“Three-term Chebyshev iteration”) in *Chebyshev Iteration.pdf* after the symbol map $\mu = \tfrac{a+b}{2}$ (center) and $\delta = \tfrac{b-a}{2}$ (half-width), using the identification $\beta_k = -1/c_k$ and $\alpha_k = \mu\,\beta_k$.

There is a special first step in Algorithm 1; using that initialization gives

$$
\beta_1^{-1} = \frac{a+b}{2} - \frac{(b-a)^2}{4(a+b)}.
$$

If instead one plugs $k=1$ into the recurrence above, one obtains
$\beta_1^{-1} = \tfrac{a+b}{2} - \tfrac{(b-a)^2}{8(a+b)}$,
which differs by a factor of 2. This difference is due to the special $k=1$ initialization in Algorithm 1; from $k\ge 2$ the two formulations are identical. -->

The primary advantage of this second-order method is that **it gives the optimal result for every step $k$** (apart from rounding errors), eliminating the need to choose the total number of iterations $p$ a priori.

### Advantages of the 2nd Order Chebyshev Iterative Method

The second-order method offers distinct advantages over its first-order counterpart:

*   **Numerical Stability:** The two-step version of the Chebyshev iterative method is **numerically stable**.
*   **Optimal Performance at Every Step:** Unlike the one-step method which optimizes the result only for a predetermined step $k$, the two-step method is designed so that, apart from rounding errors, it gives the **optimal result for every step $k$**.
*   **Elimination of A Priori Step Count:** The method eliminates the requirement to choose the number of iteration steps $k$ *a priori*.

Both the first-order (when stable) and the second-order methods share the optimal asymptotic rate of convergence, which increases only as the square root of the condition number of $A$.

## Requirements and Sensitivity to Eigenvalue Bounds

The Chebyshev iteration method has **excellent convergence properties** because the number of iterations required for a constant error reduction is proportional to $\sqrt{\kappa}$, which is a vast improvement over basic stationary methods (like Jacobi or Gauss-Seidel) where the number of iterations is typically proportional to $\kappa$.

However, this method **requires precise knowledge of the interval that contains the eigenvalues** of the iteration matrix.

1.  **Real and Positive/Negative Eigenvalues:** The basic method is derived assuming the eigenvalues are real and positive (or all negative, allowing for a transformation). If the eigenvalues are complex, the method is still applicable, provided they are contained within an ellipse in the right or left half of the complex plane, symmetric with respect to the real axis.
2.  **Sensitivity to the Lower Bound:** For the optimal convergence rate derived above, the method critically depends on the ratio $b/a$, where $a$ is the smallest positive eigenvalue.
    *   In many problems, particularly those arising from discretized PDEs, the **smallest eigenvalues ($a$) can be very close to 0** (e.g., $O(h^2)$ in many finite difference schemes).
    *   Determining this lower bound $a$ accurately can be **difficult and computationally expensive**.
    *   If the estimate for the lower bound $a$ is slightly too small, the estimated condition number $\kappa = b/a$ becomes inflated, leading to a much slower predicted rate of convergence. If the estimate for $a$ is too large (i.e., outside the true interval), the method may diverge.
    *   The convergence rate curve can be **extremely sensitive near the optimal value** of the parameters, highlighting the danger of slightly inaccurate eigenvalue estimates.

## Pros, Cons, and Applications

The method's reliance on precise spectral bounds means that despite its theoretical efficiency, it is often avoided in favor of **parameter-free methods** like the Conjugate Gradient (CG) method, which typically converges even faster in a certain norm without requiring eigenvalue estimates.

The following table summarizes the key properties, advantages (pros), and disadvantages (cons) of the 1st order (one-step) and 2nd order (two-step) Chebyshev iterative methods:

| Feature | 1st Order Chebyshev Iterative Method (One-Step) | 2nd Order Chebyshev Iterative Method (Two-Step) |
| :--- | :--- | :--- |
| **Recurrence Relation** | Uses a one-step recurrence relation: $x^{l+1} = x^l - \tau_l(Ax^l - b)$. | Uses a two-step recurrence relation: $x^{l+1} = \alpha_l x^l + (1 - \alpha_l) x^{l-1} - \beta_l r^l$. |
| **Primary Pros** | Parameters are computed from a **closed expression**—that is, the one-step method's parameters are explicit functions of $a$, $b$, and $k$. Shares the **optimal asymptotic rate of convergence** (proportional to the square root of the condition number, $O(\sqrt{\kappa})$). | **Numerically stable**. Gives the **optimal result for every step $l$**, apart from rounding errors, unlike the one-step method which is only optimized for a pre-chosen final step $k$. Eliminates the requirement to choose the number of iteration steps $k$ *a priori*. |
| **Primary Cons** | Major disadvantage is potential **numerical instability**. Instability occurs because, for certain iteration indices $l$, the iteration matrices can have a spectral radius much larger than 1. Requires the user to choose the total number of iteration steps, $k$, **a priori**. The iteration parameters $\tau_l$ must be selected and used in a specific order (large values balanced by small values) to mitigate instability. | Parameters are derived from **recursions**. Requires more memory and complexity due to the two-step dependence. |

**Shared Requirements**: Requires estimates of extreme eigenvalues for optimal performance. Convergence rate curve can be **extremely sensitive near the optimal value** of the parameters if eigenvalue estimates are inaccurate.

**Applications and Use Cases:**

1.  **Preconditioning:** The method is applicable when combined with a preconditioning matrix $M$ (such as SSOR) to solve the system $M^{-1}A\mathbf{x} = M^{-1}\mathbf{b}$, provided $M^{-1}A$ has positive eigenvalues. This combination (SSOR with Chebyshev acceleration) can substantially reduce the iteration count.
2.  **Other Specific Matrix Classes:** Chebyshev methods can be extended to systems where the eigenvalues lie in two disjoint intervals or ellipses (such as indefinite symmetric matrices after transformation).