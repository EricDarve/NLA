# Eigenvalue Computation

Eigenvalues and eigenvectors are fundamental concepts in linear algebra. They provide deep insight into the behavior of a matrix $A$ by revealing the directions $v$ in which the matrix acts as a simple scalar $\lambda$, such that $A v = \lambda v$.

This single property unlocks a wide range of powerful applications:

* **Diagonalization:** For a diagonalizable matrix, the eigendecomposition $A = X \Lambda X^{-1}$ provides a "canonical" representation of the linear transformation, simplifying its structure.
* **Matrix Operations:** This decomposition makes complex operations computationally trivial. It turns matrix powers into scalar powers ($A^k = X \Lambda^k X^{-1}$) and enables the general computation of matrix functions:

 $$f(A) = X f(\Lambda) X^{-1}.$$
 
* **Dynamic Systems:** Eigenvalues are essential for studying time-evolving systems. The solution to a system of linear differential equations, $\frac{dx}{dt} = M x$, is given by 

$$x(t) = \exp(Mt) x_0.$$

The stability and long-term behavior of this solution are dictated entirely by the eigenvalues of $M$, as the solution depends on terms of the form $\exp(\lambda_i t)$.

Given their critical importance in physics, engineering, data analysis, and countless other scientific disciplines, the ability to *compute* eigenvalues and eigenvectors is a core practical problem. The theoretical definitions, such as solving the characteristic polynomial $\det(A - \lambda I) = 0$, do not lead to stable or efficient algorithms for general matrices.

In this section, we shift our focus from the *theory* of eigenvalues to the *practice* of their numerical computation. We will develop the robust and scalable algorithms that are used in modern software to solve eigenvalue problems efficiently and accurately.

## The Challenge: No Direct Solution

Unlike other problems we have studied, such as solving $Ax = b$ using $LU$ factorization, there is no "direct" method or finite formula to compute the eigenvalues of a general $n \times n$ matrix for $n \ge 5$.

The root of this difficulty lies in a deep, fundamental connection between finding eigenvalues and finding the roots of a polynomial.

## Equivalence to Polynomial Root-Finding

We know that the eigenvalues of $A$ are the roots of the **characteristic polynomial**, $p(\lambda) = \det(A - \lambda I)$. This shows that finding eigenvalues is *at least as hard* as finding polynomial roots.

It turns out these two problems are mathematically **equivalent**. We can also prove that finding the roots of *any* polynomial $p(x)$ is equivalent to finding the eigenvalues of a specific matrix.

Consider a general monic polynomial of degree $n$:

$$
p(x) = x^n + a_{n-1} x^{n-1} + \cdots + a_1 x + a_0
$$

We can construct a special matrix $A$, known as the **companion matrix** of the polynomial $p(x)$:

$$
A =
\begin{pmatrix}
    0 & 1 & & & \\
    & 0 & 1 & & \\
    & & \ddots & \ddots & \\
    & & & 0 & 1 \\
    -a_0 & -a_1 & \cdots & -a_{n-2} & -a_{n-1}
\end{pmatrix}
$$

````{prf:theorem} Eigenvalues of the Companion Matrix
:label: thm:companion_matrix
The eigenvalues of this matrix $A$ are *exactly* the roots of the polynomial $p(x)$.
````

````{prf:proof} Let $z \in \mathbb{C}$ be an eigenvalue and $u = (u_1, \ldots, u_n)^T$ be the corresponding eigenvector, so $A u = z u$. Let's write out the first $n-1$ equations from this system:

$$
\begin{aligned}
u_2 &= z u_1 \\
u_3 &= z u_2 = z^2 u_1 \\
&\vdots \\
u_n &= z u_{n-1} = z^{n-1} u_1
\end{aligned}
$$

This defines the structure of the eigenvector. If we set $u_1 = 1$, the eigenvector is $u = (1, z, z^2, \ldots, z^{n-1})^T$.

Now, we look at the $n$-th and final equation from $A u = z u$, which corresponds to the last row of $A$:

$$
-a_0 u_1 - a_1 u_2 - \cdots - a_{n-1} u_n = z u_n
$$

Substituting $u_k = z^{k-1}$:

$$
-a_0 - a_1 z - \cdots - a_{n-1} z^{n-1} = z (z^{n-1})
$$

Finally, moving all terms to one side, we get:

$$
z^n + a_{n-1} z^{n-1} + \cdots + a_1 z + a_0 = 0
$$

This is precisely $p(z) = 0$. Therefore, $z$ is an eigenvalue of $A$ if and only if $z$ is a root of $p(x)$.
````

## Iterative Methods are a Necessity

This equivalence is profound. The **Abel-Ruffini theorem**, a landmark result from the 19th century, proves that there is **no general formula in radicals** (i.e., no equivalent of the quadratic formula) for the roots of polynomials of degree 5 or higher.

Because of this:

1.  If a direct, exact algorithm for eigenvalues existed, we could use it to find the roots of any polynomial exactly, via its companion matrix.
2.  The Abel-Ruffini theorem tells us such a root-finding algorithm is impossible.
3.  Therefore, **no direct, exact algorithm for eigenvalues can exist** for general $n \times n$ matrices.

This forces us to use a different approach. All general-purpose eigenvalue solvers are, by necessity, **iterative**. They do not compute the exact eigenvalues in a finite number of steps. Instead, they generate a sequence of approximate solutions that converge to the true eigenvalues.

Fortunately, the algorithms we will study are so effective that this is not a practical limitation. They converge extremely fast (often quadratically or cubically) and are numerically stable, allowing us to compute eigenvalues to machine accuracy. For all practical purposes, the results are indistinguishable from an "exact" solution.