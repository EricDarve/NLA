# Dot Product and Vector Norms

## Dot Product

For vectors $x, y \in \mathbb{R}^n$, the **dot product** (also called the **inner product** in this context) is defined as

$$
x^T y = x_1 y_1 + x_2 y_2 + \cdots + x_n y_n.
$$

For vectors $x, y \in \mathbb{C}^n$, we use the **conjugate transpose**:

$$
x^H y = \overline{x}_1 y_1 + \overline{x}_2 y_2 + \cdots + \overline{x}_n y_n,
$$
where $\overline{x}_i$ denotes the complex conjugate of $x_i$.

The computational cost of computing a dot product is $O(n)$.

## Orthogonality

Two vectors $x$ and $y$ are **orthogonal** if

$$
x^T y = 0
$$
(in the complex case, $x^H y = 0$).

This definition naturally extends to subspaces.

### Orthogonal Complements

For a subspace $S \subset \mathbb{R}^n$, the **orthogonal complement** $S^\perp$ is defined as:

$$
S^\perp = \{ y \in \mathbb{R}^n \;|\; \forall x \in S,\; x^T y = 0 \}.
$$

Key properties:

- $S^\perp$ is itself a subspace of $\mathbb{R}^n$.
- $S \cap S^\perp = \{0\}$, and any vector $v \in \mathbb{R}^n$ can be uniquely written as the sum of a vector in $S$ and a vector in $S^\perp$. In the notation of the previous section, $\mathbb{R}^n = S \oplus S^\perp$.
- $(S^\perp)^\perp = S$.
- The dimensions satisfy:

$$
\dim(S) + \dim(S^\perp) = n.
$$

## Vector Norms

The dot product allows us to define the **Euclidean norm** (or **2-norm**) of a vector $x$:

$$
\|x\|_2 = \sqrt{x^T x} = \left( \sum_{i=1}^n x_i^2 \right)^{1/2}
\qquad (x \in \mathbb{R}^n).
$$

For $x \in \mathbb{C}^n$ the conjugate transpose is required, since $\sum_i x_i^2$ need not be a non-negative real number:

$$
\|x\|_2 = \sqrt{x^H x} = \left( \sum_{i=1}^n |x_i|^2 \right)^{1/2}.
$$

### General Definition of a Norm

A **norm** is a function $\|\cdot\|$ mapping vectors to non-negative real numbers that satisfies:

1. **Positive definiteness:** $\|x\| = 0 \iff x = 0$.
2. **Homogeneity:** $\|\alpha x\| = |\alpha| \, \|x\|$ for all scalars $\alpha$.
3. **Triangle inequality:** $\|x + y\| \leq \|x\| + \|y\|$.

### Common Vector Norms

1. **1-norm** (Manhattan norm):

$$
\|x\|_1 = \sum_{i=1}^n |x_i|.
$$

2. **2-norm** (Euclidean norm):

$$
\|x\|_2 = \left( \sum_{i=1}^n |x_i|^2 \right)^{1/2}.
$$

3. **Infinity norm** (max norm):

$$
\|x\|_\infty = \max_{1 \le i \le n} |x_i|.
$$

4. **p-norm** (for $p \ge 1$):

$$
\|x\|_p = \left( \sum_{i=1}^n |x_i|^p \right)^{1/p}.
$$

**Convention:** In this book, $\|x\|$ means $\|x\|_2$ unless otherwise stated.

## Equivalence Between Norms

In a finite-dimensional vector space like $\mathbb{R}^n$, all norms are **equivalent**. This means that for any two vector norms $\|\cdot\|_a$ and $\|\cdot\|_b$, there exist positive constants $c_1$ and $c_2$ such that for every vector $\mathbf{x}$:

$$c_1 \|\mathbf{x}\|_b \le \|\mathbf{x}\|_a \le c_2 \|\mathbf{x}\|_b$$

This property is powerful because it guarantees that if a sequence of vectors converges in one norm, it converges in all norms. 📐

### Specific Inequalities Between 1, 2, and Infinity Norms

For any $x \in \mathbb{R}^n$:

$$
\|x\|_2 \le \|x\|_1 \le \sqrt{n} \, \|x\|_2,
$$

$$
\|x\|_\infty \le \|x\|_2 \le \sqrt{n} \, \|x\|_\infty,
$$

$$
\|x\|_\infty \le \|x\|_1 \le n \, \|x\|_\infty.
$$

These follow from the **Cauchy–Schwarz inequality** (see below) and basic properties of maxima and sums.

To illustrate why the dimension $n$ is critical, consider the vector $\mathbf{x} = [1, 1, \dots, 1]^T \in \mathbb{R}^n$.
* $\|\mathbf{x}\|_\infty = 1$
* $\|\mathbf{x}\|_2 = \sqrt{1^2 + \dots + 1^2} = \sqrt{n}$
* $\|\mathbf{x}\|_1 = 1 + \dots + 1 = n$

These values exactly match the scaling factors in the inequalities.

### General Ordering of p-Norms

For any vector $\mathbf{x}$, the value of its $p$-norm is a non-increasing function of $p$. This provides a simple and elegant ordering.

For any $p > q \ge 1$:

$$\|\mathbf{x}\|_p \le \|\mathbf{x}\|_q$$

This leads to the most frequently cited chain of inequalities:

$$\cdots \le \|\mathbf{x}\|_3 \le \|\mathbf{x}\|_2 \le \|\mathbf{x}\|_1$$

The $\infty$-norm is the limit of the $p$-norm as $p \to \infty$, making it the smallest of all $p$-norms:

$$\|\mathbf{x}\|_\infty \le \cdots \le \|\mathbf{x}\|_2 \le \|\mathbf{x}\|_1$$

## Geometric Interpretation

The **unit ball** of a norm $\|\cdot\|$ is the solid region it bounds; the **unit sphere** is that region's boundary:

$$
B = \{ x \in \mathbb{R}^n \;|\; \|x\| \le 1 \},
\qquad
\partial B = \{ x \in \mathbb{R}^n \;|\; \|x\| = 1 \}.
$$

In $\mathbb{R}^2$, the unit ball is a disk for the 2-norm, a diamond for the 1-norm, and a square for the infinity norm.

Since $\|x\|_p$ decreases with $p$, the constraint $\|x\|_p \le 1$ weakens as $p$ grows: the unit balls are nested and expand from the diamond at $p=1$ to the square at $p=\infty$.

## Why Different Norms Matter

- **Geometric Interpretation in Optimization:** The choice of norm is critical in optimization and machine learning, as it defines the geometry of the "constraint region" for your solution. In many problems, we minimize a loss function subject to a constraint that the solution vector's norm must be small ($
\|x\|_p \le C$). The optimal solution is often found where the level curves of the loss function first touch the boundary of this constraint region (the "norm ball").

- **Inducing Sparsity ($p \le 1$):**
    * The L1-norm ($\|x\|_1$) constraint region is shaped like a diamond (or hyper-diamond in higher dimensions), with sharp corners that lie on the axes.     
    * Because of these corners, the expanding level curves of the loss function are most likely to make contact at a point on an axis.
    * A point on an axis means that the other components of the vector are zero. This is why L1 regularization (used in LASSO) produces **sparse solutions**, which is useful for feature selection.

- **Encouraging Small, Non-Zero Values ($p \ge 2$):**
    * The L2-norm ($\|x\|_2$) constraint region is a disk (or ball), which is perfectly round and has no corners. The solution can occur anywhere on its boundary.
    * This norm penalizes large values heavily ($x_i^2$), so it tends to find solutions where all components are small and non-zero rather than forcing some to be exactly zero. This is the basis for Ridge Regression.
    * For norms with $p>2$, the penalty on large components is even more severe, further encouraging solutions where all entries have similar, non-zero magnitudes.

## Cauchy–Schwarz and Hölder Inequalities

**Hölder’s inequality:**  
For $x, y \in \mathbb{R}^n$ and conjugate exponents $1 \le p, q \le \infty$ satisfying $\frac{1}{p} + \frac{1}{q} = 1$ (with the convention $1/\infty = 0$),

$$
|x^T y| \le \|x\|_p \, \|y\|_q.
$$

The pair $p=1$, $q=\infty$ gives the frequently used bound $|x^T y| \le \|x\|_1 \, \|y\|_\infty$.

**Cauchy–Schwarz inequality:**  
Special case $p = q = 2$:

$$
|x^T y| \le \|x\|_2 \, \|y\|_2,
$$

with equality if and only if $x$ and $y$ are linearly dependent.

For nonzero $x$ and $y$, this bound is what makes the **angle** $\theta$ between them well defined:

$$
\cos \theta = \frac{x^T y}{\|x\|_2 \, \|y\|_2},
$$

since Cauchy–Schwarz forces $|\cos\theta| \le 1$.

**Why is This Bound Important?**

- **Proving Algorithm Convergence:** Many iterative algorithms in NLA work by generating a sequence of vectors that get progressively closer to a solution. The Cauchy-Schwarz inequality is often used to prove that the **error** at each step is decreasing, guaranteeing that the algorithm will eventually converge.
- **Geometric Insight:** The inequality has a clear geometric meaning. The dot product $x^T y$ is largest when the vectors are aligned ($\theta = 0$), smallest when they are opposed ($\theta = \pi$), and zero when they are orthogonal ($\theta = \pi/2$). This intuition is invaluable when developing new algorithms. 

## Application: Pythagorean Theorem

If $x, y \in \mathbb{R}^n$ are orthogonal ($x^T y = 0$), then:

$$
\|x + y\|_2^2 = \|x\|_2^2 + \|y\|_2^2.
$$
This follows immediately from expanding $\|x+y\|_2^2$ using the dot product definition.

**Summary Table:**

| Norm          | Formula                                              | Unit Ball Shape ($\mathbb{R}^2$) |
|---------------|------------------------------------------------------|-----------------------------------|
| 1-norm        | $\sum_{i=1}^n \lvert x_i\rvert$                      | Diamond                           |
| 2-norm        | $(\sum_{i=1}^n \lvert x_i\rvert^2)^{1/2}$          | Disk                              |
| Infinity-norm | $\max_i \lvert x_i\rvert$                            | Square                            |
| p-norm        | $(\sum_{i=1}^n \lvert x_i\rvert^p)^{1/p}$           | Smooth transition from diamond to square as $p$ increases |
