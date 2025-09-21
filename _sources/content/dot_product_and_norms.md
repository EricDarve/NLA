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
- Any vector $v \in \mathbb{R}^n$ can be uniquely written as the sum of a vector in $S$ and a vector in $S^\perp$.
- The dimensions satisfy:

$$
\dim(S) + \dim(S^\perp) = n.
$$

## Vector Norms

The dot product allows us to define the **Euclidean norm** (or **2-norm**) of a vector $x$:

$$
\|x\|_2 = \sqrt{x^T x} = \left( \sum_{i=1}^n x_i^2 \right)^{1/2}.
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
\|x\|_2 = \left( \sum_{i=1}^n x_i^2 \right)^{1/2} = \sqrt{x^T x}.
$$

3. **Infinity norm** (max norm):

$$
\|x\|_\infty = \max_{1 \le i \le n} |x_i|.
$$

4. **p-norm** (for $p > 1$):

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

The **unit ball** of a norm $\|\cdot\|$ is the set:

$$
\{ x \in \mathbb{R}^n \;|\; \|x\| = 1 \}.
$$

- For the 2-norm, the unit ball is a sphere (circle in $\mathbb{R}^2$).
- For the 1-norm, the unit ball in $\mathbb{R}^2$ is a diamond shape.
- For the infinity norm, the unit ball is a square in $\mathbb{R}^2$.

As $p \to \infty$, the $p$-norm unit ball transitions from “diamond-like” (near $p=1$) to “square-like” (as $p \to \infty$).

## Cauchy–Schwarz and Hölder Inequalities

**Hölder’s inequality:**  
For $x, y \in \mathbb{R}^n$ and $p, q > 0$ such that $\frac{1}{p} + \frac{1}{q} = 1$,

$$
|x^T y| \le \|x\|_p \, \|y\|_q.
$$

**Cauchy–Schwarz inequality:**  
Special case $p = q = 2$:

$$
|x^T y| \le \|x\|_2 \, \|y\|_2.
$$

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
| 2-norm        | $(\sum_{i=1}^n x_i^2)^{1/2}$                         | Circle                            |
| Infinity-norm | $\max_i \lvert x_i\rvert$                            | Square                            |
| p-norm        | $(\sum_{i=1}^n \lvert x_i\rvert^p)^{1/p}$           | Smooth transition from diamond to square as $p$ increases |
