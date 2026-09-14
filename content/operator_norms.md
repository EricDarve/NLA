# Operator and Matrix Norms

A **matrix norm** is a function that assigns a non-negative number to a matrix, providing a measure of its "size" or "magnitude". Norms are fundamental tools in numerical analysis for understanding the behavior of matrix operations and for quantifying errors in computations.

## The Operator Norm (Induced p-Norm)

The most common family of matrix norms are **operator norms**, also known as **induced norms**. These norms are derived directly from vector p-norms and are defined by how much a matrix can "stretch" a vector.

The operator norm is defined as the maximum possible ratio of $\|Ax\|_p$ to $\|x\|_p$ for any non-zero vector $x$. This is equivalent to finding the maximum length of the vector $Ax$ for any vector $x$ with a length of 1.

$$\| A \|_p = \sup_{x \neq 0} \frac{\|A x\|_p}{\|x\|_p} = \max_{\|x\|_p = 1} \| Ax \|_p$$

Geometrically, take $p = 2$: applying $A$ to the Euclidean unit sphere gives an ellipsoid (an ellipse in $\mathbb{R}^2$), and $\|A\|_2$ is the length of the longest vector from the origin to a point on it. For other values of $p$ the unit ball is not round, so its image is not an ellipsoid; $\|A\|_p$ is then the largest $p$-norm attained on that image.

The three most widely used operator norms are:
* **1-Norm**: The maximum absolute column sum. It is computed as $\|A\|_1 = \max_j \sum_i |a_{ij}|$.
* **$\infty$-Norm**: The maximum absolute row sum. It is computed as $\|A\|_\infty = \max_i \sum_j |a_{ij}|$.
* **2-Norm (Spectral Norm)**: This norm is equal to the largest **singular value** of the matrix, denoted $\sigma_{\max}(A)$. While it corresponds to the true geometric "maximum stretch," it is generally more computationally intensive than the 1-norm or $\infty$-norm.

Singular values are defined later in this chapter, in the section on the singular value decomposition. Until then, read $\sigma_{\max}(A)$ as the largest factor by which $A$ stretches a vector in the 2-norm.

## The Frobenius Norm

Another useful and easy-to-compute norm is the **Frobenius norm**. Unlike the operator norms, it is not induced by a vector norm. Instead, it treats the matrix as a single long vector and calculates its standard Euclidean norm (the square root of the sum of the squares of all its elements).

$$\|A\|_F = \left( \sum_{i=1}^m \sum_{j=1}^n |a_{ij}|^2 \right)^{1/2}$$

That it is not an induced norm is easy to check: every induced norm satisfies $\|I\| = 1$, whereas $\|I_n\|_F = \sqrt{n}$.

The Frobenius norm is very convenient for analyzing matrices that can be decomposed into blocks.

## Key Inequalities for Operator Norms

While many inequalities exist, they primarily fall into a few key categories that define how matrix norms behave with respect to multiplication and each other.

### Sub-multiplicative Property

This is the most important property for analyzing products of matrices. For any two matrices $A$ and $B$ whose product is defined, the operator and Frobenius norms satisfy:

$$\|AB\| \le \|A\| \|B\|$$

This inequality is crucial for bounding errors in matrix computations and for proving the convergence of iterative methods, such as in solving systems of linear equations.

### Consistency with Vector Norms

The consistency property for the operator norm connects the matrix norm back to the vector norm it's based on. For any matrix $A$ and vector $x$:

$$\|Ax\| \le \|A\| \|x\|$$

This inequality states that the norm of a matrix tells you the maximum factor by which it can stretch any vector.

## Equivalence Between Norms

Just like vector norms, all matrix norms on a finite-dimensional space are **equivalent**. This means that for any two matrix norms $\|\cdot\|_a$ and $\|\cdot\|_b$, there exist positive constants $c_1$ and $c_2$ such that:

$$c_1 \|A\|_b \le \|A\|_a \le c_2 \|A\|_b$$

The most useful specific inequalities relate the common operator norms ($1$-norm, $2$-norm, $\infty$-norm) and the Frobenius norm ($F$) for any $n \times n$ matrix $A$. 📐

* $\|A\|_2 \le \|A\|_F \le \sqrt{n}\|A\|_2$
* $\frac{1}{\sqrt{n}}\|A\|_\infty \le \|A\|_2 \le \sqrt{n}\|A\|_\infty$
* $\frac{1}{\sqrt{n}}\|A\|_1 \le \|A\|_2 \le \sqrt{n}\|A\|_1$
* $\frac{1}{n}\|A\|_\infty \le \|A\|_1 \le n\|A\|_\infty$
* $\|A\|_2 \le \sqrt{\|A\|_1 \|A\|_\infty}$, which carries no factor of $n$ and is usually the sharpest of these bounds

## Other Fundamental Properties

All matrix norms must also satisfy the standard properties of a norm:
1.  **Triangle Inequality**: $\|A+B\| \le \|A\| + \|B\|$
2.  **Absolute Homogeneity**: $\|\alpha A\| = |\alpha| \|A\|$ for any scalar $\alpha$.
3.  **Positive Definiteness**: $\|A\| \ge 0$, and $\|A\|=0$ if and only if $A=0$.

## Other Matrix Norms

While the operator and Frobenius norms are the most common, other specialized norms are valuable in various applications. We introduce two such examples: the simple but non-sub-multiplicative **max norm**, and the more general family of **Schatten norms**, which are based on a matrix's singular values.

### Max Norm

The largest element in the matrix by magnitude. Note that this norm is **not** sub-multiplicative ($\|AB\| \le \|A\|\|B\|$ does not hold).

$$\|A\|_{\text{max}} = \max_{ij} |a_{ij}|$$

Taking $A = B = \begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix}$ shows this: $\|A\|_{\text{max}} = \|B\|_{\text{max}} = 1$, while $AB = \begin{pmatrix} 2 & 2 \\ 2 & 2 \end{pmatrix}$ has $\|AB\|_{\text{max}} = 2$.

### Schatten Norms (based on singular values)

The Schatten $p$-norm is defined by applying the vector $p$-norm to the vector of the matrix's singular values, $\sigma_i$.

$$\|A\|_{S,p} = \left( \sum_{i=1}^{\min(m,n)} \sigma_i^p \right)^{1/p}$$

This family generalizes several important norms:

* **Schatten 1-Norm (Nuclear Norm or Trace Norm)**: The sum of the singular values. It's widely used in machine learning for matrix completion and rank minimization problems.

    $$\|A\|_{S,1} = \sum_{i=1}^{\min(m,n)} \sigma_i$$

* **Schatten 2-Norm**: The square root of the sum of the squares of the singular values. This is exactly the same as the **Frobenius norm**.

    $$\|A\|_{S,2} = \left( \sum_{i=1}^{\min(m,n)} \sigma_i^2 \right)^{1/2} = \|A\|_F$$

* **Schatten $\infty$-Norm**: The limit as $p \to \infty$, which is simply the largest singular value. This is exactly the same as the **spectral norm (operator 2-norm)**.

    $$\|A\|_{S,\infty} = \max_i \sigma_i = \|A\|_2$$