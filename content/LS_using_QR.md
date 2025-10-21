# Solving Least-Squares using QR Factorization

While the method of normal equations is straightforward, its numerical instability, caused by squaring the condition number, makes it unsuitable for many practical problems. A more robust and widely used approach is based on the **QR factorization** of the matrix $A$. This method is generally the standard for solving dense linear least-squares problems.

## Derivation from the Orthogonality Principle

We start from the same geometric principle as before: the residual vector $r = b - Ax$ must be orthogonal to the column space of $A$.

$$
A^T(Ax - b) = 0
$$

Now, let's introduce the full QR factorization, $A = QR$, where $Q$ is an $m \times m$ orthogonal matrix and $R$ is an $m \times n$ upper trapezoidal matrix. The key insight is that the columns of $A$ and the columns of $Q$ span the same spaces. Specifically, the first $n$ columns of $Q$ form an orthonormal basis for the column space of $A$, $\text{span}(A)$.

Therefore, the orthogonality condition can be restated as requiring the residual to be orthogonal to the basis vectors in $Q$:

$$
Q^T(Ax - b) = 0
$$

Now, we substitute $A=QR$ into this equation:

$$
Q^T(QRx - b) = 0
$$

Since $Q$ is orthogonal, we know that $Q^T Q = I$. Applying this property, we get:

$$
\begin{aligned}
(Q^T Q)Rx - Q^T b &= 0 \\
IRx - Q^T b &= 0
\end{aligned}
$$

This simplifies the least-squares problem into a much nicer linear system.


### The Transformed System

The equation becomes:

$$
Rx = Q^T b
$$

This system is an $m \times n$ upper trapezoidal system. If we partition $R$ and $Q^T b$:

$$
R = \begin{bmatrix} R_1 \\ 0 \end{bmatrix}, \quad Q^T b = \begin{bmatrix} b_1 \\ b_2 \end{bmatrix}
$$

where $R_1$ is an $n \times n$ upper-triangular matrix and $b_1$ is an $n$-vector, the system becomes:

$$
\begin{bmatrix} R_1 \\ 0 \end{bmatrix} x = \begin{bmatrix} b_1 \\ b_2 \end{bmatrix}
$$

This gives us two equations:
1.  $R_1 x = b_1$
2.  $0 \cdot x = b_2$

The solution $x$ is found by solving the first equation. The norm of the residual is simply $\|b_2\|_2$.

The solution to the least-squares problem is therefore found by solving the **square, upper-triangular system**:

```{math}
:label: eq-qr-system
R_1 x = b_1
```

This is easily solved using **back-substitution**.

## Algorithm Summary

The QR method for solving $\|Ax - b\|_2$ can be summarized in three steps:

1.  **Factorize:** Compute the QR factorization of $A$, yielding $Q$ and $R$.
2.  **Transform:** Compute the vector $Q^T b$.
3.  **Solve:** Solve the $n \times n$ upper-triangular system $R_1 x = b_1$ for $x$ using back-substitution, where $R_1$ and $b_1$ are the top $n$ rows of $R$ and $Q^T b$, respectively.

## Analysis and Advantages

### Numerical Stability

This is the **primary advantage** of the QR method. By avoiding the formation of the Gram matrix $A^T A$, we work with a system whose condition number is $\kappa_2(R_1) = \kappa_2(A)$, not $\kappa_2(A)^2$.

$$
\kappa_2(\text{Normal Equations}) = \kappa_2(A)^2 \quad \text{vs.} \quad \kappa_2(\text{QR Method}) = \kappa_2(A)
$$

This dramatic improvement in conditioning means the solution is far less sensitive to perturbations and rounding errors. It's the reason QR is the workhorse for least-squares problems. ✅

### Computational Cost

The computational cost is dominated by the QR factorization step (e.g., using Householder transformations), which takes approximately $2mn^2 - \frac{2}{3}n^3$ flops. The total cost is $O(mn^2)$, which is comparable to the normal equations method.

## Derivation using the Normal Equations

For this derivation, we first need to discuss the full vs thin QR factorization. The previous derivation used the full QR factorization. But now it will be clearer if we use the thin QR factorization.

### Full vs. Thin QR and the Invertible Factor

### 1. Full QR Factorization

For an $m \times n$ matrix $A$ (with $m \ge n$), the **full** QR factorization is:

$$
A_{m \times n} = Q_{m \times m} R_{m \times n}
$$

* $Q$ is a **square** $m \times m$ orthogonal matrix.
* $R$ is a **rectangular** $m \times n$ upper trapezoidal matrix. It has the structure:

    $$
    R = \begin{bmatrix} R_1 \\ 0 \end{bmatrix}
    $$

    where $R_1$ is a **square** $n \times n$ upper-triangular matrix, and the lower block is an $(m-n) \times n$ matrix of zeros.

In this form, $R$ is rectangular and **not invertible**. The invertible part is its submatrix $R_1$.

### 2. Reduced (or Thin) QR Factorization

The **reduced** QR factorization is a more compact version that is often more useful in practice:

$$
A_{m \times n} = Q_{1, m \times n} R_{1, n \times n}
$$

*   $Q_1$ is a **rectangular** $m \times n$ matrix with orthonormal columns. (It's the first $n$ columns of the full $Q$).
*   $R_1$ is a **square** $n \times n$ upper-triangular matrix. (It is identical to the $R_1$ submatrix from the full QR).

The key takeaway is that $R_1$ from the reduced factorization **is square and invertible** (assuming $A$ has full column rank).

## Derivation using Reduced QR

Let's re-derive the solution starting from the normal equations using the **reduced QR** factorization, $A = Q_1 R_1$.

1.  **Start with the normal equations:**

    $$
    A^T A x = A^T b
    $$

2.  **Substitute $A = Q_1 R_1$:**

    $$
    \begin{aligned}
    (Q_1 R_1)^T (Q_1 R_1) x &= (Q_1 R_1)^T b \\
    (R_1^T Q_1^T) (Q_1 R_1) x &= R_1^T Q_1^T b
    \end{aligned}
    $$

3.  **Simplify using orthogonality:**
    The columns of $Q_1$ are orthonormal, so $Q_1^T Q_1 = I_n$ (the $n \times n$ identity matrix). The equation simplifies:

    $$
    \begin{aligned}
    (R_1^T I_n R_1) x &= R_1^T Q_1^T b \\
    R_1^T R_1 x &= R_1^T Q_1^T b
    \end{aligned}
    $$

4.  **Cancel the invertible term:**
    Since $R_1$ is an $n \times n$ invertible matrix (assuming full rank), its transpose $R_1^T$ is also invertible. We can now correctly left-multiply both sides by $(R_1^T)^{-1}$:

    $$
    (R_1^T)^{-1} R_1^T R_1 x = (R_1^T)^{-1} R_1^T Q_1^T b
    $$

    This leaves us with the clean, square, upper-triangular system:

    $$
    R_1 x = Q_1^T b = b_1
    $$

    This is the same system we derived earlier using the orthogonality principle (see equation {eq}`eq-qr-system`).
