# Solving Least-Squares using QR Factorization

The **QR method** solves a least-squares problem without forming $A^TA$. It avoids squaring the condition number of the matrix used in the solve and is a standard approach for dense least-squares problems.

Throughout this page, $A\in\mathbb{R}^{m\times n}$ has full column rank, $m\ge n$, and $b\in\mathbb{R}^m$. These assumptions ensure that the least-squares minimizer is unique.

## Derivation from the Orthogonality Principle

We start from the same geometric principle as before: the residual vector $r = b - Ax$ must be orthogonal to the column space of $A$.

$$
A^T(Ax - b) = 0
$$

Introduce the full QR factorization $A=QR$, where $Q$ is an $m\times m$ orthogonal matrix and $R$ is an $m\times n$ upper trapezoidal matrix. Partition the factors as

$$
Q=\begin{bmatrix}Q_1&Q_2\end{bmatrix},
\qquad R=\begin{bmatrix}R_1\\0\end{bmatrix}.
$$

The first $n$ columns, $Q_1$, form an orthonormal basis for the column space of $A$. The remaining columns, $Q_2$, span its orthogonal complement.

The residual must be orthogonal to the columns of $Q_1$:

$$
Q_1^T(Ax-b)=0.
$$

Since $A=Q_1R_1$ and $Q_1^TQ_1=I_n$, this gives

$$
R_1x=Q_1^Tb.
$$

### The Transformed System

To identify the minimum residual as well as the solution, partition the transformed right-hand side:

$$
Q^Tb=\begin{bmatrix}b_1\\b_2\end{bmatrix},
\qquad b_1=Q_1^Tb\in\mathbb{R}^n,
\qquad b_2=Q_2^Tb\in\mathbb{R}^{m-n}.
$$

Orthogonality preserves the residual norm, so

$$
\begin{aligned}
\|Ax-b\|_2^2
&=\|Q^T(Ax-b)\|_2^2\\
&=\left\|\begin{bmatrix}R_1x-b_1\\-b_2\end{bmatrix}\right\|_2^2\\
&=\|R_1x-b_1\|_2^2+\|b_2\|_2^2.
\end{aligned}
$$

The term $\|b_2\|_2^2$ does not depend on $x$. Since $R_1$ is nonsingular, we can make the first term zero by solving the **square, upper-triangular system**:

```{math}
:label: eq-qr-system
R_1 x = b_1
```

Solve this system by **back-substitution**. The minimum residual norm is $\|b_2\|_2$, and the residual vector is $r=Q_2b_2$. We do not require $b_2=0$; a nonzero $b_2$ is precisely the part of $b$ that cannot be fitted.

## Algorithm Summary

The QR method for minimizing $\|Ax-b\|_2$ has three steps:

1.  **Factorize:** Compute Householder QR, storing the reflections and the triangular factor $R_1$.
2.  **Transform:** Apply the stored reflections to compute $Q^Tb$, without forming $Q$ explicitly.
3.  **Solve:** Solve the $n \times n$ upper-triangular system $R_1 x = b_1$ for $x$ using back-substitution, where $R_1$ and $b_1$ are the top $n$ rows of $R$ and $Q^T b$, respectively.

## Analysis and Advantages

### Numerical Stability

Householder QR has the small factorization and orthogonality errors described in the [earlier analysis](householder_reflections.md), under the arithmetic assumptions stated there. It also avoids forming the Gram matrix $A^TA$. The matrices used in the two approaches have condition numbers

$$
\kappa_2(A^TA)=\kappa_2(A)^2,
\qquad \kappa_2(R_1)=\kappa_2(A).
$$

The second identity follows because multiplying by $Q_1$ preserves lengths, so $A$ and $R_1$ have the same singular values. QR avoids the additional accuracy risk of forming $A^TA$, but does not remove the sensitivity of the original least-squares problem. A backward-stable method can still have a large forward error when the problem is sensitive to perturbations.

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
* $R$ is an $m\times n$ upper trapezoidal matrix. It has the structure:

    $$
    R = \begin{bmatrix} R_1 \\ 0 \end{bmatrix}
    $$

    where $R_1$ is a **square** $n \times n$ upper-triangular matrix, and the lower block is an $(m-n) \times n$ matrix of zeros.

When $m>n$, $R$ is rectangular and has no two-sided inverse. The square factor $R_1$ is invertible under our full-column-rank assumption. When $m=n$, $R=R_1$.

### 2. Reduced (or Thin) QR Factorization

The **reduced** QR factorization is a more compact version that is often more useful in practice:

$$
A_{m \times n} = Q_{1, m \times n} R_{1, n \times n}
$$

*   $Q_1$ is an $m\times n$ matrix with orthonormal columns, rectangular when $m>n$. It consists of the first $n$ columns of the full $Q$.
*   $R_1$ is a **square** $n \times n$ upper-triangular matrix. (It is identical to the $R_1$ submatrix from the full QR).

The key takeaway is that $R_1$ from the reduced factorization **is square and invertible** (assuming $A$ has full column rank).

## Derivation using Reduced QR

The normal equations lead to the same triangular system. Substitute $A=Q_1R_1$ into $A^TAx=A^Tb$. Since $Q_1^TQ_1=I_n$, we obtain

$$
R_1^TR_1x=R_1^TQ_1^Tb.
$$

Full column rank makes $R_1$ and its transpose invertible. Multiplying both sides by $(R_1^T)^{-1}$ gives

$$
R_1x=Q_1^Tb=b_1,
$$

which is the same system as {eq}`eq-qr-system`. This is an algebraic connection between the methods; the QR least-squares method never needs to form $A^TA$.
