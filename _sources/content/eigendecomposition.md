# Eigendecomposition

## What is an Eigenvalue?

The **eigendecomposition** is a method for breaking down a square matrix ($A$) into its fundamental constituents: its eigenvalues and eigenvectors. This breakdown reveals the core properties of the linear transformation represented by the matrix, particularly how it behaves when applied repeatedly.

For any square matrix $A$, a non-zero vector $x$ is called an **eigenvector** if applying the matrix $A$ to $x$ results only in scaling $x$ by a scalar factor $\lambda$. This scalar factor $\lambda$ is known as the **eigenvalue**.

The relationship is expressed by the key equation:

$$
Ax = \lambda x
$$

Eigenvectors represent special directions in the space that are not changed by the transformation, but are only stretched, compressed, or flipped.

## Existence of Eigenvalues

````{prf:theorem} Existence of Eigenvalues
:label: thm:eigenvalue_existence
A cornerstone of linear algebra is that every $n \times n$ matrix with complex entries has at least one eigenvalue and its corresponding eigenvector.
````

````{prf:proof} Using the Characteristic Polynomial

This is the most common proof and relies on the determinant and the Fundamental Theorem of Algebra.

1.  The defining equation for an eigenvalue is $Ax = \lambda x$, which can be rewritten as:

    $$(A - \lambda I)x = 0$$

    where $I$ is the identity matrix.

2.  For this equation to have a non-zero solution for $x$, the matrix $(A - \lambda I)$ must be **singular**. A square matrix is singular if and only if its determinant is zero.

    $$\det(A - \lambda I) = 0$$

3.  The expression $p(\lambda) = \det(A - \lambda I)$ is a polynomial in the variable $\lambda$. If $A$ is an $n \times n$ matrix, this **characteristic polynomial** has degree $n$. For example, if $A = \begin{pmatrix} a & b \\ c & d \end{pmatrix}$, then:

    $$
    \det(A - \lambda I) = \det\begin{pmatrix} a-\lambda & b \\ c & d-\lambda \end{pmatrix} = (a-\lambda)(d-\lambda) - bc = \lambda^2 - (a+d)\lambda + (ad-bc)
    $$

4.  The **Fundamental Theorem of Algebra** states that any non-constant polynomial with complex coefficients has at least one root in the complex numbers.

5.  Since the characteristic polynomial $p(\lambda)$ is a polynomial of degree $n \ge 1$, it must have at least one complex root. Any such root $\lambda_0$ satisfies $\det(A - \lambda_0 I) = 0$, which by definition makes $\lambda_0$ an eigenvalue of $A$.
````

````{prf:proof} Using Linear Dependence

This algebraic proof avoids the use of determinants.

1.  Take any non-zero vector $x \in \mathbb{C}^n$. Consider the following set of $n+1$ vectors:

    $$\{x, Ax, A^2x, \dots, A^nx\}$$

2.  Since these $n+1$ vectors reside in an $n$-dimensional space, they must be **linearly dependent**. Therefore, there exist complex scalars $c_0, c_1, \dots, c_n$, not all zero, such that:

    $$c_0 x + c_1 Ax + c_2 A^2x + \dots + c_n A^nx = 0$$

3.  This equation can be written as a polynomial in the matrix $A$ acting on the vector $x$:

    $$P(A)x = 0 \quad \text{where} \quad P(z) = c_0 + c_1 z + \dots + c_n z^n$$

4.  By the Fundamental Theorem of Algebra, we can factor the polynomial $P(z)$ in terms of its roots $\lambda_1, \dots, \lambda_k$ (where $k \le n$):

    $$P(z) = c_n(z - \lambda_1)(z - \lambda_2)\cdots(z - \lambda_k)$$

5.  Substituting the matrix $A$ back into this factored form gives:

    $$c_n(A - \lambda_1 I)(A - \lambda_2 I)\cdots(A - \lambda_k I)x = 0$$

6.  Let's read this expression from right to left. Let $v_k = x$, $v_{k-1} = (A - \lambda_k I)v_k$, $v_{k-2} = (A - \lambda_{k-1} I)v_{k-1}$, and so on. Since the final result is the zero vector, and we started with a non-zero vector $x$, there must be some step where a non-zero vector was transformed into the zero vector.

7.  That is, there must be some index $i$ and some non-zero vector $v_i$ such that $(A - \lambda_i I)v_i = 0$. This is the definition of $\lambda_i$ being an eigenvalue with corresponding eigenvector $v_i$. If no such step existed, all the matrices $(A-\lambda_i I)$ would be invertible, and their product would also be invertible, which contradicts the fact that it maps the non-zero vector $x$ to zero.
````

## The (Complex) Schur Decomposition

The existence of at least one eigenvalue leads directly to the **Schur Decomposition**, a powerful result that exists for *every* square matrix. Unlike the eigendecomposition, the Schur decomposition is guaranteed to exist regardless of whether the matrix is diagonalizable.

The Schur decomposition represents the matrix $A$ in the form:

$$
A = Q T Q^{-1}
$$

### Components of the Schur Decomposition

1.  **$T$**: This is an **upper triangular matrix**.
2.  **$Q$**: This is a **unitary matrix**. A unitary matrix means that its inverse is equal to its conjugate transpose ($Q^{-1} = Q^H$). The columns of $Q$ form an orthonormal basis for the space $\mathbb{C}^n$.

This decomposition establishes that $A$ is unitarily similar to an upper triangular matrix $T$. The use of unitary transformations ($Q$) ensures numerical stability, as these transformations preserve lengths and do not amplify errors. This makes the Schur decomposition one of the most reliable methods for computing eigenvalues in practice.

### Relationship to Eigenvalues

The **eigenvalues of $A$ are exactly the diagonal entries of the triangular matrix $T$**.

### Proof of Existence

The existence of one eigenvalue for a square matrix $A$ in the complex domain leads to the decomposition $A = Q T Q^{-1}$.

````{prf:proof} Existence of the Schur Decomposition

The existence of the Schur decomposition can be proven by induction on the size of the matrix, $n$.

* **Base Case (n=1)**: The result is trivially true, as a $1 \times 1$ matrix is already in upper triangular form.

* **Inductive Step**: Assume the decomposition exists for all matrices of size $(n-1) \times (n-1)$.
    1.  We know every $n \times n$ matrix $A$ has at least one eigenvalue, $\lambda_1$, with a corresponding unit eigenvector, $x_1$.
    2.  We can extend this vector $x_1$ to form an orthonormal basis for $\mathbb{C}^n$. Let the matrix with these basis vectors as its columns be $Q_1 = [x_1, v_2, \dots, v_n]$. By construction, $Q_1$ is unitary.
    3.  Consider the transformation $Q_1^H A Q_1$. Its first column is $Q_1^H A x_1 = Q_1^H (\lambda_1 x_1) = \lambda_1 (Q_1^H x_1)$. Since $x_1$ is the first column of $Q_1$, $Q_1^H x_1$ is the vector $(1, 0, \dots, 0)^T$.
    4.  This means the transformed matrix has a block structure:

        $$
        Q_1^H A Q_1 = \begin{pmatrix}
        \lambda_1 & \mathbf{w}^H \\
        \mathbf{0} & A_2
        \end{pmatrix}
        $$

        where $A_2$ is an $(n-1) \times (n-1)$ matrix.
    5.  By our induction hypothesis, the smaller matrix $A_2$ has its own Schur decomposition, $A_2 = Q_2 T_2 Q_2^H$.
    6.  We can now "embed" this smaller decomposition back into the larger matrix structure. Define a new unitary matrix $Q$ as $Q = Q_1 \begin{pmatrix} 1 & \mathbf{0}^T \\ \mathbf{0} & Q_2 \end{pmatrix}$.
    7.  This choice of $Q$ transforms $A$ into an upper triangular matrix, completing the induction:

        $$
        Q^H A Q = \begin{pmatrix}
        \lambda_1 & \mathbf{w}^H Q_2 \\
        \mathbf{0} & Q_2^H A_2 Q_2
        \end{pmatrix} = \begin{pmatrix}
        \lambda_1 & \mathbf{w}^H Q_2 \\
        \mathbf{0} & T_2
        \end{pmatrix}
        $$

````

### Why It's So Important for Computation

* **It's How We Actually Find Eigenvalues:** State-of-the-art numerical methods (like the famous QR algorithm) do not try to compute the eigendecomposition directly. Instead, they are designed to reliably and stably compute the **Schur decomposition**.
* **A Stable Pathway:** Once the Schur form is found, the eigenvalues can be simply read off the diagonal of $T$. If the eigenvectors are also needed, they can be calculated from $Q$ and $T$. The Schur decomposition is the practical and stable **pathway** to finding eigenvalues and eigenvectors.
* **The Power of $Q$:** The reason this approach works so well is the presence of the orthogonal matrix $Q$. Because orthogonal transformations preserve length, they don't amplify rounding errors during computation. This makes the entire process **numerically stable** and trustworthy.

## The Real Schur Decomposition

When dealing exclusively with real matrices, a complication arises if the matrix $A$ possesses complex eigenvalues. Since complex eigenvalues of real matrices always occur in conjugate pairs, it is impossible for the matrix $T$ in the decomposition $A = Q T Q^T$ to be both real and triangular, because its diagonal must contain the complex eigenvalues.

The **real Schur decomposition** provides a mechanism to keep the entire decomposition within the real numbers.

The real Schur form, often denoted $S$, is achieved through a real orthogonal matrix $Q$ (where $Q^T = Q^{-1}$) such that $A = Q S Q^T$. The resulting matrix $S$ is **block upper triangular**.

This form replaces the complex eigenvalues on the diagonal with $2 \times 2$ blocks that contain the information needed for the complex conjugate pairs.

For example, a $5 \times 5$ real matrix with one real eigenvalue ($\lambda_1$) and two pairs of complex conjugate eigenvalues ($a \pm ib$ and $c \pm id$) would have the following real Schur form:

$$
S = \begin{pmatrix} \lambda_1 & * & * & * & * \\ 0 & a & b & * & * \\ 0 & -b & a & * & * \\ 0 & 0 & 0 & c & d \\ 0 & 0 & 0 & -d & c \end{pmatrix}
$$

## Proof of Existence

````{prf:proof} Real Schur Decomposition

**Proof by Induction.** We prove the theorem by induction on the dimension $n$ of the matrix.

### Base Case (n=1)
For $n=1$, the matrix $A = [a]$ is a scalar. The decomposition is trivial: $A = [1][a][1]^T$, where $Q = [1]$ is orthogonal and $S = [a]$ is quasi-upper triangular. The theorem holds.

### Inductive Hypothesis
Assume that the Real Schur Decomposition exists for all real matrices of size $(n-1) \times (n-1)$ or smaller.

### Inductive Step
Let $A$ be an $n \times n$ real matrix. The characteristic polynomial of $A$ has real coefficients. Therefore, its roots (the eigenvalues of $A$) are either real or appear in complex conjugate pairs. We consider two cases.

#### Case 1: $A$ has a real eigenvalue
Let $\lambda_1$ be a real eigenvalue of $A$, and let $q_1$ be a corresponding real eigenvector with unit norm ($\|q_1\|_2 = 1$).

1.  We can extend this single vector $q_1$ to form a full orthonormal basis for $\mathbb{R}^n$. Let the matrix containing these basis vectors as its columns be $Q_1 = [q_1, v_2, \dots, v_n]$. By construction, $Q_1$ is a real orthogonal matrix.
2.  Now, consider the similarity transformation $Q_1^T A Q_1$. The first column of this new matrix is:

    $$
    Q_1^T A q_1 = Q_1^T (\lambda_1 q_1) = \lambda_1 (Q_1^T q_1) = \lambda_1 \begin{pmatrix} 1 \\ 0 \\ \vdots \\ 0 \end{pmatrix}
    $$

3.  This means the transformation results in a block upper triangular matrix:

    $$
    Q_1^T A Q_1 = \begin{pmatrix}
    \lambda_1 & \mathbf{w}^T \\
    \mathbf{0} & A_2
    \end{pmatrix}
    $$

    where $A_2$ is a real matrix of size $(n-1) \times (n-1)$.
4.  By the inductive hypothesis, there exists a real orthogonal matrix $Q_2$ and a quasi-upper triangular matrix $S_2$ such that $A_2 = Q_2 S_2 Q_2^T$.
5.  Define a new $n \times n$ orthogonal matrix $Q$:

    $$
    Q = Q_1 \begin{pmatrix} 1 & \mathbf{0}^T \\ \mathbf{0} & Q_2 \end{pmatrix}
    $$

    This matrix is orthogonal as it is the product of two orthogonal matrices.
6.  Finally, applying this transformation to $A$ yields the desired form:

    $$
    Q^T A Q = \begin{pmatrix} 1 & \mathbf{0}^T \\ \mathbf{0} & Q_2^T \end{pmatrix} (Q_1^T A Q_1) \begin{pmatrix} 1 & \mathbf{0}^T \\ \mathbf{0} & Q_2 \end{pmatrix} = \begin{pmatrix} \lambda_1 & \mathbf{w}^T Q_2 \\ \mathbf{0} & Q_2^T A_2 Q_2 \end{pmatrix} = \begin{pmatrix} \lambda_1 & \mathbf{w}^T Q_2 \\ \mathbf{0} & S_2 \end{pmatrix}
    $$

    This matrix, which we call $S$, is quasi-upper triangular. Thus, $A = Q S Q^T$.

#### Case 2: $A$ has no real eigenvalues
If $A$ has no real eigenvalues, it must have a pair of complex conjugate eigenvalues, $\lambda = a + ib$ and $\bar{\lambda} = a - ib$, where $b \neq 0$.

1.  Let the complex eigenvector for $\lambda$ be $z = x + iy$, where $x$ and $y$ are real vectors. From $Az = \lambda z$, we get:

    $$
    A(x+iy) = (a+ib)(x+iy) = (ax - by) + i(bx + ay)
    $$

2.  Equating the real and imaginary parts gives two real vector equations:
    * $Ax = ax - by$
    * $Ay = bx + ay$
3.  These equations show that the subspace $W = \text{span}\{x, y\}$ is a two-dimensional invariant subspace under $A$. (The vectors $x$ and $y$ must be linearly independent, otherwise $\lambda$ would be real).
4.  Use the Gram-Schmidt process (this will be covered later, please accept this result for now) to find an orthonormal basis $\{q_1, q_2\}$ for the subspace $W$. Extend this to a full orthonormal basis for $\mathbb{R}^n$ and form the real orthogonal matrix $Q_1 = [q_1, q_2, \dots, q_n]$.
5.  Because $W$ is an invariant subspace, $Aq_1$ and $Aq_2$ are linear combinations of only $q_1$ and $q_2$. This means that the transformation $Q_1^T A Q_1$ will have a block structure:

    $$
    Q_1^T A Q_1 = \begin{pmatrix}
    B & C \\
    0 & D
    \end{pmatrix}
    $$

    where $B$ is a $2 \times 2$ real matrix whose eigenvalues are the pair $\lambda, \bar{\lambda}$, and $D$ is a real matrix of size $(n-2) \times (n-2)$.
6.  By the inductive hypothesis, there is a real Schur decomposition for $D$, so $D = Q_2 S_2 Q_2^T$.
7.  Define a new orthogonal matrix $Q$:

    $$
    Q = Q_1 \begin{pmatrix} I_2 & 0 \\ 0 & Q_2 \end{pmatrix}
    $$

8.  This transformation yields the final quasi-upper triangular form:

    $$
    Q^T A Q = \begin{pmatrix} B & C Q_2 \\ 0 & S_2 \end{pmatrix}
    $$

    This matrix $S$ has a $2 \times 2$ block and a quasi-upper triangular block on its diagonal, and is therefore itself quasi-upper triangular.

Since both cases lead to the desired decomposition, the proof is complete by induction.
````

## The Eigendecomposition (Diagonalization)

The eigendecomposition is a special, highly useful case of the Schur form. It is guaranteed to exist only when the matrix $A$ is **diagonalizable**. A matrix is diagonalizable if and only if it possesses a full set of $n$ linearly independent eigenvectors.

When $A$ is diagonalizable, the upper triangular matrix $T$ becomes a diagonal matrix, commonly written as $\Lambda$. The unitary matrix $Q$ becomes $X$, a matrix whose columns are the eigenvectors of $A$.

The decomposition takes the form:

$$
A = X \Lambda X^{-1}
$$

### Components of the Eigendecomposition

1.  **$\Lambda$**: A diagonal matrix containing the eigenvalues ($\lambda_i$).
2.  **$X$**: A matrix whose columns are the $n$ linearly independent eigenvectors of $A$.

## Interpretation, Application, and Connection to Time ($A^n$)

The true utility of the eigendecomposition is revealed in the analysis of dynamical systems, which are systems that evolve over time. This conceptual link leads to the association of eigendecomposition with **"Time"**.

### Application: Calculating Powers of a Matrix

Consider a system whose state at time $k+1$, denoted $x_{k+1}$, is derived from its state at time $k$, $x_k$, by the linear rule $x_{k+1} = Ax_k$. To find the state after $k$ steps, we must calculate $x_k = A^k x_0$.

Calculating $A^k$ directly is often computationally complex. However, if $A$ is diagonalizable, the computation simplifies immensely:

```{prf:proof} Computing Matrix Powers $A^k$

If $A$ is diagonalizable such that $A = X \Lambda X^{-1}$, then calculating the $k$-th power of $A$ involves a simple simplification due to the telescoping nature of the product:

$$
A^k = (X \Lambda X^{-1})^k = (X \Lambda X^{-1})(X \Lambda X^{-1})\cdots(X \Lambda X^{-1}) = X \Lambda^k X^{-1}
$$

Since $\Lambda$ is a diagonal matrix, calculating $\Lambda^k$ only requires raising its diagonal entries (the eigenvalues) to the $k$-th power:

$$
\Lambda^k = \begin{pmatrix} \lambda_1^k & 0 & \cdots & 0 \\ 0 & \lambda_2^k & \cdots & 0 \\ \vdots & \vdots & \ddots & \vdots \\ 0 & 0 & \cdots & \lambda_n^k \end{pmatrix}
$$
```

### Interpretation: Dynamics and Stability

This simplification shows that the **long-term behavior** (the system's dynamics over "time") is governed entirely by the magnitudes of its eigenvalues, $|\lambda_i|$.

*   If the magnitude of all eigenvalues is less than 1 ($|\lambda_i| < 1$), the system is **stable**, meaning the state approaches zero as $k$ approaches infinity.
*   If the magnitude of any eigenvalue is greater than 1 ($|\lambda_i| > 1$), the system is **unstable** and will diverge.

This governing principle also applies to continuous systems described by differential equations, where the matrix exponential $e^{At}$ is computed using the eigendecomposition, showing how eigenvalues control the exponential growth or decay of the system's modes over time.