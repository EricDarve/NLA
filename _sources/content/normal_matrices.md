# Normal matrices

A normal matrix is unitarily diagonalizable because the algebraic condition of being normal is precisely what guarantees that the matrix has a complete set of orthogonal eigenvectors. This fundamental result is known as the **Spectral Theorem**.

## What is a Normal Matrix?

A square matrix $A$ is **normal** if it commutes with its conjugate transpose, $A^H$. The conjugate transpose is the result of taking the transpose of the matrix and then taking the complex conjugate of each entry.

The defining property of a normal matrix is:

$$A A^H = A^H A$$

While this definition seems purely algebraic, it has a profound geometric meaning: the linear transformation represented by $A$ does not "shear" space in a way that misaligns its eigenvectors.

### Examples of Normal Matrices

Several important families of matrices are normal, which is why this property is so useful.

* **Hermitian matrices**: These are the complex equivalent of real symmetric matrices and satisfy $A^H = A$. They are clearly normal since $A \cdot A = A \cdot A$. Their eigenvalues are always real.
* **Unitary matrices**: These matrices represent rotations and reflections in complex space and satisfy $Q^H Q = I$. Since $Q^{-1} = Q^H$, they also satisfy $Q Q^H = I$, so $Q Q^H = Q^H Q$.
* **Skew-Hermitian matrices**: These satisfy $A^H = -A$. They are normal because $A(-A) = (-A)A = -A^2$. Their eigenvalues are purely imaginary.
* **Diagonal matrices**: Any diagonal matrix is normal because diagonal matrices always commute with each other, and the conjugate transpose of a diagonal matrix is also diagonal.

## The Spectral Theorem for Normal Matrices

The **Spectral Theorem** provides the crucial link between the algebraic property of normality and the geometric property of its eigenvectors.

````{prf:theorem} Spectral Theorem
A complex square matrix $A$ is **unitarily diagonalizable** if and only if it is **normal**.
````

To be **unitarily diagonalizable** means that $A$ can be decomposed into the form:

$$A = Q \Lambda Q^H$$

This decomposition has a beautiful interpretation:
* $\Lambda$ is a **diagonal matrix** whose entries are the **eigenvalues** of $A$.
* $Q$ is a **unitary matrix** whose columns are the corresponding **eigenvectors** of $A$.
* The fact that $Q$ is unitary means its columns form an **orthonormal basis** for the vector space $\mathbb{C}^n$.

In short, the Spectral Theorem guarantees that for any normal matrix, we can find a set of perpendicular axes (the eigenvectors) along which the transformation acts simply as a stretch or compression (the eigenvalues).

````{prf:proof} Spectral Theorem

The proof of the Spectral Theorem has two parts.

*Unitarily Diagonalizable $\implies$ Normal.* This is the straightforward direction. Assume $A = Q \Lambda Q^H$. We check if it commutes with its conjugate transpose:

$$A^H = (Q \Lambda Q^H)^H = Q \Lambda^H Q^H$$

Now we compute the two products:
* $A A^H = (Q \Lambda Q^H)(Q \Lambda^H Q^H) = Q \Lambda (Q^H Q) \Lambda^H Q^H = Q (\Lambda \Lambda^H) Q^H$
* $A^H A = (Q \Lambda^H Q^H)(Q \Lambda Q^H) = Q \Lambda^H (Q^H Q) \Lambda Q^H = Q (\Lambda^H \Lambda) Q^H$

Since diagonal matrices always commute ($\Lambda \Lambda^H = \Lambda^H \Lambda$), the two results are identical. Thus, $A$ is normal.

**Normal $\implies$ Unitarily Diagonalizable.** This direction is more involved but reveals the core of the theorem.

1.  Start with the **Schur Decomposition**, which states that *any* square matrix $A$ can be written as $A = Q T Q^H$, where $Q$ is unitary and $T$ is upper triangular.
2.  Now, we use the fact that $A$ is normal. As shown in the first part, if $A$ is normal, then $T$ must also be normal, meaning $T T^H = T^H T$.
3.  The final step is to prove that an **upper triangular matrix that is also normal must be a diagonal matrix**. We can show this by comparing the diagonal entries of the product $T T^H$ and $T^H T$. For the first entry on the diagonal, $(1,1)$:
    * $(T^H T)_{11} = |T_{11}|^2$
    * $(T T^H)_{11} = |T_{11}|^2 + |T_{12}|^2 + \dots + |T_{1n}|^2$
4.  For these to be equal, the sum of the squares of the off-diagonal elements in the first row ($|T_{12}|^2 + \dots$) must be zero. This forces all those elements to be zero.
5.  By continuing this process down the diagonal, we can prove that all off-diagonal elements of $T$ must be zero.
6.  Therefore, $T$ is a diagonal matrix (we'll call it $\Lambda$), and the Schur decomposition becomes $A = Q \Lambda Q^H$. This shows that $A$ is unitarily diagonalizable.
````

## The Spectral Theorem for Real Symmetric Matrices

The Spectral Theorem for real symmetric matrices is a more specific and powerful version of the general theorem for normal matrices. It guarantees that any real symmetric matrix can be diagonalized by a real orthogonal matrix, which has a profound geometric meaning.

````{prf:theorem} Spectral Theorem for Real Symmetric Matrices
The theorem states that a real matrix $A$ is **orthogonally diagonalizable** if and only if it is **symmetric**.
````

This means that for any real symmetric matrix $A$ (where $A^T = A$), there exists a decomposition:

$$A = Q \Lambda Q^T$$

Here, the components have special real-valued properties:
* **$A$** is an $n \times n$ real symmetric matrix.
* **$Q$** is an $n \times n$ **real orthogonal matrix**. This means its columns are a set of $n$ mutually perpendicular unit vectors (an orthonormal basis), and its inverse is simply its transpose ($Q^{-1} = Q^T$). The columns of $Q$ are the eigenvectors of $A$.
* **$\Lambda$** is an $n \times n$ **real diagonal matrix**. Its diagonal entries are the eigenvalues of $A$.

## Key Properties of Symmetric Matrices

This special result emerges from two crucial properties that all real symmetric matrices possess.

### 1. Eigenvalues of a Symmetric Matrix are Always Real
A symmetric matrix, even if it has complex entries in theory, will always have purely real eigenvalues.

````{prf:proof} Real Eigenvalues

Start with the eigenvalue equation, $Ax = \lambda x$. Taking the conjugate transpose of both sides gives $x^H A^H = \bar{\lambda} x^H$. Since $A$ is real and symmetric, $A^H = A$. Thus, $x^H A = \bar{\lambda} x^H$. Right-multiplying by $x$ gives:

$$x^H A x = \bar{\lambda} x^H x = \bar{\lambda} \|x\|^2$$

Now, if we left-multiply the original equation $Ax = \lambda x$ by $x^H$, we get:

$$x^H A x = \lambda x^H x = \lambda \|x\|^2$$

Equating the two expressions shows that $\lambda \|x\|^2 = \bar{\lambda} \|x\|^2$. Since eigenvectors are non-zero, we can conclude that $\lambda = \bar{\lambda}$, meaning the eigenvalue $\lambda$ must be a real number.
````

### 2. Eigenvectors for Distinct Eigenvalues are Orthogonal
For a symmetric matrix, eigenvectors corresponding to different eigenvalues are always perpendicular to each other.

````{prf:proof} Orthogonal Eigenvectors

Let $\lambda_1 \neq \lambda_2$ be two distinct eigenvalues with corresponding eigenvectors $x_1$ and $x_2$. We start with the expression $\lambda_1(x_1^T x_2)$:

$$\lambda_1(x_1^T x_2) = (\lambda_1 x_1)^T x_2 = (Ax_1)^T x_2 = x_1^T A^T x_2$$

Since $A$ is symmetric ($A^T = A$), this becomes:

$$x_1^T A x_2 = x_1^T (\lambda_2 x_2) = \lambda_2 (x_1^T x_2)$$

So, we have $\lambda_1 (x_1^T x_2) = \lambda_2 (x_1^T x_2)$, which rearranges to $(\lambda_1 - \lambda_2) (x_1^T x_2) = 0$. Because the eigenvalues are distinct, $\lambda_1 - \lambda_2 \neq 0$, which forces the dot product $x_1^T x_2$ to be zero. Thus, the eigenvectors are orthogonal.
````

Even if eigenvalues are repeated, it is always possible to find an orthonormal basis for the corresponding eigenspace.

## Geometric Interpretation 📐

The decomposition $A = Q \Lambda Q^T$ describes the action of a symmetric matrix as a sequence of three simple geometric steps:

1.  **Rotation ($Q^T x$)**: The space is rotated by the orthogonal matrix $Q^T$ so that the standard axes align with the orthogonal eigenvectors of $A$.
2.  **Scaling ($\Lambda (Q^T x)$)**: In this new, rotated orientation, the transformation is a simple scaling along each axis. The scaling factors are the real eigenvalues on the diagonal of $\Lambda$. There is no "shear" or rotation in this step.
3.  **Rotation Back ($Q (\Lambda Q^T x)$)**: The space is rotated back to its original orientation by $Q$.

This means that any transformation by a symmetric matrix can be thought of as a pure stretch or compression along a set of perpendicular axes. The eigenvectors define these **principal axes**, and the eigenvalues define the amount of stretching along them.