# Singular Value Decomposition

The Singular Value Decomposition (SVD) is one of the most important matrix factorizations used in data science and machine learning. It is a vital tool used in methods like Principal Component Analysis (PCA) and Latent Semantic Analysis (LSA).

## Explanation of Singular Value Decomposition (SVD)

````{prf:theorem}
:label: thm:svd
The SVD expresses any $m \times n$ real matrix $A$ as the product of three specific matrices:

$$ A = U \Sigma V^T $$

*   $A$ is an $m \times n$ matrix.
*   $U$ is an orthogonal matrix ($m \times m$), whose columns are the **left singular vectors**.
*   $V$ is an orthogonal matrix ($n \times n$), whose columns are the **right singular vectors**.
*   $\Sigma$ (or $\Lambda$) is an $m \times n$ rectangular diagonal matrix. Its diagonal elements, $\sigma_i$, are real, non-negative numbers called **singular values**. These singular values are arranged in descending order: $\sigma_1 \geq \sigma_2 \geq \cdots \geq \sigma_p \geq 0$, where $p = \min(m, n)$.

````

The SVD exists for any real matrix, even if it is not a square matrix. It can be regarded as a generalization of the diagonalization of a square matrix.

## SVD and the Association with "Space" (Scales of the Matrix)

SVD is conceptually associated with **"Space"** because it describes the geometry of a transformation, unlike eigendecomposition, which is associated with "Time". SVD systematically represents the **multiple scales** involved when the matrix $A$ is applied to a vector.

Geometrically, an $m \times n$ matrix $A$ represents a linear transformation from the $n$-dimensional space $\mathbb{R}^n$ to the $m$-dimensional space $\mathbb{R}^m$. The SVD decomposes this complex transformation into three simple, sequential transformations:

1.  **Rotation/Reflection ($V^T$):** A rotation or reflection transformation of the orthogonal coordinate system in $\mathbb{R}^n$.
2.  **Scaling ($\Sigma$):** A scaling transformation of the coordinate axes by the factors $\sigma_1, \sigma_2, \ldots$. These singular values ($\sigma_i$) represent the pure scaling applied by the matrix.
3.  **Rotation/Reflection ($U$):** A final rotation or reflection transformation of the coordinate system in $\mathbb{R}^m$.

If $A$ transforms a unit ball in $\mathbb{R}^n$ into an ellipsoid, **the lengths of the axes of this ellipsoid are the singular values of $A$**.

## Connection between SVD and the Four Fundamental Spaces

The vectors composing the orthogonal matrices $U$ and $V$ provide orthonormal bases for the four fundamental subspaces associated with matrix $A$:

Let $r$ be the rank of matrix $A$, which is equal to the number of non-zero singular values ($\sigma_r > 0$, $\sigma_{r+1} = 0$, etc.).

| Subspace | Basis Vectors (Dimension $r$) | Orthonormal Bases Constituted by SVD Vectors |
| :--- | :--- | :--- |
| **Row Space** ($R(A^T)$) | Subspace of $\mathbb{R}^n$ | The first $r$ right singular vectors: $v_1, \ldots, v_r$ |
| **Null Space** ($N(A)$) | Subspace of $\mathbb{R}^n$ | The remaining $n-r$ right singular vectors: $v_{r+1}, \ldots, v_n$ |
| **Column Space** ($R(A)$) | Subspace of $\mathbb{R}^m$ | The first $r$ left singular vectors: $u_1, \ldots, u_r$ |
| **Left Null Space** ($N(A^T)$) | Subspace of $\mathbb{R}^m$ | The remaining $m-r$ left singular vectors: $u_{r+1}, \ldots, u_m$ |

The SVD bases demonstrate that the row space ($R(A^T)$) and the null space ($N(A)$) are orthogonal complements in $\mathbb{R}^n$, while the column space ($R(A)$) and the left null space ($N(A^T)$) are orthogonal complements in $\mathbb{R}^m$.

## Connection with Rank and the Rank-Nullity Theorem

The **rank** of matrix $A$ is precisely equal to $r$, the number of positive singular values $\sigma_i$.

The SVD provides orthonormal bases that explicitly recover the **Rank-Nullity Theorem**, which states that the dimension of the column space (rank) plus the dimension of the null space (nullity) equals the number of columns ($n$).
*   $\text{dim}(R(A^T)) = r$ (Rank).
*   $\text{dim}(N(A)) = n - r$ (Nullity).
*   $r + (n-r) = n$.

## Proof of the Existence of the SVD

The existence of the SVD for any $m \times n$ real matrix $A$ is guaranteed by a basic theorem (Theorem 15.1). A constructive proof can be derived by analyzing the related symmetric matrices. Another approach involves diagonalizing a related symmetric matrix:

````{prf:proof}

Let $A$ be an $m \times n$ matrix. Define the symmetric matrix $B$ as:

$$ B = \begin{bmatrix} 0 & A \\ A^T & 0 \end{bmatrix} $$

Since $B$ is symmetric, it is unitarily diagonalizable:

$$ B = Q \Lambda Q^T $$

where $\Lambda$ is a diagonal matrix of eigenvalues and $Q$ is an orthogonal matrix of eigenvectors.
We consider an eigenvector $[x; y]$ of $B$ corresponding to a non-zero eigenvalue $\lambda$. The definition of the eigenvalue/eigenvector relationship $B [x; y] = \lambda [x; y]$ yields the equations:

$$
\begin{gather}
A y = \lambda x \\
A^T x = \lambda y
\end{gather}
$$

Applying $A^T$ to the first equation and $A$ to the second, we find:

$$
\begin{gather}
A^T A y = \lambda (A^T x) = \lambda (\lambda y) = \lambda^2 y \\
A A^T x = \lambda (A y) = \lambda (\lambda x) = \lambda^2 x
\end{gather}
$$

This shows that $\lambda^2$ is an eigenvalue of both $A^T A$ and $A A^T$. We define the singular values $\Sigma$ such that the diagonal elements are the positive square roots of these eigenvalues: $\sigma_i = \sqrt{\lambda_i^2}$.

It can also be shown that $[x; -y]$ is an eigenvector corresponding to the eigenvalue $-\lambda$.

If $X$ denotes the eigenvectors of $A A^T$ and $Y$ denotes the eigenvectors of $A^T A$ (corresponding to non-zero $\lambda^2$), we can structure the matrices $Q$ and $\Lambda$ related to $B$ as:

$$ Q = \frac{1}{\sqrt{2}} \begin{bmatrix} X & X \\ Y & -Y \end{bmatrix}, \quad \Lambda = \begin{bmatrix} \Sigma & 0 \\ 0 & -\Sigma \end{bmatrix} $$

By substitution into $B = Q \Lambda Q^T$, the decomposition $A = X \Sigma Y^T$ is recovered. If we identify $U=X$ and $V=Y$, the SVD $A = U \Sigma V^T$ is established.
````

## Relation of SVD with Eigendecomposition of Symmetric Matrices

The SVD is intimately linked to the eigendecomposition of the symmetric matrices $A^T A$ and $A A^T$.

1.  **Right Singular Vectors ($V$):** The columns of $V$ (the right singular vectors) are the **eigenvectors of $A^T A$**. The diagonalization of $A^T A$ is given by:
2.  
    $$ A^T A = V \Sigma^T \Sigma V^T = V \Sigma^2 V^T $$

3.  **Left Singular Vectors ($U$):** The columns of $U$ (the left singular vectors) are the **eigenvectors of $A A^T$**. The diagonalization of $A A^T$ is given by:
4.  
    $$ A A^T = U \Sigma \Sigma^T U^T = U \Sigma^2 U^T $$

5.  **Singular Values ($\Sigma$):** The squares of the singular values ($\sigma_i^2$) are the **eigenvalues of both $A^T A$ and $A A^T$**. Thus, the singular values $\sigma_j$ are the square roots of these eigenvalues: $\sigma_j = \sqrt{\lambda_j}$.

## Thin or Truncated SVD

The full SVD, $A = U \Sigma V^T$, is sometimes referred to as the full Singular Value Decomposition. In practice, **compact** and **truncated** forms are commonly used.

1.  **Compact SVD (Thin SVD):**
    *   This form is obtained when the SVD has a rank equal to the rank $r$ of the original matrix $A$.
    *   It is represented as $A = U_r \Sigma_r V_r^T$.
    *   $\Sigma_r$ is an $r$-order diagonal matrix containing only the $r$ positive singular values. $U_r$ contains the first $r$ columns of $U$ (the orthonormal bases for $R(A)$), and $V_r$ contains the first $r$ columns of $V$ (the orthonormal bases for $R(A^T)$).
    *   Compact SVD corresponds to **lossless compression** of the data.

2.  **Truncated SVD:**
    *   This is the form typically referred to in practical applications.
    *   It is obtained by taking only the part corresponding to the largest $k$ singular values, where $k < r$ (the rank of $A$).
    *   The approximation is $A \approx U_k \Sigma_k V_k^T$.
    *   The matrix $A_k = U_k \Sigma_k V_k^T$ has a rank of $k$, which is lower than the rank of the original matrix. This decomposition provides an optimal low-rank approximation of $A$.
    *   Truncated SVD corresponds to **lossy compression**.

We discuss this result further in the next section with the Eckart-Young-Mirsky Theorem.

## The Eckart-Young-Mirsky Theorem: Optimal Low-Rank Approximation

````{prf:theorem} The Eckart-Young-Mirsky Theorem
:label: thm:eckart-young-mirsky
Let the SVD of a matrix $A$ be $A = U \Sigma V^T$. Let $A_k$ be the truncated SVD matrix of rank $k$ obtained by keeping only the $k$ largest singular values:

$$
A_k = U_k \Sigma_k V_k^T
$$

The matrix $A_k$ is the best rank-$k$ approximation of $A$ in both the spectral norm and the Frobenius norm. That is, for any matrix $B$ with $\text{rank}(B) = k$:

$$
\|A - A_k\|_2 \le \|A - B\|_2 \quad \text{and} \quad \|A - A_k\|_F \le \|A - B\|_F
$$

The approximation error is given by the first neglected singular value: $\|A - A_k\|_2 = \sigma_{k+1}$.
````

This theorem is the theoretical bedrock for using SVD in data compression, noise reduction, and machine learning. It guarantees that truncating the SVD is not just a heuristic but the mathematically optimal way to reduce the dimensionality (rank) of your data while minimizing the error in these important norms.

## Connection with Matrix Norms and the Determinant of a Matrix

### Matrix Norms (Frobenius Norm)

The SVD provides an optimal method of matrix approximation in the sense of the **Frobenius norm**. The Frobenius norm ($\|A\|_F$) is a generalization of the $L_2$ norm of a vector.

The Frobenius norm of a matrix $A$ is related directly to its singular values:

$$ \|A\|_F = \left( \sigma_1^2 + \sigma_2^2 + \cdots + \sigma_n^2 \right)^{1/2} $$

Furthermore, the truncated SVD with rank $k$ yields the optimal $k$-rank approximation $X=A_k$ in the sense of the Frobenius norm. The error of this optimal approximation is given by the neglected singular values:

$$ \|A - X\|_F = \left( \sigma_{k+1}^2 + \sigma_{k+2}^2 + \cdots + \sigma_n^2 \right)^{1/2} $$

The spectral norm (operator 2-norm) of $A$ is equal to the largest singular value, $\sigma_1(A)$:

$$ \|A\|_2 = \sigma_1(A) $$

A key conceptual link exists between the Singular Value Decomposition (SVD) and the definition of various **matrix norms**, especially those used in analyzing low-rank approximations. 

While the Frobenius norm is mentioned as the metric for optimal approximation, SVD is also crucial for defining the **Schatten $p$-norms**.

### Connection to Schatten $p$-Norms

The Schatten $p$-norm is a generalization that uses the singular values $\sigma_i$ in a manner analogous to how the $L_p$ vector norm uses vector elements.

The **Schatten $p$-norm** of a matrix $A$, denoted as $\|A\|_p$, is defined using its singular values ($\sigma_i$):

$$ \|A\|_p = \left( \sum_{i=1}^r \sigma_i^p \right)^{1/p} $$

where $r$ is the rank of the matrix $A$ (the number of non-zero singular values).

The Schatten $p$-norms are derived by treating the singular values as a vector and calculating the vector $p$-norm of that vector.

The Schatten $p$-norms generalize three highly important matrix norms derived from SVD:

1.  **Schatten 2-Norm (Frobenius Norm):**
    When $p=2$, the Schatten norm is equivalent to the **Frobenius norm**.

    $$ \|A\|_F = \|A\|_2 = \left( \sum_{i=1}^r \sigma_i^2 \right)^{1/2} $$

    The Frobenius norm is used to quantify the "square loss" in matrix approximation. The SVD provides the optimal rank-$k$ approximation regarding this norm.

2.  **Schatten 1-Norm (Nuclear Norm or Trace Norm):**
    When $p=1$, the Schatten norm is called the **Nuclear Norm** (or Trace Norm).
    
    $$ \|A\|_* = \|A\|_1 = \sum_{i=1}^r \sigma_i $$
    
    This norm is widely used in machine learning for promoting low-rank solutions (known as matrix completion or sparse PCA).

3.  **Schatten $\infty$-Norm (Spectral Norm or Operator 2-Norm):**
    When $p \to \infty$, the Schatten norm converges to the **Spectral Norm** (also known as the Operator 2-norm, or $\|A\|_2$):
    
    $$ \|A\|_2 = \sigma_1(A) $$
    
    The spectral norm is equal to the largest singular value, $\sigma_1(A)$. This norm is crucial as it measures the maximum stretching factor of the matrix transformation.

## Determinant and Singular Values

The absolute value of the determinant of a square matrix is equal to the product of its singular values.

````{prf:theorem}
:label: thm:det-svd
For any square matrix $A \in \mathbb{C}^{n \times n}$ with singular values $\sigma_1, \sigma_2, \dots, \sigma_n$, the following relationship holds:

$$
|\det(A)| = \prod_{i=1}^n \sigma_i
$$

````

````{prf:proof}

Start with the **Singular Value Decomposition (SVD)** of the matrix $A$:

$$
A = U \Sigma V^H
$$

where $U$ and $V$ are unitary matrices, and $\Sigma$ is a diagonal matrix containing the singular values $\sigma_i$ on its diagonal.

Take the determinant of both sides of the equation:

$$
\det(A) = \det(U \Sigma V^H)
$$

Using the multiplicative property of determinants, $\det(XYZ) = \det(X)\det(Y)\det(Z)$, we can separate the terms:

$$
\det(A) = \det(U) \det(\Sigma) \det(V^H)
$$

The determinant of a diagonal matrix is the product of its diagonal entries. The diagonal entries of $\Sigma$ are the singular values:

$$
\det(\Sigma) = \prod_{i=1}^n \sigma_i
$$

For any unitary matrix $Q$, its determinant has a magnitude of one, i.e., $|\det(Q)|=1$. Therefore, $|\det(U)| = 1$ and $|\det(V^H)| = 1$.

Now, take the absolute value of the entire determinant equation:

$$
|\det(A)| = |\det(U)| \cdot |\det(\Sigma)| \cdot |\det(V^H)|
$$

Substituting the known values gives:

$$
|\det(A)| = 1 \cdot \left(\prod_{i=1}^n \sigma_i\right) \cdot 1
$$

This simplifies to the final result:

$$
|\det(A)| = \prod_{i=1}^n \sigma_i
$$
````

### Geometric Interpretation 📐

This relationship has a clear geometric meaning. The **singular values** ($\sigma_i$) represent the scaling factors that the matrix applies along its principal axes (the directions of maximum stretch). The **determinant** represents the total volume scaling factor of the transformation. This theorem shows that the total change in volume is simply the product of the individual stretches along these principal, orthogonal directions.