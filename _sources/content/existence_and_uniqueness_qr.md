# Existence and Uniqueness of QR Decomposition

We have now seen three different algorithms (Householder, Givens, and Modified Gram-Schmidt) for computing the QR decomposition of a matrix $A$. A natural question arises: does this factorization always exist, and if so, is it unique?

## Existence

````{prf:theorem} Existence of QR Decomposition
:label: thm:existence_qr
For any real matrix $A \in \mathbb{R}^{m \times n}$, there exists a QR decomposition such that $A = QR$, where $Q$ is an orthogonal matrix and $R$ is an upper triangular (or upper trapezoidal) matrix.
````

The constructive nature of the algorithms we've studied is a proof of this.

* **Full QR**: For any $A$, we can find an orthogonal matrix $Q \in \mathbb{R}^{m \times m}$ and an upper trapezoidal matrix $R \in \mathbb{R}^{m \times n}$ such that $A=QR$.
* **Thin QR**: If $m \ge n$, we can find a matrix $Q_1 \in \mathbb{R}^{m \times n}$ with orthonormal columns and an upper triangular matrix $R_1 \in \mathbb{R}^{n \times n}$ such that $A=Q_1 R_1$.

The algorithms we have discussed can always be carried out, guaranteeing existence.

## Uniqueness

The uniqueness of the QR decomposition is more nuanced and depends on the properties of the matrix $A$ and any constraints we place on the matrix $R$.

### **The Full Rank Case**

````{prf:theorem} Unique Thin QR Decomposition for Full Column Rank Matrices
:label: thm:unique_thin_qr
Let's assume $A \in \mathbb{R}^{m \times n}$ with $m \ge n$ has **full column rank** (i.e., its columns are linearly independent).

In this case, the **thin QR decomposition** $A = Q_1 R_1$ is **unique** if we impose the condition that all the **diagonal elements of $R_1$ must be positive** ($r_{ii} > 0$).
````

````{prf:proof} 
Consider the product $A^T A$:

$$
A^T A = (Q_1 R_1)^T (Q_1 R_1) = R_1^T Q_1^T Q_1 R_1
$$

Since $Q_1$ has orthonormal columns, $Q_1^T Q_1 = I$ (the identity matrix). This simplifies the expression to:

$$
A^T A = R_1^T R_1
$$

The matrix $A^T A$ is symmetric and positive definite (since $A$ has full rank). The equation above is precisely the **Cholesky decomposition** of $A^T A$. The Cholesky decomposition is unique when we require the diagonal elements of the triangular factor ($R_1$ in this case) to be positive. Therefore, $R_1$ is uniquely determined.

Once $R_1$ is uniquely determined, we can find $Q_1$ from the relation $Q_1 = A R_1^{-1}$. Since $A$ has full rank and $R_1$ is nonsingular, $Q_1$ is also uniquely determined.
````

*Note on Full QR*: Even in this full-rank case, the **full QR decomposition is not unique**. The first $n$ columns of $Q$ are the unique columns of $Q_1$, but the remaining $m-n$ columns only need to form an orthonormal basis for the orthogonal complement of the column space of $A$. There are infinitely many ways to choose this basis.

### **The Rank-Deficient Case**

If the matrix $A$ is **rank-deficient** (its columns are linearly dependent), then the QR decomposition is **not unique**.

In this case, at least one of the diagonal elements of $R$ will be zero. This introduces degrees of freedom into the system $A^T A = R^T R$, and the Cholesky factorization is no longer unique. This lack of uniqueness in $R$ carries over to $Q$.