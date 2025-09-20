# Schur Decomposition

The Schur decomposition is a fundamental result in linear algebra that guarantees that any square matrix can be represented as the product of a unitary matrix and an upper triangular matrix. Unlike eigendecomposition, which only works for diagonalizable matrices, the Schur decomposition exists for *every* square matrix, making it a more general and numerically reliable tool.

## The (Complex) Schur Decomposition

````{prf:theorem} Schur Decomposition

The theorem states that for any $n \times n$ square matrix $A$ with complex entries, there exists a **unitary matrix** $Q$ and an **upper triangular matrix** $T$ such that:

$$A = Q T Q^H$$

````

This decomposition has two key components:
* **$T$**: An upper triangular matrix whose diagonal entries are the **eigenvalues** of $A$.
* **$Q$**: A unitary matrix ($Q^H Q = I$), which means its columns form an orthonormal basis for the space $\mathbb{C}^n$.

The Schur decomposition is incredibly useful because it reveals the eigenvalues of a matrix while using unitary transformations. These transformations are numerically stable (they preserve lengths and do not amplify errors), making the Schur decomposition one of the most reliable methods for computing eigenvalues in practice.

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

## The Real Schur Decomposition

A complication arises when working with real matrices. If a real matrix $A$ has complex eigenvalues (which must appear in conjugate pairs), it is impossible to make the matrix $T$ in the decomposition $A = Q T Q^T$ both real and triangular, as its diagonal would have to contain the complex eigenvalues.

The **real Schur decomposition** provides a clever workaround. 

````{prf:theorem} Real Schur Decomposition
It states that for any real square matrix $A$, there exists a **real orthogonal matrix** $Q$ ($Q^T Q = I$) such that:

$$A = Q S Q^T$$

Here, $S$ is a **real block upper triangular matrix**, also known as a **quasi-triangular matrix**.

The structure of $S$ is upper triangular, but its diagonal may contain both:
* **$1 \times 1$ blocks**: These are the real eigenvalues of $A$.
* **$2 \times 2$ blocks**: These correspond to pairs of complex conjugate eigenvalues. A block of the form $\begin{pmatrix} a & b \\ -b & a \end{pmatrix}$ represents the complex conjugate pair of eigenvalues $a \pm ib$.

````

For example, a $5 \times 5$ real matrix with one real eigenvalue $\lambda_1$ and two pairs of complex conjugate eigenvalues ($a \pm ib$ and $c \pm id$) would have a real Schur form like this:

$$
S = \begin{pmatrix}
\lambda_1 & * & * & * & * \\
0 & a & b & * & * \\
0 & -b & a & * & * \\
0 & 0 & 0 & c & d \\
0 & 0 & 0 & -d & c
\end{pmatrix}
$$

This form allows the entire decomposition to remain within the real numbers, which is often highly advantageous for computational algorithms.

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