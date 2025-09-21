# Orthogonal Matrices

An **orthogonal matrix** is a square matrix whose columns and rows are orthonormal vectors. These matrices represent rigid transformations, like rotations and reflections, that preserve lengths and angles. They are fundamental in many areas of mathematics and computer science, especially in numerical linear algebra, where they contribute to the stability and accuracy of algorithms.

## Definition and Core Properties

A real square matrix $Q \in \mathbb{R}^{n \times n}$ is **orthogonal** if its transpose is equal to its inverse:

$$
Q^T = Q^{-1}
$$

This leads to the more common definition:

$$
Q^T Q = Q Q^T = I
$$

where $I$ is the identity matrix.

This defining property implies that the columns (and rows) of $Q$ form an **orthonormal set**. This means:
1.  **Orthogonality**: The dot product of any two distinct columns (or rows) is zero.
2.  **Unit Norm**: The length (Euclidean norm) of each column (or row) is one.

A non-square matrix $Q \in \mathbb{R}^{m \times n}$ with $m > n$ is also called orthogonal if its columns are orthonormal, satisfying $Q^T Q = I_n$. However, in this case, $Q Q^T \neq I_m$.

**Unitary Matrix**: The equivalent concept for complex matrices is a **unitary matrix**, which satisfies $Q^H Q = I$, where $Q^H$ is the conjugate transpose of $Q$.

## Key Results and Mathematical Formulas

Orthogonal matrices have several important properties that make them incredibly useful.

### Isometry (Length and Angle Preservation)
Orthogonal transformations are **isometries**, meaning they preserve the Euclidean norm (length) of a vector.

$$
\|Qx\|_2^2 = (Qx)^T(Qx) = x^T Q^T Q x = x^T I x = x^T x = \|x\|_2^2
$$

Therefore, $\|Qx\|_2 = \|x\|_2$. They also preserve the dot product, and thus the angles between vectors:

$$
(Qx) \cdot (Qy) = (Qx)^T(Qy) = x^T Q^T Q y = x^T y = x \cdot y
$$

### Determinant
The determinant of an orthogonal matrix is always either **+1** or **-1**.

$$
\det(Q^T Q) = \det(Q^T)\det(Q) = (\det(Q))^2 = \det(I) = 1 \implies \det(Q) = \pm 1
$$

* $\det(Q) = +1$: The transformation is a **rotation** (a proper rotation). It preserves the orientation of the space. The set of these matrices forms the **special orthogonal group SO(n)**.
* $\det(Q) = -1$: The transformation is a **reflection** or an **improper rotation** (e.g., a reflection followed by a rotation). It reverses the orientation of the space.

## The Cartan–Dieudonné Theorem

The **Cartan–Dieudonné theorem** is a fundamental result in geometry that provides a simple, constructive way to think about orthogonal transformations. It states that any orthogonal transformation can be broken down into a series of simpler reflections.

````{prf:theorem} Cartan–Dieudonné theorem
:label: thm:cartan_dieudonne
Every orthogonal transformation in an $n$-dimensional Euclidean space ($\mathbb{R}^n$) can be described as the composition of at most **n** reflections.
````

This means for any orthogonal matrix $Q$, we can write:

$$
Q = H_1 H_2 \cdots H_k
$$

where each $H_i$ is a reflection matrix (like a Householder matrix) and $k \le n$. If $k=0$, $Q$ is the identity matrix.

### Proof by Induction

````{prf:proof}

The proof works by showing that we can use one reflection to simplify the problem to a lower-dimensional space.

* **Base Case (n=1)**: In a 1-dimensional space ($\mathbb{R}$), an orthogonal transformation is either $x \mapsto x$ (the identity, 0 reflections) or $x \mapsto -x$ (a reflection, 1 reflection). The theorem holds.

* **Inductive Step**: Assume the theorem is true for transformations in $\mathbb{R}^{n-1}$. Let $Q$ be an orthogonal transformation in $\mathbb{R}^n$.
    1.  Consider the first basis vector $e_1 = (1, 0, \dots, 0)^T$. Let $v = Qe_1$.
    2.  **Case 1**: If $v = e_1$, then $Q$ fixes the first basis vector. This means that the $(n-1)$-dimensional subspace $S$ orthogonal to $e_1$ is invariant under $Q$. By our induction hypothesis, this transformation on $S$ can be achieved with at most $n-1$ reflections. So, the theorem holds for $Q$.
    3.  **Case 2**: If $v \neq e_1$, we can find a reflection $H_1$ that maps $v$ back to $e_1$. The specific reflection is the one across the hyperplane that is orthogonal to the vector $(v - e_1)$. This gives us $H_1 v = e_1$.
    4.  Now consider the new transformation $Q' = H_1 Q$. It is also orthogonal since it's a product of two orthogonal transformations. And importantly, it fixes the vector $e_1$:
        
        $$
        Q'e_1 = (H_1 Q)e_1 = H_1(Qe_1) = H_1 v = e_1
        $$
        
    5.  Since $Q'$ fixes $e_1$, it behaves like an $(n-1)$-dimensional orthogonal transformation on the subspace orthogonal to $e_1$. By the induction hypothesis, $Q'$ can be written as the composition of at most $n-1$ reflections, say $Q' = H_2 H_3 \cdots H_k$ where $k-1 \le n-1$.
    6.  Substituting back, we get $H_1 Q = H_2 H_3 \cdots H_k$.
    7.  Since $H_1$ is a reflection, its inverse is itself ($H_1^{-1} = H_1$). We can solve for $Q$:
        
        $$
        Q = H_1^{-1} (H_2 H_3 \cdots H_k) = H_1 H_2 H_3 \cdots H_k
        $$
        
    8.  This expresses $Q$ as a product of $k$ reflections, where $k \le (n-1)+1 = n$.

Thus, the theorem is proven.
````