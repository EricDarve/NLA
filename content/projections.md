# Projections

A **projection** is a linear transformation that takes a vector and maps it onto a subspace. Intuitively, you can think of it as casting a shadow. If you shine a light on an object (a vector), the shadow it creates on a surface (a subspace) is its projection. Projections are fundamental in areas like computer graphics, statistics (for least-squares regression), and solving systems of linear equations.

All projections are represented by a square matrix $P$ that is **idempotent**, meaning that applying the projection more than once has no further effect. Algebraically, this is the defining property of a projection matrix:

$$P^2 = P$$

## Types of Projections: Orthogonal vs. Oblique

The distinction between the two main types of projections depends on the direction of the projection relative to the target subspace. 📐

### Orthogonal Projections
An **orthogonal projection** is the most common type. Here, the direction of projection is **perpendicular (orthogonal)** to the subspace you are projecting onto. This is like the shadow cast by the sun when it's directly overhead.

A key feature of an orthogonal projection is that it finds the point in the subspace that is **geometrically closest** to the original vector. This "best approximation" property is the foundation of least-squares solutions to overdetermined systems of equations.

### Oblique Projections
An **oblique projection** is one where the direction of projection is **not perpendicular** to the target subspace. This is analogous to the long, angled shadow cast by the sun when it's low in the sky. While less common than orthogonal projections, they are used in specific applications like perspective projections in computer graphics.

## Mathematical Formulas

The specific formula for a projection matrix depends on the type of projection and the subspace being projected onto.

### Orthogonal Projection Formula
To project a vector onto the subspace spanned by the columns of a matrix $A$ (where the columns of $A$ are linearly independent), the orthogonal projection matrix $P$ is given by:

$$P = A(A^T A)^{-1} A^T$$

The projection of a vector $b$ onto this subspace is then calculated as $p = Pb$.

#### Special Case: Projection onto a Line
If you are projecting onto a line spanned by a single non-zero vector $a$, the formula simplifies significantly. In this case, $A$ is just the column vector $a$.

$$P = \frac{aa^T}{a^T a}$$

The projection of a vector $b$ onto the line defined by $a$ is:

$$p = \left(\frac{a \cdot b}{\|a\|^2}\right) a$$

This formula has a clear geometric meaning: it scales the direction vector $a$ by a factor determined by the dot product of $a$ and $b$.

### Oblique Projection Formula
For an oblique projection, you must specify both the subspace to project **onto** (the range, spanned by the columns of matrix $A$) and the direction to project **along** (the null space, which can be defined by a matrix $B$). The projection matrix is:

$$P = A(B^T A)^{-1} B^T$$

Notice that if you choose $B = A$, the projection direction becomes orthogonal to the range, and this formula reduces to the orthogonal projection formula.

## Key Results and Properties

Projection matrices have several important and elegant properties.

* **Idempotence**: As mentioned, all projection matrices satisfy $P^2 = P$.
* **Symmetry**: A projection matrix $P$ represents an **orthogonal projection** if and only if it is a **symmetric matrix** ($P = P^T$). Oblique projection matrices are generally not symmetric.
* **Complementary Projection**: If $P$ is a projection matrix, then the matrix $Q = I - P$ is also a projection matrix. $Q$ is called the **complementary projection**. It projects onto the null space of $P$ along the range of $P$. For an orthogonal projection $P$, its complement $I-P$ is also an orthogonal projection onto the orthogonal complement of the original subspace.

## Orthogonal Matrices and Projections

The connection is simple and direct: if the columns of a matrix $A$ are **orthonormal**, meaning they form an **orthogonal matrix** $Q$, the formula for the orthogonal projection matrix onto the column space of $A$ simplifies dramatically.

Specifically, the projection matrix $P$ becomes $P = QQ^T$.

### Derivation and Explanation

The standard formula for an orthogonal projection matrix $P$ that projects vectors onto the column space of a matrix $A$ is:

$$P = A(A^T A)^{-1} A^T$$

An **orthogonal matrix** $Q$ (which must be a tall or square matrix with orthonormal columns) has the defining property that its columns are mutually perpendicular and have a length of one. Algebraically, this means:

$$Q^T Q = I$$

where $I$ is the identity matrix.

If we substitute our orthogonal matrix $Q$ for the general matrix $A$ in the projection formula, the $(Q^T Q)$ term simplifies to the identity matrix $I$:

$$P = Q(\underbrace{Q^T Q}_{I})^{-1} Q^T = Q(I)^{-1}Q^T = QIQ^T = QQ^T$$

So, the orthogonal projection matrix onto the subspace spanned by a set of orthonormal vectors (the columns of $Q$) is simply $P = QQ^T$.

### Properties of the Projection Matrix $P = QQ^T$

This simplified form, $P=QQ^T$, neatly demonstrates the two defining properties of an orthogonal projection matrix.

1.  **Symmetry**: A matrix is symmetric if $P^T = P$.

    $$P^T = (QQ^T)^T = (Q^T)^T Q^T = QQ^T = P$$

    The matrix is its own transpose, so it is **symmetric**.

2.  **Idempotence**: A matrix is idempotent if $P^2 = P$.
    
    $$P^2 = (QQ^T)(QQ^T) = Q(\underbrace{Q^T Q}_{I})Q^T = QIQ^T = QQ^T = P$$

    Projecting a vector twice is the same as projecting it once.

## Geometric Interpretation 📐

When we project a vector $x$ onto the column space of an orthogonal matrix $Q$, the operation $p = Px = QQ^Tx$ can be seen as a two-step process:

1.  $Q^Tx$: First, we multiply by $Q^T$. This calculates the dot products of $x$ with each orthonormal column of $Q$. In essence, it finds the coordinates of $x$ in the basis defined by the columns of $Q$.

2.  $Q(Q^Tx)$: Next, we multiply the resulting coordinate vector by $Q$. This creates a linear combination of the orthonormal basis vectors (the columns of $Q$) using those coordinates.

The result is a new vector $p$ that lives within the column space of $Q$ and is the closest possible vector in that subspace to the original vector $x$.