# The Four Fundamental Subspaces

Every $m \times n$ matrix $A$ defines a linear transformation and is associated with four fundamental vector subspaces. These spaces provide a complete geometric understanding of the matrix's behavior, revealing what happens to vectors when they are transformed by $A$. 🗺️



## 1. The Column Space (Range)

The **column space**, also known as the **range**, is the set of all possible output vectors of the transformation. It consists of all vectors that can be formed by multiplying the matrix $A$ by an input vector $x$.

$$R(A) = \{y \in \mathbb{R}^m \,|\, y = Ax \text{ for some } x \in \mathbb{R}^n\}$$

Geometrically, the column space is the **span of the column vectors** of $A$. It is a subspace of the output space, $\mathbb{R}^m$. The dimension of the column space is the **rank** of the matrix, denoted by $r$. The rank represents the number of linearly independent columns in the matrix.

## 2. The Null Space (Kernel)

The **null space**, also known as the **kernel**, is the set of all input vectors that are mapped to the zero vector by the transformation.

$$N(A) = \{x \in \mathbb{R}^n \,|\, Ax = 0\}$$

The null space is a subspace of the input space, $\mathbb{R}^n$. If the null space contains vectors other than the zero vector, it means the transformation is "many-to-one"—multiple different input vectors are "squashed" onto the same output vector (zero). The dimension of the null space is called the **nullity**.

## 3. The Row Space

The **row space** is the set of all linear combinations of the row vectors of $A$. It is equivalent to the column space of the transpose of $A$.

$$R(A^T)$$

The row space is a subspace of the input space, $\mathbb{R}^n$. A fundamental result in linear algebra is that the dimension of the row space is also equal to the **rank** ($r$) of the matrix. This means a matrix always has the same number of linearly independent rows as it does linearly independent columns.

## 4. The Left Null Space

The **left null space** is the null space of the transpose of $A$.

$$N(A^T) = \{y \in \mathbb{R}^m \,|\, A^T y = 0\}$$

It is called the "left" null space because the equation can be transposed to $y^T A = 0^T$, showing that the vectors $y$ act on $A$ from the left. The left null space is a subspace of the output space, $\mathbb{R}^m$.

## The Fundamental Theorem of Linear Algebra

This theorem connects the dimensions and relationships of the four subspaces.

### The Rank-Nullity Theorem
This theorem relates the dimensions of the column space and the null space. It states that the rank of a matrix plus its nullity is equal to the number of columns in the matrix.

$$\text{dim}(R(A)) + \text{dim}(N(A)) = n$$

Or more simply: **rank + nullity = number of columns**.

This theorem provides a beautiful intuition: the dimension of the input space ($n$) is split. Part of it is preserved and forms the output image (the rank), and the other part is collapsed to zero (the nullity).

### Orthogonal Complements
The subspaces come in orthogonal pairs:

* The **row space** and the **null space** are orthogonal complements in the input space $\mathbb{R}^n$. This means every vector in the row space is perpendicular to every vector in the null space.
* The **column space** and the **left null space** are orthogonal complements in the output space $\mathbb{R}^m$.