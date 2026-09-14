# Matrix-Matrix Multiplications

## Composition of Linear Transformations

If $T_1: \mathbb{R}^n \to \mathbb{R}^m$ has matrix $B$ and $T_2: \mathbb{R}^m \to \mathbb{R}^\ell$ has matrix $A$,  
then the composition $T_2 \circ T_1$ has matrix $AB$, where:

$$
(AB)x = A(Bx).
$$

The formula for **matrix-matrix multiplication** is:

$$
c_{ij} = \sum_{k=1}^m a_{ik} b_{kj}
$$

where $C = AB$. Column by column, this says that the $j$-th column of $AB$ is $A b_j$, where $b_j$ is the $j$-th column of $B$, which extends the column picture of the matrix-vector product.

Each of the $\ell n$ entries of $C$ is a dot product of length $m$, so forming $AB$ takes about $2 \ell m n$ flops. For square matrices this is $O(n^3)$, the baseline against which the cost of every factorization in later chapters is measured.

**Properties:**
- Defined only if the number of columns of $A$ equals the number of rows of $B$.
- **Associative:** $(AB)C = A(BC)$.
- **Not commutative:** $AB \ne BA$ in general.
- **Transpose rule:** $(AB)^T = B^T A^T$, and $(AB)^H = B^H A^H$.

---

## Identity and Inverse Matrices

The **identity matrix** $I_n$ has $1$ on the diagonal and $0$ elsewhere, and satisfies:

$$
I_n x = x
$$

for all $x \in \mathbb{R}^n$. It also acts as the identity for matrix multiplication: $A I_n = I_m A = A$ for any $m \times n$ matrix $A$.

A square matrix $A$ is **invertible** if there exists $A^{-1}$ such that:

$$
A^{-1}A = AA^{-1} = I_n.
$$

When it exists, $A^{-1}$ is unique. For invertible $A$ and $B$ of the same size, inverting a product reverses the order:

$$
(AB)^{-1} = B^{-1} A^{-1}.
$$

```{admonition} Fact
A square matrix $A$ is invertible if and only if its columns are linearly independent.
```

---

## Change of Basis

Let $\{v_1, \dots, v_n\}$ be a basis of $\mathbb{R}^n$. Any $x \in \mathbb{R}^n$ can be written as:

$$
x = \sum_{i=1}^n \alpha_i v_i
$$

where $\alpha = (\alpha_1, \dots, \alpha_n)^T$ are the **coordinates of $x$ in this basis**.

If $V$ is the matrix with $v_i$ as columns, then:

$$
x = V\alpha, \quad \alpha = V^{-1}x.
$$

Let $T: \mathbb{R}^n \to \mathbb{R}^n$, so that its matrix is square and the same basis is used for the input and the output. If $A$ is the matrix of $T$ in the standard basis and $B$ is the matrix in the basis $\{v_i\}$, then:

$$
A = V B V^{-1},
\qquad
B = V^{-1} A V.
$$

This is the **change of basis formula**. Matrices related this way are called **similar**, and the map $B \mapsto V B V^{-1}$ is a **similarity transformation**. Similar matrices represent the same linear transformation in different bases; later chapters use similarity transformations to reduce a matrix to a simpler form without changing its eigenvalues.