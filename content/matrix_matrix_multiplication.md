# Matrix-Matrix Multiplication

## Composition of Linear Transformations

If $T_1: \mathbb{R}^n \to \mathbb{R}^m$ has matrix $B$ and $T_2: \mathbb{R}^m \to \mathbb{R}^\ell$ has matrix $A$,  
then the composition $T_2 \circ T_1$ has matrix $AB$, where:

$$
(AB)x = A(Bx).
$$

The formula for **matrix-matrix multiplication** is:

$$
c_{ij} = \sum_{k} a_{ik} b_{kj}
$$

where $C = AB$.

<video controls width="640">
  <source src="../_static/CompositionScene.mp4" type="video/mp4">
</video> 

**Properties:**
- Defined only if the number of columns of $A$ equals the number of rows of $B$.
- **Associative:** $(AB)C = A(BC)$.
- **Not commutative:** $AB \ne BA$ in general.
- **Transpose rule:** $(AB)^T = B^T A^T$.

---

## Identity and Inverse Matrices

The **identity matrix** $I_n$ has $1$ on the diagonal and $0$ elsewhere, and satisfies:

$$
I_n x = x
$$

for all $x \in \mathbb{R}^n$.

A square matrix $A$ is **invertible** if there exists $A^{-1}$ such that:

$$
A^{-1}A = AA^{-1} = I_n.
$$

```{admonition} Fact
$A$ is invertible if and only if its columns are linearly independent.
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

If $A$ is the matrix of $T$ in the standard basis and $B$ is the matrix in the basis $\{v_i\}$, then:

$$
A = V B V^{-1}.
$$

This is the **change of basis formula**.