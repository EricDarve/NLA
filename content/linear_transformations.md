# Linear Transformations and Matrices

## Linear Transformations

A **linear transformation** from $\mathbb{R}^n$ to $\mathbb{R}^m$ is a map  
$T: \mathbb{R}^n \to \mathbb{R}^m$ that satisfies **linearity**:

$$
T(x + \alpha y) = T(x) + \alpha T(y)
$$

for all $x, y \in \mathbb{R}^n$ and scalars $\alpha \in \mathbb{R}$.

A linear transformation is completely determined by its action on a basis of $\mathbb{R}^n$.  
If $\{v_1, \dots, v_n\}$ is a basis of $\mathbb{R}^n$, and any vector $x$ can be written as

$$
x = \sum_{i=1}^n \alpha_i v_i
$$

then

$$
T(x) = \sum_{i=1}^n \alpha_i T(v_i).
$$

### Examples

1. $T: \mathbb{R}^2 \to \mathbb{R}^2$  
   $(1, 0) \mapsto (2, 0)$, $(0, 1) \mapsto (0, -1)$  
   Stretches the $x$-direction by a factor of $2$ and flips the $y$-direction.

2. $T: \mathbb{R}^2 \to \mathbb{R}^2$  
   $(1, 0) \mapsto \left(\frac{1}{\sqrt{2}}, \frac{1}{\sqrt{2}}\right)$,  
   $(0, 1) \mapsto \left(-\frac{1}{\sqrt{2}}, \frac{1}{\sqrt{2}}\right)$  
   Rotates the plane by $\pi/4$.

3. $T: \mathbb{R}^2 \to \mathbb{R}$  
   $(1, 0) \mapsto 1$, $(0, 1) \mapsto 0$  
   Orthogonal projection onto the $x$-axis.

4. $T: \mathbb{R} \to \mathbb{R}^2$  
   $1 \mapsto (1, 1)$  
   Maps $\mathbb{R}$ onto the line $y = x$.

---

## Linear Transformations as Matrices

A particularly useful basis for $\mathbb{R}^n$ is the **standard basis**:

$$
e_i =
\begin{pmatrix}
0 \\ \vdots \\ 1 \\ \vdots \\ 0
\end{pmatrix}
$$

with the $1$ in the $i$-th position.

Any vector $x \in \mathbb{R}^n$ can be written as:

$$
x = \sum_{i=1}^n x_i e_i.
$$

If $T(e_i) = a_i \in \mathbb{R}^m$, we can write $a_i$ in the standard basis as:

$$
a_i =
\begin{pmatrix}
a_{1i} \\ a_{2i} \\ \vdots \\ a_{mi}
\end{pmatrix}.
$$

Placing these column vectors $a_i$ side-by-side forms the **matrix representation** of $T$:

$$
A =
\begin{pmatrix}
a_{11} & a_{12} & \cdots & a_{1n} \\
a_{21} & a_{22} & \cdots & a_{2n} \\
\vdots & \vdots & \ddots & \vdots \\
a_{m1} & a_{m2} & \cdots & a_{mn}
\end{pmatrix}.
$$

---

## Matrix-Vector Product

For $x \in \mathbb{R}^n$, we have:

$$
T(x) = \sum_{i=1}^n x_i a_i.
$$

This is a **linear combination of the columns of $A$**, with coefficients given by the entries of $x$.  
The formula for the **matrix-vector product** is:

$$
Ax =
\begin{pmatrix}
\sum_{j=1}^n a_{1j} x_j \\
\sum_{j=1}^n a_{2j} x_j \\
\vdots \\
\sum_{j=1}^n a_{mj} x_j
\end{pmatrix}.
$$

---

## Transpose and Conjugate Transpose

- **Transpose** $A^T$: $(i, j)$ entry of $A^T$ is $a_{ji}$.
- **Conjugate transpose** $A^H$: Take complex conjugate of each entry of $A^T$.

Key property for real matrices:

$$
(Ax)^T y = x^T (A^T y).
$$

For complex matrices:

$$
(Ax)^H y = x^H (A^H y).
$$

---

## Special Matrices

- **Symmetric**: $A^T = A$
- **Skew-symmetric**: $A^T = -A$
- **Hermitian**: $A^H = A$
- **Skew-Hermitian**: $A^H = -A$

Vectors are treated as column vectors ($n \times 1$ matrices). Their transposes are row vectors ($1 \times n$ matrices).