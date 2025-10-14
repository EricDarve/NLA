# Givens Rotations

## Introduction

We previously discussed Householder transformations, which use reflections to introduce zeros into a vector or a column of a matrix. We'll now explore another orthogonal transformation technique: **Givens rotations**.

A Givens rotation is a rotation within a two-dimensional plane spanned by two coordinate axes. These were introduced to numerical analysts in the 1950s by Wallace Givens. While Householder transformations are like using a large mirror to reflect an entire vector at once, Givens rotations are more like a series of smaller, targeted adjustments. They rotate a vector in a specific 2D plane to zero out a single element, leaving other elements in that plane changed but preserving the vector's overall length.

Their main purpose, like Householder transformations, is to create zeros in matrices, most commonly as part of a **QR decomposition**. While generally slower for dense matrices, they have significant advantages in specific scenarios, such as parallel computing and processing sparse matrices.

## The Givens Rotation Algorithm

The core of the method is a $2 \times 2$ rotation. Suppose we have a vector $\begin{pmatrix} a \\ b \end{pmatrix}$ and we want to rotate it so that the second component becomes zero. We need to find an orthogonal matrix $G$ such that:

$$
G \begin{pmatrix} a \\ b \end{pmatrix} = \begin{pmatrix} r \\ 0 \end{pmatrix}
$$

The rotation matrix that accomplishes this is:

$$
G = \begin{pmatrix} c & s \\ -s & c \end{pmatrix}
$$

where $c = \cos(\theta)$ and $s = \sin(\theta)$ for some angle $\theta$. Applying this transformation, we get:

$$
\begin{pmatrix} c & s \\ -s & c \end{pmatrix} \begin{pmatrix} a \\ b \end{pmatrix} = \begin{pmatrix} ca + sb \\ -sa + cb \end{pmatrix} = \begin{pmatrix} r \\ 0 \end{pmatrix}
$$

For the second entry to be zero, we must have $cb = sa$. We also know that for the transformation to be a rotation, we need $c^2 + s^2 = 1$. The solution that satisfies these conditions is:

$$
c = \frac{a}{\sqrt{a^2 + b^2}}, \qquad s = \frac{b}{\sqrt{a^2 + b^2}}, \qquad r = \sqrt{a^2+b^2}
$$

Notice that the resulting vector length, $r$, is the original 2-norm of the vector $(a, b)$, as expected for an orthogonal transformation. 

:::{note}
A common convention is to use a counter-clockwise rotation matrix, which would zero out the second element using $s \leftarrow -b/r$. The choice is a matter of convention; both achieve the goal.
:::

To apply this to a larger matrix $A \in \mathbb{R}^{m \times n}$, we embed this $2 \times 2$ rotation into an $m \times m$ identity matrix. A Givens rotation that acts on rows $i$ and $j$ (where $i < j$) has the form:

$$
G(i, j, \theta) =
\begin{pmatrix}
1 & \cdots & 0 & \cdots & 0 & \cdots & 0 \\
\vdots & \ddots & \vdots & & \vdots & & \vdots \\
0 & \cdots & c & \cdots & s & \cdots & 0 \\
\vdots & & \vdots & \ddots & \vdots & & \vdots \\
0 & \cdots & -s & \cdots & c & \cdots & 0 \\
\vdots & & \vdots & & \vdots & \ddots & \vdots \\
0 & \cdots & 0 & \cdots & 0 & \cdots & 1
\end{pmatrix}
\quad
\begin{matrix}
 \\ \\ \leftarrow \text{row } i \\ \\ \leftarrow \text{row } j \\ \\ \\
\end{matrix}
$$

When you multiply a matrix $A$ from the left by $G(i, j, \theta)$, it only affects rows $i$ and $j$ of $A$. The operation rotates the row vectors in the $(i,j)$-plane. By selecting the correct $c$ and $s$ based on the entries $a_{ik}$ and $a_{jk}$ in a specific column $k$, we can introduce a zero at the $(j, k)$ position.

## Numerical Stability

A naive implementation of the formulas for $c$ and $s$ can be numerically unstable. The direct calculation of $r = \sqrt{a^2+b^2}$ can cause **overflow** if $a$ or $b$ are very large, or **underflow** and loss of precision if they are very small.

A more robust method, often implemented in libraries like LAPACK and as the `hypot` function in many languages, avoids this problem. The key is to factor out the element with the larger magnitude.

Let's assume $|a| \ge |b|$. The algorithm is as follows:
1.  If $b=0$, then $c=1, s=0$, and $r=a$.
2.  If $b \ne 0$, compute $t = b/a$. Since $|a| \ge |b|$, we know $|t| \le 1$.
3.  Calculate $r = |a| \sqrt{1+t^2}$. This avoids squaring large or small numbers.
4.  Then, $c = \frac{a}{r} = \frac{a}{|a|\sqrt{1+t^2}} = \frac{\text{sign}(a)}{\sqrt{1+t^2}}$ and $s = \frac{b}{r} = \frac{at}{|a|\sqrt{1+t^2}} = \frac{t \cdot \text{sign}(a)}{\sqrt{1+t^2}}$.

A similar calculation is performed if $|b| > |a|$ by setting $t=a/b$. This stable procedure is crucial for reliable numerical software.

## Computational Cost & Comparison with Householder

### **Dense Matrices**
To perform a QR factorization on a dense $m \times n$ matrix, we must zero out all elements below the main diagonal. This requires applying a sequence of Givens rotations.
* Applying one Givens rotation to an $m \times n$ matrix (i.e., updating two rows of length $n$) costs approximately $6n$ floating-point operations (flops).
* To triangularize the matrix, we need to eliminate roughly $\frac{1}{2} n(2m-n-1)$ entries for $m \ge n$.
* The total cost is approximately $\approx 3mn^2 - n^3$ flops. For a square $n \times n$ matrix, this is $\approx 2n^3$.

In comparison, the Householder QR algorithm costs $\approx 2mn^2 - \frac{2}{3}n^3$ flops ($\approx \frac{4}{3}n^3$ for square matrices). Householder is therefore about twice as fast for dense matrices, not only because of the lower flop count but also because its reliance on matrix-vector products (BLAS-2) is more efficient on modern computer architectures than the vector-scalar operations (BLAS-1) used in Givens rotations.

### **Sparse Matrices**
The true power of Givens rotations becomes apparent with sparse matrices. 
* **Targeted Elimination**: Unlike Householder reflections, which alter an entire column below the diagonal, a Givens rotation only affects two rows. If those rows are sparse, the computational cost is proportional to the number of non-zero elements, not the full row length $n$.
* **Upper Hessenberg Matrices**: For an $n \times n$ upper Hessenberg matrix, there is only one non-zero element to eliminate below the diagonal in each of the first $n-1$ columns. This requires only $n-1$ Givens rotations. Each rotation costs $O(n)$, leading to a total cost of **$O(n^2)$** for the QR factorization—a dramatic improvement over Householder's $O(n^3)$.
* **Tridiagonal Matrices**: The situation is even better for a tridiagonal matrix. Eliminating the single subdiagonal element in a column only affects a handful of elements. The cost for the entire QR factorization is just **$O(n)$**.

:::{note}
An **upper Hessenberg matrix** is a square matrix where all the entries below the first subdiagonal are zero. In other words, for a matrix $A$, if $i > j+1$, then the element $a_{ij} = 0$.

A **tridiagonal matrix** is a square matrix where the only non-zero entries are on the main diagonal, the first diagonal below it (the subdiagonal), and the first diagonal above it (the superdiagonal). For a matrix $A$, if $|i-j| > 1$, then $a_{ij}=0$.
:::

In summary:
| Method      | Dense QR        | Hessenberg QR | Tridiagonal QR | Key Feature                                  |
|-------------|-----------------|---------------|----------------|----------------------------------------------|
| **Householder** | $O(mn^2)$ (faster) | $O(mn^2)$     | $O(mn^2)$      | Reflects and zeros entire sub-columns at once |
| **Givens** | $O(mn^2)$ (slower) | $O(mn)$       | $O(m)$         | Rotates and zeros single elements            |

The choice between the two methods is a classic example of a performance trade-off that depends entirely on the structure of the input matrix.