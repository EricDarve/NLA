# Cholesky Factorization

For symmetric positive definite (SPD) matrices, the LU factorization simplifies to the Cholesky factorization. It is twice as fast as LU, uses half the memory, and is numerically stable without pivoting.

## Symmetric Positive Definite Matrices

````{prf:definition} Symmetric Positive Definite Matrix (SPD)
$A$ real, square matrix $A ∈ \mathbb R^{n×n}$ is SPD if it satisfies:

- Symmetry: $A^T = A$
- Positive definiteness: $x^T A x > 0$ for all nonzero $x ∈ R^n$
````
For symmetric $A$, there is an orthogonal eigendecomposition $A = Q \Lambda Q^T$ with real eigenvalues. 

If $A$ is SPD, then all eigenvalues $\lambda_i > 0$. Proof: let $q_i$ be an eigenvector associated with $\lambda_i$; then

$$
x = q_i \implies x^T A x = q_i^T (\lambda_i q_i) = \lambda_i \|q_i\|^2 > 0 \Rightarrow \lambda_i > 0.
$$

## From LU to Cholesky

We start with a hand-wavey argument for why the Cholesky factorization exists for SPD matrices. The details are in the next section.

The hand-wavey justification is based on the **uniqueness** of the LU factorization without pivoting. For a matrix where LU factorization exists without pivoting, the Doolittle factorization ($L$ has ones on the diagonal) is unique. Since $A$ is symmetric ($A=A^T$), its factorization must honor that symmetry. If we factor $A$ to get $LU$, we could also factor $A^T$ to get $U^T L^T$. But since $A = A^T$, their unique factorizations must be one and the same. This implies that the factors $L$ and $U$ can't be independent entities; they must be related through a transpose operation, leading us toward a symmetric factorization like $LL^T$.

One consequence is that in $LL^T$, the diagonal of $L$ will no longer be all ones. Instead, we will require the diagonal entries of $L$ to be positive to ensure uniqueness.

Let's now derive the Cholesky factorization from LU. For SPD matrices, we seek a lower triangular L with positive diagonal such that:

$$
A = L L^T.
$$

As for the LU factorization, we use the outer-product approach. This means that we compute $L$ one column at a time, updating the trailing submatrix at each step.

A convenient derivation uses block matrices. Let's partition

$$
A =
\begin{pmatrix}
a_{11} & \mathbf{a}_{21}^T \\
\mathbf{a}_{21} & A_{22}
\end{pmatrix}, \quad
L =
\begin{pmatrix}
l_{11} & \mathbf{0}^T \\
\mathbf{l}_{21} & L_{22}
\end{pmatrix}.
$$

Then

$$
L L^T =
\begin{pmatrix}
l_{11}^2 & l_{11}\mathbf{l}_{21}^T \\
l_{11}\mathbf{l}_{21} & \mathbf{l}_{21}\mathbf{l}_{21}^T + L_{22} L_{22}^T
\end{pmatrix}
=
\begin{pmatrix}
a_{11} & \mathbf{a}_{21}^T \\
\mathbf{a}_{21} & A_{22}
\end{pmatrix}.
$$

Equating blocks gives:

$$
\begin{aligned}
l_{11} &= \sqrt{a_{11}},\\
\mathbf{l}_{21} &= \frac{1}{l_{11}}\,\mathbf{a}_{21},\\
L_{22} L_{22}^T &= A_{22} - \mathbf{l}_{21}\mathbf{l}_{21}^T.
\end{aligned}
$$

This yields a recursive algorithm: compute the first column of L, update the trailing submatrix, and repeat on the (n−1)×(n−1) subproblem.

## The Cholesky Algorithm

Below is the complete algorithm for the Cholesky factorization. It overwrites the lower triangle of $A$ with $L$. This is the in-place outer-product version.

```python
import numpy as np

def cholesky_in_place(A: np.ndarray) -> np.ndarray:
    """
    In-place, outer-product Cholesky factorization.
    Overwrites the lower triangle of A with L such that A_original = L @ L.T.
    A: symmetric positive definite matrix (n × n). Only the lower triangle is used.
    """
    n = A.shape[0]
    for k in range(n-1):
        # Compute the diagonal element
        A[k, k] = np.sqrt(A[k, k])
        
        # Compute column k below the diagonal
        A[k+1:, k] /= A[k, k]
        # Rank-1 update of the trailing submatrix
        A[k+1:, k+1:] -= np.outer(A[k+1:, k], A[k+1:, k])

    # Last diagonal element    
    A[n-1, n-1] = np.sqrt(A[n-1, n-1])
```
## Algorithm Properties

**Computational Cost**

The cost of the Cholesky factorization is dominated by the rank-1 update of the trailing submatrix inside the main loop. In step $k$, this update operates on a matrix of size $(n-k-1) \times (n-k-1)$ and costs approximately $(n-k-1)^2$ multiplications and $(n-k-1)^2$ additions. Summing over $k$ from $0$ to $n-2$:

$$
\text{Flops} \approx \sum_{k=0}^{n-2} (n-k-1)^2 = \sum_{j=1}^{n-1} j^2 \approx \frac{n^3}{3}
$$

This is half the computational cost of LU factorization, which requires approximately $2n^3/3$ floating-point operations.

**Storage Requirements**

Since the input matrix $A$ is symmetric, only its lower (or upper) triangular part needs to be stored, which amounts to $\frac{n(n+1)}{2}$ elements. The algorithm can be implemented in-place, overwriting the lower triangle of $A$ with the Cholesky factor $L$. This means no significant additional storage is required, making the method highly memory-efficient.

**Numerical Stability**

The Cholesky factorization is numerically stable for SPD matrices without any need for pivoting. As shown in the final section, the elements of $L$ are bounded: $|l_{ik}| ≤ \sqrt{a_{ii}}$. This lack of element growth prevents the amplification of rounding errors, which is a key advantage over LU factorization where pivoting is generally necessary for stability.

## Existence and Uniqueness

The previous derivation assumed that all the pivots are positive. This is true for SPD matrices, as we now prove.

````{prf:theorem} Existence and Uniqueness of Cholesky Factorization
If $A ∈ \mathbb R^{n×n}$ is SPD, then there exists a unique lower triangular matrix $L$ with positive diagonal such that $A = LL^T$.
````

### Formal proof

Now we arrive at the formal proof. It's a beautiful argument that uses induction to show that the recursive algorithm we just derived will always succeed for any SPD matrix.

````{prf:proof}

The proof proceeds by induction on the dimension $n$ of the matrix.

**Existence**

Base Case ($n=1$):

If $A$ is a $1 \times 1$ SPD matrix, then $A = [a_{11}]$. By the definition of positive definiteness (with $x=e_1=[1]$), we have 

$$a_{11} = e_1^T A e_1 > 0.$$

We can therefore set $l_{11} = \sqrt{a_{11}}$, which is lower triangular and has a positive diagonal. The factorization $A=LL^T$ holds trivially.

Inductive Step:

Assume that every $(n-1) \times (n-1)$ SPD matrix has a unique Cholesky factorization. Now, let's consider an $n \times n$ SPD matrix $A$. We partition $A$ into a block form, separating out the first row and column:

$$
A =
\begin{pmatrix}
a_{11} & \mathbf{c}^T \\
\mathbf{c} & B
\end{pmatrix}
$$

where $a_{11}$ is a scalar, $\mathbf{c} \in \mathbb{R}^{n-1}$, and $B$ is an $(n-1) \times (n-1)$ symmetric matrix. As in the base case, $a_{11} = \mathbf{e}_1^T A \mathbf{e}_1 > 0$.

The crux of the proof is to show that the trailing submatrix after one step of factorization, the **Schur complement** of $a_{11}$, is also SPD. It is equal to:

$$S = B - \frac{1}{a_{11}} \mathbf{c}\mathbf{c}^T$$

First, $S$ is clearly symmetric since $B$ is symmetric. 

To prove $S$ is positive definite, we use the following construction. For any non-zero vector $\mathbf{y} \in \mathbb{R}^{n-1}$, we want to show that $\mathbf{y}^T S \mathbf{y} > 0$. Let's construct a special vector $\mathbf{x} \in \mathbb{R}^n$ of the form 

$$\mathbf{x} = \begin{pmatrix} x_1 \\ \mathbf{y} \end{pmatrix}.$$ 

Since $A$ is SPD, we know $\mathbf{x}^T A \mathbf{x} > 0$.

$$\mathbf{x}^T A \mathbf{x} = \begin{pmatrix} x_1 & \mathbf{y}^T \end{pmatrix} \begin{pmatrix} a_{11} & \mathbf{c}^T \\ \mathbf{c} & B \end{pmatrix} \begin{pmatrix} x_1 \\ \mathbf{y} \end{pmatrix} = a_{11} x_1^2 + 2 x_1 (\mathbf{c}^T \mathbf{y}) + \mathbf{y}^T B \mathbf{y}$$

The magic happens when we choose $x_1$ specifically to simplify this expression. To do that, we need to develop some intuition. The expression is a quadratic form in $x_1$. This suggests completing the square to see what happens when we vary $x_1$.

Factor out $a_{11}$ from the first part:

$$= a_{11} \left( x_1^2 + 2 x_1 \frac{\mathbf{c}^T \mathbf{y}}{a_{11}} \right) + \mathbf{y}^T B \mathbf{y}$$

Now, let's complete the square inside the parenthesis. We get:

$$\mathbf{x}^T A \mathbf{x} = a_{11} \left( x_1 + \frac{\mathbf{c}^T \mathbf{y}}{a_{11}} \right)^2 + \mathbf{y}^T \left( B - \frac{1}{a_{11}}\mathbf{c}\mathbf{c}^T \right) \mathbf{y}$$

This beautifully decomposes the quadratic form. Our goal is to relate $\mathbf{x}^T A \mathbf{x}$ to the Schur complement term. The most direct way to do that is to make the first term disappear. The first term, being a scaled square, is zero if and only if:

$$x_1 + \frac{\mathbf{c}^T \mathbf{y}}{a_{11}} = 0 \implies x_1 = -\frac{\mathbf{c}^T \mathbf{y}}{a_{11}}$$

With this choice of $x_1$, we have:

$$
\mathbf{x}^T A \mathbf{x} = \mathbf{y}^T \left( B - \frac{1}{a_{11}}\mathbf{c}\mathbf{c}^T \right) \mathbf{y} = \mathbf{y}^T S \mathbf{y}.
$$

Since $\mathbf{y} \neq \mathbf{0}$, our constructed vector $\mathbf{x}$ is also non-zero, which means $\mathbf{x}^T A \mathbf{x} > 0$. Therefore, we have shown that $\mathbf{y}^T S \mathbf{y} > 0$ for any non-zero $\mathbf{y}$, so $S$ is SPD.

By our inductive hypothesis, the $(n-1) \times (n-1)$ SPD matrix $S$ has a Cholesky factorization $S = L_S L_S^T$. We can now construct the Cholesky factor for $A$:

$$
L =
\begin{pmatrix}
\sqrt{a_{11}} & \mathbf{0}^T \\
\frac{1}{\sqrt{a_{11}}}\mathbf{c} & L_S
\end{pmatrix}
$$

This matrix is lower triangular with a positive diagonal. Let's verify that it works:

$$
LL^T =
\begin{pmatrix}
\sqrt{a_{11}} & \mathbf{0}^T \\
\frac{1}{\sqrt{a_{11}}}\mathbf{c} & L_S
\end{pmatrix}
\begin{pmatrix}
\sqrt{a_{11}} & \frac{1}{\sqrt{a_{11}}}\mathbf{c}^T \\
\mathbf{0} & L_S^T
\end{pmatrix}
=
\begin{pmatrix}
a_{11} & \mathbf{c}^T \\
\mathbf{c} & \frac{1}{a_{11}}\mathbf{c}\mathbf{c}^T + L_S L_S^T
\end{pmatrix}
$$

Since 

$$L_S L_S^T = S = B - \frac{1}{a_{11}}\mathbf{c}\mathbf{c}^T,$$

the bottom-right block becomes 

$$\frac{1}{a_{11}}\mathbf{c}\mathbf{c}^T + \left(B - \frac{1}{a_{11}}\mathbf{c}\mathbf{c}^T\right) = B.$$

Thus, $LL^T = A$, and the existence is proven.

**Uniqueness**

Suppose we have two Cholesky factorizations $A = L_1 L_1^T$ and $A = L_2 L_2^T$, where $L_1$ and $L_2$ are lower triangular with positive diagonals.

$$L_1 L_1^T = L_2 L_2^T$$

Since $L_2$ has a positive diagonal, its determinant is non-zero, so it is invertible. We can write:

$$L_2^{-1} L_1 = L_2^T (L_1^T)^{-1} = (L_1^{-1} L_2)^T$$

Let's analyze the matrix $M = L_2^{-1} L_1$.

- The inverse of a lower triangular matrix is lower triangular. The product of two lower triangular matrices is also lower triangular. Therefore, $M$ is **lower triangular**.
- The right-hand side, $(L_1^{-1} L_2)^T$, is the transpose of a lower triangular matrix, making it **upper triangular**.

The only way a matrix can be both lower and upper triangular is if it is a **diagonal matrix**. Let's call it $D.$ So, $M=D$. This implies $L_1 = L_2 D$. 

Substituting this back into the original equality:

$$(L_2 D)(L_2 D)^T = L_2 L_2^T \implies L_2 D D^T L_2^T = L_2 L_2^T$$

Since $D$ is diagonal, $D=D^T$, so we have 

$$L_2 D^2 L_2^T = L_2 L_2^T.$$

We can left-multiply by $L_2^{-1}$ and right-multiply by $(L_2^T)^{-1}$ to get:

$$D^2 = I$$

This means the diagonal entries of $D$ must satisfy $d_{ii}^2=1$, so $d_{ii} = \pm 1$. However, the diagonal entries of $L_1$ and $L_2$ are both positive by definition. The diagonal entries of $D$ are given by $d_{ii} = [L_1]_{ii} / [L_2]_{ii}$. Since both are positive, their ratio must be positive.

$$d_{ii} = +1$$

Therefore, $D$ must be the identity matrix $I$. From $L_1 = L_2 D$, we conclude that $L_1=L_2$, proving uniqueness.
````

## The Remarkable Stability of Cholesky Factorization

One of the most powerful and celebrated features of the Cholesky factorization is its unconditional numerical stability. For general matrices, we saw that LU factorization can be disastrously unstable without a proper pivoting strategy. For SPD matrices, the Cholesky factorization requires no pivoting at all and is guaranteed to be stable. This is a remarkable property.

The reason for this stability is beautifully simple and can be seen directly from the factorization equation itself.

From the relation $A = LL^T$, the formula for the $i$-th diagonal entry of $A$ is the dot product of the $i$-th row of $L$ with itself:

$$a_{ii} = \sum_{k=1}^{n} l_{ik} (L^T)_{ki} = \sum_{k=1}^{i} l_{ik}^2 = l_{i1}^2 + l_{i2}^2 + \dots + l_{ii}^2$$

The sum only goes up to $k=i$ because $L$ is lower triangular.

This simple equation tells us something profound. Since every term in the sum, $l_{ik}^2$, is non-negative, the sum must be greater than or equal to any of its individual terms.

$$l_{ik}^2 \le \sum_{j=1}^{i} l_{ij}^2 = a_{ii} \quad \text{for any } k \le i$$

Taking the square root gives us a tight bound on the magnitude of every element in the Cholesky factor $L$:

$$|l_{ik}| \le \sqrt{a_{ii}}$$

This is a powerful result. It guarantees that the elements of $L$ can never grow large. The magnitude of any entry in $L$ is bounded by the square root of the corresponding diagonal entry of the original matrix $A$.

In numerical analysis, an algorithm is considered stable if it does not amplify the rounding errors that are inherent in floating-point arithmetic. The absence of element growth is the hallmark of a stable factorization. See for example our previous section ({prf:ref}`thm:backward_error_lu`) on the stability of the LU decomposition.

Because the entries of $L$ remain controlled, the algorithm is not susceptible to catastrophic cancellation or error amplification.