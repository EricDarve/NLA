# Modified Gram-Schmidt

## Introduction

We previously covered the QR factorization using Householder and Givens transformations. These methods fall under the category of **orthogonal triangularization**, where we find an orthogonal matrix $Q$ and apply it to $A$ to produce an upper triangular matrix $R$, such that $Q^T A = R$.

Let's examine a different approach: the **Gram-Schmidt process**. This method performs **triangular orthogonalization**. Instead of finding $Q$ first, we find an upper triangular matrix $R$ such that the columns of $A R^{-1}$ form an orthonormal matrix $Q$. This can be rewritten as our familiar $A = QR$. The Gram-Schmidt process is particularly well-suited for matrices that are **tall and thin** (i.e., $m \gg n$).

GS has an advantage in iterative methods. This is a crucial point of comparison with Householder transformations. While Householder methods are generally more stable, the Gram-Schmidt process produces the orthogonalized vectors sequentially. The j-th vector, $q_j$, is finalized at the j-th iteration. In contrast, Householder reflections produce the full orthogonal matrix $Q$ only at the very end. This sequential nature makes Gram-Schmidt suitable for **iterative methods like the Arnoldi iteration,** which is used for finding eigenvalues of large matrices.

We will focus exclusively on the **Modified Gram-Schmidt (MGS)** algorithm, which is a numerically stable version of the classical process.

## The Modified Gram-Schmidt Algorithm

The core idea of MGS is to build the columns of $Q$ one by one. At each step $k$, we generate the vector $q_k$ by normalizing the current $k$-th column of our working matrix. Then, we immediately make all subsequent columns orthogonal to this new vector $q_k$.

### An Outer-Product Perspective

We will develop MGS from an outer-product point of view. The computational strategy is quite similar to the methods we used for LU and Cholesky factorizations.

We start by viewing the product $A=QR$ as a sum of rank-one matrices:

$$
A = \sum_{k=1}^{n} q_k r_k^T
$$

where $q_k$ are the column vectors of $Q$ and $r_k^T$ are the row vectors of $R$. Our goal is to find one component ($q_k r_k^T$) at each step, subtract it from our matrix, and then repeat the process on the remainder.

### Algorithm Derivation

Let's begin with the first column ($k=1$). The first term in our sum is $q_1 r_1^T$.

1.  **Find $q_1$ and $r_{11}$.**
    From the relationship $A=QR$, we know the first column of $A$ is given by $a_1 = r_{11} q_1$. Since $q_1$ must be a unit vector (${\|q_1\|}_2 = 1$), we can determine $r_{11}$ and $q_1$ directly by normalizing the first column of $A$:

    $$
    r_{11} = \|a_1\|_2, \qquad q_1 = \frac{a_1}{\|a_1\|_2}
    $$

2.  **Find the rest of the first row of R ($r_1^T$).**
    To find the other elements of the first row of $R$ ($r_{12}, r_{13}, \dots, r_{1n}$), we consider the equation for any column $j$ of $A$:

    $$
    a_j = r_{1j} q_1 + r_{2j} q_2 + \cdots + r_{jj} q_j
    $$

    By using the fact that the columns of $Q$ are orthogonal ($q_k^T q_i = 0$ for $k \ne i$), we can find $r_{1j}$ by left-multiplying by $q_1^T$:

    $$
    q_1^T a_j = q_1^T (r_{1j} q_1 + r_{2j} q_2 + \cdots) = r_{1j} (q_1^T q_1) + 0 + \cdots = r_{1j}
    $$

    Thus, we find each element by **projecting column $a_j$ onto our new basis vector $q_1$**. This gives us the complete first row of $R$.

3.  **Subtract the first component and repeat.**
    Now that we have the full first component, $q_1 r_1^T$, we can subtract it from $A$.

    $$
    A^{(1)} = A - q_1 r_1^T
    $$

    This update has a crucial effect: the **first column of the new matrix $A^{(1)}$ is now zero**, and every subsequent column of $A^{(1)}$ is now orthogonal to $q_1$. We can then **repeat the entire process** on $A^{(1)}$ to find $q_2$ and $r_2^T$, and so on.

### The Complete Algorithm

This leads to a simple and elegant algorithm where we iteratively build $Q$ and $R$.

For $k = 1, \dots, n$:
1.  Normalize the $k$-th column:

    $$
    r_{kk} = \|a_k\|_2, \quad
    q_k = a_k / r_{kk}
    $$

2.  Calculate the remaining elements of the $k$-th row of R:

    $$
    r_{kj} = q_k^T a_j \quad (\text{for } j=k+1, \dots, n)
    $$

3.  Update the subsequent columns of the matrix:

    $$
    a_j \leftarrow a_j - q_k r_{kj} \quad (\text{for } j=k+1, \dots, n)
    $$

## Numerical Stability

When implemented in finite-precision arithmetic, different algorithms can produce very different results due to rounding errors. The classical version of Gram-Schmidt is known to be **numerically unstable** because the computed vectors $q_k$ can quickly lose their orthogonality.

The modified version we have just discussed is significantly more stable. The crucial difference is that when we orthogonalize $a_j$ against $q_k$, the vector $a_j$ has *already* been made orthogonal to all previous vectors $q_1, \dots, q_{k-1}$. By continuously "cleaning" the remaining vectors at each step, we prevent the accumulation of orthogonality errors that plague the classical method. This makes MGS a much more reliable tool for practical computation.

## Computational Cost

Let's analyze the cost of the MGS algorithm for an $m \times n$ matrix.

* The outer loop runs $n$ times (for each column).
* Inside the loop for column $k$:
    * Normalizing column $k$ (a vector of length $m$) takes about $2m$ flops.
    * The inner loop runs for the remaining $n-k$ columns. For each of these columns $j$:
        * Calculating the dot product $q_k^T a_j$ costs $O(m)$ flops.
        * The vector update $a_j - r_{kj} q_k$ also costs $O(m)$ flops.
    * The cost of the inner loop is roughly $(n-k) \times O(m)$ flops.

Summing over the outer loop, the total cost is approximately:

$$
\sum_{k=1}^{n} (2m + (n-k) \cdot 2m) \approx \sum_{k=1}^{n} 2m(n-k+1) \approx 2m \frac{n(n+1)}{2} \approx mn^2
$$

The total computational cost is **$O(mn^2)$ flops**. This is asymptotically the same as the Householder QR algorithm.