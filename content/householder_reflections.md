# Householder Reflections

We now explore methods to solve least-squares problems. The main technique is the **QR factorization**, which decomposes a matrix $A$ into the product of an orthogonal matrix $Q$ and an upper triangular matrix $R$. The most widely used algorithm for computing the QR factorization for dense matrices is based on **Householder transformations**. This method is efficient, numerically stable, and provides a beautiful geometric interpretation of the factorization process.

## The Big Picture: Orthogonal Triangularization

The goal of QR factorization is to decompose an $m \times n$ matrix $A$ into $A = QR$, where $Q$ is an $m \times m$ orthogonal matrix and $R$ is an $m \times n$ upper triangular matrix.

We can rewrite this as $Q^T A = R$. This gives us a new way to think about the process: we are looking for an orthogonal matrix, $Q^T$, that transforms $A$ into an upper triangular matrix $R$. We won't find this $Q^T$ all at once. Instead, we'll build it as a sequence of simpler orthogonal transformations, $Q_k$, applied one after another.

$$Q^T = Q_n \cdots Q_2 Q_1$$

Each transformation $Q_k$ is strategically designed to introduce zeros below the diagonal in the $k$-th column of the matrix, without disturbing the zeros we created in previous columns. The tool we'll use to build these transformations is the Householder reflection.

## The Tool: Householder Reflections

A **Householder transformation** (or Householder reflector) is a matrix that represents a reflection across a hyperplane. 📐

Geometrically, it's a linear transformation that flips a vector space across a plane or hyperplane. The key properties for our purposes are:

1.  It is **orthogonal**. A reflection preserves the lengths of vectors and the angles between them, just in a mirrored way.
2.  It is **symmetric** ($P^T = P$) and **involutory** ($P^2 = I$).
3.  It can transform any given vector $x$ into another vector $y$, as long as they have the same Euclidean norm ($\|x\|_2 = \|y\|_2$).

For QR factorization, we will use a Householder reflection to map a column vector from our matrix onto a standard basis vector, effectively zeroing out most of its components.

## Constructing a Reflection

Let's focus on the first step: zeroing out the subdiagonal elements of the first column of $A$, which we'll call $x = a_1$. Our goal is to find a reflection matrix $P_1$ such that $P_1x$ is a multiple of the first standard basis vector, $e_1$.

$$P_1x = \begin{pmatrix} \sigma \\ 0 \\ \vdots \\ 0 \end{pmatrix} = \sigma e_1$$

Since reflections preserve length, we must have $|\sigma| = \|x\|_2$.

The reflection is performed across a hyperplane. To define this hyperplane, we only need to specify a vector $v$ that is normal to it. Geometrically, the vector normal to the reflection plane must be parallel to the difference between the original vector $x$ and its reflected image $\sigma e_1$.

We can therefore define this normal vector $v$ as the difference between the starting vector and the target vector:

$$
v = x - \sigma e_1
$$

:::{note}
The vector that **bisects the angle** between $x$ and $\sigma e_1$ is $x + \sigma e_1$. This bisecting vector lies *within* the reflection hyperplane and is orthogonal to the normal vector $v$.
:::

The Householder matrix that performs the reflection across the hyperplane orthogonal to $v$ is given by the formula:

$$
P = I - 2 \frac{vv^T}{v^T v}
$$

This is the famous **Householder transformation**. We often write it as $P = I - \beta vv^T$, where $\beta = 2/\|v\|_2^2$.

## The Householder QR Algorithm

The algorithm proceeds column by column.

1.  **Step 1:** For the first column $x = a_1$, we compute its corresponding Householder vector $v_1$ and form the matrix $Q_1 = I - \beta_1 v_1 v_1^T$. We then apply this transformation to the entire matrix $A$:
    $A^{(1)} = Q_1 A$.
    The result is a matrix where the first column is zero below the diagonal.

    $$
    A^{(1)} = \begin{pmatrix} r_{11} & r_{12} & \cdots & r_{1n} \\ 0 & & & \\ \vdots & & A' & \\ 0 & & & \end{pmatrix}
    $$

2.  **Step 2:** Now, we leave the first row and column alone and repeat the process on the smaller submatrix $A'$. We find a Householder vector $v'_2$ for the first column of $A'$. This defines a smaller reflection matrix $Q'_2$. We embed this into the full $m \times m$ identity matrix to form our second transformation, $Q_2$.

    $$
    Q_2 = \begin{pmatrix} 1 & 0 \\ 0 & Q'_2 \end{pmatrix}
    $$
    Applying this yields $A^{(2)} = Q_2 A^{(1)}$, which now has zeros below the diagonal in its first *two* columns.

3.  **Repeat:** We continue this process for $n$ columns (or $m-1$ if $m \le n$), successively creating zeros below the diagonal.

After all steps, we have our upper triangular matrix $R$:
$R = Q_n \cdots Q_2 Q_1 A$

The final orthogonal matrix is the product of all the individual reflections:
$Q = Q_1 Q_2 \cdots Q_n$

## Stability and Accuracy: A Crucial Detail

The numerical stability of the Householder method is one of its greatest strengths, but it hinges on a subtle choice in the definition of $v$.

Recall our definition: $v = x - \sigma e_1$. We have two choices for $\sigma$: $\|x\|_2$ or $-\|x\|_2$. If our vector $x$ is already close to $\sigma e_1$, then this calculation involves subtracting two nearly identical numbers. This is a classic recipe for **catastrophic cancellation**, where we lose significant precision.

**The Solution**: We must choose the sign of $\sigma$ to *avoid* this subtraction. We do this by choosing the sign of $\sigma$ to be the opposite of the sign of $x_1$, the first element of $x$. A robust formula is:

$$
\sigma = - \text{sign}(x_1) \|x\|_2
$$

With this choice, the first component of $v$ becomes $v_1 = x_1 - \sigma = x_1 + \text{sign}(x_1)\|x\|_2$, which is an addition of two numbers of the same sign. This completely avoids the cancellation and makes the algorithm remarkably stable.

Because of this property, Householder QR is **backward stable** and is the standard algorithm for solving dense least-squares problems in high-quality numerical software.

## A Note on Practical Implementation: Never Form $P$ Explicitly

A crucial point for both performance and memory is that we **never** explicitly form the Householder matrix $P = I - \beta vv^T$. Forming this $m \times m$ dense matrix would be computationally wasteful and require unnecessary memory storage. 💻

Instead, we exploit its structure to apply the transformation directly.

### Applying the Transformation Efficiently

To apply the transformation $P$ to a matrix $A$, we compute $PA$ not as a full matrix-matrix product, but as a rank-1 update.

Notice that:

$$
PA = (I - \beta vv^T)A = A - \beta v(v^T A)
$$

This expression gives us a much cheaper way to compute the result:

1.  **Compute the vector-matrix product**: First, calculate the row vector $w^T = v^T A$.
2.  **Compute the outer product**: Then, form the rank-1 matrix $vw^T$.
3.  **Update**: Finally, scale the result by $\beta$ and subtract it from $A$: $A_{new} = A - \beta vw^T$.

This procedure is far more efficient than forming $P$ and then multiplying. For an $m \times n$ matrix $A$, this update costs approximately $2mn$ floating-point operations, whereas forming $P$ and multiplying would cost $O(m^2n)$.

### Storing the Factorization

The full orthogonal matrix $Q = Q_1 Q_2 \dots Q_n$ is also typically not formed explicitly. It's an $m \times m$ dense matrix, and we often don't need it. For solving the least-squares problem $Rx = Q^T b$, we only need to compute the product $Q^T b$.

We can do this efficiently by applying the transformations sequentially:

$$
Q^T b = (Q_n \dots Q_2 Q_1)b
$$

This is just a sequence of cheap Householder updates applied to the vector $b$.

To do this, we only need to store the essential information for each transformation $Q_k$: the **Householder vector $v_k$** and the scalar $\beta_k$. Conveniently, the vector $v_k$ for each step can be stored in the column of $A$ that it is designed to zero out. Since we are creating zeros below the diagonal, this lower-triangular part of $A$ becomes free real estate for storing the Householder vectors. The $\beta_k$ values can be stored in a separate small array.

## Computational Cost

Let’s analyze the computational cost, measured in floating-point operations (flops). Here a "flop" typically refers to one floating-point operation (an addition, subtraction, or multiplication). Note that it is also common to count Mult+Add as one operation, but we will count them separately here.

### Cost per Step (Step $k$)

At step $k$, we are applying a Householder transformation to an $(m-k+1) \times (n-k+1)$ submatrix. Let $m' = m-k+1$ and $n' = n-k+1$.

The update $A' \leftarrow (I - \beta \mathbf{v}\mathbf{v}^T) A'$ is computed as:

1.  **$\mathbf{w}^T = \mathbf{v}^T A'$** (vector-matrix product): This requires $m'n'$ multiplications and $m'(n'-1) \approx m'n'$ additions.
    * *Cost: $\approx 2m'n' \text{ flops}$*

2.  **$A' \leftarrow A' - (\beta \mathbf{v})\mathbf{w}^T$** (outer product update): This requires $m'n'$ multiplications (to form $\beta \mathbf{v}\mathbf{w}^T$) and $m'n'$ subtractions.
    * *Cost: $\approx 2m'n' \text{ flops}$*

The total cost at step $k$ is the sum of these, which is approximately $\mathbf{4(m-k+1)(n-k+1)}$ **flops**.

### Total Cost

To find the total cost, we sum this from $k=1$ to $n$:

$$\text{Total Cost } \approx \sum_{k=1}^{n} 4(m-k+1)(n-k+1) \text{ flops}$$

We can approximate this sum with an integral. We get:

$$\text{Cost}(A=QR) \approx 2 \times \left( mn^2 - \frac{1}{3}n^3 \right) = \mathbf{2mn^2 - \frac{2}{3}n^3} \text{ flops}$$

**Two important special cases emerge from this formula:**

1.  **Square Matrix ($m=n$):** The cost is $\approx 2n^3 - \frac{2}{3}n^3 = \mathbf{\frac{4}{3}n^3}$ flops. This is about twice the cost of performing an LU factorization, which is $\frac{2}{3}n^3$ flops.

2.  **Tall and Skinny Matrix ($m \gg n$):** The $2mn^2$ term dominates, and the cost is approximately $\mathbf{2mn^2}$ flops.

### Cost of Using the Factors

Once we have the factorization, we use the stored Householder vectors to apply $Q$ or $Q^T$ as needed.

**Applying $Q^T$ to a vector (e.g., forming $Q^T b$).** This is the most critical operation for solving the least-squares problem. We need to compute $y = Q^T b$. We do this by applying the transformations sequentially: $y = Q_n \dots Q_2 Q_1 b$.

-   Applying $Q_1$ to $b$ costs $\approx 2m$ flops.
-   Applying $Q_2$ to the result costs $\approx 2(m-1)$ flops.
-   ...and so on, down to $Q_n$, which costs $\approx 2(m-n+1)$ flops.

The total cost is the sum of an arithmetic series:

$$
\text{Cost}(Q^T b) \approx \sum_{k=1}^{n} 2(m-k+1) \approx 2mn - n^2 \text{ flops}
$$

For $m \ge n$, this is an $O(mn)$ operation. This is significantly cheaper than the factorization itself, which is $O(mn^2)$.

### Summary of Costs

Here is a quick reference for an $m \times n$ matrix with $m \ge n$:

| Operation                                     | Leading Term Flop Count           | Big O Notation      |
| --------------------------------------------- | --------------------------------- | ------------------- |
| **Factorization ($A \rightarrow Q, R$)** | $2mn^2 - \frac{2}{3}n^3$          | $O(mn^2)$           |
| **Apply $Q^T$ to vector $b$ ($Q^T b$)** | $2mn - n^2$                       | $O(mn)$             |