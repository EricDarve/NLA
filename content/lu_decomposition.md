# The LU Decomposition Algorithm

## Solving Triangular Systems: The Building Block

The strategy of solving $Ax = b$ by factoring $A$ into $LU$ is only useful if solving the resulting triangular systems, $Ly = b$ and $Ux = y$, is significantly easier than solving the original problem. Fortunately, it is. The unique structure of triangular matrices—with zeros on one side of the diagonal—allows us to solve for the unknown variables one by one in a straightforward process.

### Lower Triangular Systems: Forward Substitution

Consider a lower triangular system $Lx = b$. Let's write it out for a small $3 \times 3$ case to see the structure:

$$
\begin{pmatrix}
l_{11} & 0 & 0 \\
l_{21} & l_{22} & 0 \\
l_{31} & l_{32} & l_{33}
\end{pmatrix}
\begin{pmatrix}
x_1 \\
x_2 \\
x_3
\end{pmatrix}
=
\begin{pmatrix}
b_1 \\
b_2 \\
b_3
\end{pmatrix}
$$

Written as a system of equations, this is:
1.  $l_{11}x_1 = b_1$
2.  $l_{21}x_1 + l_{22}x_2 = b_2$
3.  $l_{31}x_1 + l_{32}x_2 + l_{33}x_3 = b_3$

The solution process unfolds naturally:

-   From the first equation, we can immediately solve for $x_1$, as it's the only unknown: $x_1 = b_1 / l_{11}$.
-   Now that we know $x_1$, we can substitute it into the second equation, which now only has one unknown, $x_2$. We can solve for it: $x_2 = (b_2 - l_{21}x_1) / l_{22}$.
-   Finally, knowing $x_1$ and $x_2$, we substitute them into the third equation to find $x_3$.

This sequential process is called **forward substitution** because we solve for the variables in the forward order: $x_1, x_2, \dots, x_n$. The general formula for $x_i$, assuming we have already computed $x_1, \dots, x_{i-1}$, is:

$$x_i = \frac{1}{l_{ii}} \Big( b_i - \sum_{j=1}^{i-1} l_{ij} x_j \Big)$$

This process is well-defined as long as all diagonal entries $l_{ii}$ are non-zero, which is guaranteed if $L$ is invertible.

```python
import numpy as np

def forward_substitution(L: np.ndarray, b: np.ndarray) -> np.ndarray:
    """
    Solves the lower triangular system Lx = b using forward substitution.
    """
    n = L.shape[0]
    x = np.zeros_like(b)
    for i in range(n):
        x[i] = (b[i] - L[i, :i] @ x[:i]) / L[i, i]
    return x
```


### Upper Triangular Systems: Backward Substitution

The exact same logic applies to upper triangular systems of the form $Ux = b$, but the solution is found in the reverse order. Let's look at the $3 \times 3$ case:

$$
\begin{pmatrix}
u_{11} & u_{12} & u_{13} \\
0 & u_{22} & u_{23} \\
0 & 0 & u_{33}
\end{pmatrix}
\begin{pmatrix}
x_1 \\
x_2 \\
x_3
\end{pmatrix}
=
\begin{pmatrix}
b_1 \\
b_2 \\
b_3
\end{pmatrix}
$$

The system of equations is:
1.  $u_{11}x_1 + u_{12}x_2 + u_{13}x_3 = b_1$
2.  $u_{22}x_2 + u_{23}x_3 = b_2$
3.  $u_{33}x_3 = b_3$

Here, we start from the bottom and work our way up:
-   The last equation has only one unknown, $x_3$, which we can solve for immediately: $x_3 = b_3 / u_{33}$.
-   Knowing $x_3$, we can plug it into the second-to-last equation to solve for $x_2$.
-   And so on, until we find $x_1$.

This process is called **backward substitution**. The general formula for computing $x_i$, assuming we have already found $x_{n}, x_{n-1}, \dots, x_{i+1}$, is:

$$x_i = \frac{1}{u_{ii}} \left( b_i - \sum_{j=i+1}^{n} u_{ij} x_j \right)$$


### Computational Cost 💰

Let's analyze the cost of forward substitution. To compute each $x_i$, the formula requires:
-   $i-1$ multiplications ($l_{ij} x_j$)
-   $i-1$ subtractions
-   $1$ division

In numerical linear algebra, we often count a fused multiplication and addition/subtraction as a single **floating-point operation**, or **flop**. Thus, computing $x_i$ takes approximately $i-1$ flops for the sum and one flop for the division, for a total of $i$ flops. To find the total cost for solving the entire system, we sum this over all $i$:

$$\text{Total Flops} = \sum_{i=1}^{n} (i-1) \approx \sum_{i=1}^{n} i = \frac{n(n+1)}{2} = \frac{1}{2}n^2 + O(n)$$
*Note: A more precise count yields approximately $n^2$ flops.*

The cost for backward substitution is identical. The key takeaway is that solving a triangular system of size $n$ costs approximately $n^2$ floating-point operations. We say the complexity is **$O(n^2)$**.

This is remarkably efficient. As we will see, general methods for solving $Ax=b$ (like the LU factorization itself) cost $O(n^3)$ flops. For a large matrix, the difference between $n^2$ and $n^3$ is enormous. The cost of the two triangular solves is so low that it's considered negligible compared to the cost of the initial factorization. This is what makes factorization-based methods so powerful.

Once you have the LU factorization of a matrix $A$, solving the linear system $Ax=b$ becomes a straightforward and efficient two-step process involving only triangular solves.

This method transforms one difficult problem into two simple ones.

## The Two-Step Solution Process

The core idea is to substitute the factorization $A=LU$ into the original equation and strategically group the terms.

1.  **Start with the original system:**

    $$Ax = b$$

2.  **Substitute the factorization:**
   
    $$(LU)x = b$$

3.  **Group the terms.** Using the associative property of matrix multiplication, we can write:
   
    $$L(Ux) = b$$

4.  **Introduce an intermediate vector.** Let's define a temporary vector $z$ such that:
   
    $$z = Ux$$

    Substituting $z$ into the equation from step 3 gives us our first problem:

    $$Lz = b$$

This decouples the original system into two manageable triangular systems:

**Step 1: Solve for $z$ using Forward Substitution**

First, we solve the lower triangular system $Lz = b$ for the intermediate vector $z$. As we've seen, this is computationally inexpensive.

**Step 2: Solve for $x$ using Backward Substitution**

Once we have computed $z$, we solve the upper triangular system $Ux = z$ for our final solution vector $x$. This is also computationally inexpensive.

This two-stage process—forward substitution followed by backward substitution—is the standard method for solving a linear system once its LU factorization is known.

### Computational Cost 💰

The efficiency of this approach is its main advantage. Assuming the LU factorization is already available:

* **Cost of Step 1 (Forward Substitution):** Solving $Lz=b$ takes **$n^2$** flops.
* **Cost of Step 2 (Backward Substitution):** Solving $Ux=z$ also takes **$n^2$** flops.

The **total computational cost** to solve the system is the sum of these two steps: $n^2 + n^2 = \mathbf{2n^2}$ flops. Therefore, the overall complexity is $O(n^2)$.

This provides the core motivation for LU factorization. While finding the factors $L$ and $U$ is an expensive $O(n^3)$ operation, once you have them, you can solve for any right-hand side $b$ very quickly. This is a massive advantage in applications where the same matrix $A$ must be used with many different $b$ vectors.


## Two Views of Matrix Multiplication

The product of two matrices can be viewed as a sum of **outer products**. This perspective is a powerful tool for developing and understanding matrix factorization algorithms.

Let's consider the product of two $n \times n$ matrices, $A = BC$.

**1. The Inner Product View (The Standard Method)**

You're likely most familiar with the "inner product" or "dot product" view. To find the entry $a_{ij}$, you take the dot product of the **$i$-th row** of $B$ with the **$j$-th column** of $C$.

$$a_{ij} = (\text{row } i \text{ of } B) \cdot (\text{column } j \text{ of } C) = \sum_{k=1}^n b_{ik}c_{kj}$$

Here, we compute the final matrix $A$ one scalar entry at a time.

**2. The Outer Product View**

The outer product view reframes the entire calculation. Instead of a sum of scalars, we see the matrix $A$ as a **sum of matrices**. Specifically, it's the sum of the outer products of the columns of $B$ with the corresponding rows of $C$.

$$A = \sum_{k=1}^n (\text{column } k \text{ of } B) (\text{row } k \text{ of } C) = \sum_{k=1}^n b_{:,k} c_{k,:}$$

Each term in this sum, $b_{:,k} c_{k,:}$, is an outer product between a column vector (size $n \times 1$) and a row vector (size $1 \times n$). The result of each outer product is a full $n \times n$ matrix, often called a **rank-one matrix**. The final matrix $A$ is constructed by adding these rank-one matrices together.

### A Concrete Example

Let's see this in action for a simple $2 \times 2$ case where $A = BC$.

$$A = \underbrace{\begin{pmatrix} b_{11} \\ b_{21} \end{pmatrix}}_{b_{:,1}} \underbrace{\begin{pmatrix} c_{11} & c_{12} \end{pmatrix}}_{c_{1,:}} + \underbrace{\begin{pmatrix} b_{12} \\ b_{22} \end{pmatrix}}_{b_{:,2}} \underbrace{\begin{pmatrix} c_{21} & c_{22} \end{pmatrix}}_{c_{2,:}}$$

First, compute the two rank-one matrices:

$$\begin{pmatrix} b_{11}c_{11} & b_{11}c_{12} \\ b_{21}c_{11} & b_{21}c_{12} \end{pmatrix} + \begin{pmatrix} b_{12}c_{21} & b_{12}c_{22} \\ b_{22}c_{21} & b_{22}c_{22} \end{pmatrix}$$

Then, sum them to get the final result:

$$A = \begin{pmatrix} b_{11}c_{11} + b_{12}c_{21} & b_{11}c_{12} + b_{12}c_{22} \\ b_{21}c_{11} + b_{22}c_{21} & b_{21}c_{12} + b_{22}c_{22} \end{pmatrix}$$

As you can see, each entry matches the result from the standard inner product definition.

### Why is This View Important? 🤔

While it might seem more complex, the outer product perspective is crucial for algorithm design. It shows us how a matrix can be built up iteratively. Many factorization algorithms, including LU, are based on the idea of "peeling off" or subtracting these rank-one components from the original matrix one at a time to reveal its underlying structure. We will use this exact idea to derive the LU factorization algorithm.


## Deriving the Algorithm via Outer Products

The LU factorization algorithm elegantly computes the factors $L$ and $U$ by systematically eliminating entries in the matrix $A$. The outer product perspective provides a clear way to understand this process as a sequence of **rank-one updates**.

The core idea is to build $L$ and $U$ iteratively. We start with the outer product formulation of the factorization:

$$A = \sum_{k=1}^n l_{:,k} u_{k,:} = l_{:,1} u_{1,:} + l_{:,2} u_{2,:} + \dots + l_{:,n} u_{n,:}$$

This expresses the matrix $A$ as a sum of rank-one matrices. Our goal is to determine one pair of vectors—a column of $L$ and a row of $U$—at each step.

### Step 1: Determining $l_{:,1}$ and $u_{1,:}$

Let's isolate the first term ($k=1$). Due to the triangular structures of $L$ and $U$:

* The first column of $U$ is $(u_{11}, 0, \dots, 0)^T$.
* The first row of $L$ is $(l_{11}, 0, \dots, 0)$.

This structure simplifies the first column and first row of the product $LU$:

* **First Column:** $a_{:,1} = (LU)_{:,1} = L u_{:,1} = l_{:,1}u_{11}$.
* **First Row:** $a_{1,:} = (LU)_{1,:} = l_{1,:} U = l_{11}u_{1,:}$.

We now have two equations, but more unknowns than constraints. To get a unique solution, we must impose a condition. The standard convention is to require the diagonal entries of $L$ to be 1. This is known as a **Doolittle factorization**.

Setting $l_{11} = 1$:

1.  From $a_{1,:} = l_{11}u_{1,:} = 1 \cdot u_{1,:}$, we immediately get the first row of $U$:
   
    $$u_{1,:} = a_{1,:}$$
2.  From $u_{1,:}$, we know that its first element is $u_{11} = a_{11}$. Substituting this into the column equation $a_{:,1} = l_{:,1} u_{11}$, we can solve for the first column of $L$:
   
    $$l_{:,1} = \frac{a_{:,1}}{u_{11}} = \frac{a_{:,1}}{a_{11}}$$

This step completely determines the first column of $L$ and the first row of $U$. This is only possible if our pivot element $a_{11} \neq 0$, an assumption we will revisit later.

### Step 2: The Rank-One Update

We have now successfully "peeled off" the first term of the outer product sum. We can define an updated matrix, $A^{(1)}$, which represents the remainder of the sum:

$$A^{(1)} = A - l_{:,1} u_{1,:} = \sum_{k=2}^n l_{:,k} u_{k,:}$$

By construction, the first row and first column of the matrix $l_{:,1}u_{1,:}$ are identical to the first row and column of $A$. Therefore, the first row and column of $A^{(1)}$ are zero.

$$
A^{(1)} = 
\begin{pmatrix}
0 & 0 & \dots & 0 \\
0 & & & \\
\vdots & & A' & \\
0 & & & 
\end{pmatrix}
$$

The problem now reduces to finding the LU factorization of the smaller $(n-1) \times (n-1)$ submatrix $A'$ in the lower-right corner.

### The General Algorithm

We can repeat this process iteratively. At each step $k$ (from $1$ to $n-1$), we:

1.  **Identify** the $k$-th column of $L$ and $k$-th row of $U$ using the current updated matrix $A^{(k-1)}$.
2.  **Perform** a rank-one update to compute the next matrix: $A^{(k)} = A^{(k-1)} - l_{:,k} u_{k,:}$.

This process continues until all columns of $L$ and rows of $U$ have been determined.

### Computational Cost 💰

Let's analyze the cost. At each step $k$, the main work is the rank-one update of the lower-right submatrix of size $(n-k) \times (n-k)$. This update requires approximately $2(n-k)^2$ floating-point operations.

To find the total cost, we sum this over all steps:

$$\text{Total Flops} \approx \sum_{k=1}^{n-1} 2(n-k)^2 = 2 \sum_{j=1}^{n-1} j^2 \approx 2 \frac{(n-1)^3}{3} \approx \frac{2}{3}n^3$$

The computational cost of computing the LU factorization is **$O(n^3)$**. This is significantly more expensive than the $O(n^2)$ triangular solves that follow, which is why we separate the factorization and solving stages.


## The Algorithm Step-by-Step

The algorithm constructs the columns of $L$ and the rows of $U$ iteratively, from $k=1$ to $n$. Let's walk through the process for a general step $k$, assuming the first $k-1$ steps are complete. At this point, the first $k-1$ columns of $L$ and rows of $U$ are finalized.

1.  **Determine Row $k$ of $U$**: The $k$-th row of $U$ is simply the $k$-th row of the *current*, modified matrix $A$. The elements $u_{kj}$ for $j < k$ are zero because $U$ is upper triangular.
   
    $$u_{k, k:n} = a_{k, k:n}$$

2.  **Determine Column $k$ of $L$**: The $k$-th column of $L$ is found by taking the $k$-th column of the current matrix $A$ and scaling it by the diagonal element $a_{kk}$, which we call the **pivot**. By convention, we set $l_{kk}=1$.
   
    $$l_{k:n, k} = \frac{a_{k:n, k}}{a_{kk}}$$

3.  **Update the Submatrix (Schur Complement)**: This is the core of the algorithm. We form the outer product of the just-computed vectors ($l_{:,k}$ and $u_{k,:}$) and subtract it from the matrix $A$. This update effectively zeroes out the $k$-th row and column's influence on the rest of the matrix, leaving a smaller problem to solve in the next iteration.
   
    $$A \leftarrow A - l_{:,k} u_{k,:}$$
    Specifically, this update only affects the submatrix to the lower right of the pivot, from row and column $k+1$ to $n$.

This loop continues until all columns of $L$ and rows of $U$ are determined. In a typical implementation, the computed values for $L$ (below the diagonal) and $U$ (on and above the diagonal) are stored directly in the matrix $A$ to save space.

### In-Place LU Factorization (no pivoting)

The factors are stored in A: the strict lower triangle holds L (with unit diagonal), and the upper triangle holds U.

```python
import numpy as np
def lu_inplace(A: np.ndarray) -> np.ndarray:
    """
    Performs in-place LU factorization (Doolittle, no pivoting) and overwrites the input matrix A.
    On exit:
      - L is in the strict lower triangle with implicit unit diagonal
      - U is in the upper triangle (including diagonal)
    Assumes:
      - All pivots are nonzero
      - No pivoting is performed (numerically unstable for some matrices)
    """
    n = A.shape[0]
    for k in range(n-1):
        # Update the k-th column of L
        A[k+1:n, k] /= A[k, k]
        # Rank-one update of the trailing submatrix
        A[k+1:n, k+1:n] -= np.outer(A[k+1:n, k], A[k, k+1:n])
```

### The Pivot Problem: When Things Go Wrong 🚧

The algorithm has a critical weak point: the division by the pivot element $a_{kk}$ in step 2.

The algorithm **fails** if at any step $k$, the pivot element $a_{kk}$ of the *current* matrix is zero. This would require division by zero, bringing the entire process to a halt.