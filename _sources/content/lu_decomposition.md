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