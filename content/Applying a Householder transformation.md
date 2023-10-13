We never form the matrix:
$$
P = I - \beta vv^T
$$
with $\beta = \frac{2}{v^Tv}.$

Instead, say we are given a vector $y$ and want to compute $Py$. We can follow these steps:

1. $v^T y$
2. $y = \beta (v^T y) v.$

**Cost:**

Time: $O(n)$; space $O(n)$.

For a matrix, if we want to compute $PA$, the cost goes up to $O(n^2)$.

[[Householder transformation]]