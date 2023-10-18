Say we want to apply a Householder transformation to a vector $y$ or a matrix $A$, e.g. $Py$ or $PA$. In this case, we do not form the matrix:
$$
P = I - \beta vv^T
$$
with $\beta = \frac{2}{v^Tv},$ and then compute the product. This is very inefficient.

Instead, say we are given a vector $y$ and want to compute $Py$. We can follow these steps:

1. $v^T y$
2. $y = \beta (v^T y) v.$

**Cost:**

Time: $O(n)$; space $O(n)$.

For a matrix, if we want to compute $PA$, the cost goes up to $O(n^2)$.

[[Householder transformation]]