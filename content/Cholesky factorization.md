- Applies only to [[Symmetric Positive Definite Matrices]].
- [[Triangular factorization|Classical LU]]: $A = LU$. 
- Since $A$ is SPD, we can hope for $A = L L^T$. 
- This is called the **Cholesky** factorization. 
- We will prove later on that $L$ exists for any SPD matrix $A$.

To derive an algorithm to compute $L$, we use the same method as before.
$$
A = \sum_{k=1}^n l_{,k} \, l_{,k}^T
$$
where $L$ is lower triangular. Using the same approach as [[Triangular factorization|LU]], we find that $l_{11} = \sqrt{a_{11}}$ and
$$
l_{,1} = a_{,1} / \sqrt{a_{11}}
$$

**General process:**

Loop for $k = 1$ to $n$
- $l_{,k} = a_{,k} / \sqrt{a_{kk}}$
- $A \leftarrow A - l_{,k} * l_{,k}^T$

See [[LU algorithm]]

Comments:
- Simple algorithm
- Requires half the storage and flops compared to LU
- Requires the pivots $a_{kk}$ to be positive. We will [[Existence of the Cholesky factorization|prove]] later that the pivots are always positive. So this is not an issue.

- [[Row pivoting|Pivoting]] is not required for Cholesky. 
- The algorithm always completes and is very accurate.

[[Triangular factorization]], [[LU algorithm]], [[Symmetric Positive Definite Matrices]]