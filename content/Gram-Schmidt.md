Goal: $A = QR$

Summary of the different approaches:
1. Householder/Givens: Find $Q$ such that $Q^T A = R$.
2. Gram-Schmidt: Find $R$ such that $A R^{-1} = Q$. **This is preferred when $A$ is tall and thin.**

Regular QR factorization with a square matrix:

![[Gram-Schmidt 2023-10-15 17.54.00.excalidraw.svg]]

QR factorization of a thin matrix $A$. The matrix $Q$ is square. $R$ is thin. $R$ has 0 below the diagonal. This is a transformation typically obtained using Householder or Givens transformations.

![[Drawing 2023-10-15 18.35.05.excalidraw.svg]]

QR factorization of a thin matrix. $Q$ is thin and has the same size as $A$. $R$ is upper triangular and square. This is a factorization that can be obtained using Gram-Schmidt.

![[Drawing 2023-10-15 18.36.19.excalidraw.svg]]

The computational strategy is similar to LU and Cholesky.

We start with the product $QR$ in outer-form:
$$
A = \sum_k q_{,k} \, r_{k,}
$$
The first column is simple. We just need it to be of norm 1. Define
$$
\begin{gather}
r_{11} = \|a_{,1}\|_2 \\
q_{,1} = \frac{a_{,1}}{\|a_{,1}\|_2}
\end{gather}
$$
How do we get $r_{1,}$?

Consider column $j$ of $A$:

$a_{,j} = r_{1j} q_{,1} + \cdots + r_{jj} q_{,j}$

Use the fact that the $q_{,k}$ are orthogonal.

$r_{1j} = q_{,1}^T \, a_{,j}$

We are projecting $a_{,j}$ onto $q_{,1}$.

Then subtract:

$A \leftarrow A - q_{,1} \, r_{1,}$

The first column of $A$ is now zero. Repeat to get all the other columns.

Summary: $A$ is $m \times n$.

Loop: $k$ from 1 to $n$
- $r_{kk} = \|a_{,k}\|_2$
- $q_{,k} = \|a_{,k}\|_2^{-1} \, a_{,k}$
- $r_{kj} = q_{,k}^T \, a_{,j}$, $j > k$.
- $A \leftarrow A - q_{,k} \, r_{k,}$

Each iteration requires $O(n^2)$ flops. So the total computational cost is $O(n^3)$ flops.

If matrix $A$ is $m \times n$, the cost is $O(mn^2)$.

[[QR factorization]], [[QR using Householder transformations]], [[QR using Givens transformations]]
