Goal: $A = QR$

Summary of the different approaches:
1. Householder/Givens: Find $Q$ such that $Q^T A = R$. Orthogonal triangularization.
2. Gram-Schmidt: Find $R$ such that $A R^{-1} = Q$. **This is preferred when $A$ is tall and thin.** This is a triangular orthogonalization process.

## Square matrix case

Regular QR factorization with a square matrix:

![[Gram-Schmidt 2023-10-15 17.54.00.excalidraw.svg]]

## Thin matrix case

**QR factorization of a thin matrix $A$.** The matrix $Q$ is square. $R$ is thin. $R$ has 0 below the diagonal. This is a transformation typically obtained using Householder or Givens transformations.

![[Drawing 2023-10-15 18.35.05.excalidraw.svg]]

**Thin QR factorization.** $Q$ is thin and has the same size as $A$. $R$ is upper triangular and square. This is a factorization that can be obtained using Gram-Schmidt.

![[Drawing 2023-10-15 18.36.19.excalidraw.svg]]

## Steps in the algorithm

The computational strategy is similar to LU and Cholesky.

We start with the product $QR$ in outer-form:
$$
A = \sum_k q_{,k} \, r_{k,}
$$
The first column is simple. We just need it to be of norm 1. Define
$$
r_{11} = \|a_{,1}\|_2, \qquad
q_{,1} = \frac{a_{,1}}{\|a_{,1}\|_2}
$$
How do we get $r_{1,}$?

Consider column $j$ of $A$:

$$
a_{,j} = r_{1j} \, q_{,1} + \cdots + r_{jj} \, q_{,j}
$$
Use the fact that the $q_{,k}$ are orthogonal.

$r_{1j} = q_{,1}^T \, a_{,j}$

We are projecting $a_{,j}$ onto $q_{,1}$.

Then subtract:

$A \leftarrow A - q_{,1} \, r_{1,}$

The first column of $A$ is now zero. Repeat to get all the other columns.

## Complete algorithm

Assume $A$ is $m \times n$.

Loop: $k$ from 1 to $n$
- $r_{kk} = \|a_{,k}\|_2$
- $q_{,k} = \|a_{,k}\|_2^{-1} \, a_{,k}$
- $r_{kj} = q_{,k}^T \, a_{,j}$, $j > k$.
- $A \leftarrow A - q_{,k} \, r_{k,}$

Each iteration requires $O(mn)$ flops. So the total computational cost is $O(m n^2)$ flops.

[[QR factorization]], [[QR using Householder transformations]], [[QR using Givens transformations]]
