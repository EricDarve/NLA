Remember [[Upper Hessenberg form for the QR iteration|one of the steps]] in computing eigenvalues
$$
Q^T A Q = H
$$
This was the step where $A$ is transformed to an [[Upper Hessenberg form for the QR iteration|upper Hessenberg matrix]] to reduce the computational cost of the [[QR iteration for upper Hessenberg matrices|QR iteration algorithm.]]

Let's rewrite the equality as
$$
AQ = QH
$$
where $H$ is upper Hessenberg.

Upon inspection, once $q_1$ is chosen, the other vectors in $Q$ are fixed:
$$
\begin{gather}
(q_1, \dots, q_j) \to q_{j+1} \\
A q_j = \sum_{k=1}^{j+1} h_{kj} \, q_k
\end{gather}
$$
This algorithm is very similar to the [[Gram-Schmidt]] algorithm for the [[QR factorization]] but it uses $Aq_j$.

We get column $j$ of $h$
$$
h_{kj} = q_k^T A q_j, \quad k \le j
$$
Let's [[Gram-Schmidt|project]] $A q_j$ on $q_k$:
$$
h_{j+1,j} \, q_{j+1} = A q_j - \sum_{k=1}^j h_{kj} \, q_k
$$
$h_{j+1,j}$ is the norm of the right-hand-side vector; once we have $h_{j+1,j}$, we get:
$$
 \, q_{j+1} = h_{j+1,j}^{-1} \, \Big( A q_j - \sum_{k=1}^j h_{kj} \, q_k \Big)
$$

### Iteration

Start from $q_1$. Then:

1. Compute $A q_j$. This requires a [[Motivation of iterative methods for eigenvalue computation|sparse matrix-vector product.]]
2. $h_{kj} = q_k^T A q_j$, $1 \le k \le j$.
3. $h_{j+1,j} \, q_{j+1} = A q_j - \sum_{k=1}^j h_{kj} \, q_k$

- If you continue this iterative process to the end, you get a [[Motivation of iterative methods for eigenvalue computation|dense matrix]] $H$.
- The sparsity advantage is lost.
- But what happens if we stop before the end?
- By using the top left $k \times k$ block of $H$ we will be able to approximately compute eigenvalues and solve linear systems. As long as $k$ remains small, this can be computationally very fast.

### Details of iterative process

![[2022-10-28-15-38-28.png]]

If we stop at step $k$, we get the following equation: ^d63511
$$
Q_k H_k + h_{k+1,k} \, q_{k+1} e_k^T = A Q_k
$$
and
$$
Q_k H_k \approx A Q_k
$$
if $h_{k+1,k}$ is small. This suggests using the following matrix $H_k$
$$
H_k = Q_k^T A Q_k
$$
instead of $A$. 

From that point, iterative algorithms replace a calculation using $A$ with a calculation using $H_k$. $H_k$ is a small matrix, and "direct" methods are efficient and fast.

The main application cases are the conjugate gradients method (which we will cover later) and the Arnoldi process.
