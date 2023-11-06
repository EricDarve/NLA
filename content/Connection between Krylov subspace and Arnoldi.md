Using the definition of the [[Krylov subspace]], we form a matrix whose columns are the vectors in the Krylov subspace
$$
K_k = [q_1, Aq_1, A^2 q_1, \ldots, A^{k-1} q_1]
$$

We prove that $K_k = Q_k R_k$, where $Q_k$ is the orthogonal basis from [[Arnoldi process|Arnoldi]]. See also [[Key idea of iterative methods for eigenvalue computation]].

To prove this, we show that $Q^T K_k$ is upper triangular. Consider
$$
Q^T A^{j-1} q_1 = Q^T Q H^{j-1} Q^T q_1
= H^{j-1}e_1
$$
We define matrix $R$ as:
$$
[R_k]_{ij} = H^{j-1}(i,1)
$$
We find that $H^{j-1}(i,1) = 0$ for $j < i$. So matrix $R_k$ is upper triangular.

$Q_k$ is an orthogonal basis of the [[Krylov subspace]] ${\mathcal K}(A,q_1,k).$ 

So $Q_k R_k$ is the [[QR factorization]] of $K_k$.