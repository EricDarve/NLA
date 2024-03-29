The [[GMRES algorithm]] becomes more efficient when $A$ is symmetric. This leads to the MINRES algorithm. The initial steps of MINRES are similar to those of GMRES.

In MINRES, we also aim to minimize $\|r^{(k)}\|_2$. The primary distinction between [[GMRES]] and MINRES is analogous to the difference between the [[Arnoldi process|Arnoldi]] and [[Lanczos process|Lanczos]] processes.

Recall that, following [[Lanczos process|Lanczos,]] $A Q = Q T$ where $T$ is tri-diagonal.

We get
$$
A Q_k = Q_{k+1} Q_{k+1}^T A Q_k = Q_{k+1} \underline{T}_k
$$
$\underline{T}_k$ has size $(k+1) \times k$ and is "symmetric" tri-diagonal.

![[MINRES 2023-12-04 08.40.09.excalidraw.svg]]

[[GMRES algorithm|As before,]] the quantity we want to minimize becomes
$$
\| b - A Q_k y^{(k)} \|_2 = \| \beta_0 Q_{k+1} e_1 - Q_{k+1} \underline{T}_k y^{(k)} \|_2 
= \| \beta_0 e_1 - \underline{T}_k y^{(k)} \|_2
$$

[[GMRES algorithm|As before,]] we use Givens rotations
$$
G_k^T \cdots G_1^T \underline{T}_k = 
\begin{pmatrix}
R_k \\ 0 
\end{pmatrix}
$$ 
Since $\underline{T}_k$ is tri-diagonal, $R_k$ is upper triangular with only two upper diagonals.

As in [[GMRES algorithm|GMRES]], we then solve the least-squares problem using the $G_k$s and $R_k$.

Because of the tri-diagonal structure of $\underline{T}_k$ and the structure of $R_k$ with only two upper diagonals, we can efficiently transition from step $k-1$ to $k$ in MINRES with a [[Space and time costs of CG and GMRES|time cost]] per iteration of $O({\rm nnz} + n)$ instead of $O({\rm nnz} + kn)$. 

This is the same cost as [[Space and time costs of CG and GMRES|CG]] but for general symmetric matrices instead of [[Symmetric Positive Definite Matrices|SPD.]]

The efficient implementation of this algorithm is the MINRES algorithm of [Paige and Saunders](https://web.stanford.edu/group/SOL/software/minres/).

When $A$ is [[Symmetric Positive Definite Matrices|symmetric positive definite,]] we get an error bound similar to [[Convergence of the Conjugate Gradients|CG]]
$$
\|r^{(k)}\|_2 \le 2 \, \| r^{(0)} \|_2 \; \Big( \frac{\sqrt{\kappa}-1}{\sqrt{\kappa}+1} \Big)^k
$$