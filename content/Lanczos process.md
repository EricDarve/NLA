- Much stronger results are available for symmetric matrices.
- The algorithm is called Lanczos.
- Lanczos: $T = Q^T A Q$
- Symmetric + upper Hessenberg = tri-diagonal

$$
T =
\begin{pmatrix}
\alpha_1 & \beta_1 & \cdots & & 0 \\
\beta_1 & \alpha_2 & \cdots \\
\vdots & & \ddots & & \vdots \\
& & & \ddots & \beta_{n-1} \\
0 & & \cdots & \beta_{n-1} & \alpha_n
\end{pmatrix}
$$

As a result, the steps in Lanczos algorithm are simplified compared to [[Key idea of iterative methods for eigenvalue computation|Arnoldi:]]

1. $\alpha_k = \boldsymbol q_{k}^T A \boldsymbol q_k$
2. $\boldsymbol r_k = A \boldsymbol q_k - \beta_{k-1} \boldsymbol q_{k-1} - \alpha_k \boldsymbol q_k$; $\beta_{k-1}$ is from the previous iteration.
3. $\beta_k \, \boldsymbol q_{k+1} = \boldsymbol r_k$

See [[Algorithm for the Arnoldi process]], [[Arnoldi process]].