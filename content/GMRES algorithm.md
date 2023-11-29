These are the computational steps in this algorithm.

1. [[Arnoldi process|Arnoldi step:]] we compute column $k$ of $\underline{H}_k$ and column $k+1$ of $Q_{k+1}$.
2. We update the [[QR factorization|QR decomposition]] of $\underline{H}_k$ with the new column $k$. This requires a new [[QR using Givens transformations|Givens transformation:]]
$$
G_k^T \cdots G_1^T \underline{H}_k = 
\begin{pmatrix}
R_k \\ 0 
\end{pmatrix}
$$
The [[GMRES least-squares problem|least-squares]] problem becomes
$$
\left\| G_k^T \cdots G_1^T (\beta_0 e_1) - 
\begin{pmatrix}
R_k \\ 0 
\end{pmatrix} 
y^{(k)} \right\|_2 
$$
Let's compute
$$
\beta_0 G_k^T \cdots G_1^T e_1 = \begin{pmatrix} p_k \\ \rho_k \end{pmatrix}
$$
Then we have to solve $R_k \, y^{(k)} = p_k.$ The norm of the residual is 
$$
\| \beta_0 e_1 - \underline{H}_k y^{(k)} \|_2 = \rho_k.
$$

3. The [[GMRES]] approximate solution is finally expressed as
$$
x^{(k)} = Q_k \, y_k
$$