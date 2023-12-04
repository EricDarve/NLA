How to solve triangular systems.

![[Pasted image 20230915170125.png|500]]

$Lx = b$ with $L$ lower triangular. Algebraic equation:
$$
\sum_{j=1}^i l_{ij} x_j = b_i
$$
Solve for $x_i$:
$$
\begin{gather}
l_{ii} x_i = b_i - \sum_{j=1}^{i-1} l_{ij} x_j \\
x_i = l_{ii}^{-1} \, (b_i - \sum_{j=1}^{i-1} l_{ij} x_j)
\end{gather}
$$
We can solve iteratively starting from $x_1$, then $x_2$, $x_3$, ..., $x_n$.

This is possible because $L$ is lower triangular.

This process is fast and accurate.

The same process applies to upper triangular matrices: $Ux = b$.
$$
\begin{gather}
\sum_{j=i}^n u_{ij} x_j = b_i \\
x_i = u_{ii}^{-1} \, (b_i - \sum_{j=i+1}^n u_{ij} x_j)
\end{gather}
$$
We solve iteratively in the reverse order: $x_n$, $x_{n-1}$, ..., $x_1$.

The computational cost of solving a triangular system is $O(n^2)$.