Let's assume we have computed a Schur decomposition $A = Q T Q^H$. How can we obtain the eigenvectors?

Let's focus on the evecs of $T$. For a triangular matrix, finding evecs is relatively straightforward.
$$
T = \begin{pmatrix} T_{11} & T_{12} & T_{13} \\
0 & \lambda & T_{23} \\
0 & 0 & T_{33}
\end{pmatrix}
$$
Assume we have an eigenvalue/eigenvector pair:
$$
Tv = \lambda v
$$
This is the linear system we need to solve to find $v$:
$$
T - \lambda I = \begin{pmatrix} T_{11} - \lambda I & T_{12} & T_{13} \\
0 & 0 & T_{23} \\
0 & 0 & T_{33} - \lambda I 
\end{pmatrix}
$$
and $(T - \lambda I) v = 0.$

Do a back-substitution starting with $v_3$:
$$
(T_{33} - \lambda I) v_3 = 0.
$$
So $v_3 = 0$ assuming that $\lambda$ is not an eval of $T_{33}.$

The next equation is:
$$
0 \; v_2 + T_{23} \; 0 = 0.
$$
So $v_2$ can be chosen arbitrarily. Let's choose $v_2 = 1$.

The last equation is
$$
(T_{11} - \lambda I) v_1 + T_{12} \, 1 + T_{13} \, 0 = 0.
$$
There is a unique solution for $v_1$. The final solution is summarized as:
$$
v = \begin{pmatrix}
- (T_{11}-\lambda I)^{-1} T_{12} \\ 1 \\ 0
\end{pmatrix}
$$
The computational cost is $O(n^2)$ per eigenvector because we only need to solve linear systems with triangular matrices.