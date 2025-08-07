[[GMRES]] problem: $\min_y \| b - A Q_k y \|_2$

This is a [[Least-squares problems|least-squares problem]] but the solution is simple.

Recall that from [[Arnoldi process|Arnoldi:]] 
$$
AQ_k = Q_{k+1} Q_{k+1}^T AQ_k
$$
These are the steps to solve the least-squares problem:
$$
\begin{align}
\| b - A Q_k y \|_2
& = 
\| Q_{k+1} Q_{k+1}^T b - Q_{k+1} Q_{k+1}^T A Q_k y \|_2 \\[.3em]
& = 
\| Q_{k+1}^T b - Q_{k+1}^T A Q_k y \|_2
\end{align}
$$
From the definition of the [[Krylov subspace]] and the fact that $Q_{k+1}$ is orthogonal, we get: 
$$
Q_{k+1}^T b = \|b\|_2 \, e_1.
$$
We are already familiar with this matrix
$$
Q_{k+1}^T A Q_k \overset{def}{=} \underline{H}_k
$$

![[GMRES least-squares problem 2023-11-29 10.49.27.excalidraw.svg|400]]

This is the top left block in the [[Upper Hessenberg form for the QR iteration|upper Hessenberg form]] $H = Q^T A Q.$ Its size of $(k+1) \times k$.

The final [[Least-squares problems|least-squares problem]] we need to solve is
$$
\min_y \big\| \; \|b\|_2 \, e_1 - \underline{H}_k \, y  \; \big\|_2
$$