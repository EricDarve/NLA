The [[QR iteration|key step]] becomes computing the [[QR factorization]] of $T_k$ and applying $U_k$ to the right of $R_k$.

This can be made more efficient if $A$ is in [[QR using Givens transformations|upper Hessenberg form]]:

![[2022-10-22-16-30-53.png]]

This transformation is of the form:
$$
Q^T A Q = H
$$
Note how $Q$ is applied to the left and right. This matrix is unitarily similar to $A$. The reason why this transformation is relatively easy to compute is the fact that $H$ is upper Hessenberg, instead of triangular as in the [[Schur decomposition]].

![[2022-10-22-16-31-41.png]]

![[2022-10-22-16-31-47.png]]

![[2022-10-22-16-31-54.png]]

We can repeat this process until $A$ is in upper Hessenberg form. This process is similar to the QR factorization using [[QR using Householder transformations|Householder transformations]]. But note again that $Q$ is applied to the left and right so this is still a different process.

- Each orthogonal transformation costs $O(n^2).$ 
- The [[QR using Householder transformations|total computational cost]] is $O(n^3)$. 
- This is very fast.
- Once this is done, the [[QR iteration]] becomes much faster.