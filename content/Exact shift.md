Let’s see [[QR iteration with shift|what happens]] when $\mu$ is [[Accelerating convergence using a shift|equal]] to an eigenvalue.

The [[Algorithm for QR iteration with shift|first step]] is:
$$
U_{k+1} R_{k+1} = T_k - \mu I
$$

The matrix is singular.

Let's assume that $T_k$ is unreduced; that is, there is no zero on the sub-diagonal. $t_{i+1,i} \neq 0$.

In the [[Algorithm for QR iteration with shift|second step]], the last row of $R_{k+1}$ is 0:
$$
T_{k+1} = R_{k+1} U_{k+1} + \mu I
$$
The [[QR iteration with shift|last row]] is exactly $[0, \dots, 0, \mu]$. Matrix  $T_{k+1}$ is exactly in block upper triangular form. We recover the exact eigenvalue as expected.

In the next step, we can [[Deflation in the QR iteration|deflate]] and work with a smaller matrix.