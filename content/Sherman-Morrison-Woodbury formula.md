Assume that $A$ is invertible. If we perturb $A$ slightly, we can easily compute $A^{-1}$.

What do we mean by perturbing slightly?

Let's consider two vectors $u$ and $v$ and define the matrix
$$
B = u v^T
$$
that is $b_{ij} = u_i v_j$. Then, SMW can be used to compute the inverse of $A + B$.
$$
(A + B)^{-1}
= A^{-1} - A^{-1} u (1 + v^T A^{-1} u)^{-1} v^T A^{-1}
$$
The key observation is that $1 + v^T A^{-1} u$ is a scalar. If we have $A^{-1}$, we can compute $x = (A + u v^T)^{-1} b$ following these steps:
1. $A^{-1} b$: a vector
2. $v^T (A^{-1} b)$: a dot product. We get out a scalar.
3. Divide by the scalar $1 + v^T A^{-1} u$: $\alpha = (1 + v^T A^{-1} u)^{-1} \; v^T A^{-1} b$.
4. Multiply $u$ by that scalar: $\alpha u$.
5. Solve again using $A^{-1}$: $A^{-1} (\alpha u)$.
6. Calculate the final result $x$.

[[Invertible matrix]]