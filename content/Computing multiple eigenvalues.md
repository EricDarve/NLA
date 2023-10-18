This is great but how can we obtain the other eigenvalues as well?

We will now assume that
$$
|\lambda_1| > |\lambda_2| > \cdots > |\lambda_n|
$$
We will see a modification later to address situations when this is not true. A shift of the eigenvalues will be required.

Recall that
$$
A^k = \lambda_1^k \, \boldsymbol x_1 \boldsymbol y_1^T + \cdots + \lambda_n^k \, \boldsymbol x_n \boldsymbol y_n^T
$$
As we multiply by $A$, the vector `qk` converges to `x1`. What happens if we apply a projection orthogonal to `x1`?

That is, let's run the power iteration for $PA$ with 
$$
P = I - x_1 x_1^H
$$
with $\|x_1\|_2 = 1$.

We use the Schur decomposition of $A$.
$$
PA = P Q T Q^H = \sum_{k \ge 2} q_{,k} \, [TQ^H]_{k,}
= \sum_{k \ge 2} q_{,k} \, t_{k,} \, Q^H
$$
Define $\tilde{T}$ such that its first row is 0
$$
\tilde{T}_{1,} = 0
$$
and 
$$
\tilde{T}_{i,} = T_{i,}
$$
when $i \ge 2$. Then:
$$
PA = \sum_{k \ge 1} q_{,k} \, \tilde{t}_{k,} \, Q^H = Q \, \tilde{T} \, Q^H
$$
This is the Schur decomposition of $PA$. This shows that the largest eigenvalue of $PA$ is $\lambda_2$. Its eigenvector is $q_2$.

Therefore, the power iteration applied to $PA$ yields $\lambda_2$ and $q_2$.

We can generalize this process and get all the eigenvalues. To get $\lambda_{i+1}$ we need to apply the power iteration to:
$$
P_{ \{q_1,q_2,\dots, q_i\}^\perp } A
$$
The largest eigenvalue is now $\lambda_{i+1}$ and its eigenvector is $q_{i+1}$.

This is a nice algorithm, but it is fairly complicated. You need to calculate $q_1$, then $q_2$, $q_3$, etc. We will see later on a faster and simpler way to perform this calculation.