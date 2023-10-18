This approach goes back to the fundamental motivation behind eigenvalue:
$$
A^k = X \Lambda^k X^{-1} = X \Lambda^k Y^T
$$
We can write this product using the [[Outer form of matrix-matrix product|outer form]] of matrix-matrix products:
$$
A^k = X \Lambda^k Y^T = \lambda_1^k \, \boldsymbol x_{,1} \boldsymbol y_{,1}^T + \cdots + \lambda_n^k \, \boldsymbol x_{,n} \boldsymbol y_{,n}^T
$$
Assume that
$$
|\lambda_1| > |\lambda_2| \ge \cdots \ge |\lambda_n|
$$
Then: 
$$
A^k \approx \lambda_1^k \, \boldsymbol x_1 \, \boldsymbol y_1^T
$$
The range of $A^k$ is nearly $\text{span}\{\boldsymbol x_1\}.$

Pick a random vector $\boldsymbol q$ with $\boldsymbol y_1^T \boldsymbol q \neq 0$ and 
$$
A^k \, \boldsymbol q \sim \lambda_1^k \, (\boldsymbol y_1^T \boldsymbol q) \, \boldsymbol x_1
$$
We have that $A^k \boldsymbol{q}$ simply converges to the span of $\boldsymbol{x}_{1}$.

Algorithm pseudo-code:
```julia
while not_converged:
    zk = A * qk
    ev = dot(zk, qk)
    qk = zk / norm(zk)
```
- `ev` converges to the eigenvalue. 
- `qk` does not necessarily converge, but `span(qk)` converges to `span(x1)`.

We can make a more precise statement. Denote by
$$
\lambda_1 = \rho \: e^{i \theta}
$$for a general complex eigenvalue. $\boldsymbol{x}_1$ is a unit eigenvector associated with $\lambda_1$.

Then, because of the step `qk = zk / norm(zk)`, we get
$$
\boldsymbol{q}_k \to e^{ik \theta} \boldsymbol{x}_1
$$
So $e^{-i k \theta} \boldsymbol q_k$ converges to $\boldsymbol{x}_{1}$.

[[Eigenvalues]], [[Eigenvalues cannot be computed exactly]], [[Why eigenvalues?]]