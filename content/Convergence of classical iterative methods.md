[[Splitting methods|Recall]] the key iterative step in splitting methods:
$$
M x^{(k+1)} = b + N x^{(k)}
$$
The exact solution to $Ax=b$ satisfies
$$
M x = b + N x
$$
since $A = M-N$. Let's define the error as:
$$
e^{(k)} \overset{def}{=} x^{(k)} - x
$$
From the equations above, we get:
$$
\begin{gather}
M e^{(k+1)} = N e^{(k)} \\
e^{(k)} = M^{-1} N \, e^{(k-1)} = (M^{-1} N)^k \, e^{(0)}
\end{gather}
$$
This leads to the following theorem.

### Theorem
These iterative schemes converge for any initial guess $x^{(0)}$ and $b$ if and only if $\rho(M^{-1} N) < 1$.

$\rho$ denotes the largest eigenvalue in magnitude for a matrix.

### Proof
Denote by $G = M^{-1}N$.

$\rho(G) \ge 1$ implies no convergence when $e^{(0)} =$ eigenvector associated with the eigenvalue that is greater than 1.

Assume that $\rho(G) < 1.$ If $G$ is diagonalizable, then
$$
G^k e_0 = X \Lambda^k X^{-1} e_0
$$
and $\Lambda^k \to 0.$

When $G$ is not diagonalizable, the argument is a little more complex. Consider $Z$ a non-singular matrix, then we can define the following vector norm:
$$
\|x\|_Z = \|Zx\|_\infty
$$
We can define $\|G\|_Z$ as the induced operator norm, and we can prove that
$$
\|G\|_Z = \|Z G Z^{-1}\|_\infty
$$
Consider now the Schur decomposition of $G$:
$$
G = Q T Q^H
$$
We define the following diagonal scaling matrix $D$:
$$
d_{ii} = \delta^{-i}
$$
for some $\delta > 0$. We find that 
$$
[DTD^{-1}]_{ij} = t_{ij} \delta^{j-i}, \quad j \ge i,
$$
and 0 otherwise. Therefore:
$$
\| DTD^{-1} \|_\infty \le \max_i \{ 
|t_{ii}| + n \max_{j>i} \{ |t_{ij}| \delta^{j-i} \} 
\}
$$
We can choose $\delta$ sufficiently small such that
$$
\| DTD^{-1} \|_\infty \le \rho(G) + \epsilon < 1
$$
since the eigenvalues of $G$ are equal to $t_{ii}$. Define $Z = DQ^H$. We have proved that
$$
\| G \|_Z = \|Z G Z^{-1}\|_\infty \le \rho(G) + \varepsilon < 1
$$
Then
$$
\| G^k \|_Z \le \| G \|_Z^k \to 0
$$
since $\| G \|_Z < 1$.

$\square$

- $\rho(M^{-1} N) < 1$ is required for convergence.
- In practice though, numerical convergence may be slow.
- These methods are often used in combination with other techniques (e.g., multigrid) as accelerators.
- Their main advantage is that they can be very cheap to apply and therefore very fast.