The simplest method is the Jacobi iteration.
$$
\begin{gather}
A = D - L - U, \quad M = D \\
D x^{(k+1)} = b + (L+U) x^{(k)}
\end{gather}
$$
- We can check that it follows the general framework of [[Splitting methods|splitting methods.]]
- Very easy to implement
- Fast iteration

### Code

```julia
for k = 1:maxiter
    xnew = zeros(n)
    for i = 1:n
        for j = 1:n
            if j != i
                xnew[i] += A[i,j] * x[j]
            end
        end
        xnew[i] = (b[i] - xnew[i]) / A[i,i]
    end
    x = xnew
end
```

### Convergence of Jacobi

If $A$ is strictly row diagonally dominant
$$
|a_{ii}| > \sum_{j \neq i} |a_{ij}|
$$
then Jacobi [[Convergence of classical iterative methods|converges]]. We have a similar result if $A$ is strictly column diagonally dominant.

- Jacobi is very easy to implement on parallel computers because $D$ is diagonal (or block diagonal).
- Solution step is parallel, i.e., entries of $x^{(k+1)}$ can be computed concurrently using independent computing cores:
$$
D x^{(k+1)} = b + (L+U) x^{(k)}
$$

**Proof**

The iteration matrix in Jacobi is
$$
G = D^{-1} (L+U)
$$
Assume that $A$ is row-diagonally dominant. Then from the definition of the operator $\infty$ norm, we find that
$$
\|G\|_\infty < 1
$$
This proves the convergence of Jacobi since
$$
\|G^k\|_\infty \le \|G\|_\infty^k \to 0
$$
Assume now that $A$ is column-diagonally dominant. Then:
$$
G^k = D^{-1} [ (L+U) D^{-1} ]^k D
$$
But
$$
\|(L+U) D^{-1}\|_1 < 1
$$
So
$$
\Big \| \big((L+U) D^{-1} \big)^k \Big \|_1 \le \|(L+U) D^{-1}\|_1^k \to 0
$$
and $\|G^k\|_1 \to 0$.

There is also a proof that uses the Gershgorin Circle Theorem.

**Gershgorin Circle Theorem.** For any $n \times n$ complex matrix $A$, each eigenvalue of $A$ lies within at least one of the Gershgorin discs in the complex plane. The $i$-th Gershgorin disc is centered at $a_{ii}$ with radius $R_i = \sum_{j \neq i} |a_{ij}|$. Mathematically, this is expressed as:
$$
| \lambda - a_{ii} | \leq R_i 
$$
where $\lambda$ represents an eigenvalue of $A$.

**Convergence proof using Gershgorin's Theorem.** $G = D^{-1} (L+U)$ is the iteration matrix in Jacobi. For the Jacobi method to converge, the spectral radius $\rho(G)$ must be less than 1.

**Applying Gershgorin's Theorem.** The diagonal elements of $G$ are zero, and the off-diagonal elements are $-\frac{a_{ij}}{a_{ii}}$. The sum of the absolute values of the off-diagonal elements in the $i$-th row is:
$$
\sum_{j \neq i} \left| -\frac{a_{ij}}{a_{ii}} \right| = \frac{1}{|a_{ii}|} \sum_{j \neq i} |a_{ij}|
$$
If $A$ is strictly diagonally dominant, then:
$$
|a_{ii}| > \sum_{j \neq i} |a_{ij}|
$$
This implies:
$$
\frac{1}{|a_{ii}|} \sum_{j \neq i} |a_{ij}| < 1
$$
Therefore, each Gershgorin disc for $G$ is centered at zero with a radius less than 1. Consequently, all eigenvalues of $G$ lie within the unit circle.

Since all eigenvalues of $G$ have magnitudes less than 1, the spectral radius $\rho(G) < 1$. This ensures that the Jacobi iteration converges for any initial guess $x^{(0)}$. 

The same proof applies for column diagonally dominant matrices.