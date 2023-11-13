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