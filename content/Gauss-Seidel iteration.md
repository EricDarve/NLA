This is a better method compared to [[Jacobi iteration|Jacobi]] but at the cost of additional complexity.
$$
\begin{gather}
A = D - L - U \\
(D - L) x^{(k+1)} = b + U x^{(k)}
\end{gather}
$$
Again this method belongs to the class of [[Splitting methods|splitting methods.]]

### Algorithm

```julia
for k = 1:maxiter
    for i = 1:n
        sigma = 0.
        for j = 1:n
            if j != i
                sigma += A[i,j] * x[j]
            end
        end
        x[i] = (b[i] - sigma)/A[i,i]
    end
end
```

Gauss-Seidel uses the most recent estimate of `x[i]`, whereas [[Jacobi iteration|Jacobi]] uses the previous `x` to compute `xnew`.

## [[Convergence of classical iterative methods|Convergence]] of Gauss-Seidel

**Theorem 1.** If $A$ is strictly row diagonally dominant
$$
|a_{ii}| > \sum_{j \neq i} |a_{ij}|
$$
then both Jacobi and Gauss-Seidel always converge.

An analogous result holds if $A$ is strictly column diagonally dominant (i.e., $A^T$ is strictly row diagonally dominant).

**Theorem 2.** If $A$ is a [[Symmetric Positive Definite Matrices|symmetric positive definite matrix,]] then Gauss-Seidel always converges.

- Gauss-Seidel typically converges faster than Jacobi. 
- In many cases, the number of iterations with GS is about half that of Jacobi. See details in the textbook.
- However, GS is less parallel than Jacobi because of the additional data dependency.