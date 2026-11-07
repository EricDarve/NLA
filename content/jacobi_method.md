# The Jacobi Iteration Method

Jacobi iteration updates each component of an approximate solution using the previous values of all the other components. This makes it a simple example of a stationary splitting method for solving $Ax=b$.

One update of the full vector is called an **iteration** or a **sweep**. As before, $x^{(k)}$ denotes the approximation after $k$ such updates. We assume $A$ is real and nonsingular, with nonzero diagonal entries.

## Defining the Splitting and Iteration

Recall that in a splitting method, we decompose the coefficient matrix $A$ into $A = M - N$, where linear systems with $M$ are inexpensive to solve. The iterative step is defined by $M x^{(k+1)} = b + N x^{(k)}$.

For the Jacobi method, we define the splitting based on the diagonal, strictly lower triangular, and strictly upper triangular parts of $A$:

$$A = D - L - U$$

where $D = \text{diag}(A)$ is the diagonal matrix of $A$, $-L$ is the strictly lower triangular part, and $-U$ is the strictly upper triangular part.

The Jacobi method uses $M = D$ and $N = L + U$. Substituting these into the general splitting framework yields the Jacobi iteration:

$$D x^{(k+1)} = b + (L+U) x^{(k)} $$

Since $D$ is diagonal, solving with $D$ requires only componentwise division. The iteration can be written as

$$x^{(k+1)} = D^{-1} (L+U) x^{(k)} + D^{-1} b $$

Define the Jacobi iteration matrix $R_J=D^{-1}(L+U)$ and the constant vector $c_J=D^{-1}b$. The stationary form is

$$x^{(k+1)} = R_J x^{(k)} + c_J $$

Expanding row $i$ and dividing by $a_{ii}$ gives the component update. Every term on the right uses the previous approximation:

$$x_i^{(k+1)} = \frac{1}{a_{ii}} \left( b_i - \sum_{j \neq i} a_{ij} x_j^{(k)} \right) $$

## Python Code Example

The code uses the residual correction form of the splitting iteration. With $M=D$, each component update is

$$
x_i^{(k+1)}=x_i^{(k)}+\frac{r_i^{(k)}}{a_{ii}},
\qquad r^{(k)}=b-Ax^{(k)}.
$$

Dividing the residual componentwise by the diagonal gives the same Jacobi update without forming $D$ or $D^{-1}$. The implementation uses real NumPy arrays.

We cannot check the error $x^{(k)}-x$ directly because the exact solution is unknown. Instead, we check the residual relative to the scale of the system, using the [backward-error measure](lu_pivoting.md) introduced for linear solves:

$$
\eta(x^{(k)})=
\frac{\|b-Ax^{(k)}\|_\infty}
{\|A\|_\infty\|x^{(k)}\|_\infty+\|b\|_\infty}.
$$

A small change between iterates alone does not guarantee a small residual. A small $\eta$ measures how accurately the equations are satisfied; the solution error also depends on the conditioning of $A$. In the code, `x` is the current approximation and `tol` is the requested upper bound on $\eta$.

```python
import numpy as np

def jacobi(A, b, x0=None, tol=1e-10, max_iter=1000):
    """Return (x, updates), stopping when the scaled residual is <= tol.

    A and b are real arrays. Raise RuntimeError if the iteration limit
    is reached without satisfying the stopping criterion.
    """
    A = np.asarray(A, dtype=float)
    b = np.asarray(b, dtype=float)
    if b.ndim != 1 or A.shape != (b.size, b.size):
        raise ValueError("A must be square and b a matching vector.")
    if not np.isfinite(tol) or tol <= 0:
        raise ValueError("tol must be positive and finite.")
    if not isinstance(max_iter, (int, np.integer)) or max_iter < 0:
        raise ValueError("max_iter must be a nonnegative integer.")

    diagonal = np.diag(A)
    if np.any(diagonal == 0):
        raise ValueError("Jacobi requires nonzero diagonal entries.")

    x = np.zeros_like(b) if x0 is None else np.array(x0, dtype=float, copy=True)
    if x.shape != b.shape:
        raise ValueError("x0 must have the same shape as b.")

    norm_A = np.linalg.norm(A, ord=np.inf)
    norm_b = np.linalg.norm(b, ord=np.inf)
    for k in range(max_iter + 1):
        r = b - A @ x
        scale = norm_A * np.linalg.norm(x, ord=np.inf) + norm_b
        if np.linalg.norm(r, ord=np.inf) <= tol * scale:
            return x, k
        if k < max_iter:
            # Every correction uses the same old vector x.
            x = x + r / diagonal

    raise RuntimeError("Jacobi did not converge within max_iter updates.")
```

The returned count is the number of updates performed; it is zero if the initial guess already satisfies the tolerance. If $b=0$ and $x=0$, the residual is exactly zero and the code accepts it without dividing by a zero scale.

A sparse implementation can use the same correction formula with a sparse matrix-vector product. The example above uses dense arrays for simplicity.

## Implementation and Parallelism

All entries of $x^{(k+1)}$ depend only on $x^{(k)}$. They can therefore be computed independently once the old vector is available. Keep the old values until every new component has been computed; using a newly updated component would change the method.

This is why Jacobi is also called **simultaneous iteration**. It is well suited to parallel computation. Gauss-Seidel, developed next, uses new values immediately and introduces dependencies between component updates.

## Convergence Criteria

As with any splitting method, the Jacobi iteration converges for every initial guess if and only if $\rho(R_J)<1$.

We can often check this condition without computing the eigenvalues. The key bound is $\rho(R_J)\le\|R_J\|_\infty$, where the infinity norm is the largest absolute row sum.

To see why, recall the [operator norm inequality](operator_norms.md): the operator norm bounds how much the matrix can stretch any vector. Apply it to an eigenvector $v\ne0$ with $R_Jv=\lambda v$:

$$
|\lambda|\,\|v\|_\infty
=\|\lambda v\|_\infty
=\|R_Jv\|_\infty
\le\|R_J\|_\infty\,\|v\|_\infty.
$$

Since $\|v\|_\infty>0$, division gives $|\lambda|\le\|R_J\|_\infty$. This holds for every eigenvalue, so their largest magnitude, $\rho(R_J)$, satisfies the same bound. Thus, if every absolute row sum is less than one, all eigenvalues lie inside the unit circle and Jacobi converges.

The same argument works for any induced norm, including the 1-norm used below, which is the largest absolute column sum.

### Diagonal Dominance Condition

```{prf:theorem} Convergence of the Jacobi Method
:label: thm:jacobi_convergence

The Jacobi iteration is guaranteed to converge if the matrix $A$ is **strictly row diagonally dominant**:

$$|a_{ii}| > \sum_{j \neq i} |a_{ij}| \text{ for all } i $$

The same conclusion holds if $A$ is **strictly column diagonally dominant**, meaning $|a_{jj}|>\sum_{i\ne j}|a_{ij}|$ for every column $j$.
```

:::{prf:proof}
For strict row diagonal dominance, the iteration matrix satisfies

$$
\rho(R_J)\le\|R_J\|_\infty
=\max_i\frac{\sum_{j\ne i}|a_{ij}|}{|a_{ii}|}<1.
$$

Multiplication by $D^{-1}$ on the left scales rows, so strict column dominance does not directly bound $\|R_J\|_1$. Instead, use the similar matrix

$$
B=DR_JD^{-1}=(L+U)D^{-1}.
$$

Its column sums give

$$
\rho(R_J)=\rho(B)\le\|B\|_1
=\max_j\frac{\sum_{i\ne j}|a_{ij}|}{|a_{jj}|}<1.
$$

The splitting convergence theorem therefore applies in both cases.
:::

### An Alternative Proof Using Gershgorin Discs

Gershgorin's theorem describes eigenvalue locations using row sums. For a matrix $B=(b_{ij})$, every eigenvalue $\lambda$ lies in at least one disc

$$
|\lambda-b_{ii}|\le\sum_{j\ne i}|b_{ij}|.
$$

To see why, take an eigenvector $v$ and choose $i$ so that $|v_i|$ is its largest component magnitude. Row $i$ of $Bv=\lambda v$ gives

$$
|\lambda-b_{ii}|\,|v_i|
=\left|\sum_{j\ne i}b_{ij}v_j\right|
\le\sum_{j\ne i}|b_{ij}|\,|v_j|
\le\left(\sum_{j\ne i}|b_{ij}|\right)|v_i|.
$$

Dividing by $|v_i|>0$ proves the claim. Each disc is centered at a diagonal entry, with radius equal to the absolute sum of the other entries in that row.

For $R_J$, the diagonal entries are zero and the radii are $\sum_{j\ne i}|a_{ij}|/|a_{ii}|$. Under strict row diagonal dominance, these radii are less than one. Every disc therefore lies inside the open unit circle, giving $\rho(R_J)<1$.

Strict diagonal dominance is sufficient, but is not necessary: other matrices also give convergent Jacobi iterations.

Symmetric positive definiteness alone does not guarantee Jacobi convergence. For example, a $3\times3$ matrix with diagonal entries $1$ and all off-diagonal entries $3/4$ has positive eigenvalues $5/2,1/4,1/4$. But $D=I$, so $R_J=I-A$ has eigenvalue $-3/2$ and the iteration can diverge. Gauss-Seidel, developed next, does converge for every SPD matrix.

## Practical Considerations

Each sparse Jacobi update requires a matrix-vector product and componentwise division, so its cost is proportional to the number of nonzero entries in $A$. Low cost and parallelism are its main attractions, but convergence can be slow. The next section develops [Gauss-Seidel iteration](gauss_seidel_method.md), which uses new component values as soon as they are available.
