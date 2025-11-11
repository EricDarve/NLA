# The Jacobi Iteration Method

The Jacobi iteration, proposed by Carl Gustav Jacob Jacobi in 1845, is the **simplest** of the stationary iterative methods for solving large sparse linear systems $A\mathbf{x} = \mathbf{b}$. It serves as an excellent foundational example of a splitting method.

## Defining the Splitting and Iteration

Recall that in a splitting method, we decompose the coefficient matrix $A$ into $A = M - N$, where $M$ is easily invertible. The iterative step is defined by $M \mathbf{x}^{(k+1)} = \mathbf{b} + N \mathbf{x}^{(k)}$.

For the Jacobi method, we define the splitting based on the diagonal, strictly lower triangular, and strictly upper triangular parts of $A$:

$$A = D - L - U$$

where $D = \text{diag}(A)$ is the diagonal matrix of $A$, $-L$ is the strictly lower triangular part, and $-U$ is the strictly upper triangular part.

The Jacobi method uses $M = D$ and $N = L + U$. Substituting these into the general splitting framework yields the Jacobi iteration:

$$D \mathbf{x}^{(k+1)} = \mathbf{b} + (L+U) \mathbf{x}^{(k)} \quad \text{}$$

Since $D$ is a diagonal matrix, $D^{-1}$ is trivial to compute, and the iteration can be explicitly written as a fixed-point iteration:

$$\mathbf{x}^{(k+1)} = D^{-1} (L+U) \mathbf{x}^{(k)} + D^{-1} \mathbf{b} \quad \text{}$$

We define the Jacobi iteration matrix (or tensor) as $G = R_J = D^{-1}(L+U)$ and the constant vector $\mathbf{c}_J = D^{-1}\mathbf{b}$. The iterative scheme then takes the standard stationary form:

$$\mathbf{x}^{(k+1)} = G \mathbf{x}^{(k)} + \mathbf{c}_J \quad \text{}$$

In component form, this updates each entry $x_i$ using the previous iteration's values for all other components:

$$x_i^{(k+1)} = \frac{1}{a_{ii}} \left( b_i - \sum_{j \neq i} a_{ij} x_j^{(k)} \right) \quad \text{}$$

Note that this method requires the diagonal entries $a_{ii}$ to be nonzero.

## Python Code Example

Here is a simple implementation of the Jacobi method in Python:

```python
import numpy as np

def jacobi(A, b, x0=None, tol=1e-10, max_iter=1000):
    """
    Solve the linear system Ax = b using the Jacobi iteration method.
    
    Parameters:
    -----------
    A : ndarray
        Coefficient matrix (n x n)
    b : ndarray
        Right-hand side vector (n,)
    x0 : ndarray, optional
        Initial guess (n,). If None, uses zero vector.
    tol : float, optional
        Tolerance for convergence (default: 1e-10)
    max_iter : int, optional
        Maximum number of iterations (default: 1000)
    
    Returns:
    --------
    x : ndarray
        Approximate solution vector
    iterations : int
        Number of iterations performed
    """
    n = len(b)
    x = np.zeros(n) if x0 is None else x0.copy()
    x_new = np.zeros(n)
    
    # Extract diagonal D and off-diagonal (L + U)
    D = np.diag(np.diag(A))
    LU = A - D
    D_inv = np.diag(1.0 / np.diag(A))
    
    for k in range(max_iter):
        # x^(k+1) = D^(-1) * (b + (L+U) * x^(k))
        x_new = D_inv @ (b - LU @ x)
        
        # Check convergence
        if np.linalg.norm(x_new - x, ord=np.inf) < tol:
            return x_new, k + 1
        
        x = x_new.copy()
    
    print(f"Warning: Maximum iterations ({max_iter}) reached without convergence")
    return x, max_iter
```

## Implementation and Parallelism

A key attraction of the Jacobi method is its **ease of implementation and high parallelism**.

Because the matrix $D$ is diagonal (or block diagonal), the computation of all components of $\mathbf{x}^{(k+1)}$ can be performed **independently and concurrently**. This characteristic has led to the Jacobi method also being referred to as **simultaneous iteration**.

For massively parallel machines, the Jacobi method is appealing because the operation $D \mathbf{x}^{(k+1)} = \mathbf{b} + (L+U) \mathbf{x}^{(k)}$ is highly parallelized. This contrasts with methods like Gauss-Seidel, where component updates are sequential due to the dependency on newly computed values within the current iteration.

## Convergence Criteria

As with any splitting method, the Jacobi iteration converges for any initial guess $\mathbf{x}^{(0)}$ if and only if the **spectral radius** of the iteration matrix $\rho(G)$ is less than one.

We rely on the matrix properties to guarantee this condition:

### Diagonal Dominance Condition

```{prf:theorem} Convergence of the Jacobi Method
:label: thm:jacobi_convergence

The Jacobi iteration is guaranteed to converge if the matrix $A$ is **strictly row diagonally dominant**:

$$|a_{ii}| > \sum_{j \neq i} |a_{ij}| \text{ for all } i \quad \text{}$$

A similar result holds if $A$ is **strictly column diagonally dominant**.
```

```{prf:proof}

If $A$ is strictly row diagonally dominant, we can prove convergence using the operator $\infty$-norm:
1. The iteration matrix is $G = D^{-1}(L+U)$.
2. The $\infty$-norm of $G$ is given by $\|G\|_{\infty} = \max_i \sum_{j \neq i} |a_{ij}| / |a_{ii}|$.
3. Strict row diagonal dominance implies $\|G\|_{\infty} < 1$.
4. Since the spectral radius $\rho(G)$ is bounded by any matrix norm, $\rho(G) \le \|G\|_{\infty}$, convergence is assured because $\rho(G) < 1$.

For column diagonally dominant matrices, convergence can be shown similarly using the matrix 1-norm, which confirms that $\|G^k\|_1 \to 0$.

**Convergence Proof using the Gershgorin Circle Theorem**

The Gershgorin Circle Theorem provides an alternative way to establish convergence under diagonal dominance:

1. Any eigenvalue $\lambda$ of an $n \times n$ matrix $G$ lies within at least one of the Gershgorin discs. The $i$-th disc is centered at $g_{ii}$ with radius $R_i = \sum_{j \neq i} |g_{ij}|$.
2. For the Jacobi iteration matrix $G = D^{-1}(L+U)$, the diagonal entries are $g_{ii} = 0$.
3. The radius of the $i$-th disc for $G$ is $R_i = \sum_{j \neq i} |a_{ij}/a_{ii}|$.
4. If $A$ is strictly diagonally dominant, $R_i < 1$ for all $i$.
5. Since the center of every disc is 0 and all discs have radii less than 1, all eigenvalues $\lambda$ of $G$ must satisfy $|\lambda| < 1$.
6. This ensures that the spectral radius $\rho(G) < 1$, and thus the Jacobi iteration converges for any initial guess $\mathbf{x}^{(0)}$.
```

## Practical Considerations

The Jacobi method provides a crucial entry point into iterative solvers, demonstrating the fundamental trade-off between simplicity/parallelism and convergence speed.