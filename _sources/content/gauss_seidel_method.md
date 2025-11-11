# The Gauss-Seidel Iteration Method

The Gauss-Seidel (GS) method, also known as the Seidel process or successive iteration method, is a stationary iterative algorithm designed to solve $A\mathbf{x} = \mathbf{b}$. Although commonly attributed to Gauss and Seidel, historical accounts note that Gauss nowhere mentioned it, and Seidel (1874) mentioned the process but recommended against its use.

The Gauss-Seidel method improves upon the Jacobi iteration by using updated solution components as soon as they become available within the current iteration. This characteristic makes it a **successive iteration** method, in contrast to the Jacobi method, which performs updates simultaneously.

## Splitting and Iteration

Like the Jacobi method, Gauss-Seidel is a splitting method, where the coefficient matrix $A$ is decomposed as $A = M - N$. We use the standard decomposition of $A$ into its diagonal ($D$), strictly lower triangular ($L$), and strictly upper triangular ($U$) parts, such that $A = D - L - U$.

For the (forward) Gauss-Seidel method, the splitting matrices are chosen as:

$$M = D - L \quad \text{and} \quad N = U$$

Substituting this splitting into the general stationary iteration form $M \mathbf{x}^{(k+1)} = \mathbf{b} + N \mathbf{x}^{(k)}$ yields the defining equation for the Gauss-Seidel iteration:

$$(D - L) \mathbf{x}^{(k+1)} = \mathbf{b} + U \mathbf{x}^{(k)}$$

Since $D-L$ is a lower triangular matrix, computing $\mathbf{x}^{(k+1)}$ involves solving a triangular system. The Gauss-Seidel iteration matrix is $G = R_{GS} = (D-L)^{-1}U$.

In component form, the update for the $i$-th component $x_i$ uses the newly computed components $x_j^{(k+1)}$ for $j < i$ and the previously calculated components $x_j^{(k)}$ for $j > i$:

$$x_i^{(k+1)} = \frac{1}{a_{ii}} \left( b_i - \sum_{j=1}^{i-1} a_{ij} x_j^{(k+1)} - \sum_{j=i+1}^{n} a_{ij} x_j^{(k)} \right)$$

## Algorithm and Implementation

A key difference between Gauss-Seidel and Jacobi iteration lies in the dependencies: because the components of $\mathbf{x}^{(k+1)}$ must be computed sequentially (or successively), the Gauss-Seidel method is **less parallel** than the Jacobi method. Conversely, Gauss-Seidel typically converges **faster than Jacobi, often requiring about half the number of iterations in many cases.**

## Python Code Example

Here is a simple implementation of the Gauss-Seidel method in Python:

```python
import numpy as np

def gauss_seidel(A, b, x0=None, tol=1e-10, max_iter=1000):
    """
    Solve the linear system Ax = b using the Gauss-Seidel iteration method.
    
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
    
    for k in range(max_iter):
        x_old = x.copy()
        
        # Update each component using the latest available values
        for i in range(n):
            sigma = np.dot(A[i, :i], x[:i]) + np.dot(A[i, i+1:], x[i+1:])
            x[i] = (b[i] - sigma) / A[i, i]
        
        # Check convergence
        if np.linalg.norm(x - x_old, ord=np.inf) < tol:
            return x, k + 1
    
    print(f"Warning: Maximum iterations ({max_iter}) reached without convergence")
    return x, max_iter
```

## Conditions for Convergence

The convergence of the Gauss-Seidel iteration is determined by the spectral radius of the iteration matrix $G = (D-L)^{-1}U$. The scheme converges for any initial guess $\mathbf{x}^{(0)}$ if and only if $\rho(G) < 1$.

We present sufficient conditions that guarantee convergence:

```{prf:theorem} Convergence of Gauss-Seidel Method
:label: thm:gauss_seidel_convergence

The associated Gauss-Seidel iteration converges for any initial guess $\mathbf{x}^{(0)}$ if $A$ satisfies one of the following conditions:

1.  **Strict Diagonal Dominance:** $A$ is **strictly row diagonally dominant** (or strictly column diagonally dominant), i.e.,

    $$|a_{ii}| > \sum_{j \neq i} |a_{ij}| \quad \text{for all } i$$

2.  **Symmetric Positive Definite (SPD):** $A$ is **Symmetric Positive Definite**.
```

```{prf:proof}

We prove the convergence under the condition of strict row diagonal dominance (SRDD).

Let $G = (D-L)^{-1}U$ be the Gauss-Seidel iteration matrix. By the splitting method convergence theorem, we need to show that if $A$ is SRDD, then $\rho(G) < 1$. This is equivalent to showing that for any eigenvalue $\lambda$ of $G$, we have $|\lambda| < 1$.

1.  **Assume Contradiction:** Assume there exists an eigenvalue $\lambda$ such that $|\lambda| \ge 1$. Let $\mathbf{x} \neq \mathbf{0}$ be the corresponding eigenvector:

    $$G \mathbf{x} = \lambda \mathbf{x}$$

2.  **Derive Singular Matrix:** Substitute the splitting definition $G = (D-L)^{-1}U$:

    $$\begin{aligned}
    (D-L)^{-1}U \mathbf{x} &= \lambda \mathbf{x} \\
    U \mathbf{x} &= \lambda (D - L) \mathbf{x}
    \end{aligned}$$

    Rearranging this gives:

    $$(\lambda D - \lambda L - U) \mathbf{x} = \mathbf{0}$$

    Since $\mathbf{x} \neq \mathbf{0}$, the matrix $A(\lambda) = \lambda D - \lambda L - U$ must be **singular**.

3.  **Establish Diagonal Dominance of $A(\lambda)$:** We now show that if $A$ is SRDD and $|\lambda| \ge 1$, then $A(\lambda)$ must also be SRDD, implying it is nonsingular (which leads to a contradiction).
    *   The entries of $A(\lambda)$ are:
        *   Diagonal: $A(\lambda)_{ii} = \lambda a_{ii}$
        *   Lower triangular (for $j < i$): $A(\lambda)_{ij} = -\lambda a_{ij}$
        *   Upper triangular (for $j > i$): $A(\lambda)_{ij} = -a_{ij}$
    *   Consider the $i$-th row of $A(\lambda)$. The magnitude of the diagonal element is $|A(\lambda)_{ii}| = |\lambda| |a_{ii}|$.
    *   The sum of the magnitudes of the off-diagonal elements in the $i$-th row is:

        $$\begin{aligned}
        \sum_{j \neq i} |A(\lambda)_{ij}| &= \sum_{j < i} |-\lambda a_{ij}| + \sum_{j > i} |-a_{ij}| \\
        &= |\lambda| \sum_{j < i} |a_{ij}| + \sum_{j > i} |a_{ij}|
        \end{aligned}$$

    *   Since $A$ is SRDD, we know that $|a_{ii}| > \sum_{j \neq i} |a_{ij}| = \sum_{j < i} |a_{ij}| + \sum_{j > i} |a_{ij}|$.
    *   Since $|\lambda| \ge 1$, and $|a_{ij}| \ge 0$:

        $$|\lambda| \sum_{j < i} |a_{ij}| \ge \sum_{j < i} |a_{ij}|$$

    *   Therefore:

        $$\sum_{j \neq i} |A(\lambda)_{ij}| \le |\lambda| \sum_{j < i} |a_{ij}| + |\lambda| \sum_{j > i} |a_{ij}| = |\lambda| \left( \sum_{j \neq i} |a_{ij}| \right)$$

    *   Since $A$ is strictly dominant, $\sum_{j \neq i} |a_{ij}| < |a_{ii}|$.
    *   Thus, $\sum_{j \neq i} |A(\lambda)_{ij}| < |\lambda| |a_{ii}| = |A(\lambda)_{ii}|$.

**Conclusion:** $A(\lambda)$ is strictly row diagonally dominant. A strictly diagonally dominant matrix is nonsingular. This contradicts the necessity that $A(\lambda)$ must be singular if $|\lambda| \ge 1$. Therefore, our initial assumption must be false, and we conclude that **$|\lambda| < 1$**. The Gauss-Seidel iteration converges.
```

## Conceptual Analogy

The Gauss-Seidel method is like having a team of specialized construction workers updating a skyscraper sequentially: the foundation workers (lower indices) immediately use freshly poured concrete from their upstream colleagues to speed up their work, which generally results in the building (the solution) reaching stability much faster than if every team waited for every other team to finish their previous task (as in Jacobi). However, the critical path (the dependency on the previous worker's job) makes it harder to use 1,000 workers simultaneously compared to the independent setup of Jacobi.