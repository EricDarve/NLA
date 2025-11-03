# The Arnoldi Process

Previously, we introduced the core idea of the Arnoldi iteration: a numerically stable process to build an orthonormal basis $Q_k = [q_1, \dots, q_k]$ for the Krylov subspace $\mathcal{K}_k(A, q_1)$. This process simultaneously generates a $k \times k$ upper Hessenberg matrix $H_k = Q_k^T A Q_k$.

Now, we will formalize this algorithm and show precisely why the eigenvalues of $H_k$ are good approximations for the eigenvalues of $A$.

## The Arnoldi Algorithm

The algorithm is an iterative application of the Modified Gram-Schmidt process. Starting with a chosen unit vector $q_1$, we generate $q_2, q_3, \dots, q_k$ one by one.

**Algorithm: The Arnoldi Iteration**

* **Input:** A matrix $A$, a number of steps $k$, and a starting vector $q_1$ with $\|q_1\|_2 = 1$.
* **Output:** $Q_k = [q_1, \dots, q_k]$ with orthonormal columns and a $k \times k$ upper Hessenberg matrix $H_k$ with entries $h_{ij}$.

```python
def arnoldi_iteration(A, q1, k):
    """
    The Arnoldi Iteration
    
    Parameters:
    -----------
    A : matrix or linear operator
        The matrix to approximate eigenvalues of
    q1 : array
        Starting vector with ||q1||_2 = 1
    k : int
        Number of iterations
    
    Returns:
    --------
    Q_k : matrix
        Orthonormal basis vectors [q_1, ..., q_k]
    H_k : matrix
        k x k upper Hessenberg matrix
    """
    n = len(q1)
    Q = [q1]  # List to store q_j vectors
    H = np.zeros((k, k))  # Upper Hessenberg matrix
    
    for j in range(k):
        # Compute matrix-vector product
        v = A @ Q[j]
        
        # Modified Gram-Schmidt orthogonalization
        for i in range(j + 1):
            H[i, j] = Q[i].T @ v  # Compute projection coefficient
            v = v - H[i, j] * Q[i]  # Subtract projection
        
        # Compute next basis vector
        if j < k - 1:
            h_next = np.linalg.norm(v)
            
            if h_next == 0:
                # Lucky breakdown: invariant subspace found
                print("Lucky breakdown at iteration", j + 1)
                return np.column_stack(Q), H[:j+1, :j+1]
            
            H[j + 1, j] = h_next
            Q.append(v / h_next)
    
    return np.column_stack(Q), H
```

This algorithm produces the $k \times k$ matrix $H_k$ and the basis $Q_k$ that satisfy the **Arnoldi relation** we saw previously:

$$
A Q_k = Q_k H_k + h_{k+1,k} q_{k+1} e_k^T
$$

The computational cost is dominated by the $k$ sparse matrix-vector products and the $O(k^2 n)$ work of the inner loop. As long as $k \ll n$, this is far more efficient than any $O(n^3)$ direct method.

## Ritz Values and Ritz Vectors

Our strategy is to use the eigenvalues of the small, dense matrix $H_k$ to approximate the eigenvalues of $A$. Let's formalize this.

Let $H_k$ have an eigendecomposition $H_k = Y_k \Lambda_k Y_k^{-1}$, where $\Lambda_k = \text{diag}(\theta_1, \dots, \theta_k)$ is the diagonal matrix of eigenvalues and $Y_k = [y_1, \dots, y_k]$ contains the corresponding eigenvectors.

* **Ritz Values:** The eigenvalues $\theta_j$ of $H_k$ are called the **Ritz values**. They serve as our approximations to the eigenvalues of $A$.
* **Ritz Vectors:** The eigenvectors $y_j$ of $H_k$ are $k$-dimensional. We can "lift" them back into the $n$-dimensional space to get our approximate eigenvectors of $A$, which are called the **Ritz vectors**:

$$
u_j = Q_k y_j
$$

Note that $u_j$ is a linear combination of our basis vectors $q_1, \dots, q_k$, so $u_j \in \mathcal{K}_k$.

## Why is this a good approximation?

Let's see just how "approximate" these eigenpairs are. We can derive an exact expression for the residual $A u_j - \theta_j u_j$.

We start with the Arnoldi relation:

$$
A Q_k = Q_k H_k + h_{k+1,k} q_{k+1} e_k^T
$$

Now, multiply from the right by a single eigenvector $y_j$ of $H_k$:

$$
A Q_k y_j = Q_k H_k y_j + h_{k+1,k} q_{k+1} (e_k^T y_j)
$$

Substitute our definitions. On the left, $Q_k y_j = u_j$. On the right, $H_k y_j = \theta_j y_j$.

$$
\begin{aligned}
A u_j &= Q_k (\theta_j y_j) + h_{k+1,k} q_{k+1} (e_k^T y_j) \\
&= \theta_j (Q_k y_j) + h_{k+1,k} q_{k+1} (e_k^T y_j) \\
&= \theta_j u_j + h_{k+1,k} q_{k+1} (y_j)_k
\end{aligned}
$$

where $(y_j)_k = e_k^T y_j$ is simply the $k$-th (i.e., last) component of the small eigenvector $y_j$.

Rearranging this gives us the **residual vector**:

$$
A u_j - \theta_j u_j = h_{k+1,k} (y_j)_k q_{k+1}
$$

This is a remarkable and powerful result. It gives us a cheap and easy way to check the quality of our approximation. By taking the 2-norm (and since $\|q_{k+1}\|_2 = 1$):

$$
\| A u_j - \theta_j u_j \|_2 = |h_{k+1,k}| \cdot |(y_j)_k|
$$

This tells us that the Ritz pair $(\theta_j, u_j)$ is a good approximation for an eigenpair of $A$ if:

1.  The quantity $h_{k+1,k}$ (the "breakdown" term) is small, OR
2.  The last component $(y_j)_k$ of the small eigenvector $y_j$ is small.

In practice, the Arnoldi iteration is run, and the eigenvalues of $H_k$ are computed. We then compute the residual norm for each Ritz pair using this formula. If the residual is below a desired tolerance, we declare that Ritz pair "converged."

```{note}
**Lucky breakdown:** If at some step $j$, we find $h_{j+1, j} = 0$, what does this mean? The Arnoldi relation becomes $A Q_j = Q_j H_j$. This implies that the subspace $\mathcal{K}_j = \text{span}(q_1, \dots, q_j)$ is an **invariant subspace** of $A$. In this case, the $j$ Ritz values (the eigenvalues of $H_j$) are *exactly* eigenvalues of $A$. This is rare in practice but is a welcome event when it occurs.
```