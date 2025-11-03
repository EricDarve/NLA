# The Lanczos Algorithm: Arnoldi for Symmetric Matrices

We now turn our attention to the important special case where the matrix $A \in \mathbb{R}^{n \times n}$ is symmetric. The Arnoldi iteration, while general, does not exploit this structure. When symmetry is present, the algorithm simplifies dramatically, leading to profound computational savings. This specialized algorithm is the **Lanczos algorithm**, one of the most influential algorithms in numerical linear algebra.

## From Hessenberg to Tridiagonal

Recall the fundamental relation from the Arnoldi iteration after $k$ steps:

$$
A Q_k = Q_k H_k + h_{k+1,k} \boldsymbol{q}_{k+1} \boldsymbol{e}_k^T
$$

where $Q_k = [\boldsymbol{q}_1 | \cdots | \boldsymbol{q}_k]$ has orthonormal columns, and $H_k$ is a $k \times k$ upper Hessenberg matrix. By pre-multiplying by $Q_k^T$ and using the orthogonality of the columns of $Q_{k+1}$, we arrived at the projection:

$$
H_k = Q_k^T A Q_k
$$

Now, let us impose the condition that $A$ is symmetric ($A = A^T$). The projected matrix $H_k$ must inherit this symmetry:

$$
H_k^T = (Q_k^T A Q_k)^T = Q_k^T A^T (Q_k^T)^T = Q_k^T A Q_k = H_k
$$

So, $H_k$ must be a symmetric matrix.

It is a remarkable fact that a matrix cannot be both upper Hessenberg and symmetric without having a very specific structure: it must be **tridiagonal**.

An upper Hessenberg matrix has zeros below its first subdiagonal ($h_{i,j} = 0$ for $i > j+1$). A symmetric matrix requires that $h_{i,j} = h_{j,i}$ for all $i,j$. If we combine these properties, any non-zero entry in the upper triangle, say $h_{j, i}$ with $i > j+1$, would imply a corresponding non-zero entry $h_{i,j}$ in the lower triangle, which would violate the Hessenberg structure. Thus, the only non-zero entries can be on the main diagonal, the first superdiagonal, and the first subdiagonal.

We therefore rename the matrix $T_k = H_k$ to emphasize its tridiagonal nature. It is standard to denote the diagonal entries by $\alpha_j$ and the off-diagonal entries by $\beta_j$. Since $T_k$ is symmetric, we have $T_{j, j-1} = T_{j-1, j}$.

In the Lanczos recurrence, we will set the diagonal entries $\alpha_k = h_{k,k}$ and the new subdiagonal entry $\beta_k = h_{k+1,k}$. The matrix $T_k$ is then:

$$
T_k = Q_k^T A Q_k =
\begin{pmatrix}
\alpha_1 & \beta_1 & & & \\
\beta_1 & \alpha_2 & \beta_2 & & \\
& \beta_2 & \ddots & \ddots & \\
& & \ddots & \alpha_{k-1} & \beta_{k-1} \\
& & & \beta_{k-1} & \alpha_k
\end{pmatrix}
$$

## The Three-Term Recurrence

This tridiagonal structure leads to a dramatic simplification of the Arnoldi recurrence. The core step in Arnoldi is the modified Gram-Schmidt process used to compute $\boldsymbol{q}_{k+1}$:

$$
h_{k+1, k} \boldsymbol{q}_{k+1} = A \boldsymbol{q}_k - \sum_{j=1}^{k} h_{j,k} \boldsymbol{q}_j
$$

Since $H_k$ (now $T_k$) is tridiagonal, for a given column $k$, the only non-zero entries $h_{j,k}$ are $h_{k-1, k}$ and $h_{k,k}$. Using our new notation, these are $\beta_{k-1}$ and $\alpha_k$, respectively. All other $h_{j,k}$ for $j < k-1$ are zero.

The summation therefore collapses from $k$ terms to just two, yielding the famous **Lanczos three-term recurrence**:

$$
\beta_k \boldsymbol{q}_{k+1} = A \boldsymbol{q}_k - \alpha_k \boldsymbol{q}_k - \beta_{k-1} \boldsymbol{q}_{k-1}
$$

### Another Derivation of the Three-Term Recurrence

Another way to derive the three-term recurrence is to start from the Arnoldi relation and use symmetry directly. We start with:

$$
A q_k = \sum_{j=1}^{k+1} h_{j,k} q_j
$$

Take a dot product with $q_i$ for $i < k-1$:

$$
q_i^T A q_k = \sum_{j=1}^{k+1} h_{j,k} q_i^T q_j = h_{i,k}
$$

Using symmetry of $A$:

$$
q_i^T A q_k = (A q_i)^T q_k
$$

However $A q_i \in \text{span}\{q_1, \dots, q_{i+1}\}$, so $q_k$ is orthogonal to $A q_i$ for $i < k-1$. Thus $q_i^T A q_k = 0$ for $i < k-1$, which implies $h_{i,k} = 0$ for $i < k-1$. This leads us again to the three-term recurrence.

## Lanczos Algorithm

Let's assemble this into the algorithm. The Python implementation of the Lanczos iteration is as follows:

```python
import numpy as np

def lanczos(A, b, m):
    """
    Lanczos iteration for symmetric matrix A.
    
    Parameters:
    A: n x n symmetric matrix
    b: starting vector (will be normalized)
    m: number of iterations
    
    Returns:
    T: m x m tridiagonal matrix
    Q: n x m matrix with orthonormal columns (optional, for Ritz vectors)
    """
    n = A.shape[0]
    
    # Initialize
    q = b / np.linalg.norm(b)
    q_prev = np.zeros(n)
    beta_prev = 0.0
    
    alpha = np.zeros(m)
    beta = np.zeros(m)
    Q = np.zeros((n, m))
    
    for k in range(m):
        Q[:, k] = q
        
        # Matrix-vector product
        v = A @ q
        
        # Compute diagonal entry
        alpha[k] = q.T @ v
        
        # Three-term recurrence
        r = v - alpha[k] * q - beta_prev * q_prev
        
        # Compute off-diagonal entry
        beta_k = np.linalg.norm(r)
        
        # Check for invariant subspace
        if beta_k < 1e-12:
            alpha = alpha[:k+1]
            beta = beta[:k]
            Q = Q[:, :k+1]
            break
        
        # Update for next iteration
        q_prev = q
        q = r / beta_k
        beta_prev = beta_k
        if k < m - 1:
            beta[k] = beta_k
    
    # Construct tridiagonal matrix T
    T = np.diag(alpha) + np.diag(beta[:-1] if len(beta) > 1 else [], 1) + \
        np.diag(beta[:-1] if len(beta) > 1 else [], -1)
    
    return T, Q
```

The eigenvalues of the tridiagonal matrix $T_m$ are the **Ritz values**, which approximate the eigenvalues of $A$. The corresponding eigenvectors of $T_m$ can be used to compute the associated **Ritz vectors**, which approximate the eigenvectors of $A$.

## Lanczos vs. Arnoldi: A Comparison

The transition from Arnoldi to Lanczos is not merely an aesthetic simplification; it represents one of the most significant performance gains in numerical methods. The cost of running for $k$ iterations is dominated by two components: matrix-vector products and the vector orthogonalization steps. For a sparse matrix $A$ where the number of non-zero entries, $\text{nnz}(A)$, is on the order of $n$, the costs are summarized below.

| Cost Component | Arnoldi Iteration (General $A$) | Lanczos Iteration (Symmetric $A$) |
| :--- | :--- | :--- |
| **Matrix-Vector Products** <br/> (Total for $k$ iterations) | $O(k \cdot \text{nnz}(A))$ | $O(k \cdot \text{nnz}(A))$ |
| **Vector Operations** <br/> (Total for $k$ iterations) | $O(k^2 n)$ <br/> *Orthogonalization against all $j-1$ previous vectors at iteration $j$* | $O(kn)$ <br/> *Operations involve only the 2 previous vectors at each iteration* |
| **Storage Requirement** | $O(kn)$ <br/> *Must store the entire basis $Q_k$* | $O(n)$ <br/> *Only needs to store the two most recent vectors* |
| **Total Asymptotic Cost** <br/> (for sparse $A$) | **$O(k^2 n)$** | **$O(kn)$** |

As the table shows, the cost of Lanczos scales linearly with the number of iterations, whereas Arnoldi scales quadratically. This makes the **Lanczos algorithm dramatically more efficient**, allowing for many more iterations to be run to achieve convergence, which is essential when dealing with the massive matrices found in modern scientific computing.

## A Practical Caveat: Loss of Orthogonality

In theory, the Lanczos vectors $\boldsymbol{q}_j$ are perfectly orthogonal. However, in finite-precision arithmetic, rounding errors accumulate, and the orthogonality among the vectors degrades over time. Specifically, $\boldsymbol{q}_k$ tends to lose orthogonality with the first few vectors $\boldsymbol{q}_j$ whose corresponding Ritz vector components are converging.

This loss of orthogonality is a famous and subtle feature of the algorithm. It manifests as the appearance of multiple copies ("ghosts") of already-converged eigenvalues in the spectrum of $T_m$. While this can be problematic, it can also be managed. Techniques like *full reorthogonalization* (which makes the algorithm as expensive as Arnoldi) or more practical *selective reorthogonalization* can be employed to maintain orthogonality where it matters most. For many applications, particularly finding the largest or smallest eigenvalues, the unmodified Lanczos algorithm is remarkably effective despite this issue.