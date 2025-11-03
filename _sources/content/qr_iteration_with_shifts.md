# The Shifted QR Iteration Algorithm

This approach combines the **shift** (to dramatically speed up convergence) with **deflation** (to reduce the problem size) to efficiently find all eigenvalues of a matrix.

The core of the process is an iterative loop. Instead of decomposing the matrix $A$ (or $T_k$) directly, we "shift" it by an approximate eigenvalue $\mu$ *before* the QR factorization.

Each step $k$ of the iteration works as follows:

1.  **Start with your matrix:** Let $T_k$ be the matrix at the current step (initially $T_0 = A$, which is assumed to be in **upper Hessenberg** form).

2.  **Choose a Shift:** Select a shift $\mu$. A very effective choice is the bottom-right element of the current matrix: $\mu = T_k[n,n]$.

3.  **Factor the Shifted Matrix:** Compute the QR factorization of the *shifted* matrix:

    $$T_k - \mu I = U_k R_k$$

4.  **Recombine and "Unshift":** Create the next matrix, $T_{k+1}$, by multiplying $R_k$ and $U_k$ in the reverse order and then adding the shift $\mu$ back:

    $$T_{k+1} = R_k U_k + \mu I$$

The process is a **similarity transformation** ($T_{k+1} = U_k^H T_k U_k$), which means all eigenvalues are preserved at every step. The algorithm is simply converging to the **Schur form**, an upper triangular matrix with the eigenvalues on the diagonal.

## Shifted QR Iteration Step-by-Step

This iteration is implemented in the following Python code:

```python
import numpy as np

# Assume A is the initial upper Hessenberg matrix
Tk = A.copy()
n = A.shape[0]

# This loop would contain logic to check for convergence and perform deflation
while not_converged:
    # Choose the shift (e.g., the bottom-right element)
    mu = Tk[n-1, n-1]
    
    # Create an identity matrix of the same size
    I = np.eye(n)
    
    # Perform QR factorization on the shifted matrix
    Qk, Rk = np.linalg.qr(Tk - mu * I)
    
    # Recombine and add the shift back to get the next iterate
    Tk = Rk @ Qk + mu * I
```

## Similarity Transformation Verification

We can verify that this algorithm preserves eigenvalues because the transformation from $T_k$ to $T_{k+1}$ is a **unitary similarity transformation**. Let's walk through the derivation in detail.

Our starting point is the two core equations of a single iteration:

1.  **The QR factorization of the shifted matrix:**

    $$T_k - \mu I = U_{k} R_{k}$$

2.  **The update step to find the next matrix in the sequence:**

    $$T_{k+1} = R_{k} U_{k} + \mu I$$

Our goal is to show that these two steps are equivalent to the similarity transformation $T_{k+1} = U_{k}^H T_k U_{k}$.

**Step 1: Isolate $R_{k}$**

We begin with the first equation. Since $U_{k}$ is a unitary matrix from the QR decomposition, its inverse is its conjugate transpose, $U_{k}^H$. We can left-multiply both sides by $U_{k}^H$ to isolate $R_{k}$:

$$
U_{k}^H (T_k - \mu I) = U_{k}^H (U_{k} R_{k})
$$

Since $U_{k}^H U_{k} = I$, the equation simplifies to:

$$
R_{k} = U_{k}^H (T_{k} - \mu I)
$$

**Step 2: Substitute $R_{k}$ into the Update Equation**

Now, we take this expression for $R_{k}$ and substitute it into the second equation, $T_{k+1} = R_{k} U_{k} + \mu I$:

$$
T_{k+1} = \left[ U_{k}^H (T_k - \mu I) \right] U_{k} + \mu I
$$

**Step 3: Expand and Simplify**

Let's distribute the terms inside the parentheses. The term $U_{k}$ on the right multiplies both $T_k$ and $\mu I$:

$$
T_{k+1} = (U_{k}^H T_k U_{k}) - (U_{k}^H (\mu I) U_{k}) + \mu I
$$

Because $\mu$ is a scalar, it can be moved outside the matrix multiplication:

$$
T_{k+1} = U_{k}^H T_k U_{k} - \mu (U_{k}^H I U_{k}) + \mu I
$$

Again, since $U_{k}^H U_{k} = I$, the equation becomes:

$$
T_{k+1} = U_{k}^H T_k U_{k} - \mu I + \mu I = U_{k}^H T_k U_{k}
$$

This proves that $T_{k+1}$ is unitarily similar to $T_k$. Since this relationship holds for every step of the iteration, the final converged matrix $T$ is unitarily similar to the original matrix $A$. Similarity transformations **preserve eigenvalues**, which is why this algorithm correctly computes the eigenvalues of $A$ while converging to the upper-triangular Schur form.

## Convergence and Deflation

This is where the shift and deflation work together. The shift is designed to dramatically speed up convergence, which in turn allows us to deflate the matrix and reduce the problem size.

### How the Shift Accelerates Convergence

The convergence of the *basic* (un-shifted) QR iteration is governed by the ratios of the eigenvalues (e.g., $|\lambda_{i+1} / \lambda_i|^k$). The shifted QR iteration is mathematically equivalent to applying the basic QR iteration to the matrix $T_k - \mu I$. Therefore, its convergence rate is governed by the ratios of the *shifted* eigenvalues.

The effectiveness of the algorithm is in *how* the shift $\mu$ is chosen. We select $\mu = T_k[n,n]$, which, as the iteration progresses, becomes an increasingly accurate approximation of an eigenvalue, $\lambda_n$.

This choice has a powerful effect on the convergence rate. The rate at which the sub-diagonal element $T_k[n, n-1]$ goes to zero is determined by the ratio:

$$
\left| \frac{\lambda_{n} - \mu}{\lambda_{n-1} - \mu} \right|
$$

Because our shift $\mu = T_k[n,n]$ is chosen to be very close to $\lambda_n$, the numerator $|\lambda_n - \mu|$ becomes extremely small. This makes the entire ratio very close to zero, leading to exceptionally fast convergence. In practice, this convergence is often **quadratic** (and **cubic** for symmetric matrices), meaning only one or two iterations are needed per eigenvalue.

This abstract convergence has a direct, observable effect on the matrix $T_k$: the sub-diagonal element $T_k[n, n-1]$ is rapidly driven to zero.

### How Deflation is Used

This convergence is the signal to "deflate."

1.  **Identify Convergence:** We check the sub-diagonal element in the last row, $T_k[n, n-1]$. When this value is small enough (e.g., below a machine-precision threshold), we declare that the $n$-th eigenvalue has converged.

2.  **Record the Eigenvalue:** The converged value $T_k[n,n]$ is now one of the eigenvalues of the original matrix $A$. We record it.

3.  **Deflate the Matrix:** When $T_k[n,n-1] \approx 0$, the matrix has the approximate form:

$$
T_k \approx \begin{pmatrix}
A_{11} & A_{12} \\
0 & \lambda_n
\end{pmatrix}
$$

The eigenvalues of $T_k$ are the eigenvalues of the $(n-1) \times (n-1)$ block $A_{11}$ combined with the eigenvalue $\lambda_n$ we just found.

**Repeat the Process:** We can now "deflate" the problem by discarding the last row and column. The *entire* shifted QR iteration process (choose shift, factor, recombine) is repeated, but only on the **smaller $(n-1) \times (n-1)$ matrix $A_{11}$**.

This cycle repeats: the algorithm runs the shifted QR iteration on the smaller matrix until *its* last row converges, it records another eigenvalue, and it deflates again. This continues until the matrix is reduced to a $1 \times 1$ matrix, and all eigenvalues have been found.