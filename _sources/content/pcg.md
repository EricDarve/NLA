# Preconditioning the Conjugate Gradient (PCG) Algorithm

Our previous discussion on the convergence of the Conjugate Gradient (CG) method highlighted its dependence on the condition number $\kappa(A)$. To accelerate convergence, we must apply a preconditioner, but we face a critical constraint: **we must maintain the Symmetric Positive Definite (SPD) property** of the system matrix.

A general-purpose left preconditioner $MA$ or right preconditioner $AM$ will destroy the symmetry of $A$. We must therefore use a **symmetric preconditioner**. 

```{note}
Below we will denote $M$ instead of $M^{-1}$ (as in the previous section) for convenience and clarity.
```

## The Transformed System

Let us assume we have an SPD preconditioner $M$. Since $M$ is SPD, it admits a factorization $M = C^T C$ (e.g., a Cholesky factorization), where $C$ is a non-singular matrix.

Following the symmetric preconditioning strategy, we will solve an *equivalent* system. Instead of $Ax=b$, we start by multiplying by $C$:

$$
C A x = C b
$$

Now, we introduce a change of variables, $x = C^T y$. Substituting this into our equation gives:

$$
(C A C^T) y = C b
$$

We have now defined a new, transformed linear system $A'y = b'$, where:

* **$A' = C A C^T$**
* **$b' = C b$**
* **$x = C^T y$**

The new matrix $A'$ is SPD (as $A$ is SPD and $C$ is non-singular), so we can apply the standard CG algorithm to $A'y = b'$. The new system $A'y=b'$ is designed to have a much smaller condition number than $A$, i.e., $\kappa(A') \ll \kappa(A)$.

## A "Direct" (but Cumbersome) Implementation

If we apply the standard CG algorithm (as derived in our previous chapter) directly to $A'y=b'$, we get the following:

1.  **Initialize:**
    * $y^{(0)}$ (e.g., $y^{(0)} = 0$)
    * $s^{(0)} = b' - A'y^{(0)} = Cb - (CAC^T)y^{(0)}$
    * $p^{(1)} = s^{(0)}$
2.  **Iterate for $k=1, 2, \dots$**
    * $q^{(k)} = A' p^{(k)} = (C A C^T) p^{(k)}$
    * $\mu_k = \frac{\| s^{(k-1)} \|_2^2}{(p^{(k)})^T q^{(k)}}$
    * $y^{(k)} = y^{(k-1)} + \mu_k p^{(k)}$
    * $s^{(k)} = s^{(k-1)} - \mu_k q^{(k)}$
    * $\tau_k = \frac{\| s^{(k)} \|_2^2}{\| s^{(k-1)} \|_2^2}$
    * $p^{(k+1)} = s^{(k)} + \tau_k p^{(k)}$

This algorithm is mathematically correct, but it is **cumbersome and impractical**.

The primary issue is the explicit appearance of $C$ and $C^T$. This implementation would require us to:

1.  First, compute the factor $C$ from $M$.
2.  In every single iteration, compute $q^{(k)}$ by performing a matrix-vector product with $C^T$, then with $A$, and then with $C$.

This is computationally expensive and complex. Furthermore, we almost never have $C$ and $C^T$ explicitly. A preconditioner $M$ is typically provided as a "black box," i.e., a function that performs a solve and gives $z = M r$ for a given $r$.

## The PCG Simplification

The **Preconditioned Conjugate Gradient (PCG)** algorithm is the solution to this problem. PCG is *algebraically equivalent* to the "cumbersome" algorithm above, but it is derived by performing a careful change of variables back to $x$ and the original residual $r = b - Ax$.

By substituting $x = C^T y$, $p_{pcg} = C^T p_{cg}$, and $s_{cg} = C r_{pcg}$ (where $cg$ denotes the "direct" algorithm's vectors and $pcg$ denotes the final algorithm's vectors), a significant simplification occurs. All references to $C$ and $C^T$ cancel out, and they are replaced by a single operation involving $M$.

This simplified algorithm avoids $C$ entirely, and its only preconditioning-related step is the solution of a single linear system and computing $z = M r$, which is exactly what a preconditioner is designed to do.

## The Krylov Subspace Connection

We provide a brief explanation of *why* the "cumbersome" factor $C$ can be algebraically eliminated, allowing us to build a practical algorithm that only uses the preconditioner $M$. Following this, we will go through the step by step derivation.

Let's break down the logic:

1.  **The "Cumbersome" Space:** If we literally ran CG on the transformed system $A'y = b'$, where $A' = CAC^T$ and $b' = Cb$, the algorithm would build its solution for $y$ from the Krylov subspace:

    $$\mathcal{K}_k(A', b') = \mathcal{K}_k(CAC^T, Cb)$$

    This is the space $\text{span} \{ b', A'b', (A')^2 b', \dots \}$.

2.  **The "Target" Space:** We don't ultimately care about $y$; we want $x$. The relationship between them is $x = C^T y$. Therefore, the space containing our actual solution $x$ is found by multiplying every vector in the $y$-space by $C^T$:

    $$\text{Space for } x = C^T \, \mathcal{K}_k(CAC^T, Cb)$$

3.  **The "Aha!" Moment (The Algebra):** This is where the magic happens. We define our preconditioner as $M = C^T C$. (This requires $M$ to be SPD, as we noted). Now, let's look at the vectors that form this new "target" space for $x$:

    * **The 0-th vector:**
        $C^T (b') = C^T(Cb) = (C^T C)b = \mathbf{Mb}$

    * **The 1-st vector:**
        $C^T (A'b') = C^T(CAC^T)(Cb) = (C^T C) A (C^T C) b = \mathbf{(MA)Mb}$

    * **The 2-nd vector:**
        $C^T ((A')^2 b') = C^T(CAC^T)(CAC^T)(Cb) = (C^T C) A (C^T C) A (C^T C) b = \mathbf{(MA)^2 Mb}$

4.  **The Conclusion:** The "target" space for $x$ is, in fact, a different (and more convenient) Krylov subspace:

    $$C^T \, \mathcal{K}_k(CAC^T, Cb) = \text{span} \{ Mb, (MA)Mb, (MA)^2 Mb, \dots \} = \mathcal{K}_k(MA, Mb)$$

This is a critical insight. It proves that although we *conceptually* started with a transformed system involving $C$, the solution $x$ we are looking for lives in a Krylov subspace generated *only by $M$ and $A$*.

The PCG algorithm is a clever implementation that finds the best solution in $\mathcal{K}_k(MA, Mb)$ without ever needing to know $C$ or $C^T$. It only needs to be able to compute matrix-vector products with $A$ and, in its "preconditioning" step, with $M$.

## Derivation of the Practical PCG Algorithm

We have established that applying the standard CG algorithm to the transformed system $A'y = b'$, where $A' = CAC^T$ and $b' = Cb$, is mathematically sound but operationally cumbersome.

We will now derive the practical **Preconditioned Conjugate Gradient (PCG)** algorithm. We do this by starting with the "cumbersome" CG algorithm (which iterates on $y, s, p$) and performing a careful change of variables to create an equivalent algorithm that works directly with our original solution $x$ and residual $r = b - Ax$.

## The "Dictionary" of Variables

Our derivation relies on a "dictionary" that relates the vectors from the cumbersome algorithm to the new, practical vectors.

* **Solution:** $x^{(k)} = C^T y^{(k)}$
* **Residual:** We use the original residual $r^{(k)} = b - Ax^{(k)}$. The relationship to the cumbersome residual $s^{(k)} = Cb - A'y^{(k)}$ is:
    
    $$
    \begin{aligned}
    s^{(k)} &= Cb - (CAC^T) y^{(k)} \\
        &= Cb - CA(C^T y^{(k)}) \\
        &= C\big(b - A x^{(k)}\big) \\
        &= C r^{(k)}
    \end{aligned}
    $$

* **Search Direction:** We define a new "practical" search direction $p_{pcg}^{(k)}$:
    
    $$p_{pcg}^{(k)} = C^T p^{(k)}$$

* **Preconditioner:** Recall that $M = C^T C$.
* **Work Vector:** We will define a new work vector $z^{(k)}$ to handle the preconditioning step.

We will now transform the cumbersome algorithm, line by line.

### 1. Solution Update

The cumbersome solution update is:

$$
y^{(k)} = y^{(k-1)} + \mu_k \, p^{(k)}
$$

To find our practical update for $x$, we simply multiply the entire equation by $C^T$:

$$
C^T y^{(k)} = C^T y^{(k-1)} + \mu_k \, (C^T p^{(k)})
$$

Using our dictionary ($x^{(k)} = C^T y^{(k)}$ and $p_{pcg}^{(k)} = C^T p^{(k)}$), this immediately simplifies to:

$$
\mathbf{x^{(k)} = x^{(k-1)} + \mu_k \, p_{pcg}^{(k)}}
$$

### 2. Residual Update

The cumbersome residual update is:

$$
s^{(k)} = s^{(k-1)} - \mu_k \, q^{(k)}
$$

We first need to transform $q^{(k)} = A' p^{(k)} = (CAC^T)p^{(k)}$. Let's define an intermediate vector $q_{pcg}^{(k)}$:

$$
q_{pcg}^{(k)} = A (C^T p^{(k)}) = A p_{pcg}^{(k)}
$$

With this definition, we can write $q^{(k)}$ as:

$$
q^{(k)} = C (A p_{pcg}^{(k)}) = C q_{pcg}^{(k)}
$$

Now, substitute the $s$ and $q$ relations into the cumbersome update:

$$
C r^{(k)} = C r^{(k-1)} - \mu_k \, (C q_{pcg}^{(k)})
$$

Since $C$ is non-singular, we can multiply from the left by $C^{-1}$ to cancel it from every term:

$$
\mathbf{r^{(k)} = r^{(k-1)} - \mu_k \, q_{pcg}^{(k)}}
$$

### 3. Search Direction Update

The cumbersome search direction update is:

$$
p^{(k+1)} = s^{(k)} + \tau_k \, p^{(k)}
$$

To find our practical update for $p_{pcg}$, we multiply the entire equation by $C^T$:

$$
C^T p^{(k+1)} = C^T s^{(k)} + \tau_k \, (C^T p^{(k)})
$$

Using our dictionary, this becomes:

$$
p_{pcg}^{(k+1)} = (C^T s^{(k)}) + \tau_k \, p_{pcg}^{(k)}
$$

Now let's simplify the $C^T s^{(k)}$ term using $s^{(k)} = C r^{(k)}$:

$$
C^T s^{(k)} = C^T (C r^{(k)}) = (C^T C) r^{(k)} = M r^{(k)}
$$

To make the algorithm efficient, we introduce a new **work vector $z^{(k)} = M r^{(k)}$**. This is the only place where the preconditioner $M$ is applied.
Substituting this in, we get the final update rule:

$$
\mathbf{p_{pcg}^{(k+1)} = z^{(k)} + \tau_k \, p_{pcg}^{(k)}}
$$

### 4. Step-Length Scalars ($\mu_k$ and $\tau_k$)

The final step is to re-write the scalars $\mu_k$ and $\tau_k$ using only our practical vectors.

**The $\tau_k$ Scalar (Ratio of Residual Norms):**

$$
\tau_k = \frac{\| s^{(k)} \|_2^2}{\| s^{(k-1)} \|_2^2}
$$

Let's analyze the numerator. We can use our new $z$ vector:

$$
\| s^{(k)} \|_2^2 = (s^{(k)})^T s^{(k)} = (C r^{(k)})^T (C r^{(k)})
= (r^{(k)})^T C^T C r^{(k)} = (r^{(k)})^T (M r^{(k)}) = \mathbf{[r^{(k)}]^T z^{(k)}}
$$

This is a great simplification. $C$ has vanished, replaced by $M$, and the norm calculation is transformed into a simple dot product.

The full scalar is therefore:

$$
\mathbf{\tau_k = \frac{[r^{(k)}]^T z^{(k)}}{[r^{(k-1)}]^T z^{(k-1)}}}
$$

**The $\mu_k$ Scalar (Step Length):**

$$
\mu_k = \frac{ \| s^{(k-1)} \|_2^2 }{(p^{(k)})^T q^{(k)}}
$$

The numerator is simply the one we just found, for step $k-1$:

$$
\| s^{(k-1)} \|_2^2 = \mathbf{[r^{(k-1)}]^T z^{(k-1)}}
$$

The denominator is:

$$
(p^{(k)})^T q^{(k)} = (p^{(k)})^T (C q_{pcg}^{(k)}) = (C^T p^{(k)})^T q_{pcg}^{(k)} = \mathbf{(p_{pcg}^{(k)})^T q_{pcg}^{(k)}}
$$

Combining these gives the final, practical expression:

$$
\mathbf{\mu_k = \frac{[r^{(k-1)}]^T z^{(k-1)}}{(p_{pcg}^{(k)})^T q_{pcg}^{(k)}}}
$$

## Summary: The Final PCG Algorithm

By collecting all the practical, simplified steps (highlighted in bold), we arrive at the complete Preconditioned Conjugate Gradient algorithm. Note how it no longer contains $C$, $C^T$, $y$, $s$, or $p$. It is a self-contained algorithm operating on $x, r, p_{pcg},$ and $z$.

```{important}
### Algorithm: Preconditioned Conjugate Gradient (PCG)

**Goal:** Solve $Ax=b$, given an SPD matrix $A$ and an SPD preconditioner $M$.

1.  **Initialize:**
    * Choose $x^{(0)}$ (e.g., $x^{(0)} = 0$)
    * $r^{(0)} = b - Ax^{(0)}$
    * $z^{(0)} = M r^{(0)}$
    * $p_{pcg}^{(1)} = z^{(0)}$
    * $\rho_0 = [r^{(0)}]^T z^{(0)}$

2.  **Iterate for $k=1, 2, \dots$ until convergence:**
    * $q_{pcg}^{(k)} = A p_{pcg}^{(k)}$
    * $\mu_k = \frac{\rho_{k-1}}{(p_{pcg}^{(k)})^T q_{pcg}^{(k)}}$
    * $x^{(k)} = x^{(k-1)} + \mu_k p_{pcg}^{(k)}$
    * $r^{(k)} = r^{(k-1)} - \mu_k q_{pcg}^{(k)}$

    * *Check for convergence (e.g., $\|r^{(k)}\|_2 < \epsilon$)*

    * $z^{(k)} = M r^{(k)}$
    * $\rho_k = [r^{(k)}]^T z^{(k)}$
    * $\tau_k = \frac{\rho_k}{\rho_{k-1}}$
    * $p_{pcg}^{(k+1)} = z^{(k)} + \tau_k p_{pcg}^{(k)}$

3.  **Return $x^{(k)}$**
```

## Python Implementation

Below is a Python implementation of the Preconditioned Conjugate Gradient (PCG) algorithm.

This implementation is written as a self-contained function. It assumes that `A` and `M` are provided as matrices or as functions that can perform a matrix-vector product.

```python
import numpy as np

def pcg(A, b, M, tol=1e-10, max_iter=None):
    """
    Solves the system Ax=b using the Preconditioned Conjugate Gradient (PCG) algorithm.

    This implementation follows our derivation, where the
    preconditioned system is A' = CAC^T and the preconditioner matrix
    is defined as M = C^T C.

    The key preconditioning step is z = M*r, which is a matrix-vector
    product, not a linear solve. This implies M is an approximation
    of A^{-1}.

    Args:
        A (np.ndarray or func): The system matrix A (n x n).
                                 Can be a dense/sparse matrix or a function
                                 that computes A @ v.
        b (np.ndarray): The right-hand side vector (n x 1).
        M (np.ndarray or func): The preconditioner matrix M (n x n).
                                 Can be a dense/sparse matrix or a function
                                 that computes M @ v.
        tol (float, optional): The convergence tolerance. The loop stops
                               when the M-norm of the residual,
                               sqrt(r.T @ z), is less than tol.
        max_iter (int, optional): The maximum number of iterations.
                                 Defaults to n*10 if not set.

    Returns:
        np.ndarray: The solution vector x (n x 1).
    """
    
    # --- Setup ---
    # Determine how to apply A and M (as matrix or function)
    if callable(A):
        matvec_A = A
    else:
        matvec_A = lambda v: A @ v
        
    if callable(M):
        matvec_M = M
    else:
        matvec_M = lambda v: M @ v

    n = len(b)
    if max_iter is None:
        max_iter = n * 10  # Set a reasonable default max iterations

    # --- Initialization (k=0) ---
    x = np.zeros(n)
    r = b.copy()          # r_0 = b - A*x_0 (since x_0=0)
    z = matvec_M(r)       # z_0 = M*r_0
    p = z.copy()          # p_1 = z_0
    
    # rho_0 = r_0^T * z_0 = ||r_0||_M^2
    rho = r @ z
    
    # --- Iteration (k=1, 2, ...) ---
    for k in range(max_iter):
        
        # M-norm of residual. (r^T * M * r)^(1/2)
        norm_r_M = np.sqrt(rho)
        if norm_r_M < tol:
            print(f"PCG converged in {k} iterations (M-norm of residual < {tol}).")
            return x

        q = matvec_A(p)     # q_k = A*p_k
        
        # mu_k = (r_{k-1}^T * z_{k-1}) / (p_k^T * q_k)
        mu = rho / (p @ q)
        
        x += mu * p         # x_k = x_{k-1} + mu_k * p_k
        r -= mu * q         # r_k = r_{k-1} - mu_k * q_k
        z = matvec_M(r)     # z_k = M*r_k
        
        rho_old = rho
        rho = r @ z         # rho_k = r_k^T * z_k
        
        tau = rho / rho_old  # tau_k = rho_k / rho_{k-1}
        p = z + tau * p      # p_{k+1} = z_k + tau_k * p_k

    print(f"PCG failed to converge after {max_iter} iterations.")
    return x
```