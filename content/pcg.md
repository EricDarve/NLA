# Preconditioned Conjugate Gradient

The [preconditioning survey](preconditioning.md) described how an inexpensive approximate solve can improve a Krylov method. For CG, the transformation must preserve symmetry and positive definiteness. **Preconditioned conjugate gradient (PCG)** implements that transformation without forming the transformed matrix or the preconditioner's factors.

We assume that $A$ and $M$ are real, symmetric positive definite (SPD), and that $M$ is fixed. As in the preceding section, **$M$ approximates $A$, and applying the preconditioner means solving $Mz=r$**. We keep the CG notation $x_k$, $r_k=b-Ax_k$, $p_k$, $\mu_k$, and $\tau_k$.

## The Transformed System

Write $M=CC^T$, with $C$ nonsingular. The symmetric transformation from the preceding section is

$$
\widetilde A=C^{-1}AC^{-T},
\qquad \widetilde b=C^{-1}b,
\qquad y=C^Tx,
$$

so $\widetilde A y=\widetilde b$. Here $C^{-T}=(C^T)^{-1}$. For every nonzero vector $v$,

$$
v^T\widetilde A v=(C^{-T}v)^TA(C^{-T}v)>0,
$$

and $\widetilde A$ is symmetric. Ordinary CG therefore applies to this transformed SPD system. A useful preconditioner improves its eigenvalue distribution or condition number. The relevant condition number is $\kappa_2(\widetilde A)$, not generally $\kappa_2(M^{-1}A)$.

The factors explain the construction. The practical algorithm below needs only products with $A$ and solves with $M$.

## Returning to the Original Variables

Start transformed CG from $y_0=C^Tx_0$. Write $s_k=\widetilde b-\widetilde A y_k$ for its residual, $\widetilde p_k$ for its search direction, and $\widetilde w_k=\widetilde A\widetilde p_k$. Tildes distinguish these directions and products from those in the original coordinates.

The solution and residual satisfy

$$
\begin{aligned}
x_k&=C^{-T}y_k,\\
s_k&=C^{-1}b-C^{-1}AC^{-T}y_k\\
&=C^{-1}(b-Ax_k)=C^{-1}r_k.
\end{aligned}
$$

Map the search direction back in the same way as the solution:

$$
p_k=C^{-T}\widetilde p_k,
\qquad w_k=Ap_k,
\qquad \widetilde w_k=C^{-1}w_k.
$$

Multiplying the transformed solution update by $C^{-T}$ and the transformed residual update by $C$ gives

$$
\begin{aligned}
x_k&=x_{k-1}+\mu_kp_k,\\
r_k&=r_{k-1}-\mu_kw_k.
\end{aligned}
$$

These have exactly the same form as ordinary CG. The change appears in the direction update. Transforming $\widetilde p_{k+1}=s_k+\tau_k\widetilde p_k$ gives

$$
\begin{aligned}
p_{k+1}
&=C^{-T}s_k+\tau_kp_k\\
&=C^{-T}C^{-1}r_k+\tau_kp_k\\
&=M^{-1}r_k+\tau_kp_k.
\end{aligned}
$$

Define the **preconditioned residual** $z_k$ by the solve

$$
\boxed{Mz_k=r_k.}
$$

The new search direction is therefore $p_{k+1}=z_k+\tau_kp_k$, starting from $p_1=z_0$.

### The Scalar Coefficients

Ordinary CG on the transformed system uses $\|s_k\|_2^2$. In original variables,

$$
\|s_k\|_2^2
=(C^{-1}r_k)^T(C^{-1}r_k)
=r_k^TM^{-1}r_k=r_k^Tz_k.
$$

The step-length denominator also transforms directly:

$$
\widetilde p_k^T\widetilde w_k
=\widetilde p_k^TC^{-1}w_k
=p_k^Tw_k.
$$

Thus the coefficients are

$$
\boxed{
\mu_k=\frac{r_{k-1}^Tz_{k-1}}{p_k^Tw_k},
\qquad
\tau_k=\frac{r_k^Tz_k}{r_{k-1}^Tz_{k-1}}.
}
$$

Only the residual inner products change from the unpreconditioned CG formulas. Since $M^{-1}$ is SPD, $r_k^Tz_k>0$ for a nonzero residual. Likewise, $p_k^TAp_k>0$ for a nonzero direction.

## The Search Space and Optimality

Transformed CG searches in $y_0+\mathcal K_k(\widetilde A,s_0)$. To identify the corresponding correction space for $x$, map its generating vectors by $C^{-T}$. The first two become

$$
\begin{aligned}
C^{-T}s_0&=M^{-1}r_0,\\
C^{-T}\widetilde A s_0&=M^{-1}AM^{-1}r_0.
\end{aligned}
$$

Continuing this pattern gives

$$
\boxed{
x_k\in x_0+\mathcal K_k(M^{-1}A,M^{-1}r_0).
}
$$

The transformed energy error equals the original energy error. If $y_*=C^Tx$ is the transformed exact solution, then

$$
\|y_*-y_k\|_{\widetilde A}^2
=(x-x_k)^TC\widetilde A C^T(x-x_k)
=\|x-x_k\|_A^2.
$$

Consequently, PCG minimizes the original energy error over its preconditioned search space. The standard CG bound applies with $\kappa_2(\widetilde A)$. We obtain a new search space while retaining the objective that made CG useful.

The orthogonality relations transform as well. In exact arithmetic, the search directions remain A-orthogonal, while residuals are orthogonal in the $M^{-1}$ inner product:

$$
p_i^TAp_j=0,
\qquad r_i^TM^{-1}r_j=0
\quad(i\ne j).
$$

The original residuals need not be Euclidean-orthogonal, unlike in unpreconditioned CG.

## The PCG Algorithm

Choose $x_0$ and compute $r_0=b-Ax_0$. If it is already sufficiently small, stop. Otherwise solve $Mz_0=r_0$ and set $p_1=z_0$. For $k=1,2,\ldots$, perform

$$
\begin{aligned}
w_k&=Ap_k,\\
\mu_k&=\frac{r_{k-1}^Tz_{k-1}}{p_k^Tw_k},\\
x_k&=x_{k-1}+\mu_kp_k,\\
r_k&=r_{k-1}-\mu_kw_k.
\end{aligned}
$$

Check for convergence before generating another direction. If more work is needed, solve $Mz_k=r_k$ and set

$$
\tau_k=\frac{r_k^Tz_k}{r_{k-1}^Tz_{k-1}},
\qquad p_{k+1}=z_k+\tau_kp_k.
$$

Each update uses one matrix-vector product and one preconditioner solve, along with inner products and vector updates. Only a fixed number of vectors is stored. These are the standard [preconditioned CG operations](https://netlib.org/linalg/html_templates/node54.html); the factors $C$ do not appear in the implementation.

## Python Implementation

Here `A` is a matrix or a function computing $Av$, and `solve_M(r)` returns the solution of $Mz=r$. It must represent a fixed linear SPD operation, as assumed in the derivation. For example, a diagonal preconditioner uses componentwise division, not multiplication by the diagonal.

The stopping target is $\|b-Ax_k\|_2\le\mathrm{tol}\,\|r_0\|_2$, relative to the initial residual, as in the GMRES example. This avoids making the tolerance depend on an arbitrary scaling of $M$. The code monitors the recursively updated residual, then explicitly verifies the residual before returning. That final check costs an additional product with $A$.

```python
import numpy as np


def pcg(A, b, solve_M, x0=None, tol=1e-10, max_iter=None):
    """Return x when the original residual meets the relative tolerance.

    A is real SPD. solve_M(r) solves M z = r for a fixed SPD M.
    Raise RuntimeError if convergence cannot be verified.
    """
    b = np.asarray(b, dtype=float)
    if b.ndim != 1 or b.size == 0:
        raise ValueError("b must be a nonempty vector.")
    if not np.isfinite(tol) or tol <= 0:
        raise ValueError("tol must be positive and finite.")
    if max_iter is None:
        max_iter = 10 * b.size
    if not isinstance(max_iter, (int, np.integer)) or max_iter < 0:
        raise ValueError("max_iter must be a nonnegative integer.")

    matvec = A if callable(A) else lambda v: A @ v
    x = np.zeros_like(b) if x0 is None else np.array(x0, dtype=float, copy=True)
    if x.shape != b.shape:
        raise ValueError("x0 must have the same shape as b.")
    r = b - matvec(x)
    norm_r0 = np.linalg.norm(r)
    target = tol * norm_r0
    if norm_r0 <= target:
        return x
    z = np.asarray(solve_M(r), dtype=float)
    rz = r @ z
    if not np.isfinite(rz) or rz <= 0:
        raise ValueError("The preconditioner must be positive definite.")
    p = z.copy()

    for k in range(max_iter):
        w = matvec(p)
        denominator = p @ w
        if not np.isfinite(denominator) or denominator <= 0:
            raise ValueError("A must be positive definite.")
        mu = rz / denominator
        x += mu * p
        r -= mu * w

        if np.linalg.norm(r) <= target:
            if np.linalg.norm(b - matvec(x)) <= target:
                return x
            raise RuntimeError("The explicit residual did not meet the tolerance.")

        z = np.asarray(solve_M(r), dtype=float)
        rz_new = r @ z
        if not np.isfinite(rz_new) or rz_new <= 0:
            raise ValueError("The preconditioner must be positive definite.")
        tau = rz_new / rz
        p = z + tau * p
        rz = rz_new

    raise RuntimeError("PCG did not converge within max_iter updates.")
```

For a Jacobi preconditioner $M=\operatorname{diag}(A)$, supply `solve_M = lambda r: r / A.diagonal()`. More advanced preconditioners replace this function with triangular, block, or multilevel solves from the preceding survey. The interface stays the same.

PCG is ordinary CG in transformed coordinates, expressed using the original variables. Its exact-arithmetic convergence and optimality follow from that connection. In computation, the quality and cost of the preconditioner, loss of orthogonality, and the residual stopping check determine how useful that connection is in practice.
