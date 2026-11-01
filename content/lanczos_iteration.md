# The Lanczos Process

For a symmetric matrix, Arnoldi's projection has additional structure: the Hessenberg matrix is also symmetric, so it is tridiagonal. This reduces the orthogonalization to a **three-term recurrence**. Each new direction depends only on the current basis vector and its predecessor.

We derive this simplification for $A=A^T\in\mathbb R^{n\times n}$ in exact arithmetic. For complex Hermitian matrices, replace transpose by conjugate transpose. The implementation below supports both cases.

## Symmetry Makes the Projection Tridiagonal

The Arnoldi basis satisfies $Q_k^TQ_k=I_k$, and its projected matrix is $H_k=Q_k^TAQ_k$. Symmetry gives

$$
H_k^T=Q_k^TA^TQ_k=Q_k^TAQ_k=H_k.
$$

The Hessenberg structure places zeros below the first subdiagonal. Symmetry reflects those zeros above the first superdiagonal. We therefore write $T_k=H_k$ and obtain

$$
T_k=
\begin{pmatrix}
\alpha_1&\beta_1&&\\
\beta_1&\alpha_2&\ddots&\\
&\ddots&\ddots&\beta_{k-1}\\
&&\beta_{k-1}&\alpha_k
\end{pmatrix}.
$$

Here $\alpha_j=h_{jj}$ and $\beta_j=h_{j+1,j}$ are the Arnoldi coefficients. The change from $H_k$ to $T_k$ emphasizes that the small projected matrix is tridiagonal.

## Why Only Two Components Must Be Removed

Start with a unit vector $q_1$ and set $q_0=0$, $\beta_0=0$. The Arnoldi column relation now has just three terms:

$$
\boxed{
Aq_j=\beta_{j-1}q_{j-1}
+\alpha_jq_j+\beta_jq_{j+1}.
}
$$

Only the components along $q_{j-1}$ and $q_j$ need to be removed before normalizing the new vector. The components along older basis vectors are already zero.

To see why, take $i\le j-2$. Symmetry allows us to move $A$ to the other side of the inner product:

$$
q_i^TAq_j=(Aq_i)^Tq_j=0.
$$

The last equality holds because $Aq_i$ lies in the span of $q_1,\ldots,q_{i+1}$, all perpendicular to $q_j$. This is what makes the short recurrence possible: orthogonality to the older directions follows from symmetry.

At step $j$, compute

$$
\begin{aligned}
w&=Aq_j-\beta_{j-1}q_{j-1},\\
\alpha_j&=q_j^Tw,\\
w&\leftarrow w-\alpha_jq_j,\\
\beta_j&=\|w\|_2,\qquad
q_{j+1}=\frac{w}{\beta_j}\quad(\beta_j>0).
\end{aligned}
$$

The first subtraction removes the known component along the preceding direction. The second removes the component along the current direction. What remains gives the next basis vector. If $\beta_j=0$, the space is invariant and the process stops, just as in Arnoldi.

## The Lanczos Relation

Collecting the first $k$ column equations gives

$$
\boxed{
AQ_k=Q_kT_k+f_ke_k^T,
\qquad Q_k^Tf_k=0.
}
$$

This is the Arnoldi relation specialized to the symmetric case. When $\beta_k>0$,

$$
f_k=\beta_kq_{k+1},\qquad \|f_k\|_2=\beta_k.
$$

The coefficient $\beta_k$ is **not an entry of $T_k$**. It couples the current space to the next direction; $T_k$ contains only $\beta_1,\ldots,\beta_{k-1}$. We retain $f_k$ to check the Ritz residuals, as explained in [Residuals and Eigenvector Accuracy](sec:lanczos-residuals). At exact breakdown, $f_k=0$ and no $q_{k+1}$ is needed.

## Basic Implementation

The following code implements the three-term recurrence without reorthogonalization. It returns $Q_k$, $T_k$, and the final remainder $f_k$, matching the output convention of the Arnoldi implementation. The maximum number of steps is $k$; early termination returns smaller matrices.

```python
import numpy as np

def lanczos_iteration(A, q1, k, breakdown_tol=None):
    """Basic Lanczos recurrence for symmetric or Hermitian A.

    Return Q, T, f with A @ Q approximately Q @ T + outer(f, e_last).
    This implementation stores Q and does not reorthogonalize.
    """
    q = np.asarray(q1)
    if q.ndim != 1 or A.shape != (q.size, q.size):
        raise ValueError("A must be square and q1 a matching vector.")
    if not isinstance(k, (int, np.integer)) or not 1 <= k <= q.size:
        raise ValueError("k must be an integer between 1 and n.")

    dtype = np.result_type(A.dtype, q.dtype, np.float64)
    q = q.astype(dtype, copy=True)
    qnorm = np.linalg.norm(q)
    if qnorm == 0:
        raise ValueError("The starting vector must be nonzero.")
    q /= qnorm

    if breakdown_tol is None:
        breakdown_tol = 10 * np.finfo(q.real.dtype).eps
    if breakdown_tol < 0:
        raise ValueError("breakdown_tol must be nonnegative.")

    q_prev = np.zeros_like(q)
    beta_prev = 0.0
    Q, alpha, beta = [], [], []

    for j in range(k):
        Q.append(q)
        Aq = A @ q
        w = Aq - beta_prev * q_prev
        a = float(np.vdot(q, w).real)
        w = w - a * q
        alpha.append(a)
        b = np.linalg.norm(w)

        # b is the final coupling, outside T, if we stop here.
        if b <= breakdown_tol * np.linalg.norm(Aq) or j + 1 == k:
            T = np.diag(alpha)
            if beta:
                T += np.diag(beta, 1) + np.diag(beta, -1)
            return np.column_stack(Q), T, w

        beta.append(b)
        q_prev, q = q, w / b
        beta_prev = b
```

The final step computes the remainder even when no further basis vector is requested. The threshold is a relative test for a small new direction, not a tolerance on eigenvalue error.

## What the Short Recurrence Saves

Without reorthogonalization, each step requires one matrix-vector product and $O(n)$ additional vector work. Building $T_k$ therefore takes $k$ products and $O(nk)$ additional work, compared with Arnoldi's $O(nk^2)$ orthogonalization work.

The recurrence itself needs only a few length-$n$ working vectors and the tridiagonal coefficients. However, **storing the basis $Q_k$ costs $O(nk)$ memory**. The implementation above does store it, so it can form Ritz vectors. The smaller $O(n)$ working-memory requirement applies when the basis is not retained. These costs concern the construction; solving the small eigenvalue problem adds further work.

## Orthogonality in Floating-Point Arithmetic

The orthogonality argument used symmetry and exact arithmetic. Rounding errors can introduce components along older directions, and the short recurrence does not explicitly remove them. The computed basis can therefore lose orthogonality.

A common consequence is the appearance of extra copies of already-converged eigenvalues among the eigenvalues of $T_k$. These are often called *ghost Ritz values*. They do not reliably indicate a true eigenvalue multiplicity.

Full reorthogonalization removes components along every stored basis vector, restoring an orthogonalization cost of $O(nk^2)$. Selective reorthogonalization instead focuses on converged Ritz directions. This controls the tendency to rediscover eigenvectors already found; [the Netlib discussion of reorthogonalization](https://netlib.org/utk/people/JackDongarra/etemplates/node108.html) explains the mechanism.

The basic implementation above illustrates the short recurrence, but long floating-point runs may require these additional orthogonalization steps. When assessing an eigenpair, an explicitly computed residual provides a useful check.

The [next section](lanczos_convergence.md) explains the second benefit of symmetry: the largest and smallest Ritz values have variational characterizations that turn accurate directions in the Krylov space into eigenvalue error bounds.
