Is it possible to vary the [[Preconditioned Conjugate Gradients algorithm|preconditioner]] at every step? For example, as we proceed through the iterations, we may be able to improve our guess for $M$. We may want to use this varying preconditioner $M^{(k)}$ in PCG.

The [[Preconditioned Conjugate Gradients algorithm|previous derivation]] becomes invalid because $M^{(k)}$ changes. This leads to a slowdown in the convergence.

The convergence may be improved, nonetheless, if we use the following modified update for $\underline p^{(k+1)}$:
$$
  \underline p^{(k+1)} = z^{(k)} + \frac{[r^{(k)}]^T (z^{(k)}-z^{(k-1)})}{[r^{(k-1)}]^T z^{(k-1)}} \, \underline p^{(k)}
$$
In regular [[All the orthogonality relations in CG|CG]] and [[Preconditioned Conjugate Gradients algorithm|PCG,]] we have $(r^{(k)})^T z^{(k-1)} = 0,$ but this is no longer true with a varying preconditioner.

The fPCG algorithm can be shown to be locally optimal. That is, it does not converge slower than the steepest descent method, which is locally optimal.