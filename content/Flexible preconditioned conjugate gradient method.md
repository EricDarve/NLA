Is it possible to vary the [[Preconditioned Conjugate Gradients algorithm|preconditioner]] at every step? For example, as we proceed through the iterations, we may be able to improve our guess for $M$. We may want to use this varying preconditioner $M^{(k)}$ in PCG.

The [[Preconditioned Conjugate Gradients algorithm|previous derivation]] becomes invalid because $M^{(k)}$ changes. This leads to a slowdown in the convergence.

The convergence may be improved, nonetheless, if we use the following modified update for $\underline p^{(k+1)}$:
$$
  \underline p^{(k+1)} = z^{(k)} + \frac{[r^{(k)}]^T (z^{(k)}-z^{(k-1)})}{[r^{(k-1)}]^T z^{(k-1)}} \, \underline p^{(k)}
$$
In our previous derivations for [[All the orthogonality relations in CG|CG]] and [[Preconditioned Conjugate Gradients algorithm|PCG,]] the numerator could be simplified, but this is no longer possible with a varying preconditioner. For example, we show below that, in [[Preconditioned Conjugate Gradients algorithm|PCG,]] we have $(r^{(k)})^T z^{(k-1)} = 0,$ so that the ratio is, in PCG:
$$
\tau_k = \frac{[r^{(k)}]^T z^{(k)}}{[r^{(k-1)}]^T z^{(k-1)}}, \quad
\underline p^{(k+1)} = z^{(k)} + \tau_k \, \underline p^{(k)}
$$
Here are some additional details about this result. [[Computationally efficient search directions|Previously in CG,]] when computing $\tau_k,$ we had this numerator:
$$
[r^{(k)}]^T (r^{(k)}-r^{(k-1)})
=
[r^{(k)}]^T r^{(k)}
= \| r^{(k)} \|_2^2
$$
since $[r^{(k)}]^T r^{(k-1)} = 0$.

In PCG, we have a similar relation. Using the [[Preconditioning the Conjugate Gradients algorithm|previous notations]] for the PCG residuals $s^{(k)},$ we get
$$
[s^{(k)}]^T s^{(k-1)} = 0
$$
Since in [[Preconditioned Conjugate Gradients algorithm|PCG,]] we have $s^{(k)} = C r^{(k)}$, we get
$$
[s^{(k)}]^T s^{(k-1)} = [r^{(k)}]^T C^TC r^{(k-1)} = [r^{(k)}]^T M r^{(k-1)} = [r^{(k)}]^T z^{(k-1)} = 0
$$
However, as pointed out above, the relation $[r^{(k)}]^T z^{(k-1)} = 0$ no longer holds exactly in the flexible PCG algorithm since $M^{(k)}$ keeps changing.

The fPCG algorithm can be shown to be locally optimal. That is, it does not converge slower than the steepest descent method, which is locally optimal.