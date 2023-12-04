We follow some of the [[Preconditioning|strategies outlined earlier]] and consider a symmetric preconditioning method. Assume we have a matrix $C$ such that $CAC^T$ has good convergence properties for CG, e.g., a smaller condition number.

$CAC^T$ is [[Symmetric Positive Definite Matrices|SPD]] if $C$ is **non-singular.** This is a key property to ensure that [[Conjugate Gradients algorithm|CG]] applies.

We now consider solving the modified linear system using [[Conjugate Gradients algorithm|CG]]: 
$$
CAC^T y = Cb, \quad x = C^T y.
$$

Here is a possible implementation that follows directly the [[Conjugate Gradients algorithm|CG algorithm:]]
$$
\begin{aligned}
  q^{(k)} & = CAC^T p^{(k)} \\[.3em]
  \mu_k & = \frac{ \| s^{(k-1)} \|_2^2 }{(p^{(k)})^T q^{(k)}} \\[.3em]
  y^{(k)} & = y^{(k-1)} + \mu_k p^{(k)} \\[.3em]
  s^{(k)} & = s^{(k-1)} - \mu_k q^{(k)} \\[.3em]
  \tau_k & = \frac{\| s^{(k)} \|_2^2}{\| s^{(k-1)} \|_2^2} \\[.3em]
  p^{(k+1)} & = s^{(k)} + \tau_k p^{(k)}
\end{aligned}
$$
This implementation is cumbersome because of the need to multiply by $C$ and $C^T$. Typically, preconditioners are of the form $M = C^T C$, and we consider $MA$ as the new matrix. In this case, we are only required to be able to multiply by $M$, not $C$ and $C^T$.

The preconditioned CG algorithm allows getting rid of $C$ in favor of $M$. This is a significant simplification in the implementation.

PCG is equivalent to the algorithm written above but does the algebra slightly differently to avoid the appearance of $C$.

## Krylov subspace

The reason why $C$ can be avoided can be traced back to how the [[Krylov subspace]] is constructed.

When using $CAC^T$, in the [[Krylov subspace]], we use powers of $CAC^T$ and we get:

$$(CAC^T)^k = CA (C^TCA)^{k-1} C^T = CA (MA)^{k-1} C^T$$

We can make this even cleaner with

$$C^T (CAC^T)^k \, Cb = C^TCA (MA)^{k-1} C^TCb = (MA)^k Mb$$

Only $M$ appears!

In terms of [[Krylov subspace]], we can write:

$$\mathcal K(CAC^T,Cb,k) = \text{span} \{Cb, CAC^T\, Cb, (CAC^T)^2 Cb, \dots, (CAC^T)^{k-1} Cb\}$$

Multiply by $C^T$ and we just get:

$$C^T \, \mathcal K(CAC^T,Cb,k) = \text{span} \{Mb, MA \, Mb, (MA)^2 Mb, \dots, (MA)^{k-1} Mb\}$$

This seems promising. We can construct a [[Krylov subspace]] with $M$ only.

Note that since we require $M = C^TC$, with $C$ non-singular, the preconditioner $M$ must be [[[Symmetric Positive Definite Matrices|SPD]].