# Backward Stability of Householder and Givens QR

This is perhaps the most important theoretical guarantee for these methods. While rounding errors mean we don't compute the exact $Q$ and $R$ for a given $A$, we can prove something stronger.

````{prf:theorem} Backward Stability of QR Factorization
:label: thm:backward_stability_qr
When QR factorization is computed using a sequence of Householder or Givens transformations in floating-point arithmetic, the computed factors $\hat{Q}$ and $\hat{R}$ are the *exact* orthogonal and upper-triangular factors for a slightly perturbed matrix $A + \delta A$. The size of the perturbation is small, bounded by:

$$
\frac{\|\delta A\|}{\|A\|} = O(\epsilon_{\text{machine}})
$$

where $\epsilon_{\text{machine}}$ is the machine precision.
````

**Why it's important:** This result provides a rock-solid guarantee. It tells us that the algorithm gives the right answer for a problem that is very close to the one we started with. For most real-world problems where the input data $A$ has some inherent uncertainty, this is a fantastic result. It's the gold standard for numerical algorithms.