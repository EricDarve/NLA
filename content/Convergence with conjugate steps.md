[[Conjugate Gradients Version 1|Recall]] that the cost function in CG is quadratic: 
$$
\| x - x^{(k)} \|_A^2 = (x-x^{(k)})^T A (x-x^{(k)}).
$$
After [[Orthogonality relations in CG|scaling]] by $A^{1/2}$, the level curves of the cost function become simple circles. At each step, we make an incremental correction to the solution that is [[Dot product|orthogonal]] to all the previous steps (with the $A^{1/2}$ scaling). This is illustrated in the figure below.

![[2022-11-16-11-09-11.png|500]]

**Corollary.** We obtain the exact solution in at most $n$ steps, where $n$ is the size of matrix $A$.