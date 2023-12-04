[[Conjugate Gradients algorithm|CG:]] applies when the matrix is symmetric positive definite.

GMRES: applies to any matrix, but there is a higher computational cost. This is similar to the difference between the [[Lanczos process]] and the [[Arnoldi process]].

GMRES is simpler than CG because fewer optimizations are possible.

“Just solve a least-squares problem.”

## Norm to measure error

The key in GMRES is to change the [[Conjugate Gradients Version 1|cost function]] we are optimizing:
$$
\| x - x^{(k)} \|_?
$$

[[Conjugate Gradients Version 1|CG choice:]] $\| z \|_A = \sqrt{z^T A z}$

GMRES choice: $\| z \|_{A^TA} = \sqrt{z^T A^TA z}$

This norm makes more sense once we apply it:
$$
\begin{align}
\| x - x^{(k)} \|_{A^TA}^2 & = (x - Q_k y)^T A^TA (x - Q_k y) \\[.5em]
& = \| b - A Q_k y \|_2^2 = \| r^{(k)} \|_2^2
\end{align}
$$

- We can now explain the acronym! 
- GMRES = Generalized Minimal RESidual method.