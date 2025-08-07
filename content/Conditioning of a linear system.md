Let's connect the forward and backward errors when solving a linear system.

Assume that $Ax = b$ and $(A + E) \tilde{x} = b + e$ where $\tilde{x}$ is the solution produced by some numerical algorithm. 

$E$ and $e$ are the backward errors.

The forward error $\tilde{x} - x$ can be bounded using:
$$
\frac{\color{green}{\| \tilde{x} - x \|}}{\|x\|} \le
\frac{\color{blue}{\kappa(A)}}{1 - \kappa(A) \frac{\|E\|}{\|A\|}}
\left( \frac{\color{red}{\|E\|}}{\|A\|} + \frac{\color{red}{\|e\|}}{\|b\|} \right)
$$
- $x$ stands for the output now, when it was an input before.
- $\kappa(A) = \| A \|_2 \; \| A^{-1} \|_2$ is the **condition number of matrix $A$.**
- This result assumes that $\kappa(A) \frac{\|E\|}{\|A\|}$ is small compared to 1.

Illustration of the condition number of matrix

![[Pasted image 20231001143903.png]]

The condition number $\kappa(A)$ cannot be changed. Only the backward errors $\|E\|$ and $\|e\|$ can be controlled. **A good algorithm has a backward error on the order of $u$.**

**Stable algorithm:** an algorithm where $\|E\|$ and $\|e\|$ can be controlled. Typically, stable algorithms achieve $\|E\| \in O(u)$ and $\|e\| \in O(u)$ where $u$ is the [[Floating point arithmetic and unit roundoff error|unit roundoff error]].

[[Stability of the LU factorization]], [[Sensitivity analysis]], [[Forward and backward error]], [[Floating point arithmetic and unit roundoff error]]