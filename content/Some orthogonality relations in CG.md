The name [[Conjugate Gradients Version 1|conjugate gradients]] comes from the fact that the search directions are [[Dot product|orthogonal]] to each other with the right metrics. A special orthogonality relation needs to be used to prove this.

**Definition.** We say that two non-zero vectors $u$ and $v$ are $A$-conjugate with respect to $A$ if
$$
u^T A v = 0
$$

Denote by 
$$
\Delta x^{(k)} = x^{(k+1)} - x^{(k)}.
$$

**Theorem.** The solution increments $\Delta x^{(l)}$ and $\Delta x^{(k)},$ $l \ne k,$ are $A$-conjugate.

**Proof.**
$$
(Q_k^T A Q_k) y_k = Q_k^T b \; \Rightarrow \; Q_k^T (b - A x^{(k)}) = 0.
$$

This is a key result: the residual $b - A x^{(k)}$ is orthogonal to the [[Krylov subspace|Krylov subspace]] ${\mathcal K}(A,b,k)$. ^9bf4bf

Denote by $r^{(k)} = b - A x^{(k)}.$ Then we have that
$$
\begin{gather}
Q_k^T r^{(k)} = 0
\; \text{ and } \;
Q_k^T r^{(k+1)} = 0 \\[.3em]
\Delta x^{(k)} = x^{(k+1)} - x^{(k)} \\[.3em]
A\Delta x^{(k)} = r^{(k)} - r^{(k+1)} \\[.3em]
\Rightarrow \quad Q_k^T A \Delta x^{(k)} = 0
\end{gather}
$$
But: 
$$
\begin{gather}
\Delta x^{(l)} = x^{(l+1)} - x^{(l)} \in Q_{l+1} \\
\Rightarrow \quad (\Delta x^{(l)})^T A \Delta x^{(k)} = 0, \; l < k
\end{gather}
$$
$\square$

The steps in CG can be visualized as shown below. If we multiply the vectors by $A^{1/2}$ then each step is [[Dot product|orthogonal]] to all the previous ones. This looks like a street map of Manhattan!

![[2022-11-16-11-07-38.png|600]]
