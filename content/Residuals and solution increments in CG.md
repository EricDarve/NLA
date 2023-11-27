We state and prove a key theorem to build the CG algorithm.

**Theorem.** As long as we have not converged, then
$$
\begin{aligned}
\text{span} \left\{ r^{(0)}, \ldots, r^{(k)} \right\} 
& = \text{span} \left\{ \Delta x^{(0)}, \ldots, \Delta x^{(k)} \right\} \\[.3em]
& = \mathcal K(A, b, k+1)
\end{aligned}
$$
where
$$
r^{(k)} = b - A x^{(k)}, \quad \Delta x^{(k)} = x^{(k+1)} - x^{(k)}
$$

### Proof

**Step 1**
[[Krylov methods for sparse systems|Recall]] that we search a solution in the Krylov subspace:
$$
\begin{gather}
x^{(k)} \in \mathcal K(A, b, k)\\[.3em]
\Delta x^{(k)} = x^{(k+1)} - x^{(k)}
\end{gather}
$$
This implies that
$$
\text{span} \left\{ \Delta x^{(0)}, \ldots, \Delta x^{(k)} \right\} \subset \mathcal K(A, b, k+1)
$$
The steps $\Delta x^{(j)}$ are [[Orthogonality relations in CG|conjugate]] so they must be linearly independent. Since the dimensions of the spaces match, we have
$$
\text{span} \left\{ \Delta x^{(0)}, \ldots, \Delta x^{(k)} \right\} = \mathcal K(A, b, k+1)
$$

**Step 2**
$$
\begin{gather}
x^{(k)} \in \mathcal K(A, b, k)\\[.3em]
A x^{(k)} \in \mathcal K(A, b, k+1)
\end{gather}
$$
by definition of the [[Krylov subspace]]. This implies that:
$$
\Rightarrow \quad r^{(k)} = b - A x^{(k)} \in \mathcal K(A, b, k+1)
$$
So
$$
\text{span} \left\{ r^{(0)}, \ldots, r^{(k)} \right\}
\subset
\mathcal K(A, b, k+1)
$$
Let's assume that the dimension of $\text{span} \left\{ r^{(0)}, \ldots, r^{(k)} \right\}$ is less than $k+1$.

Then there is a $j$ such that $r^{(j)} \in \mathcal K(A, b, j).$ [[Orthogonality relations in CG#^9bf4bf|We also have]] that $r^{(j)} \perp \mathcal K(A, b, j).$ So $r^{(j)} = 0,$ which is a contradiction. So
$$
\text{span} \left\{ r^{(0)}, \ldots, r^{(k)} \right\}
=
\mathcal K(A, b, k+1)
$$
We have proved that
$$
\text{span} \left\{ r^{(0)}, \ldots, r^{(k)} \right\}
=
\text{span} \left\{ \Delta x^{(0)}, \ldots, \Delta x^{(k)} \right\}
$$
unless CG has converged to the [[Convergence with conjugate steps|exact solution.]]

$\square$