[[Computing multiple eigenvalues|This process]] can be simplified into a single kind of iteration: the orthogonal iteration.
$$
Z = A Q_k, \qquad
Q_{k+1} R_{k+1} = Z
$$
Column $i+1$ of $Q_{k+1}$ converges to column $i+1$ of $P_{ \{q_1,q_2,\dots, q_i\}^\perp } A Q_k$.

In fact, $Q_k \to Q$.

Note: in the [[Schur decomposition]], $Q$ is not unique:

$Q T Q^H$

We can multiply $Q$ by any diagonal matrix $Q e^{i \Theta}$ with $\Theta$ a diagonal real matrix. This does not change the diagonal of $T$ or its upper triangular structure. It does change the strict upper diagonal entries, though.

So $Q_k \to Q$ up to some diagonal matrix $e^{i \Theta_k}$. Or more simply
$$
\text{span}{\{[Q_k]_{,j}\}} \to \text{span}{\{q_{,j}\}}
$$
The [[Angle between subspaces|angle]] between these subspaces goes to 0.