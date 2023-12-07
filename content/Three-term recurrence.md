We derive an efficient [[Lanczos process|three-term relation]] for the search directions $p^{(l)}$.

Using the [[All the orthogonality relations in CG|orthogonality conditions,]] we can derive an efficient three-term recurrent relation.

[[CG search directions|Recall that]]
$$
r^{(k)} = p^{(k+1)} + \sum_{l=1}^k \alpha_l \, p^{(l)}
$$

The search directions $p^{(l)}$ are [[All the orthogonality relations in CG|conjugate.]] So this is, in fact, a type of QR decomposition of the matrix
$$
[r^{(0)}, \dots, r^{(k)}]
$$
using $A$-orthogonality. Note that the search directions $p^{(l)}$ do not have norm 1. They are just conjugate.

Recall the algorithms to compute the [[Gram-Schmidt|QR factorization]] $A = QR$ and the [[Lanczos process|Lanczos process,]] $AQ = QT$. In both cases, we found the columns of $Q$ using a recurrence and orthogonality. The same ideas can be applied here. Consider
$$
r^{(k)} = p^{(k+1)} + \sum_{l=1}^k \alpha_l \, p^{(l)}
$$
and multiply to the left by $[p^{(r)}]^T A$, we get
$$
[p^{(r)}]^T A r^{(k)} = [p^{(r)}]^T A p^{(k+1)} + \sum_{l=1}^k \alpha_l \, [p^{(r)}]^T A p^{(l)}
$$
From the [[All the orthogonality relations in CG|orthogonality relations,]] we have that for $r < k$:
$$
[p^{(r)}]^T A r^{(k)} = 0
$$ 
and 
$$
[p^{(r)}]^T A p^{(l)} = 0, \quad r \neq l
$$
So for $r < k$, we get
$$
0 = \alpha_r \, [p^{(r)}]^T A p^{(r)}
$$
and $[p^{(r)}]^T A p^{(r)} \neq 0$ since $A$ is SPD. So $\alpha_r = 0$, $r < k$. We are only left with the term $l=k$:
$$
r^{(k)} = p^{(k+1)} + \alpha_k p^{(k)}
$$
There is only one term left. Let us denote for convenience:
$$
r^{(k)} = p^{(k+1)} - \tau_k p^{(k)}
$$
Multiply by $A$ and take a dot product with $p^{(k)}$:
$$
(p^{(k)})^T A r^{(k)} = 0 - \tau_k (p^{(k)})^T A p^{(k)} 
$$

**Theorem.** The search directions in CG satisfy:
$$
\begin{gather}
p^{(k+1)} = r^{(k)} + \tau_k p^{(k)} \\[.3em]
\tau_k = - \frac{ (p^{(k)})^T A r^{(k)} }{ (p^{(k)})^T A p^{(k)} }
\end{gather}
$$
$\tau_k$ is the negative of the $A$-projection of $r^{(k)}$ onto $p^{(k)}$.

We find the new search direction by taking a vector in $K_{k+1}$ and making it conjugate to the previous $p^{(k)}$.

$p^{(k+1)}$ is a linear combination of the residual and the previous search direction.