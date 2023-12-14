We derive all the [[Some orthogonality relations in CG|orthogonality relations]] in CG that you need to know. This is the foundation to build the final steps of CG.

[[CG search directions|Recall]] our definition of the search directions and residual vectors:
$$
\Delta x^{(k)} = \mu_{k+1} p^{(k+1)}, \qquad
r^{(k)} = b - A x^{(k)}.
$$
From these definitions, we can derive all the results below. Make sure you understand the reasoning behind each result.

**Property 1.**
$$
x^{(k)} \in \mathcal K_k.
$$
This is because of the basic construction of iterative methods using the [[Krylov methods for sparse systems|Krylov subspace.]]

**Property 2.**
$$
\Delta x^{(k)} \in \mathcal K_{k+1}
$$ 
and
$$
p^{(k+1)} \in \mathcal K_{k+1}.
$$
This is true because of the definition of the [[CG search directions|search directions]] above.

**Property 3.**
$$
r^{(k)} \in \mathcal K_{k+1}.
$$
This is true because of the definition of the [[Residuals and solution increments in CG|residual vectors]] and the [[Krylov subspace]].

**Property 4.** If $y \in K_k$, then $Ay \in \mathcal K_{k+1}.$

This is true because of the definition of the [[Krylov subspace]].

Note that these orthogonality relations are not a “new thing.” It’s a direct consequence of the [[Conjugate Gradients Version 1|basic starting point]] of CG
$$
\min_{y} \| x - Q_k y \|_A
$$
and the [[Krylov subspace]]. We have seen [[Some orthogonality relations in CG#^9bf4bf|previously]] that this implies that $r^{(k)}$ is $\perp$ to $\mathcal K_k.$

**Property 5.**
$$
r^{(k)} \perp \mathcal K_k.
$$
See the [[Some orthogonality relations in CG#^9bf4bf|previous section.]]

**Property 6.**
$$
r^{(k)} \perp r^{(l)}, \quad k \neq l.
$$
This is a consequence of the results listed above (Properties 3 and 5).

**Property 7.**
$$
r^{(k)} \perp \mathcal p^{(l)}, \quad l \le k.
$$
This is a consequence of Properties 2 and 5.

**Property 8.** $p^{(k)} \perp A p^{(l)},$ $k \neq l.$ So $p^{(k)}$ and $p^{(l)}$ are [[Some orthogonality relations in CG|conjugate.]]

This is a [[Some orthogonality relations in CG|consequence of the conjugacy]] of the $\Delta x^{(k)}$ and $\Delta x^{(k)}$, along with 
$$
\Delta x^{(k)} = \mu_{k+1} \: p^{(k+1)}.
$$

**Property 9.** $r^{(k)} \perp A p^{(l)},$ $l \le k-1;$ $r^{(k)}$ and $p^{(l)}$ are conjugate.

**Proof.** Since $r^{(k)} \perp Q_k$, we have $r^{(k)} \perp A Q_{k-1}.$ Since $p^{(l)} \in \mathcal K_l,$
$$
r^{(k)} \perp Ap^{(l)}, \quad l \le k-1.
$$
$\square$