We [[Three-term recurrence|know]] how to compute $p^{(k+1)}$. We can now calculate the optimal step sizes $\mu_{k+1}$ using the [[All the orthogonality relations in CG|orthogonality relations.]]

[[CG search directions|Recall]] the basic definitions:
$$
\begin{gather}
\Delta x^{(k)} = \mu_{k+1} p^{(k+1)}\\[.3em]
x^{(k+1)} - x^{(k)} = \mu_{k+1} p^{(k+1)}
\end{gather}
$$
Multiply by $A$
$$
A x^{(k+1)} - A x^{(k)} = \mu_{k+1} A p^{(k+1)}
$$
We use the definition of the [[Residuals and solution increments in CG|residual:]] $r^{(k)} = b - A x^{(k)}$. So we get:
$$
-r^{(k+1)} + r^{(k)} = \mu_{k+1} A p^{(k+1)}
$$
[[All the orthogonality relations in CG|Recall]] that $r^{(k+1)} \perp \mathcal p^{(k+1)}$. Let's multiply by $(p^{(k+1)})^T$ to the left:
$$
\begin{gather}
0 + (p^{(k+1)})^T r^{(k)} = \mu_{k+1} (p^{(k+1)})^T A p^{(k+1)}\\[.3em]
\mu_{k+1} = \frac{(p^{(k+1)})^T r^{(k)}}{(p^{(k+1)})^T A p^{(k+1)}}
\end{gather}
$$
We can simplify it a bit more using our [[Three-term recurrence|three-term]] recurrence:
$$
p^{(k+1)} = r^{(k)} + \tau_k \, p^{(k)}, \quad \text{and} \quad
r^{(k)} \perp p^{(k)}.
$$
Take a dot product with $r^{(k)}$:
$$
[r^{(k)}]^T p^{(k+1)} = \|r^{(k)}\|_2^2 + \tau_k \, 0
$$
We have proved that:

**Theorem**. The optimal step-size in CG is given by:
$$
\begin{gather}
\Delta x^{(k)} = x^{(k+1)} - x^{(k)} = \mu_{k+1} \, p^{(k+1)}\\[.3em]
\mu_{k+1} = \frac{\|r^{(k)}\|_2^2}{(p^{(k+1)})^T A p^{(k+1)}}
\end{gather}
$$