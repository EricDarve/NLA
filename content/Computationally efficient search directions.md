We have derived all the important relations in CG. We need one more step to make computing the search directions a little bit computationally more efficient.

Recall the basic [[Three-term recurrence|three-term]] recurrence for the search directions:
$$
\begin{gather}
p^{(k+1)} = r^{(k)} + \tau_k p^{(k)}\\[.3em]
\tau_k = - \frac{ (p^{(k)})^T A r^{(k)} }{ (p^{(k)})^T A p^{(k)} }\\[.5em]
\tau_k = 
-\frac{(p^{(k)})^T A r^{(k)}}{(p^{(k)})^T A p^{(k)}}
=
-\frac{(Ap^{(k)})^T r^{(k)}}{(p^{(k)})^T A p^{(k)}}
\end{gather}
$$

Let's substitute $Ap^{(k)}$. From our [[CG search directions|previous]] derivations:
$$
\begin{gather}
x^{(k)} - x^{(k-1)} = \mu_k p^{(k)}\\[.3em]
A x^{(k)} - A x^{(k-1)} = \mu_k A p^{(k)}\\[.3em]
-r^{(k)} + r^{(k-1)} = \mu_k A p^{(k)}\\[.3em]
-A p^{(k)} = \frac{r^{(k)} - r^{(k-1)}}{\mu_k}
\end{gather}
$$
Substitute in $\tau_k$:
$$
\tau_k = 
\frac{(r^{(k)} - r^{(k-1)})^T r^{(k)}}{\mu_k (p^{(k)})^T A p^{(k)}}
$$
[[All the orthogonality relations in CG|Recall]] that $r^{(k)} \perp r^{(k-1)}:$
$$
\tau_k = 
\frac{\| r^{(k)} \|_2^2}{\mu_k \: (p^{(k)})^T A p^{(k)}}
$$
[[Optimal step size|Use the equation]] for $\mu_k$:
$$
\mu_k = \frac{\|r^{(k-1)}\|_2^2}{(p^{(k)})^T A p^{(k)}}
$$
So:
$$
\mu_k \: (p^{(k)})^T A p^{(k)} = \|r^{(k-1)}\|_2^2
$$
Substitute back in $\tau_k$:
$$
\tau_k = 
\frac{\| r^{(k)} \|_2^2}{\|r^{(k-1)}\|_2^2}
$$
We have proved that:

**Theorem.** The optimal search directions in CG are given by:
$$
\begin{gather}
p^{(k+1)} = r^{(k)} + \tau_k \, p^{(k)}\\[.5em]
\tau_k = 
\frac{\| r^{(k)} \|_2^2}{\|r^{(k-1)}\|_2^2}
\end{gather}
$$