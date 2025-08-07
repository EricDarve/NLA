We [[Residuals and solution increments in CG|previously]] showed that
$$
\text{span} \left\{ r^{(0)}, \ldots, r^{(k)} \right\} 
= \text{span} \left\{ \Delta x^{(0)}, \ldots, \Delta x^{(k)} \right\}
$$
So we can always write that:
$$
\begin{gather}
\Delta x^{(k)} = \mu_{k+1} \: p^{(k+1)} \\[.3em]
r^{(k)} = p^{(k+1)} + \sum_{l=1}^k \alpha_l \: p^{(l)}
\end{gather}
$$
We will later see how we can compute $\mu_{k+1}$ and $\alpha_l$ using the right orthogonality properties.