We have now derived all the equations we need to formulate the Conjugate Gradients algorithm. Each step of CG is described below.

Compute the [[Optimal step size|step-size:]]

$$\mu_k = \frac{ \| r^{(k-1)} \|_2^2 }{(p^{(k)})^T A p^{(k)}}$$

The [[Some orthogonality relations in CG|solution increments]] are conjugate with all the previous steps:
$$
x^{(k)} = x^{(k-1)} + \mu_k \, p^{(k)}
$$
The residual is [[All the orthogonality relations in CG|orthogonal]] to all the previous residuals:
$$
r^{(k)} = r^{(k-1)} - \mu_k A p^{(k)}
$$
We use our [[Computationally efficient search directions|computationally efficient]] equation for the [[CG search directions|search directions:]]
$$
\tau_k = \frac{\| r^{(k)} \|_2^2}{\| r^{(k-1)} \|_2^2}
$$
The [[Three-term recurrence|search directions]] are conjugate with each other:
$$
p^{(k+1)} = r^{(k)} + \tau_k p^{(k)}
$$
Finally, we increase the iteration index:
$$
k \leftarrow k+1
$$