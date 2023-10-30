- Let's assume $A$ is [[Upper Hessenberg form for the QR iteration|upper Hessenberg]]. 
- Let us choose an approximate [[Eigenvalues|eigenvalue]] $\mu$. 
- For example, $\mu = [T_k]_{nn}$. 
- This should approximate the smallest [[Eigenvalues|eigenvalue]] of $A$. 

Then:
$$
|\lambda_n - \mu| \ll |\lambda_{n-1} - \mu|
$$
We had mentioned this idea [[Accelerating convergence using a shift|previously.]]

For our [[Convergence of the orthogonal iteration|previous theoretical result,]] we have that 
$$
\text{span}(Q_k[:,1:n-1]) \to \text{span}(Q[:,1:n-1])
$$
with rate 
$$
\left| \frac{\lambda_{n} - \mu}{\lambda_{n-1} - \mu} \right|^k
$$
But since the matrices are orthogonal, we also get convergence of the last column, which is orthogonal to the subspace spanned by the previous $n-1$ columns. So
$$
\text{span}(Q_k[:,n]) \to \text{span}(Q[:,n])
$$
Consider [[QR iteration|again:]]
$$
T_k = Q_k^H A Q_k, \; \text{and} \;
T = Q^H A Q.
$$
The last row of $T_k$ converges to the last row of $T$: $[0, \dots, 0, \lambda_{n}]$.

![[QR iteration with shift 2023-10-30 11.49.00.excalidraw]]

We can now use [[Deflation in the QR iteration|deflation]], and work with a smaller $n-1 \times n-1$ matrix.

By repeating this process, we make $A$ smaller and smaller and on the way, we obtain all the eigenvalues of $A$.