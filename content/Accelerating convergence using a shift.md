We are [[Convergence of the orthogonal iteration|now in a position]] to understand how convergence can be accelerated. Assume we use $A - \lambda I$ instead of $A$. Then the eigenvalues are shifted by $\lambda$: $\lambda_i - \lambda$. Geometrically, shifting is like temporarily moving the origin of our coordinate system to accelerate the convergence of specific eigenvalues.

This can be used to accelerate convergence.

For example, assume that $\lambda \approx \lambda_n$ and we shift. We get:
$$
|\lambda - \lambda_n| \ll |\lambda - \lambda_{n-1}|.
$$
The last eigenvalue will converge very rapidly. Shifting allows accelerating convergence. We will see that, if the shifting is done correctly, we get a quadratic convergence of the eigenvalue!

This method works well. 

Note that the focus is now on the [[Computing eigenvectors using the Schur decomposition|eigenvalues]] rather than $Q$. We need to efficiently compute $T_k$. It turns out that there is a simple algorithm to do that.