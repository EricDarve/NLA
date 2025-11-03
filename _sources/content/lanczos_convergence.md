# Convergence of the Lanczos Process

The benefits of moving from Arnoldi to Lanczos for symmetric matrices extend beyond computational efficiency. The convergence behavior of the Ritz values becomes far more structured and predictable. For non-symmetric matrices, the Ritz values from Arnoldi can appear anywhere in the complex plane and may not exhibit clear convergence patterns. For symmetric matrices, the Lanczos algorithm provides powerful guarantees, particularly for the extremal eigenvalues.

## Monotonic Convergence of Extremal Ritz Values

Let's consider the largest eigenvalue, $\lambda_1$, of a symmetric matrix $A$. It can be characterized by the Rayleigh quotient:

$$ \lambda_1 = \max_{\boldsymbol{x} \neq \boldsymbol{0}} \frac{\boldsymbol{x}^T A \boldsymbol{x}}{\boldsymbol{x}^T \boldsymbol{x}}
$$

The largest eigenvalue of the $k$-th Lanczos matrix $T_k$, which we denote $\lambda_1(T_k)$, is the largest Ritz value. By definition of the Ritz procedure, this value is the maximum of the Rayleigh quotient of $A$ over the Krylov subspace $\mathcal{K}_k$:

$$
\lambda_1(T_k) = \max_{\boldsymbol{z} \in \mathcal{K}_k, \boldsymbol{z} \neq \boldsymbol{0}} \frac{\boldsymbol{z}^T A \boldsymbol{z}}{\boldsymbol{z}^T \boldsymbol{z}}
$$

Since the maximization for $\lambda_1(T_k)$ is performed over a subspace of the full space $\mathbb{R}^n$, it is an immediate consequence that the largest Ritz value can never exceed the true largest eigenvalue:

$$
\lambda_1(T_k) \le \lambda_1
$$

Furthermore, since the Krylov subspaces are nested ($\mathcal{K}_k \subset \mathcal{K}_{k+1}$), the space of available test vectors grows with each iteration. This implies that the maximum of the Rayleigh quotient can only increase (or stay the same):

$$
\lambda_1(T_k) \le \lambda_1(T_{k+1})
$$

This establishes a fundamental property of the Lanczos algorithm: the approximations to the largest eigenvalue, $\lambda_1(T_k)$, form a **monotonically increasing sequence that is bounded above by $\lambda_1$**. Likewise, the approximations to the smallest eigenvalue, $\lambda_n(T_k)$, form a monotonically decreasing sequence bounded below by $\lambda_n$. This structured, monotonic convergence is a significant improvement over the often erratic behavior of Ritz values in the general Arnoldi case.

## A Quantitative Convergence Rate

We can achieve a much more precise, quantitative understanding of the convergence rate by again invoking the connection between Krylov subspaces and matrix polynomials. Recall that any vector in the Krylov subspace $\mathcal{K}_k(\boldsymbol{q}_1, A)$ can be written as $\boldsymbol{z} = p(A)\boldsymbol{q}_1$ for some polynomial $p$ of degree less than $k$. The largest Ritz value can therefore be expressed as an optimization over this class of polynomials:

$$
\lambda_1(T_k) = \max_{p \in \mathcal{P}_{k-1}} \frac{ (p(A)\boldsymbol{q}_1)^T A (p(A)\boldsymbol{q}_1) }{ (p(A)\boldsymbol{q}_1)^T (p(A)\boldsymbol{q}_1) } = \max_{p \in \mathcal{P}_{k-1}} \frac{ \boldsymbol{q}_1^T p(A)^T A p(A) \boldsymbol{q}_1 }{ \boldsymbol{q}_1^T p(A)^2 \boldsymbol{q}_1 }
$$

where $\mathcal{P}_{k-1}$ is the set of polynomials of degree at most $k-1$.

This expression reveals that the Lanczos algorithm implicitly finds the polynomial that maximizes this ratio. By making a clever choice of a specific polynomial—the **Chebyshev polynomial**—one can derive a powerful bound on the error of the largest Ritz value. The full derivation is detailed, but the celebrated result is as follows:

$$
\lambda_1 - \lambda_1(T_k) \le (\lambda_1 - \lambda_n) \left( \frac{\tan \phi_1}{T_{k-1}(1 + 2\rho_1)} \right)^2
$$

Let's dissect this important formula:

* $T_{k-1}(x)$ is the Chebyshev polynomial of the first kind of degree $k-1$.
* $\phi_1$ is the angle between the starting vector $\boldsymbol{q}_1$ and the true eigenvector $\boldsymbol{v}_1$ corresponding to $\lambda_1$. This term reminds us that convergence depends on the starting vector not being pathologically orthogonal to the eigenvector we seek.
* The crucial term is the **gap ratio**, $\rho_1$, defined by the eigenvalues of $A$:

    $$
    \rho_1 = \frac{\lambda_1 - \lambda_2}{\lambda_2 - \lambda_n}
    $$

This ratio quantifies how well-separated the largest eigenvalue $\lambda_1$ is from the rest of the spectrum, relative to the width of the rest of the spectrum.

The power of this bound comes from the behavior of Chebyshev polynomials. For an argument $x > 1$, $T_{k-1}(x)$ grows exponentially with $k$. In our case, the argument is $1 + 2\rho_1$, which is always greater than 1. The growth is approximately:

$$
T_{k-1}(1 + 2\rho_1) \approx \frac{1}{2} (4 \rho_1)^{k-1}
$$

when $\rho_1 \gg 1$.

Since this term appears squared in the denominator of the error bound, the error $\lambda_1 - \lambda_1(T_k)$ converges to zero at a remarkable exponential rate. The speed of this convergence is governed by the gap ratio $\rho_1$. If $\lambda_1$ is well-separated from the other eigenvalues, $\rho_1$ is large, and convergence is exceptionally fast.

## Convergence of Interior Eigenvalues

The remarkable convergence rate we established for the largest eigenvalue, $\lambda_1$, is a cornerstone of the Lanczos method. The same theory applies to the smallest eigenvalue, $\lambda_n$. By analyzing the convergence of Lanczos on the matrix $-A$, whose largest eigenvalue is $-\lambda_n$, we find that the smallest Ritz value, $\lambda_k(T_k)$, converges monotonically downward to $\lambda_n$ with a rate governed by the gap ratio at the low end of the spectrum, $\rho_n = (\lambda_{n-1} - \lambda_n) / (\lambda_1 - \lambda_{n-1})$.

This raises a natural question: what about the "inner" eigenvalues, $\lambda_i$ for $1 < i < n$?

The good news is that the Ritz values $\lambda_i(T_k)$ do converge to the true eigenvalues $\lambda_i$. The challenging news is that they do so more slowly, and the convergence is no longer guaranteed to be monotonic. The Lanczos algorithm, at its core, is a procedure that excels at approximating the function $f(x)=x$ with a polynomial on the spectrum of $A$. Polynomials are notoriously good at approximating functions near the ends of an interval, but less effective in the middle. This intuition is borne out by the theory.

A similar, though more complex, bound exists for the convergence of the $i$-th Ritz value, $\lambda_i(T_k)$. It is given by:

$$
0 \le \lambda_i - \lambda_i(T_k) \le (\lambda_i - \lambda_n) \left( \frac{\kappa_i \tan \phi}{T_{k-i}(1 + 2\rho_i)} \right)^2
$$

This formula reveals why the inner eigenvalues converge more slowly. Let's analyze its components in comparison to the bound for $\lambda_1$.

1. **Reduced Polynomial Degree:** The argument of the Chebyshev polynomial is now $T_{k-i}$, not $T_{k-1}$. This is the most critical change. It implies that the algorithm effectively uses $k-i$ iterations to find the $i$-th eigenvalue. Intuitively, the Lanczos process "spends" its first degrees of freedom pinning down the larger, outer eigenvalues before it can effectively resolve the inner ones. With a lower-degree polynomial, the exponential growth is significantly stunted, leading to slower convergence.

2. **The Gap Ratio $\rho_i$:** The gap ratio for an inner eigenvalue $\lambda_i$ is typically defined relative to its nearest neighbor and the far end of the spectrum:

    $$
    \rho_i = \frac{\lambda_i - \lambda_{i+1}}{\lambda_{i+1} - \lambda_n}
    $$
    
    For many matrices, the gaps between interior eigenvalues ($\lambda_i - \lambda_{i+1}$) are much smaller than the gaps at the spectral edges ($\lambda_1 - \lambda_2$). A smaller gap ratio $\rho_i$ means the argument to the Chebyshev polynomial is smaller, further reducing its growth and slowing convergence.

3. **The "Interference" Factor $\kappa_i$:** The new term $\kappa_i$ is a product that depends on the Ritz values that are larger than the one we are targeting:

    $$
    \kappa_i = \prod_{j=1}^{i-1} \frac{\lambda_j(T_k) - \lambda_n}{\lambda_j(T_k) - \lambda_i}
    $$
    
    As the outer Ritz values $\lambda_j(T_k)$ converge to their true counterparts $\lambda_j$, this term approaches a constant. However, it quantifies the "cost" of having to approximate the larger eigenvalues first. Each factor in the product is greater than 1, so $\kappa_i$ acts as a penalty term that grows with $i$, further increasing the error bound.

In summary, the Lanczos algorithm exhibits a distinct **convergence hierarchy**:

1. **Extremal eigenvalues ($\lambda_1, \lambda_n$) converge first and fastest.** Their convergence rate is determined by the relative gaps at the edges of the spectrum.
2. **Interior eigenvalues converge later and more slowly.** The convergence for $\lambda_i$ only truly begins after the larger eigenvalues $\lambda_1, \dots, \lambda_{i-1}$ are reasonably well-approximated.

This behavior is a fundamental trade-off. While Lanczos is exceptionally powerful for finding a few extremal eigenvalues of a large symmetric matrix, it is not the ideal tool if one needs to find all eigenvalues, or specifically target an eigenvalue deep inside the spectrum without modification. For such tasks, more advanced techniques like shift-and-invert Lanczos are required.

## A Stronger Theory for a Special Case

This brings us to the final, crucial comparison between Lanczos and Arnoldi.

* **For symmetric matrices (Lanczos)**, the eigenvalues lie on the real line. The convergence of extremal eigenvalues is monotonic and is governed by the *gaps* between eigenvalues. We have sharp, predictive, and powerful convergence bounds like the one above that tell us *how fast* we can expect to converge based on the intrinsic properties (the spectrum) of the matrix $A$.

* **For non-symmetric matrices (Arnoldi)**, the eigenvalues lie in the complex plane. There is no simple notion of "monotonic" convergence, and the Ritz values can wander unpredictably before settling. While convergence theory exists, it is far more complex and less predictive. The bounds often depend not just on the eigenvalues, but on the geometry of the *field of values* and the *pseudospectrum* of the matrix, which are much harder to analyze. The problem of finding the best polynomial approximation on a 2D domain in the complex plane is fundamentally harder than the 1D version for symmetric matrices.

In summary, the Lanczos algorithm is not just faster per-iteration than Arnoldi; it is backed by a much stronger and more predictive convergence theory that guarantees the rapid and monotonic convergence of the extremal Ritz values. This combination of speed and theoretical elegance makes it an indispensable tool for symmetric eigenvalue problems.