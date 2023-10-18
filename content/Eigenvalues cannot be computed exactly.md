- Computing eigenvalues is as hard as computing the roots of polynomials.
- There is no general method to exactly compute the roots of polynomials of degree 5 or higher. 
- This theorem is known as the Abel–Ruffini theorem. To be precise, this theorem posits that there is no solution in radicals for general polynomial equations of degree five or higher with arbitrary coefficients​. In other words, there isn't a formula akin to the quadratic formula for these polynomials, meaning a formula expressed in terms of radicals (like square roots and cube roots). When it comes to solving these higher-degree polynomials, root finding consists of finding approximations of the roots, rather than exact solutions in radicals​.

We can prove that computing eigenvalues is as difficult as computing the roots of polynomials. These two problems are equivalent.

This is the proof. Define the following polynomial:
$$
p(x) = x^n + a_{n-1} x^{n-1} + \cdots + a_0
$$
Finding the roots of this polynomial is equivalent to the following eigenvalue problem:
$$
A = 
\begin{pmatrix}
    0 & 1 \\
    & 0 & 1 \\
    & & 0 \\
    & \vdots \\
    &&&& 0 & 1 \\
    -a_0 & -a_1 & -a_2 & \cdots & -a_{n-2} & -a_{n-1}
\end{pmatrix}
$$
Indeed consider $Au = z \, u$ , $z \in \mathbb C$. Then, from the form of the matrix, we must have:
$$
u = (1,z,z^2,\ldots,z^{n-1})
$$
You can check that
$$
zu = (z,z^2,\ldots,z^n)
$$
This shows that
$$
p(z) = 0 \Leftrightarrow Au = z \, u
$$

- If we had an exact algorithm for eigenvalues, we would have an algorithm for roots.
- Therefore, such an algorithm cannot exist.
- By necessity, the methods are approximate.
- But they converge very fast and the error is close to machine accuracy.
- So they are indistinguishable from a direct “exact” method.

[[Eigenvalues]], [[Computing eigenvalues]]