Consider the simplified 2x2 [[QR iteration with shift|example:]]
$$
A = \begin{pmatrix}
    w & x \\ \varepsilon & z
\end{pmatrix}
$$
We perform a [[QR iteration with shift|shift]] using $z$. After rescaling, we can reduce the problem to
$$
A = \begin{pmatrix}
    1 & x \\ \varepsilon & 0
\end{pmatrix}
$$
assuming that $w \neq z$.

Let's perform one step of the [[QR iteration]]. QR factorization + RQ multiplication:
$$
\begin{pmatrix}
\frac{(1 + \varepsilon^2 + x \varepsilon)}{1 + \varepsilon^2} & \frac{x - \varepsilon (1 + \varepsilon^2)}{1 + \varepsilon^2} \\[6pt]
- \frac{x\varepsilon^2}{1 + \varepsilon^2} & - \frac{x \varepsilon}{1 + \varepsilon^2}
\end{pmatrix}
$$
If we assume $|\varepsilon| \ll 1$, then
$$
\Big\vert - \frac{x\varepsilon^2}{1 + \varepsilon^2} \Big\vert
\sim \vert x \varepsilon^2 \vert
$$

**[[Symmetric and unsymmetric QR iteration|Unsymmetric]] case**

- (2,1) entry: $\varepsilon \to O(\varepsilon^2)$.
- **Quadratic convergence**

**[[Symmetric and unsymmetric QR iteration|Symmetric]] case**

$$
\begin{gather}
x = \varepsilon \\
A = \begin{pmatrix}
    1 & \varepsilon \\ \varepsilon & 0
\end{pmatrix} \\
\Big\vert - \frac{x\varepsilon^2}{1 + \varepsilon^2} \Big\vert
\sim \vert \varepsilon^3 \vert \\
\varepsilon \to O(\varepsilon^3)
\end{gather}
$$

We get a **cubic convergence.**