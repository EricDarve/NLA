A simple formula can be used to estimate the error resulting from floating point arithmetic.

**Fundamental rule:** a floating-point operation must approximate the corresponding real number arithmetic operation by rounding any result that is not a floating-point number to the nearest floating-point number.

In short: `a fl(op) b = fl(a op b)`, where `op` = `+,*,-,/`.

Therefore: a fl(op) b = a op b + $\epsilon$ (a op b), where $|\epsilon| \le u$, and:
$$
u = \frac{1}{2} \; \times \; \text{(distance between 1 and the next largest floating point)}
$$

- $u$ is called the **unit roundoff.**
- $u \approx 10^{-7.2}$ in single precision, and $u \approx 10^{-15.9}$ in double precision.

[[Floating point arithmetic]], [[Floating point arithmetic is different from regular arithmetic]]