It is important to understand how computers make small errors during a computation. As a result, some algorithms may yield completely incorrect errors because these errors may be allowed to grow out of control. This happens for **unstable** algorithms.

What is a floating point number?

Floating point numbers with base 2 have the form 
$$
\pm (1 + \sum_{i=1}^{p-1} d_i 2^{-i} )\; 2^e,
$$
where $e$ is called the exponent, $p$ is the precision, and 
$$
1 + \sum_{i=1}^{p-1} d_i 2^{-i}
$$
is the significand.

Example of floating-point number: 3.140625.

It is positive, so the sign is $+$.

It's between $2$ (aka, $2^1$) and $4$ (aka, $2^2$). So, the exponent $e$ is equal to $1$.

The significand is 1.5703125. Decomposition of the significand:
$$
1.5703125 = \frac{1}{2} + \frac{1}{16} + \frac{1}{128} = 2^{-1} + 2^{-4} + 2^{-7}
$$

Bits for significand: $d_1 = 1$, $d_4 = 1$, $d_7 = 1$.

[[Stability of the LU factorization]]