It is defined as the forward error divided by the backward error:
$$
\frac{\|f(x) - \tilde{f}(x)\|}{\|x - \tilde{x}\|}
= \frac{\|f(x) - f(\tilde{x})\|}{\|x - \tilde{x}\|} = {\rm sensitivity}
$$

$f$: exact map

$\tilde{f}$: numerical solution

$\tilde{x}$: defined as the input such that $\tilde{f}(x) = f(\tilde{x})$.

![[Pasted image 20231001143824.png]]

**Conditioning = relative sensitivity.** 

This quantity is dimensionless.
$$
\frac{\|f(x) - \tilde{f}(x)\|}{\|f(x)\|} \,
\frac{\| x \|}{\|x - \tilde{x}\|}
$$
The conditioning can be defined for any algorithm or map $f$.

[[Forward and backward error]], [[Stability of the LU factorization]]