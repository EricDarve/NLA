It is defined as the forward error divided by backward error:
$$
\frac{\|f(x) - \tilde{f}(x)\|}{\|x - \tilde{x}\|}
= \frac{\|f(x) - f(\tilde{x})\|}{\|x - \tilde{x}\|} = {\rm sensitivity}
$$

![[Pasted image 20231001143824.png]]

**Conditioning = relative sensitivity.** This quantity is dimensionless.
$$
\frac{\|f(x) - f(\tilde{x})\|}{\|f(x)\|} \,
\frac{\| x \|}{\|x - \tilde{x}\|}
$$
The conditioning can be defined for any algorithm of map $f$.

[[Forward and backward error]], [[Stability of the LU factorization]]