### Measuring the size of a matrix

This is the main norm we will use. It can be derived from the [[Vector norms|vector norms]] and extended to matrices
$$
\| A \|_p = \sup_{x \neq 0} \frac{\|A x\|_p}{\|x\|_p} = \max_{\|x\|_p = 1} \| Ax \|_p
$$
Example: $A = \begin{pmatrix} 2 & 1 \\ 0 & 1 \end{pmatrix}$

2-norm / 1-norm / $\infty$-norm:
![[Pasted image 20230914160537.png]]

Other norms that are useful.

Frobenius norm: it is very convenient to analyze matrices that can be decomposed into blocks
$$
\|A\|_F = \Big( \sum_{ij} a_{ij}^2 \Big)^{1/2}
$$

### Useful inequalities for proofs
$$
\|A x \|_p \le \|A\|_p \|x\|_p
$$
This is a consequence of the definition of the norm. Exercise: prove this result using the definition.

For Frobenius and $p$-norms:
$$
\| A B \| \le \|A\| \; \|B\|
$$
This is called a sub-multiplicative norm.

[[Vectors and matrices]], [[Vector norms]], [[Matrix-vector and matrix-matrix product]]