- Key tool for proofs and formal derivations
- Essential to describe and understand the properties of objects in LA
- Used to convert a vector to a single scalar number, its "size."
- The most common norm is the 2-norm:
$$
\| x \|_2 = \Big( \sum_{i=1}^n (x_i)^2 \Big)^{1/2}
= \sqrt{x^T x}
$$

[[Dot product]] and the 2-norm are connected through:
$$
x^T y = \| x \|_2 \, \| y \|_2 \cos \theta
$$
So the dot product can be used to measure the angle between two vectors.

If two vectors are orthogonal, their dot product is equal to 0.

A unit vector is a vector with norm 1: $\|x\|_2 = 1$.

In data science and linear algebra, it's common to use different norms. They differ by the weight they assign to the components of a vector.

- 1-norm: all scales contribute equally
$$
\| x \|_1 = \sum_{i=1}^n |x_i|
$$
- 2-norm: emphasizes the larger entries
- Largest entry:
$$
  \| x \|_\infty = \max_{1 \le i \le n} |x_i|
$$
- Intermediate; as $p \to \infty$, the norm becomes closer to the max-norm
$$
  \| x \|_p = \Big( \sum_{i=1}^n |x_i|^p \Big)^{1/p}
$$

**Balls in different norms**
![[20220923_155344_image.png]]
[[Vectors and matrices]]