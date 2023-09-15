### Vector subspace
Vector space: a set of elements that can be added together and multiplied by scalars (for this class, in $\mathbb R$ or $\mathbb C$). [[Vectors and matrices]]

Example: $\mathbb R^n$.

A subspace of $\mathbb R^n$ is a vector space that is a subset of $\mathbb R^n$.

Example: take $k$ vectors $x_1$, ..., $x_k$, you can check that all vectors of the form
$$
a_1 x_1 + \cdots + a_k x_k
$$
form a subspace. We will denote this subspace as $S = \text{span}\{ x_1, \dots, x_k \}$.

### Linear independence
We say that $x_1$, ..., $x_k$ are linearly independent if and only if any vector $x \in S$ has a unique decomposition as
$$
x = a_1 x_1 + \cdots + a_k x_k
$$
In particular, if 
$$
0 = a_1 x_1 + \cdots + a_k x_k
$$
then $a_1 = \cdots = a_k = 0$.

If the vectors are linearly independent, $k$ is the dimension of $S$.

Linear independence will be important when solving linear systems. It will guarantee that the solution is unique.

### Direct sum
If $U$ and $V$ are subspaces, then $U+V$ is a subspace. We say that $W=U \oplus V$ is the direct sum of $U$ and $V$ if $U \cap V = \{0\}$. The direct sum means that if a vector is decomposed into its $U$ and $V$ components, this decomposition is unique.

Example: verify that if $x_1$, ..., $x_k$ are linearly independent, then
$$
S = \text{span}\{x_1\} \oplus \cdots \oplus \text{span}\{x_k\}
$$
