A matrix $Q$ is orthogonal if
$$
Q^T Q = I
$$
This means that the columns of $Q$ are orthogonal to each other.

Unitary: same for complex matrices: $Q^H Q = I$.

Orthogonal matrices are very important in numerical linear algebra. They often guarantee that an algorithm will be accurate.

Orthogonal matrices are isometries, that is, they conserve the 2-norm of a vector:
$$
\|Qx\|_2 = \|x\|_2
$$
If $Q$ is square, then its inverse is $Q^T$. In that case, $Q$ can be decomposed into a sequence of rotations and reflections.

_Fun fact:_ Cartan–Dieudonné theorem. Every orthogonal transformation in $\mathbb C^n$ can be described as the composition of at most $n$ reflections.

The matrix $QQ^T$ is not equal to $I$ if $Q$ is rectangular (i.e., $m > n$). 

The matrix $P=QQ^T$ represents an orthogonal projection onto $R(Q)$ along $N(Q^T)$.

![[Pasted image 20230914170201.png]]

[[Dot product]], [[Vector norms]], [[Operator and matrix norms]], [[Projection]]