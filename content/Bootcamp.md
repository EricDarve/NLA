- [[Vectors and matrices]]
	- Introduces notations
- [[Matrix block operations]]	
    - Block operations will be very useful in many proofs and algorithms.
- [[Subspace and linear independence]]
	- Core concept in linear algebra
	- Important when solving linear systems
- [[Dot product]]
	- Shows up in many places.
	- Example: vector norms, matrix-vector and matrix-matrix products.
	- Used to define orthogonality
- [[Vector norms]]
	- How to measure things
	- Key to calculating errors in numerical methods
- [[Projection]]
	- Using the dot product and norm for projection
- [[Pythagorean theorem]]
	- How to simply calculate the length of a vector given its decomposition into orthogonal subspaces
	- This is key to computing the norm of a vector in certain situations.
- [[Cauchy-Schwarz]]
	- Key in proofs to derive upper bounds on error
	- Considered one of the most important and widely used inequalities in mathematics
- [[Matrix-vector and matrix-matrix product]]
	- Can be either viewed algebraically or interpreted as an operator
	- Matrix-vector: applying a linear operator to transform a vector
	- Matrix-matrix: operator composition
- [[Invertible matrix]]
	- A requirement to solve a linear system and obtain a unique solution
- [[Sherman-Morrison-Woodbury formula]]
	- How to solve a linear system when we make a small perturbation
- [[Operator and matrix norms]]
	- Measuring the size of operators
	- Key when deriving error bounds and for proof.
- [[The four fundamental spaces]]
	- Understanding the structure of linear operators
	- How they transform the input vector and map subspaces
- [[Orthogonal matrix and projector]]
	- Orthogonal matrices will be key because they act as isometries
	- Useful for building algorithms with low error
	- Key in many matrix decompositions or factorizations
- [[Eigenvalues]]
	- Key to analyzing the powers of a matrix, $A^k$ , and solving differential equations $\dot{x} = Ax$.
	- This is important for time evolution and repeated applications of an operator.
	- Long-term evolution of a dynamical system.
	- Not useful to understand what happens when applying the operator once.
- [[Diagonalizable matrices]]
	- A special matrix with a basis of eigenvectors.
- [[Determinant]]
	- How a matrix changes the volume of a subspace.
- [[Trace]]
	- Connect matrix with vector field
	- Trace = divergence of vector field = flux through a unit square
- [[Unitarily diagonalizable matrices]]
	- The simplest and most accurate case
	- Diagonalizable + orthogonal matrices!
- [[Hermitian and symmetric matrices]]
	- An important special case of [[Unitarily diagonalizable matrices|unitarily diagonalizable matrices]].
- [[Schur decomposition]]
	- This will be key to computing eigenvalues.
	- The fact that it uses orthogonal matrices will be important to ensure the accuracy of the algorithm.
	- **Exists for all square matrices.**
- [[Singular value decomposition]]
	- Key to understanding how a matrix scales and transforms space
	- A linear operator always transforms the unit ball to an $n$-ellipsoid.
	- Connection with the [[The four fundamental spaces|four fundamental spaces]].
	- Connection with [[Eigenvalues|eigenvalues]] and [[Determinant|determinant]] of a matrix.
	- Key when solving least-squares problems: $\min_x \|Ax-b\|_2$.
	- Relation to [[Operator and matrix norms|operator and matrix norms]].