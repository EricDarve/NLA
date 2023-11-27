- [[Krylov methods for sparse systems]]
	- This provides a basic introduction to all the methods covered in this chapter.
	- As [[Brief introduction to Conjugate Gradients|alluded to previously,]] we search for the solution in the Krylov subspace that satisfies some appropriate optimality condition.
	- The optimality condition depends on the properties of the matrix: [[Symmetric Positive Definite Matrices|symmetric positive definite,]] symmetric, and general matrix.
- [[Conjugate Gradients Version 1]]
	- This is a short and basic derivation of the algorithm.
	- This method produces the correct sequence of iterates $x_k$.
	- However, its space (memory) and time (number of flops) costs are not optimal.
	- This will be refined in the sections below.
- [[Some orthogonality relations in CG]]
	- We illustrate the first of several orthogonality relations satisfied by the CG vector sequences.
	- The orthogonality relation is special. It is defined using a scaling by $A^{1/2}$.
	- In this section, we prove that the solution increments $\Delta x^{(k)}$ are $A$-conjugate to each other (i.e., orthogonal with $A^{1/2}$ scaling).
- [[Convergence with conjugate steps]]
	- Based on the derivation from the [[Some orthogonality relations in CG|previous section]], we can prove that CV converges in at most $n$ steps, where $n$ is the size of $A$.
- [[Residuals and solution increments in CG]]
	- We now start deriving a series of results to help us turn the [[Conjugate Gradients Version 1|current version]] of CG into an efficient algorithm.
	- We prove that the subspace spanned by the residuals is the same as the subspace spanned by the solution increments.
- [[CG search directions]]
	- Using the [[Residuals and solution increments in CG|previous section,]] we establish that the residuals can be written as linear combinations of the increments.
	- We will see later on how the coefficients in the linear combination can be easily computed using orthogonality relations.
- [[All the orthogonality relations in CG]]
	- We derive all the orthogonality relations in CG that you need to know. 
	- This is the foundation to build the final steps of CG in the sections below.
- [[Three-term recurrence]]
	- Using the [[CG search directions|previous result]] and the [[All the orthogonality relations in CG|orthogonality relations,]] we derive a three-term recurrence for the search directions $p^{(l)}$.
- [[Optimal step size]]
	- Using our [[All the orthogonality relations in CG|orthogonality relations,]] we derive a computationally efficient formula for the step sizes $\mu_{k+1}$.
- [[Computationally efficient search directions]]
	- We need one more step to make computing the [[CG search directions|search directions]] $p^{(k)}$ efficient.
- [[Conjugate Gradients algorithm]]
	- Based on all the previous equations, we describe each step in the CG algorithm.
- [[Conjugate Gradients code]]
	- We give the Julia code for CG based on the [[Conjugate Gradients algorithm|algorithm]] we derived.
- [[Convergence of the Conjugate Gradients]]
	- We give an estimate for the convergence of CG.
	- In some sense, this convergence can be considered to be optimal.