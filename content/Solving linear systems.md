How to solve $Ax = b$? One of the most important computational tasks in NLA.

Notations | Description
--- | ---
$A$ | Matrix
$n$ | Size of the matrix
$A[1: k, 1: k]$ | Top left $k \times k$ block of $A$
$a_{ij}$ | entry $(i,j)$
$a_{,j}$ | column $j$
$a_{i,}$ | row $i$
det | determinant

- This section covers the use of the LU factorization to solve linear systems like $Ax = b$. 
- This method is fast and nearly optimal in terms of floating point operations. 
- However, it suffers from stability and accuracy issues in some cases.

- [[Solving triangular systems]]
	- This is the starting point for efficient solution of linear systems.
	- Solving triangular systems is easy and computationally efficient.
- [[Solving linear systems using LU]]
	- How the LU factorization allows us to solve linear systems
- [[Outer form of matrix-matrix product]]
	- How can we get $L$ and $U$ starting from $A$?
	- The starting point is the outer form of the product $LU$.
- [[Triangular factorization]]
	- Solving triangular systems allows us to solve general systems of equations, provided the matrix is factorized as $A = LU$ using triangular factors.
	- Computing $L$ and $U$ can be done iteratively column by column.
	- This can be done using the outer form for the matrix-matrix product.
- [[LU algorithm]]
	- This is the step-by-step algorithm.
- [[Existence of LU]]
	- The LU factorization suffers from stability issues.
	- This can lead to inaccurate solutions in some cases, or the algorithm may even break down.
	- Under what condition does the LU factorization exist? What are the situations where the algorithm breaks down?
- [[LU and determinant]]
	- The determinant seems very complicated to calculate.
	- But using LU, we can get the determinant very easily.
	- Uses the fact that $\det U = \prod_{i=1}^n u_{ii}$ for triangular matrices.
- [[Stability of the LU factorization]]
	- The existence result describes what happens when a pivot is 0.
	- What about a very small pivot? What can we expect in that case?
- [[Floating point arithmetic]]
	- A consequence of executing algorithms on computers
	- Small errors can be magnified by unstable algorithms and lead to wrong answers
- [[Floating point arithmetic is different from regular arithmetic]]
	- Understanding this difference is important to understand how large errors can creep into a calculation
- [[Unit roundoff error]]
	- How to model and estimate roundoff errors
	- This is important to estimate errors and provide accuracy bounds on calculations
- [[Forward and backward error]]
	- These are the main concepts for analyzing numerical errors in algorithms.
	- This method can be used to prove that an algorithm is stable, that is, small perturbations in the input lead to small perturbations in the output.
	- The opposite is an unstable algorithm in which errors cannot be controlled.
- [[Sensitivity analysis]]
	- This is the concept that connects the forward and backward error estimates.
- [[Conditioning of a linear system]]
	- Application of the concept of sensitivity and conditioning to the problem of solving a linear system.
- [[Backward error analysis for LU]]
	- Apply previous concepts to the LU factorization algorithm
	- The current LU algorithm without pivoting is, in fact, backward **unstable.**
	- This is consistent with our previous observation in [[Stability of the LU factorization]].
- [[Row pivoting]]
	- The simplest and most efficient method to make the LU factorization backward stable.
	- This is the most common implementation of the LU factorization.