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
- [[LU and determinant]]
	- The determinant seems very complicated to calculate.
	- But using LU, we can get the determinant very easily.
	- Uses the fact that $\det U = \prod_{i=1}^n u_{ii}$ for triangular matrices.
- [[Existence of LU]]
	- The LU factorization suffers from stability issues.
	- This can lead to inaccurate solutions in some cases, or the algorithm may even break down.
	- Under what condition does the LU factorization exist? What are the situations where the algorithm breaks down?
- [[Stability of the LU factorization]]
	- The existence result describes what happens when a pivot is 0.
	- What about a very small pivot? What can we expect in that case?
- [[Floating point numbers]]
	- A consequence of storing numbers and executing algorithms on computers
- [[Floating point arithmetic and unit roundoff error]]
	- How to model and estimate roundoff errors
	- This is important to estimate errors and provide bounds on the accuracy of calculations.
	- Small errors are magnified by unstable algorithms and lead to wrong answers
- [[Floating point arithmetic is different from regular arithmetic]]
	- Understanding this difference is important to understand how large errors can creep into a calculation
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
	- This is consistent with our previous observation regarding the [[Stability of the LU factorization|stability of the LU factorization]].
- [[Row pivoting]]
	- The simplest and most efficient method to make the LU factorization backward stable.
	- This is the most common implementation of the LU factorization.
- [[Symmetric Positive Definite Matrices]]
	- Many algorithms (such as [[Cholesky factorization|Cholesky]]) require that the matrix is SPD.
	- These matrices satisfy very strong properties and as a result, very fast and accurate algorithms exist for these matrices.
	- This is a very important class of matrices in NLA.
- [[Cholesky factorization]]
	- The Cholesky factorization applies to any SPD matrix.
	- This is a triangular factorization of the form $A=LL^T$ where $L$ is lower triangular.
	- This factorization is faster and requires less memory than LU.
- [[Existence of the Cholesky factorization]]
	- We prove that the Cholesky factorization exists and is unique.
	- [[Row pivoting|Pivoting]] is not required.
- [[Stability of the Cholesky factorization]]
	- The algorithm is always stable even without any pivoting.