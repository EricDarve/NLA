The QR decomposition is the decomposition of $A$ into $A = QR$ where $Q$ is orthogonal and $R$ is upper triangular.

We will assume that $r_{ii}$ is non-negative.

Why QR?

1. $Q$ is orthogonal and conserves the 2-norm.
2. $R$ is triangular and easy to solve.

Example of application:
- Least-squares problem
- Find $x^*$ such that $x^* = \min_x \| Ax - b \|_2$. 
- We can use $A = QR$ for this. 
- The solution satisfies $Rx^* = Q^T b$.