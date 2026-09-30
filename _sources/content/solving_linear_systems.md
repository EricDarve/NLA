# Solving Linear Systems

Solving a system of linear equations, written in matrix form as $Ax = b$, is a central problem in numerical linear algebra. Such systems arise throughout science and engineering, from simulating airflow over a wing and analyzing electrical circuits to training machine learning models and pricing financial derivatives.

If $A$ is an invertible square matrix, the unique solution is $x = A^{-1}b$. This formula suggests computing the inverse of $A$ and multiplying it by $b$. In practice, forming the inverse is generally more expensive and can be less accurate than solving the system directly. Our first principle is therefore: **avoid computing a matrix inverse explicitly to solve a linear system.**

The central strategy of this chapter is **factorization**: expressing $A$ as a product of simpler matrices. We will focus on the **LU factorization**, where we write $A = LU$, with $L$ lower triangular and $U$ upper triangular.

Triangular systems are inexpensive to solve. Once we have factored $A$, we can solve $Ax = b$ in two steps:

1.  Solve $Ly = b$ for $y$ (using *forward substitution*).
2.  Solve $Ux = y$ for $x$ (using *backward substitution*).

The same factors can be reused to solve systems with different right-hand sides.

To use this strategy reliably, we need to address three questions:

1.  **Existence and uniqueness:** Can every invertible matrix $A$ be factored as $A = LU$? The answer is no. We will examine when the factorization exists and what makes it unique.

2.  **Numerical stability:** Computers use finite-precision floating-point arithmetic, so calculations can introduce *roundoff errors*. These errors can accumulate or be amplified enough to make a computed solution inaccurate. How can we assess their effect?

3.  **Efficiency:** Large systems require methods that use time and memory efficiently as the matrix size grows.

Introducing **row pivoting** resolves the existence problem for invertible matrices and leads to the more general $PA = LU$ factorization. Pivoting also plays a central role in controlling roundoff error.

To assess accuracy, we will study floating-point arithmetic and develop a framework for **error analysis**. We will distinguish between *forward error*, which measures how close a computed solution is to the exact solution, and *backward error*, which measures how much the problem must change for the computed solution to be exact. The **condition number** of a matrix describes sensitivity to small perturbations and helps explain when small backward errors imply small forward errors.

We will begin with triangular systems, then introduce LU factorization, examine how it can fail, and develop a pivoting strategy. Next, we will use floating-point arithmetic and backward error analysis to understand why LU factorization with pivoting is usually stable in practice and where its guarantees have limits. Finally, we will study **Cholesky factorization**, an efficient and stable method for symmetric positive definite matrices.

By the end of the chapter, you should be able to solve linear systems efficiently and assess the accuracy and limitations of the computed solutions.
