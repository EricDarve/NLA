# Solving Linear Systems

Welcome to the heart of numerical linear algebra. The problem of solving a system of linear equations, written in matrix form as $Ax = b$, is arguably one of the most fundamental and frequently encountered computational problems in all of science and engineering. From simulating the airflow over a wing and analyzing electrical circuits to training machine learning models and pricing financial derivatives, the need to solve for an unknown vector $x$ is everywhere.

At first glance, the problem might seem simple. If $A$ is an invertible square matrix, a unique solution exists, and we can write it formally as $x = A^{-1}b$. So, why not just compute the inverse of $A$ and multiply it by $b$? While mathematically sound, this approach is a computational pitfall. Calculating a matrix inverse is significantly more expensive and, more importantly, often less numerically accurate than solving the system directly. Our first core principle is this: **we almost never compute a matrix inverse explicitly.**

So, what's the alternative? The central strategy of this chapter—and a recurring theme in numerical analysis—is **factorization**. The idea is to decompose the often dense and complicated matrix $A$ into a product of simpler matrices. The star of this chapter is the **LU factorization**, where we write $A = LU$, with $L$ being a lower triangular matrix and $U$ being an upper triangular one.

Why is this so helpful? Because solving triangular systems is computationally cheap and easy. By factoring $A$, we transform the single difficult problem $Ax = b$ into two simple ones:

1.  Solve $Ly = b$ for $y$ (using *forward substitution*).
2.  Solve $Ux = y$ for $x$ (using *backward substitution*).

This two-step process is dramatically faster and more stable than the inverse-based approach.

However, this elegant strategy brings its own set of challenges that form the core of our study:

1.  **Existence and Uniqueness:** Can any invertible matrix $A$ be factored as $A=LU$? As we will see, the answer is no. This leads us to the first practical hurdle.

2.  **Numerical Stability:** This is the most critical challenge. Computers do not work with real numbers; they use finite-precision floating-point arithmetic. This means that every calculation introduces a tiny *roundoff error*. While individual errors are small, they can accumulate catastrophically, rendering our computed solution meaningless. How can we trust our results?

3.  **Efficiency:** For the massive systems that arise in practice, the speed of our algorithm is paramount. We need methods that scale efficiently as the size of the matrix grows.

To navigate these challenges, we will develop a powerful set of tools and concepts. We'll see that the problem of existence can be solved by introducing **row pivoting**, which leads to the more general and robust $PA = LU$ factorization. Pivoting also turns out to be the key to ensuring numerical stability. 

To understand *why* it works, we must dive into the world of floating-point arithmetic and develop a rigorous framework for **error analysis**. We will distinguish between *forward error* (how close is our answer to the true answer?) and *backward error* (is our answer the exact solution to a slightly perturbed problem?). The concept of the **condition number** of a matrix will emerge as a crucial indicator of how sensitive a problem is to small perturbations, telling us when we can expect accurate solutions and when the problem is simply too ill-behaved.

This chapter will guide you through this landscape step-by-step. We will begin with the simplest case—solving triangular systems. We will then introduce the LU factorization and its algorithm, investigate its potential for failure, and develop the pivoting strategy that saves it. Following that, we'll take a necessary detour to formalize floating-point arithmetic and backward error analysis. Armed with these tools, we will prove why the LU factorization with pivoting is, in fact, numerically stable in practice. Finally, we will explore a special case for symmetric positive definite matrices, where the even faster and more stable **Cholesky factorization** can be used.

By the end of this chapter, you will not only know *how* to solve a linear system effectively but also understand *why* the method works, how to analyze its limitations, and how to trust the results you get from the computer.