# Least Squares Problems

In this chapter, we're diving into one of the most fundamental and widely applied problems in all of computational science and engineering: the **linear least-squares problem**. This problem arises whenever we have more data than parameters to describe it—a situation that's incredibly common in data fitting, machine learning, and statistical modeling.

## The Problem: Overdetermined Systems

Consider the familiar equation $Ax = b$. We've spent some time on cases where $A$ is a square, invertible matrix, which guarantees a unique solution. But what happens when $A$ is a "tall and skinny" matrix, meaning it's an $m \times n$ matrix with more rows than columns ($m > n$)? This is called an **overdetermined system**.

In this scenario, the vector $b$ generally does not lie in the column space of $A$. As a result, there is **no exact solution** $x$ that satisfies $Ax = b$. The system is inconsistent.

So, what do we do? We find the next best thing. If we can't find an $x$ that makes the residual vector $r = b - Ax$ exactly zero, we find an $x$ that makes the residual as small as possible. We measure the "size" of the residual using its 2-norm, which leads to the formal definition of the least-squares problem:

Find $x^*$ such that 

$$\| Ax^* - b \|_2 = \min_{x \in \mathbb{R}^n} \| Ax - b \|_2.$$

Geometrically, this means we are seeking the vector $Ax^*$ in the column space of $A$ that is closest to $b$. This vector is, in fact, the orthogonal projection of $b$ onto the column space of $A$.

## Three Paths to a Solution

In this section of the course, we will explore several powerful methods for solving this problem, each with its own strengths and trade-offs. We will cover three main approaches:

1.  **The Normal Equations**: Our first method is a classic one derived from calculus. By setting the gradient of the squared residual, $\|Ax - b\|_2^2$, to zero, we arrive at the square $n \times n$ linear system $(A^T A) x = A^T b$. These are the **normal equations**. While elegant and conceptually simple, this approach can suffer from numerical instability. The issue is that the condition number of $A^T A$ is the square of the condition number of $A$, which can amplify rounding errors in our computations.

2.  **The QR Factorization**: To develop a more numerically robust method, we will turn to orthogonal factorization. The core idea is to decompose our matrix $A$ into the product of an **orthogonal matrix $Q$** and an **upper triangular matrix $R$**, so that $A = QR$. This is an incredibly powerful approach. Since orthogonal matrices preserve the 2-norm, the least-squares problem is transformed into solving a simple and well-behaved triangular system $Rx = Q^T b$. We will learn three distinct algorithms to compute this factorization:

    * **Householder Reflections**: The stable, efficient workhorse algorithm for QR decomposition of dense matrices.
    * **Givens Rotations**: A tool for selectively introducing zeros, which is particularly useful for sparse matrices or parallel computing architectures.
    * **Gram-Schmidt Orthogonalization**: A conceptually straightforward method that builds the orthogonal basis step-by-step. It is most suitable for very thin matrices.

3.  **The Singular Value Decomposition (SVD)**: Our final approach is the most powerful and general of all. What if the columns of $A$ are not linearly independent, making the matrix **rank-deficient**? In this case, there are infinitely many solutions to the least-squares problem. The SVD provides a definitive answer. It allows us to solve the least-squares problem for *any* matrix $A$, regardless of its shape or rank. This will lead us to the concept of the **pseudo-inverse**, which gives us the unique solution $x^*$ that not only minimizes the residual but also has the smallest possible 2-norm.

Mastering these techniques will provide you with essential tools used in countless modern applications, from GPS navigation and image processing to the foundations of data science and machine learning.