# Classical Iterative Methods

For a large sparse linear system $Ax=b$, computing a factorization may require too much work or storage. Iterative methods offer another approach: start with a guess and improve it through a sequence of inexpensive corrections.

This chapter introduces simple methods that help explain how iteration works, why it can succeed or fail, and how it can be accelerated.

## Why Use an Iterative Method?

A sparse matrix contains relatively few nonzero entries. Such matrices arise naturally when each unknown interacts with only a few others. For example, a differential equation discretized on a grid often couples each grid point to its immediate neighbors.

Working with these local interactions can be inexpensive. Factorization, however, can introduce additional nonzero entries into the factors. This **fill-in** increases the work and memory required, sometimes substantially.

An iterative method can work directly with the sparse matrix and a few vectors. Its attraction is that each correction can be much cheaper than a factorization. The total cost then depends on how many corrections are needed to obtain an accurate enough solution. Sparse direct methods remain useful; the choice depends on the problem and the accuracy required.

## The Basic Idea: Repeated Inexpensive Corrections

Begin with an approximate solution and check how well it satisfies the equations. The mismatch is called the **residual**. Use that mismatch to compute a correction, update the approximation, and repeat.

The key is to make the correction affordable. A **splitting method** separates the matrix into a part that gives an easy linear solve and a remainder handled using the current approximation. Different choices lead to different iterative methods.

When the same update rule is used throughout, the method is called **stationary**. These methods are simple to understand and implement, but cheap updates alone do not guarantee a useful solver. The approximations may converge slowly or fail to converge. We therefore need to understand both whether a method converges and how quickly it does so.

## What This Chapter Develops

The [next section](splitting_methods_definition.md) makes the splitting idea precise and establishes a convergence criterion. The following sections develop the main methods:

- **[Jacobi](jacobi_method.md)** updates every component using the previous approximation. The updates are independent, making the method well suited to parallel computation.
- **[Gauss-Seidel](gauss_seidel_method.md)** uses newly computed values as soon as they become available. This changes how information passes through the system and introduces dependencies between updates.
- **[Successive over-relaxation (SOR)](sor_method.md)** adjusts the size of each Gauss-Seidel correction. A suitable choice can accelerate convergence.
- **[Chebyshev acceleration](chebyshev_acceleration.md)** chooses changing coefficients using bounds for the eigenvalues. It connects the polynomial ideas from the eigenvalue chapter to solving linear systems.

On a first reading, focus on how each method computes its corrections and what guarantees convergence. The arguments marked optional give further justification, but can be skipped without losing the main progression.

## Their Role in Larger Solvers

Classical methods can be slow when used on their own, but their inexpensive corrections remain useful inside more powerful algorithms.

One role is **preconditioning**: a simple approximate solve helps another iterative method make faster progress. We return to this idea in the later [preconditioning section](preconditioning.md).

Another role is **smoothing** in multigrid. For suitable grid-based problems, classical iterations can quickly reduce error that oscillates between neighboring points. Coarser grids then help remove error that varies more slowly. This division of work can make a simple iteration valuable even when it would be slow as a complete solver.
