# Summary

We can think of these three methods as different tools for the same job: making a matrix triangular.

## Householder Transformations (The Sledgehammer)

Geometric Intuition: Each Householder transformation is a perfect reflection across a plane or hyperplane. It's like placing a mirror in space to move an entire vector (and the sub-column of the matrix) onto a specific axis in one single, powerful step.

The Go-To for Stability: It's generally considered the most numerically stable of the three methods for QR factorization. When in doubt with a dense matrix, Householder is a robust and efficient choice.

Efficiency: It's the fastest method for dense matrices because it operates on entire columns at once, which maps well to modern computer memory and processing (using efficient BLAS-2 operations).

## Givens Rotations (The Scalpel)

Surgical Precision: A Givens rotation is more delicate. It operates on a 2D plane within the larger space, affecting only two rows at a time to zero out a single, specific element. It's the perfect tool when you don't want to disturb the whole matrix.

Great for Parallel Computing: Because Givens rotations can operate on non-overlapping pairs of rows independently (e.g., rotating rows 1-2 and 3-4 simultaneously), the algorithm is much easier to parallelize than Householder.

Connection to 3D Graphics: The 3D rotation matrices for roll, pitch, and yaw are just specific examples of Givens rotations. Composing these rotations can generate any orientation in 3D space, a fundamental concept in graphics and aerospace engineering.

## Gram-Schmidt Process (The Bricklayer)

Sequential Construction: Unlike the other two, which modify the existing matrix, Gram-Schmidt builds the orthonormal basis (Q) one vector at a time. This "brick-by-brick" approach is its most important feature.

Essential for Iterative Methods: This sequential nature is critical for advanced iterative algorithms like the Arnoldi iteration, which is used to find eigenvalues of very large, sparse matrices. Householder and Givens can't be used in the same way because they only produce the full orthogonal matrix at the end.

A Lesson in Stability: The contrast between the unstable Classical Gram-Schmidt and the stable Modified Gram-Schmidt is a classic case study in numerical analysis, showing how a small change in the order of operations can dramatically impact the result due to rounding errors.