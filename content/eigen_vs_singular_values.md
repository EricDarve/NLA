# Eigenvalues and Singular Values

The fundamental difference between eigenvalues and singular values is this: **eigenvalues** tell you about the long-term behavior of a system ($A^n$), while **singular values** tell you about the immediate, one-time amplification or "stretching" of a matrix. A matrix can have small, stable eigenvalues but produce enormous transient growth, a phenomenon captured by its singular values.


## The Setup: An Illustrative Example

Let's construct a matrix $A$ whose eigenvalues are well-behaved, but whose singular values indicate a potential for massive amplification.

1.  **Start with two nearly-aligned vectors**, where $\epsilon$ is a small positive number:

    $$
    u = \begin{bmatrix} 1 \\ \epsilon \end{bmatrix}, \quad v = \begin{bmatrix} 1 \\ -\epsilon \end{bmatrix}
    $$

    We use these to form the matrix $X = [u, v]$. As $\epsilon \to 0$, the columns become nearly parallel, making $X$ very **ill-conditioned**.

    $$
    X = \begin{bmatrix} 1 & 1 \\ \epsilon & -\epsilon \end{bmatrix}
    $$

2.  **Define a simple diagonal matrix** $D$:

    $$
    D = \begin{bmatrix} 1 & 0 \\ 0 & -1 \end{bmatrix}
    $$

3.  **Construct our matrix** $A$ using a similarity transform: $A = XDX^{-1}$. Let's compute this explicitly. The inverse of $X$ is:

    $$
    X^{-1} = \frac{1}{\det(X)} \begin{bmatrix} -\epsilon & -1 \\ -\epsilon & 1 \end{bmatrix} = \frac{1}{-2\epsilon} \begin{bmatrix} -\epsilon & -1 \\ -\epsilon & 1 \end{bmatrix} = \begin{bmatrix} 1/2 & 1/(2\epsilon) \\ 1/2 & -1/(2\epsilon) \end{bmatrix}
    $$

    Notice that as $\epsilon \to 0$, the entries of $X^{-1}$ blow up. Now, we find $A$:

    $$
    \begin{gathered}
    A = \begin{bmatrix} 1 & 1 \\ \epsilon & -\epsilon \end{bmatrix} \begin{bmatrix} 1 & 0 \\ 0 & -1 \end{bmatrix} \begin{bmatrix} 1/2 & 1/(2\epsilon) \\ 1/2 & -1/(2\epsilon) \end{bmatrix} = \begin{bmatrix} 1 & -1 \\ \epsilon & \epsilon \end{bmatrix} \begin{bmatrix} 1/2 & 1/(2\epsilon) \\ 1/2 & -1/(2\epsilon) \end{bmatrix} \\[1em]
    A = \begin{bmatrix} 0 & 1/\epsilon \\ \epsilon & 0 \end{bmatrix}
    \end{gathered}
    $$

Now we have a simple matrix $A$ that depends on our small parameter $\epsilon$. Let's analyze its eigenvalues and singular values.

## Eigenvalue Analysis 🔬

By construction ($A = XDX^{-1}$), the eigenvalues of $A$ are the diagonal entries of $D$.

$$\lambda_1 = 1, \quad \lambda_2 = -1$$

The eigenvalues are perfectly stable and have a magnitude of 1, regardless of how small $\epsilon$ gets. From an eigenvalue perspective, this matrix looks completely harmless. If we consider the discrete dynamical system $x_{k+1} = Ax_k$, the magnitudes of the eigenvalues ($|\lambda_i|=1$) suggest that the system will not blow up over time. In fact, $A^2=I$, so the system is **stable and periodic.**

## Singular Value Analysis 🧐

The singular values ($\sigma_i$) are the square roots of the eigenvalues of $A^T A$.

$$A^T A = \begin{bmatrix} 0 & \epsilon \\ 1/\epsilon & 0 \end{bmatrix} \begin{bmatrix} 0 & 1/\epsilon \\ \epsilon & 0 \end{bmatrix} = \begin{bmatrix} \epsilon^2 & 0 \\ 0 & 1/\epsilon^2 \end{bmatrix}$$

The eigenvalues of this diagonal matrix are clearly $\epsilon^2$ and $1/\epsilon^2$. The singular values of $A$ are their square roots:

$$\sigma_1 = \sqrt{1/\epsilon^2} = 1/\epsilon, \quad \sigma_2 = \sqrt{\epsilon^2} = \epsilon$$

As $\epsilon \to 0$, the largest singular value **$\sigma_1 \to \infty$**.

The largest singular value is the 2-norm of the matrix, $\|A\|_2 = \sigma_{max}$. So, even though the eigenvalues are 1 and -1, the norm of the matrix is huge! This tells us that applying the matrix $A$ just *once* can stretch a vector by an enormous factor of $1/\epsilon$.


## The Takeaway: Normality is Key

So, why the dramatic difference? The answer lies in the **non-normality** of matrix $A$.

* A matrix is **normal** if $A^T A = AA^T$. For normal matrices, the singular values are simply the absolute values of the eigenvalues ($\sigma_i = |\lambda_i|$), and their eigenvectors are orthogonal.
* Our matrix $A$ is **not normal**:

    $$
    A^T A = \begin{bmatrix} \epsilon^2 & 0 \\ 0 & 1/\epsilon^2 \end{bmatrix} \quad \neq \quad AA^T = \begin{bmatrix} 1/\epsilon^2 & 0 \\ 0 & \epsilon^2 \end{bmatrix}
    $$

This non-normality has a critical geometric consequence: the eigenvectors of $A$ (the columns of $X$) are **not orthogonal**. As $\epsilon \to 0$, they become nearly parallel. In contrast, the singular vectors of $A$ form an orthogonal basis.

**In summary:**

* **Eigenvalues** describe the behavior of $A$ with respect to its (possibly non-orthogonal) eigenvectors. They reveal long-term, asymptotic behavior but can hide short-term transient effects.
* **Singular values** describe the behavior of $A$ with respect to an optimal, orthonormal basis. They reveal the maximum possible amplification the matrix can produce in a single application, defining its norm and its potential for creating large transient growth.

This example starkly illustrates that for non-normal matrices, the eigenvalues alone give an incomplete and potentially misleading picture of the matrix's behavior. You must also consider the singular values to understand the full story.