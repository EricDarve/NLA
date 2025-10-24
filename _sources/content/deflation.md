# The Method of Deflation

In this section, we will start using the Schur decomposition to compute multiple eigenvalues and eigenvectors of a matrix $A$. We will build on the power method introduced in the previous section. Recall that the Schur decomposition of $A$ is:

$$
A = Q T Q^H,
$$

where $Q$ is unitary and $T$ is upper triangular. The diagonal entries of
$T$ are the eigenvalues of $A$, and the columns of $Q$ are the corresponding Schur vectors.

Assume the power iteration has converged, yielding the dominant eigenvalue $\lambda_1$ and a corresponding Schur vector $\boldsymbol{q}_1$ (which is the same as the eigenvector $\boldsymbol{x}_1$ of $A$). The goal is now to find the second largest eigenvalue, $\lambda_2$, and its corresponding Schur vector, $\boldsymbol{q}_2$. The key idea to find $\lambda_2$ is to "deflate" the matrix $A$, removing the influence of $\lambda_1$ and $\boldsymbol{q}_1$. We can then apply the power method to this new, deflated matrix, which will now have $\lambda_2$ as its dominant eigenvalue.

## Construct the Deflation Projector

We define an orthogonal projector $P$ that maps any vector onto the subspace *orthogonal* to $\boldsymbol{q}_1$. Assuming we have normalized $\boldsymbol{q}_1$ such that $\|\boldsymbol{q}_1\|_2 = 1$, this projector is:

$$
P = I - \boldsymbol{q}_1 \boldsymbol{q}_1^H
$$

This matrix has a simple action:

* For any vector $\boldsymbol{v}$ parallel to $\boldsymbol{q}_1$ (i.e., $\boldsymbol{v} = c \boldsymbol{q}_1$), $P \boldsymbol{v} = \boldsymbol{0}$.
* For any vector $\boldsymbol{w}$ orthogonal to $\boldsymbol{q}_1$ (i.e., $\boldsymbol{q}_1^H \boldsymbol{w} = 0$), $P \boldsymbol{w} = \boldsymbol{w}$.

## Define the Deflated Matrix

We create a new matrix, $M$, by applying this projector:

$$M = P A$$

Our strategy is to run the power iteration on $M$. To understand why this works, we must find the eigenvalues of $M$.

## Analyze the Eigenvalues of $M$

We will use the Schur decomposition of $A$, $A = Q T Q^H$, where 

$$
Q = [\boldsymbol{q}_1, \boldsymbol{q}_2, \dots, \boldsymbol{q}_n]
$$ 

is unitary and $T$ is upper triangular. The diagonal of $T$ contains the eigenvalues $\{\lambda_1, \lambda_2, \dots, \lambda_n\}$.

Let's look at the matrix $M$ in the Schur basis by computing the similarity transformation $Q^H M Q$:

$$
Q^H M Q = Q^H (P A) Q = (Q^H P Q) (Q^H A Q)
$$

Let's analyze the two parts of this product:

1.  **$Q^H A Q$**: This is, by definition, the upper triangular matrix $T$.

2.  **$Q^H P Q$**: Let's substitute $P = I - \boldsymbol{q}_1 \boldsymbol{q}_1^H$:

$$
Q^H P Q = Q^H (I - \boldsymbol{q}_1 \boldsymbol{q}_1^H) Q = Q^H I Q - (Q^H \boldsymbol{q}_1) (\boldsymbol{q}_1^H Q)
$$

* $Q^H I Q = Q^H Q = I$.
* $Q^H \boldsymbol{q}_1$ is the first column of $Q^H Q$, which is $\boldsymbol{e}_1 = (1, 0, \dots, 0)^T$.
* $\boldsymbol{q}_1^H Q$ is the first row of $Q^H Q$, which is $\boldsymbol{e}_1^T = (1, 0, \dots, 0)$.
* Therefore, $(Q^H \boldsymbol{q}_1) (\boldsymbol{q}_1^H Q) = \boldsymbol{e}_1 \boldsymbol{e}_1^T = \text{diag}(1, 0, \dots, 0)$.

This gives us $Q^H P Q = I - \text{diag}(1, 0, \dots, 0) = \text{diag}(0, 1, \dots, 1)$.

Now, we can compute the full transformation:

$$
Q^H M Q = \underbrace{\text{diag}(0, 1, \dots, 1)}_{\tilde{P}} \underbrace{(T)}_{\text{Schur form}}
$$

Multiplying the upper-triangular matrix $T$ by this diagonal matrix $\tilde{P}$ simply **zeros out the entire first row of $T$**. Let's call this new matrix $\tilde{T}$:

$$
\tilde{T} = \tilde{P} T =
\begin{pmatrix}
0 & 0 & \cdots & 0 \\
0 & \lambda_2 & t_{23} & \cdots \\
\vdots & & \ddots & \\
0 & \cdots & & \lambda_n
\end{pmatrix}
$$

The matrix $M$ is similar to $\tilde{T}$ (since $M = Q \tilde{T} Q^H$), which means they share the same eigenvalues. The eigenvalues of $\tilde{T}$ are its diagonal entries:

$$
\text{Eigenvalues}(M) = \{0, \lambda_2, \lambda_3, \dots, \lambda_n\}
$$

## Apply the Power Method

We now apply the power iteration to the matrix $M = PA$.

* The eigenvalues of $M$ are $\{0, \lambda_2, \dots, \lambda_n\}$.
* We assumed that $|\lambda_2| > |\lambda_3| \ge \dots \ge |\lambda_n|$.
* Therefore, the **strictly dominant eigenvalue of $M$ is $\lambda_2$**.

The power method applied to $M$ will converge to its dominant eigenvalue, $\lambda_2$, and its corresponding eigenvector, $\boldsymbol{q}_2$.

We can prove $\boldsymbol{q}_2$ is the eigenvector: 

$$
M \boldsymbol{q}_2 = PA \boldsymbol{q}_2 = P(T_{12}\boldsymbol{q}_1 + \lambda_2 \boldsymbol{q}_2) = \lambda_2 \boldsymbol{q}_2
$$

This process, called **deflation**, can be repeated. After finding $\boldsymbol{q}_2$, we can form $P_2 = I - \boldsymbol{q}_1 \boldsymbol{q}_1^H - \boldsymbol{q}_2 \boldsymbol{q}_2^H$ and apply the power method to $P_2 A$ to find $\lambda_3$, and so on. This generalization is described next.

## Deflation: A General Step for $\lambda_{i+1}$

This method of deflation can be applied sequentially to find all the eigenvalues. Let's assume we have already found the first $i$ Schur vectors, $\boldsymbol{q}_1, \dots, \boldsymbol{q}_i$. Our goal is to find $\lambda_{i+1}$.

### Construct the General Deflation Projector

First, we define a projector $P_i$ that maps any vector onto the subspace *orthogonal* to the entire subspace spanned by our known vectors, $\text{span}\{\boldsymbol{q}_1, \dots, \boldsymbol{q}_i\}$.

Let $Q_i = [\boldsymbol{q}_1, \dots, \boldsymbol{q}_i]$ be the $n \times i$ matrix with these vectors as its columns. Since these are orthonormal Schur vectors, $Q_i^H Q_i = I_i$ (the $i \times i$ identity).

The projector *onto* this subspace is $Q_i Q_i^H$. The projector *orthogonal* to this subspace is therefore:

$$
P_i = I - Q_i Q_i^H
$$

### Define the Deflated Matrix

As before, we create a new deflated matrix, $M_i$, by applying this projector to $A$:

$$
M_i = P_i A
$$

We will now run the power iteration on $M_i$.

### Analyze the Eigenvalues of $M_i$

We use the same strategy: find the eigenvalues of $M_i$ by performing a similarity transformation with the full Schur basis $Q = [\boldsymbol{q}_1, \dots, \boldsymbol{q}_n]$.

$$
Q^H M_i Q = Q^H (P_i A) Q = (Q^H P_i Q) (Q^H A Q)
$$

Again, let's analyze the two parts:

1.  **$Q^H A Q$**: This is just the upper triangular matrix $T$.
2.  **$Q^H P_i Q$**: We substitute the definition of $P_i$:

$$
Q^H P_i Q = Q^H (I - Q_i Q_i^H) Q = Q^H Q - (Q^H Q_i) (Q_i^H Q)
$$

* $Q^H Q = I$ (the $n \times n$ identity).
* $Q^H Q_i = Q^H [\boldsymbol{q}_1, \dots, \boldsymbol{q}_i]$ is the first $i$ columns of $Q^H Q = I$. This is an $n \times i$ matrix, $\begin{pmatrix} I_i \\ 0 \end{pmatrix}$.
* $Q_i^H Q = [\boldsymbol{q}_1, \dots, \boldsymbol{q}_i]^H Q$ is the first $i$ rows of $Q^H Q = I$. This is an $i \times n$ matrix, $\begin{pmatrix} I_i & 0 \end{pmatrix}$.
* Their product is:

$$
(Q^H Q_i) (Q_i^H Q) = \begin{pmatrix} I_i \\ 0 \end{pmatrix} \begin{pmatrix} I_i & 0 \end{pmatrix} = \begin{pmatrix} I_i & 0 \\ 0 & 0 \end{pmatrix}
$$

This is an $n \times n$ block-diagonal matrix, which is $\text{diag}(\underbrace{1, \dots, 1}_{i \text{ times}}, 0, \dots, 0)$.

This gives us 

$$
Q^H P_i Q = I - \text{diag}(1, \dots, 1, 0, \dots, 0) = \text{diag}(\underbrace{0, \dots, 0}_{i \text{ times}}, 1, \dots, 1).
$$

Let's call this projector $\tilde{P}_i$.

### Apply the Power Method

Now we compute the full transformation:

$$
Q^H M_i Q = \tilde{P}_i T = \text{diag}(0, \dots, 0, 1, \dots, 1) \cdot T
$$

Multiplying $T$ by this diagonal matrix $\tilde{P}_i$ **zeros out the first $i$ rows of $T$**. The resulting matrix, $\tilde{T}_i$, looks like this:

$$
\tilde{T}_i =
\begin{pmatrix}
0 & \cdots & 0 & 0 & \cdots & 0 \\
\vdots & \ddots & \vdots & \vdots & & \vdots \\
0 & \cdots & 0 & 0 & \cdots & 0 \\
0 & \cdots & 0 & \lambda_{i+1} & t_{i+1, i+2} & \cdots \\
\vdots & & \vdots & & \ddots & \\
0 & \cdots & 0 & 0 & \cdots & \lambda_n
\end{pmatrix}
$$

The matrix $M_i$ is similar to $\tilde{T}_i$, so their eigenvalues are the same. The eigenvalues of $\tilde{T}_i$ are its diagonal entries:

$$
\text{Eigenvalues}(M_i) = \{\underbrace{0, \dots, 0}_{i \text{ times}}, \lambda_{i+1}, \lambda_{i+2}, \dots, \lambda_n\}
$$

Assuming we have a strict separation, $|\lambda_{i+1}| > |\lambda_{i+2}|$, the **strictly dominant eigenvalue of $M_i$ is $\lambda_{i+1}$**.

Therefore, applying the power iteration to $M_i = P_i A$ will cause the iterates to converge to $\lambda_{i+1}$. The corresponding eigenvector can be shown to be $\boldsymbol{q}_{i+1}$.