We briefly outline the key idea for Conjugate Gradients based on the [[Key idea of iterative methods for eigenvalue computation|previous recurrence]] we derived for $H_k$ and $Q_k$.

The key idea is to replace
$$
Ax = b
$$
by the solution of a smaller linear system. Denote by $x_k = Q_k y.$ Then, we set
$$
Q_k^T A x_k = Q_k^T A Q_k \, y = H_k \, y = Q_k^T b
$$
From 
$$
y = H_k^{-1} \, Q_k^T b,
$$
we recover $x_k$ which approximates the exact solution $x$. The Conjugate Gradients will expand on this idea and make it computationally very efficient.