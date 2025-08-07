We have identified three methods to solve the [[Least-squares problems|least-squares problem]].

[[Method of normal equation|Method 1]]: normal equation for tall skinny matrix $A$:
$$
x = (A^TA)^{-1} \, A^T \, b
$$
[[Least-squares solution using QR|Method 2]]: recommended method using [[QR factorization|QR]]:
$$
\begin{gather}
A = QR \\
x = R^{-1} \, Q^T \, b
\end{gather}
$$
[[Least-squares solution using SVD|Method 3]]: for rank deficient $A$, we use the [[Singular value decomposition|SVD]]:
$$
\begin{gather}
A = U \Sigma V^T \\
x = V \, \Sigma^{-1} \, U^T \, b
\end{gather}
$$