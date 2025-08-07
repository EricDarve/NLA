In some of the subsequent discussions, it is convenient to discuss how a subspace may converge to another subspace.

Consider two subspaces $U$ and $V$. We can define two orthogonal projects $P_U$ and $P_V$ unto $U$ and $V$. Then the generalized angle or distance between these subspaces is defined by
$$
\|P_U - P_V\|_2
$$
We can make sense of this definition on a simple example. Consider $x$ and $y$ 2 unit vectors. What is the distance between span$\{x\}$ and span$\{y\}$?

Denote by $\theta \ge 0$ the angle between $x$ and $y$, and denote by  $c=\cos(\theta)$, $s=\sin(\theta)$. Then, there exists a unit vector $z$ orthogonal to $x$ in span$\{x,y\}$ such that
$$
y = cx + sz
$$
Then
$$
\begin{align}
x x^T - y y^T
=
\begin{pmatrix}
x & z
\end{pmatrix}
\begin{pmatrix}
1-c^2 & -cs \\
-cs & -s^2
\end{pmatrix}
\begin{pmatrix}
x^T \\ z^T
\end{pmatrix} \\
= s
\begin{pmatrix}
x & z
\end{pmatrix}
\begin{pmatrix}
s & -c \\
-c & -s
\end{pmatrix}
\begin{pmatrix}
x^T \\ z^T
\end{pmatrix}
\end{align}
$$
Since the matrix
$$
\begin{pmatrix}
x & z
\end{pmatrix}
$$
is orthogonal, we have that:
$$
\|x x^T - y y^T\|_2 = s
\; \Big\| \!
\begin{pmatrix}
s & -c \\
-c & -s
\end{pmatrix}
\Big\|_2
= s
$$
since
$$
\begin{pmatrix}
s & -c \\
-c & -s
\end{pmatrix}
$$
is orthogonal.

So the distance between these subspaces is equal to $s=\sin(\theta)$. That distance is 0 if the subspaces are equal. The maximum distance is 1.