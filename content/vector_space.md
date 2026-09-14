# Vector spaces

We will primarily work with vectors in $\mathbb{R}^n$ and occasionally $\mathbb{C}^n$.

An element of $\mathbb{R}^n$ is a column vector

$$
x = \begin{pmatrix}
x_1\\
x_2\\
\vdots\\
x_n
\end{pmatrix},
\quad x_i \in \mathbb{R}.
$$

Throughout this page, a subscript on $x$ denotes a **component** of a single vector, while a list of distinct vectors is written $v_1, \dots, v_k$.

## Vector addition

Given $x,y \in \mathbb{R}^n$,

$$
\begin{pmatrix}
x_1\\
x_2\\
\vdots\\
x_n
\end{pmatrix}
+
\begin{pmatrix}
y_1\\
y_2\\
\vdots\\
y_n
\end{pmatrix}
=
\begin{pmatrix}
x_1+y_1\\
x_2+y_2\\
\vdots\\
x_n+y_n
\end{pmatrix}.
$$

## Scalar multiplication

For $\alpha \in \mathbb{R}$ and $x \in \mathbb{R}^n$,

$$
\alpha
\begin{pmatrix}
x_1\\
x_2\\
\vdots\\
x_n
\end{pmatrix}
=
\begin{pmatrix}
\alpha x_1\\
\alpha x_2\\
\vdots\\
\alpha x_n
\end{pmatrix}.
$$

Vectors in $\mathbb{C}^n$ are defined the same way, with scalars $\alpha \in \mathbb{C}$.

## Definition: Vector space

A **vector space** $V$ over a field $\mathbb{F}$ (such as $\mathbb{R}$ or $\mathbb{C}$) is a set equipped with two operations,

$$
+ \;:\; V \times V \to V,
\qquad
\cdot \;:\; \mathbb{F} \times V \to V,
$$

written $(x,y)\mapsto x+y$ and $(\alpha,x)\mapsto \alpha x$. Both operations take their values in $V$, so $V$ is **closed** under addition and scalar multiplication. The operations must satisfy, for all $x,y,z\in V$ and $\alpha,\beta\in\mathbb{F}$:

1. $x+(y+z)=(x+y)+z$  (Associativity of addition)
2. $x+y=y+x$  (Commutativity of addition)
3. There exists $0\in V$ with $x+0=x$  (Additive identity)
4. For each $x$ there exists $-x$ with $x+(-x)=0$  (Additive inverse)
5. $\alpha(\beta x)=(\alpha\beta)x$  (Compatibility with field multiplication)
6. $1x=x$  (Multiplicative identity of the field acts as identity on vectors)
7. $\alpha(x+y)=\alpha x+\alpha y$  (Distributivity of scalar multiplication over vector addition)
8. $(\alpha+\beta)x=\alpha x+\beta x$  (Distributivity of scalar multiplication over field addition)

## Subspaces

A **subspace** $S \subset \mathbb{R}^n$ is a nonempty set closed under linear combinations. Equivalently, for any $x,y \in S$ and $\alpha,\beta \in \mathbb{R}$,

$$
\alpha x + \beta y \in S.
$$

Nonemptiness is essential: choosing $\alpha=\beta=0$ for any $x \in S$ shows that every subspace contains the zero vector.

Every subspace of $\mathbb{R}^n$ is itself a vector space under the same operations. Two subspaces always exist: the **trivial subspace** $\{0\}$, and $\mathbb{R}^n$ itself.

## Linear combinations and span

Given vectors $v_1,\dots,v_k \in \mathbb{R}^n$, a **linear combination** has the form

$$
\alpha_1 v_1 + \alpha_2 v_2 + \cdots + \alpha_k v_k,
\quad \alpha_1,\dots,\alpha_k \in \mathbb{R}.
$$

The **span** of $\{v_1,\dots,v_k\}$ is the set of all linear combinations:

$$
\operatorname{span}\{v_1,\dots,v_k\}
=
\left\{
\sum_{i=1}^k \alpha_i v_i \;:\; \alpha_1,\dots,\alpha_k \in \mathbb{R}
\right\}.
$$

The span is always a subspace of $\mathbb{R}^n$. By convention, the span of the empty set is the trivial subspace, $\operatorname{span} \emptyset = \{0\}$.

## Linear independence

Vectors $v_1,\dots,v_k$ are **linearly independent** if the only solution to the homogeneous combination equaling zero is the trivial one:

$$
\sum_{i=1}^k \alpha_i v_i = 0
\quad \Rightarrow \quad
\alpha_1=\cdots=\alpha_k=0.
$$

If there exists a nontrivial choice of coefficients yielding zero, the vectors are **linearly dependent**.

## Bases and dimension

A set of vectors $v_1,\dots,v_k$ is a **basis** for a subspace $S$ if:

1. $v_1,\dots,v_k$ are linearly independent, and
2. $\operatorname{span}\{v_1,\dots,v_k\} = S$.

Every subspace of $\mathbb{R}^n$ has a basis. While a subspace can have many different bases, every basis of $S$ has the same number of vectors, and that common number is the **dimension** of $S$, written

$$
\dim(S) = k.
$$

Dimension is therefore well defined. A few consequences are worth recording:

- $\dim(\{0\}) = 0$: the trivial subspace has the empty set as a basis.
- If $S \subset \mathbb{R}^n$ is a subspace, then $\dim(S) \le n$, with equality only for $S = \mathbb{R}^n$.
- If $\dim(S) = k$, then any $k$ linearly independent vectors of $S$ already form a basis of $S$, and any set spanning $S$ has at least $k$ vectors.

**Example (a plane in $\mathbb{R}^3$).** Let

$$
a_1=\begin{pmatrix}1\\0\\0\end{pmatrix},\quad
a_2=\begin{pmatrix}1\\1\\0\end{pmatrix}.
$$

Then $a_1$ and $a_2$ are linearly independent, and

$$
\operatorname{span}\{a_1,a_2\}
=
\left\{
\begin{pmatrix}
\alpha\\
\beta\\
0
\end{pmatrix}
:\; \alpha,\beta\in\mathbb{R}
\right\},
$$

which is the $x_1$-$x_2$ plane in $\mathbb{R}^3$. Therefore $\dim(\operatorname{span}\{a_1,a_2\})=2$.

## Sums of subspaces

Given subspaces $U, W \subset \mathbb{R}^n$, their **sum** is the set of all vectors obtained by adding an element of each:

$$
U + W = \{\, u + w \;:\; u \in U,\; w \in W \,\}.
$$

The sum $U+W$ is again a subspace of $\mathbb{R}^n$; it is the smallest subspace containing both $U$ and $W$. Its dimension satisfies

$$
\dim(U+W) = \dim(U) + \dim(W) - \dim(U \cap W).
$$

## Direct sum

The sum of $U$ and $W$ is a **direct sum**, written $U \oplus W$, if every vector of $U+W$ has a **unique** decomposition

$$
x = u + w, \qquad u \in U,\; w \in W.
$$

This uniqueness is the defining property. For **two** subspaces it is equivalent to the simple test

$$
U \cap W = \{0\}:
$$

if $u+w = u'+w'$ with $u,u' \in U$ and $w,w' \in W$, then $u-u' = w'-w$ belongs to $U \cap W$, so both sides vanish precisely when the intersection is trivial.

When the sum is direct, dimensions add:

$$
\dim(U \oplus W) = \dim(U) + \dim(W),
$$

which is the case $\dim(U \cap W)=0$ of the formula in the previous section.

### More than two subspaces

For subspaces $U_1,\dots,U_k$, the sum $U_1+\cdots+U_k$ is a **direct sum**, written $U_1 \oplus \cdots \oplus U_k$, if the only way to write the zero vector as

$$
0 = u_1 + \cdots + u_k, \qquad u_i \in U_i,
$$

is the trivial one, $u_1=\cdots=u_k=0$. This is equivalent to uniqueness of the decomposition of every vector in the sum, and also to the condition

$$
U_i \cap \sum_{j \neq i} U_j = \{0\}
\qquad \text{for each } i=1,\dots,k.
$$

For a direct sum of $k$ subspaces, dimensions again add:

$$
\dim(U_1 \oplus \cdots \oplus U_k) = \sum_{i=1}^k \dim(U_i).
$$

**Exercise.** Verify that if $v_1,\dots,v_k$ are linearly independent, then the space $S$ spanned by these vectors is the direct sum of their individual spans:

$$
S = \operatorname{span}\{v_1\} \oplus \cdots \oplus \operatorname{span}\{v_k\}.
$$

Use the zero-vector criterion above: a decomposition $0 = u_1 + \cdots + u_k$ with $u_i = \alpha_i v_i \in \operatorname{span}\{v_i\}$ is exactly a linear dependence among $v_1,\dots,v_k$, so independence forces every $u_i$ to vanish.

A central example of a direct sum appears in the next section: every subspace $S \subset \mathbb{R}^n$ satisfies $\mathbb{R}^n = S \oplus S^\perp$, where $S^\perp$ is the orthogonal complement of $S$.
