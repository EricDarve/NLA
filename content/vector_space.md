# Vector space

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

### Vector addition

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

### Scalar multiplication

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

---

### Definition: Vector space

A **vector space** $V$ over a field $\mathbb{F}$ (such as $\mathbb{R}$ or $\mathbb{C}$) is a set equipped with addition $(x,y)\mapsto x+y$ and scalar multiplication $(\alpha,x)\mapsto \alpha x$ satisfying, for all $x,y,z\in V$ and $\alpha,\beta\in\mathbb{F}$:

1. $x+(y+z)=(x+y)+z$  (Associativity of addition)
2. $x+y=y+x$  (Commutativity of addition)
3. There exists $0\in V$ with $x+0=x$  (Additive identity)
4. For each $x$ there exists $-x$ with $x+(-x)=0$  (Additive inverse)
5. $\alpha(\beta x)=(\alpha\beta)x$  (Associativity of scalar multiplication / compatibility with field multiplication)
6. $1x=x$  (Multiplicative identity of the field acts as identity on vectors)
7. $\alpha(x+y)=\alpha x+\alpha y$  (Left distributivity of scalar over vector addition)
8. $(\alpha+\beta)x=\alpha x+\beta x$  (Right distributivity of scalar addition over scalar multiplication)

---

## Subspaces

A **subspace** $S \subset \mathbb{R}^n$ is a nonempty set closed under linear combinations. Equivalently, for any $x,y \in S$ and $\alpha,\beta \in \mathbb{R}$,

$$
\alpha x + \beta y \in S.
$$

Every subspace of $\mathbb{R}^n$ is itself a vector space under the same operations.

---

## Linear combinations and span

Given vectors $x_1,\dots,x_k \in \mathbb{R}^n$, a **linear combination** has the form

$$
\alpha_1 x_1 + \alpha_2 x_2 + \cdots + \alpha_k x_k,
\quad \alpha_1,\dots,\alpha_k \in \mathbb{R}.
$$

The **span** of $\{x_1,\dots,x_k\}$ is the set of all linear combinations:

$$
\operatorname{span}\{x_1,\dots,x_k\}
=
\left\{
\sum_{i=1}^k \alpha_i x_i \;:\; \alpha_1,\dots,\alpha_k \in \mathbb{R}
\right\}.
$$

The span is always a subspace of $\mathbb{R}^n$.

---

## Linear independence

Vectors $x_1,\dots,x_k$ are **linearly independent** if the only solution to the homogeneous combination equaling zero is the trivial one:

$$
\sum_{i=1}^k \alpha_i x_i = 0
\quad \Rightarrow \quad
\alpha_1=\cdots=\alpha_k=0.
$$

If there exists a nontrivial choice of coefficients yielding zero, the vectors are **linearly dependent**.

---

## Bases and dimension

A set of vectors $x_1,\dots,x_k$ is a **basis** for a subspace $S$ if:

1. $x_1,\dots,x_k$ are linearly independent, and
2. $\operatorname{span}\{x_1,\dots,x_k\} = S$.

If $x_1,\dots,x_k$ form a basis for $S$, then $k$ is the **dimension** of $S$, written

$$
\dim(S) = k.
$$

While a subspace can have many different bases, each basis has the same number of vectors. Hence dimension is well defined.

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

--

## Direct sum

If $U$ and $V$ are subspaces, then $U+V$ is a subspace. We say that $W=U \oplus V$ is the direct sum of $U$ and $V$ if $U \cap V = \{0\}$. The direct sum means that if a vector is decomposed into its $U$ and $V$ components, this decomposition is unique.

Example: verify that if $x_1$, ..., $x_k$ are linearly independent, then

$$
S = \text{span}\{x_1\} \oplus \cdots \oplus \text{span}\{x_k\}
$$
