[[Vectors and matrices]], [[Dot product]]

A matrix can be defined as a linear map from $\mathbb R^n \to \mathbb R^m$. This is the operator view of a matrix: 
$$
A: x \to Ax
$$

Algebraically, we write
$$
y = Ax
$$
with
$$
y_i = \sum_{j=1}^n a_{ij} \, x_j
$$
Computational cost: $O(n)$.

Using the operator interpretation, we can define the product of two matrices $C=AB$ as the result of composing $A$ with $B$: 
$$
Cx = (AB)x = A(Bx)
$$
This works when $B: \mathbb R^n \to \mathbb R^p$ and $A: \mathbb R^p \to \mathbb R^m$.

Algebraically, we have
$$
c_{ij} = \sum_{k=1}^p a_{ik} \, b_{kj}
$$
For the product to be defined, the number of columns of $A$, $p$, must be equal to the number of rows of $B$.

Computational cost: $O(mnp)$.