[[QR factorization]]: $A = QR$

We can rewrite this equation as:
$$
Q^T A = R
$$
Interpretation: apply an orthogonal transformation to make the matrix triangular = orthogonal triangularization.

Start with 1st column.
$$
a_{,1} \rightarrow \| a_{,1} \|_2 \, e_1
$$

- We need to apply an orthogonal transformation that zeros out all entries except the first.
- Orthogonal transformations are either rotations or reflections.
- It turns out that reflections are easier to apply than rotations.

Orthogonal transformations are either rotations or reflections:

![[Pasted image 20231013125732.png]]

We want to find a reflection that sends the first column $a_{,1}$ to the $e_1$ axis:

![[Pasted image 20231013125829.png]]

These are the steps:

![[Pasted image 20231013125916.png]]

From this figure, we can derive the key definitions:
- Vector for projection: $v = x - \|x\|_2 e_1$.
- Reflection matrix: $P = I - \beta vv^T$ with $\beta = \frac{2}{v^T v}.$

This completely defines the Householder transformation $P$.

Properties:
- We can check that $P^T P = I$. In fact: $P=P^T$ and $P^2 = I$.
- We can check that $Px = \|x\|_2 e_1$.

[[QR factorization]]