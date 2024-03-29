![[QR using Householder transformations 2023-10-18 11.37.36.excalidraw.svg]]

QR factorization process using Householder transformations:

1. Start from $A$
2. Apply a series of Householder reflections: $Q_{n-1}^T \cdots Q_1^T A = R$ to create zeros below the diagonal. $Q_i^T$ creates zeros below $a_{ii}$.
3. QR is obtained implicitly in the form: $A = Q_1 \cdots Q_{n-1} \, R.$

The cost of applying each $Q_i^T$ is $O(n^2)$. So the total cost is $O(n^3)$.

If matrix $A$ is $m \times n$, the cost is $O(mn^2)$.

[[QR factorization]], [[Householder transformation]], [[Applying a Householder transformation]]