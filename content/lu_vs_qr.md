# LU vs. QR Decomposition

The main difference is that LU decomposition is a specialized tool for solving square linear systems, while QR decomposition is a more robust and general method primarily used for solving overdetermined least-squares problems.

## Key Differences

| Feature | LU Decomposition | QR Decomposition |
| :--- | :--- | :--- |
| **Primary Use Case** | Solving square systems: $Ax=b$ | Solving least-squares problems: $Ax \approx b$ |
| **Applicability** | Requires a **square, invertible** matrix $A$. | Works for **any** $m \times n$ matrix $A$ with full column rank. |
| **Numerical Stability**| Can be unstable without pivoting. **LU with partial pivoting** is generally stable in practice but can fail in rare cases. | **Inherently more stable** because it uses orthogonal matrices, which do not amplify rounding errors. This is its biggest advantage. |
| **Computational Cost** | Faster for square systems. Costs $\approx \frac{2}{3}n^3$ flops. | Slower for square systems. Costs $\approx \frac{4}{3}n^3$ flops (about twice as much as LU). |

## Which Method Should You Choose? 🤔

Your choice depends entirely on the problem you are trying to solve.

* **Choose LU decomposition (with partial pivoting) when:**
    You need to solve a **square linear system** $Ax=b$ and **speed is a priority**. It is the faster, standard method for this specific task and is stable enough for most applications.

* **Choose QR decomposition when:**
    You are solving a **least-squares problem** (e.g., from a regression or data-fitting task where you have more equations than unknowns). QR is the standard, most numerically reliable method for this. You should also prefer QR if you have a square system but suspect it is **very ill-conditioned** and you need maximum numerical stability, even at the cost of speed.