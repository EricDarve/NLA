In the QR iteration, we focus on computing:
$$
T_k = Q_k^H A Q_k
$$
where $Q_k$ is the sequence produced by [[Orthogonal iteration|the method of orthogonal iteration]].

Let's consider $T_k$ at the next iteration:
$$
T_{k+1} = Q_{k+1}^H A Q_{k+1}
$$
Can we calculate $T_{k+1}$ from $T_k$ without computing $Q_k$?

**Step 1:**

$Q_{k+1} R_{k+1} = A Q_k$ from [[Orthogonal iteration|orthogonal iteration]].

$$
T_k = Q_k^H (A Q_k) = Q_k^H Q_{k+1} R_{k+1} \overset{def}{=} U_{k+1} R_{k+1}
$$
where we define $U_{k+1} \overset{def}{=} Q_k^H Q_{k+1}$.

We can interpret $Q_k^H Q_{k+1}$ as being a small rotation. It corresponds to an increment in $Q_k$ as $k$ increases and $Q_k$ converges to the [[Schur decomposition|Schur orthogonal transformation.]]

$U_{k+1}$ and $R_{k+1}$ can be computed using a [[QR factorization|QR decomposition]] without the help of the [[Orthogonal iteration|orthogonal iteration]]. $Q_k$ or $Q_{k+1}$ are not needed.

**Step 2:**

This is the definition for step $k+1$:
$$
T_{k+1} = Q_{k+1}^H A Q_{k+1}
$$
and recall that $T_k = Q_k^H A Q_k$. This leads to:
$$
T_{k+1} = Q_{k+1}^H Q_k T_k Q_k^H Q_{k+1}.
$$ 
Note $Q_k$ must be square for this to be true. Recall the definition from the previous step:
$$
U_{k+1} = Q_k^H Q_{k+1}
$$
We find that $T_k$ is unitarily similar to $T_{k+1}$:
$$
T_{k+1} = U_{k+1}^H T_k U_{k+1}
$$
But recall from the step above that
$$
T_k = U_{k+1} R_{k+1}
$$
So $U_{k+1}^H T_k = R_{k+1},$ and
$$
T_{k+1} = R_{k+1} U_{k+1}
$$

This is the basis for the method of QR iteration.

### Summary of the key steps

- Step 1: $T_k = U_{k+1} R_{k+1}$: QR decomposition.
- Step 2: $T_{k+1} = R_{k+1} U_{k+1}$: apply the orthogonal transformation to the right of $R_{k+1}$.

We just switch the order of the terms!

We can compute $T_{k+1}$ from $T_k$ without computing any of the $Q_k$.

Although these steps are different from [[Orthogonal iteration|orthogonal iteration]], this is still the same iteration. For example, the sequence of $R_k$ matrices and the convergence is the same for both methods.

We can also connect the sequence of $U_k$ matrices to the [[Orthogonal iteration|orthogonal iteration]]. Since
$$
U_{k+1} = Q_k^H Q_{k+1}
$$
we can prove by induction that
$$
Q_k = U_1 \cdots U_k
$$
if we start from $Q_0 = I$. The product of the $U_i$ in the [[QR iteration]] is equal to $Q_k$ in the [[Orthogonal iteration|orthogonal iteration]], and again converges to the orthogonal matrix $Q$ from the [[Schur decomposition]].

### QR iteration algorithm

```julia
Tk = A
while not_converged:
    Uk, Rk = qr(Tk)
    Tk = Rk * Uk
```

Here are the first few iterations:

- $Q_1 R_1 = A$
- $Q_2 R_2 = R_1 Q_1$
- $Q_3 R_3 = R_2 Q_2$
- $Q_4 R_4 = R_3 Q_3$
- ... 