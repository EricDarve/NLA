Let's consider each step of the [[Preconditioning the Conjugate Gradients algorithm|naive CG algorithm]] with a preconditioner $C$. This is the step where we multiply by the matrix:
$$
\begin{aligned}
  q^{(k)} & = CAC^T p^{(k)} \\[.3em]
  \mu_k & = \frac{ \| s^{(k-1)} \|_2^2 }{(p^{(k)})^T q^{(k)}}
\end{aligned}
$$
This suggests using the following definition:
$$
\underline{p}^{(k)} = C^T p^{(k)}.
$$
The following recurrence will work
$$
\underline q^{(k)} = AC^T p^{(k)} = A \underline p^{(k)}.
$$
Note that $C$ is missing on the left.

Then we have
$$
\mu_k = \frac{ \| s^{(k-1)} \|_2^2 }{(\underline p^{(k)})^T \underline q^{(k)}}
$$
with $q^{(k)} = C \underline q^{(k)}$ from our definition of $\underline q^{(k)}$.

### Solution update step

The [[Preconditioning the Conjugate Gradients algorithm|solution update step]] for $y^{(k)}$ is:
$$
y^{(k)} = y^{(k-1)} + \mu_k \, p^{(k)}
$$
[[Preconditioning the Conjugate Gradients algorithm|Recall that]] $x^{(k)} = C^T y^{(k)}.$ Let's multiply the equation above by $C^T$
$$
C^T y^{(k)} = C^T y^{(k-1)} + \mu_k \, C^T p^{(k)}
$$
This is just
$$
x^{(k)} = x^{(k-1)} + \mu_k \, \underline p^{(k)}
$$

### Residual

Let's focus on the [[Preconditioning the Conjugate Gradients algorithm|residual update step.]] The definition of the residual for CG with symmetric preconditioning is
$$
s^{(k)} = Cb - CAC^Ty^{(k)} = Cb - CAx^{(k)}
$$
[[Residuals and solution increments in CG|Recall that]]
$$
r^{(k)} = b - Ax^{(k)}
$$
So: 
$$
s^{(k)} = C r^{(k)}.
$$

Let's see how we can [[Preconditioning the Conjugate Gradients algorithm|update the residual.]]
$$
s^{(k)} = s^{(k-1)} - \mu_k \, C \underline q^{(k)}
$$
Let's rewrite this in terms of the residual $r^{(k)}$:
$$
C r^{(k)} = C r^{(k-1)} - \mu_k \, C \underline q^{(k)}
$$
We multiply by $C^{-1}$ to the left to simplify and get our final expression:
$$
r^{(k)} = r^{(k-1)} - \mu_k \, \underline q^{(k)}
$$

### Search direction

The search direction [[Preconditioning the Conjugate Gradients algorithm|update]] is:
$$
p^{(k+1)} = s^{(k)} + \tau_k \, p^{(k)}
$$
With our definition $\underline{p}^{(k)} = C^T p^{(k)},$ we get:
$$
C^{-T} \underline p^{(k+1)} = s^{(k)} + \tau_k \, C^{-T} \underline p^{(k)}
$$
Multiply to the left by $C^T$:
$$
\underline p^{(k+1)} = C^T s^{(k)} + \tau_k \, \underline p^{(k)}
$$
Let's introduce the residual $r^{(k)}$. Recall that:
$$
s^{(k)} = C r^{(k)}
$$
So 
$$
C^T s^{(k)} = C^T C r^{(k)} = M r^{(k)}.
$$
We can substitute in the previous equation:
$$
\underline p^{(k+1)} = M r^{(k)} + \tau_k \, \underline p^{(k)}
$$
In PCG, we use a new temporary vector $z^{(k)} = M r^{(k)}$ to reduce the time cost. We finally get:
$$
\underline p^{(k+1)} = z^{(k)} + \tau_k \, \underline p^{(k)}
$$

### Norm of residual

We are almost done. The last quantity we need is $\| s^{(k)} \|_2.$ We again use the relation for the residual $s^{(k)} = C r^{(k)}.$ So
$$
\| s^{(k)} \|_2^2
= (r^{(k)})^T C^T C r^{(k)} = (r^{(k)})^T M r^{(k)}
$$
As desired, $C$ has disappeared in favor of $M$!

We can even use our temporary work vector $z^{(k)} = M r^{(k)}$ to simplify the expression further:
$$
\| s^{(k)} \|_2^2 = [r^{(k)}]^T z^{(k)}
$$

This is the last step in the Preconditioned CG algorithm!

## Summary of PCG

All the steps above are summarized to form the PCG algorithm:
$$
\begin{aligned}
  \underline q^{(k)} & = A \underline p^{(k)} \\[.3em]
  \mu_k & = \frac{ [r^{(k-1)}]^T z^{(k-1)} }{(\underline p^{(k)})^T \underline q^{(k)}} \\[.3em]
  x^{(k)} & = x^{(k-1)} + \mu_k \underline p^{(k)} \\[.3em]
  r^{(k)} & = r^{(k-1)} - \mu_k \underline q^{(k)} \\[.3em]
  z^{(k)} & = M r^{(k)} \\[.3em]
  \tau_k & = \frac{[r^{(k)}]^T z^{(k)}}{[r^{(k-1)}]^T z^{(k-1)}} \\[.3em]
  \underline p^{(k+1)} & = z^{(k)} + \tau_k \underline p^{(k)}
\end{aligned}
$$
We see that this is only a very small modification of the [[Conjugate Gradients algorithm|original]] CG algorithm!