If we know an interval containing the eigenvalues, it is possible to design a scheme with near “optimal” convergence.

Let's start with:
$$
\begin{gather}
A = M - N \\
x^{(k+1)} = M^{-1} (b + N x^{(k)}) \\
x^{(k+1)} = M^{-1} (b + (M - A) x^{(k)}) \\
\end{gather}
$$
This implies that
$$
x^{(k+1)} = x^{(k)} + M^{-1} (b - A x^{(k)})
$$
This is the residual iteration form.

Let's avoid repeating $M^{-1}$ in all the equations and let's do the substitution:
$$
M^{-1} b \rightarrow b, \quad M^{-1} A \rightarrow A.
$$
Then the iteration is:
$$
x^{(k+1)} = x^{(k)} + b - A x^{(k)}
$$
This substitution can also be interpreted as solving $M^{-1} A x = M^{-1} b$ instead of $Ax = b$.

Let's apply our idea of relaxation:
$$
x^{(k+1)} = x^{(k)} + \omega_k (b - A x^{(k)})
$$
Let's use $x = x + \omega_k \, (b - A x)$, and subtract:
$$
\begin{gather}
e^{(k)} = x^{(k)} - x \\
e^{(k+1)} = (I - \omega_k \, A) \; e^{(k)}
\end{gather}
$$
We get:
$$
\begin{gather}
e^{(k+1)} = (I - \omega_k \, A) \cdots (I - \omega_0 A) \; e^{(0)} \\[.5em]
e^{(k)} = p_k(A) \, e^{(0)}
\end{gather}
$$
where $p_k$ is a polynomial of order $k$ with $p_k(0) = 1$.

Note: there is a variant where we use
$$
x^{(k+1)} = x^{(k)} + \omega_k (b - A x^{(k)}) - \alpha_k (x^{(k)} - x^{(k-1)})
$$
Then:
$$
e^{(k+1)} = ((1-\alpha_k) I - \omega_k \, A) \; e^{(k)} + \alpha_k e^{(k-1)}
$$
In both cases: $e^{(k)} = p_k(A) \, e^{(0)}$ with $p_k(0) = 1.$

The polynomial $p_k$ needs to be designed to provide maximum error reduction. Estimating $p_k$ requires having an estimate of an interval containing the eigenvalues of $A.$

Assume that the eigenvalues of $A$ are in the interval $[m,M]$.

Then define the following polynomial:
$$
p_k(x) = \frac{T_k(z(x))}{T_k(z(0))}
$$
where we used a shift-and-scale transformation of $x$:
$$
z(x) = \frac{m + M - 2 x}{M - m}
$$
We have $q_k(0) = 1$. Moreover when $m \le x \le M$:
$$
\begin{gather}
2m \le 2x \le 2M \\[.3em]
-2M \le -2x \le -2m \\[.3em]
m-M \le m+M-2x \le M-m \\[.3em]
-1 \le z(x) \le 1
\end{gather}
$$
So $|T_k(z(x))| \le 1$ and
$$
\begin{gather}
|p_k(x)| \le \frac{1}{|T_k(z(0))|} \\[.5em]
z(0) = \frac{M + m}{M - m} = 1 + \frac{2m}{M - m} > 1
\end{gather}
$$
So as before $|T_k(z(0))| \gg 1$ and $|p_k(x)| \ll 1$ for $m \le x \le M.$

We can prove the following bound:
$$
\begin{gather}
\| e^{(k)} \|_2 \le  \| p_k(A) \|_2 \, \| e^{(0)} \|_2 \\[1em]
\| p_k(A) \|_2 \le \|X\|_2 \|X^{-1}\|_2 \|\Lambda\|_2
\le \kappa(X) \max_{x \in [m,M]} |p_k(x)|
\end{gather}
$$
We can expect a rapid convergence. But this algorithm requires a good estimation of $m$ and $M$.

We can derive a more precise bound for 
$$
\max_{x \in [m,M]} |p_k(x)|
$$

We have for $|z| > 1$:
$$
T_n(z) = \frac{1}{2} \Big[ (z + \sqrt{z^2 - 1})^n 
+ (z - \sqrt{z^2 - 1})^n \Big]
$$
where
$$
z=\frac{M + m}{M - m}
$$
Denote by $u$
$$
u = \frac{m}{M}
$$
Then using $a^2 - b^2 = (a-b)(a+b)$:
$$
\begin{align}
z + \sqrt{z^2 - 1}
& = \frac{1 + u}{1 - u} + \sqrt{ \Big(\frac{1 + u}{1 -u}\Big)^2 - 1} \\
& = \frac{1 + u + 2 \sqrt{u}}{1 - u}
= \frac{(1 + \sqrt{u})^2}{(1 - \sqrt{u})(1 + \sqrt{u})} \\
& = \frac{1 + \sqrt{u}}{1 - \sqrt{u}}
\end{align}
$$
Moreover:
$$
z - \sqrt{z^2 - 1} 
= \frac{1}{z + \sqrt{z^2 - 1}}
$$
Finally:
$$
\max_{x \in [m,M]} |p_k(x)| \le \frac{1}{|T_n(z)|} = \frac{2}{\tau^{-n} + \tau^n}
$$
with 
$$
\tau = \frac{1 - \sqrt{u}}{1 + \sqrt{u}}
$$
This expression can be simplified. Assume that the interval containing the eigenvalues is small: $m <M$ and $m \to M$. Then $u \to 1$ and $\tau \to 0$. We get:
$$
\max_{x \in [m,M]} |p_k(x)| \le \frac{1}{|T_n(z)|} \le 2 \tau^n 
\le 2 \Big( 1-\sqrt{\frac{m}{M}} \Big)^n
$$
The convergence is very fast. The figure below shows the polynomial with $m=0.2$ and $M=2.2$.

![[2022-11-09-15-59-27.png|600]]