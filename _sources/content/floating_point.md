---
jupytext:
  formats: md:myst
  text_representation:
    extension: .md
    format_name: myst
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Floating-Point Numbers

Computers store numbers with a fixed number of significant digits. Every time a result needs more digits, it must be rounded. To understand the effect of this rounding, begin by comparing the sizes of the quantities in the calculation.

**The main danger is cancellation: adding quantities of opposite signs, or subtracting nearly equal quantities, to obtain a much smaller result. Small errors in the quantities being combined can then be large relative to the result.**

You do not usually need a separate error variable for every operation to recognize this problem. Estimate the magnitudes of the intermediate quantities, identify where they cancel, and compare the rounding errors at that scale with the answer you want. The examples below develop this approach.

## Representation and Rounding

Floating-point notation is scientific notation with a fixed number of significant digits. In a decimal system with four significant digits, for example,

$$
\operatorname{fl}(3.14159265)=3.142,
\qquad
\operatorname{fl}(314159.265)=3.142\times 10^5.
$$

Here $\operatorname{fl}(x)$ denotes the stored approximation to $x$. The larger number has a larger absolute rounding error, but the relative errors are the same. For a nonzero exact value $x$ and an approximation $\widehat{x}$, these two measures are

$$
\text{absolute error}=|\widehat{x}-x|,
\qquad
\text{relative error}=\frac{|\widehat{x}-x|}{|x|}.
$$

When $x=0$, relative error is undefined; we use absolute error instead.

### Binary Floating-Point Numbers

Computers usually use base 2. A *normalized* binary floating-point number has the form

$$
\pm\left(1+\sum_{i=1}^{p-1}d_i2^{-i}\right)2^e,
\qquad d_i\in\{0,1\}.
$$

The sign determines whether the number is positive or negative, the exponent $e$ sets its scale, and the $p$ binary digits of the significand determine its precision. The leading digit is always 1 for a normalized number, so it need not be stored explicitly.

For example,

$$
3.140625
=\left(1+2^{-1}+2^{-4}+2^{-7}\right)2^1
$$

is represented exactly when at least eight significant binary digits are available. In contrast, $0.1$ has an infinite binary expansion and must be rounded.

Standard double precision (the IEEE 754 binary64 format) uses $p=53$. Within the interval $[2^e,2^{e+1})$, adjacent numbers are separated by $2^{e-52}$. Thus the absolute spacing grows with the numbers, while the relative spacing stays roughly constant. Small numbers do not generally have more significant digits than large ones.

The exponent also has a finite range. Zero has a separate representation, and *subnormal* numbers extend the range toward zero with fewer significant digits. A result that is too large causes **overflow**; a result below the normal range may lose relative precision through **underflow**. In the analysis below, we assume these range limits are not reached.

### The Error in One Operation

With rounding to nearest, a basic arithmetic operation on stored operands satisfies

$$
\operatorname{fl}(a\mathbin{\mathrm{op}}b)
=(a\mathbin{\mathrm{op}}b)(1+\delta),
\qquad |\delta|\le u,
$$

where $\mathrm{op}\in\{+,-,\times,/\}$, division requires $b\ne0$, and $u$ is the **unit roundoff**. For binary64,

$$
u=2^{-53}\approx 1.11\times10^{-16}.
$$

The gap between 1 and the next larger number is $2^{-52}=2u$; this gap is often called *machine epsilon* in software libraries.

The useful interpretation is that rounding a quantity of magnitude $M$ introduces an absolute error at most about $uM$. **This error is measured at the scale of the quantity being rounded.** It need not be small compared with a much smaller result obtained later.

The formula describes an operation on the stored operands. It does not account for errors already present in those operands. Keeping this distinction in mind resolves the apparent paradox that accurately rounded operations can produce an inaccurate final answer.

For instance, in ordinary Python floating-point arithmetic,

```{code-cell} ipython3
print(0.1 + 0.2)
print(0.1 + 0.2 == 0.3)
```

prints `0.30000000000000004` and `False`. This is a small representation and rounding effect, not a large loss of accuracy. When comparing approximate results, choose a tolerance appropriate to their scale and the accuracy required; exact equality is appropriate only when exact agreement is intended.

## Cancellation: Small Compared with What?

Suppose we want the difference of the exact decimal numbers

$$
x=1.2344,\qquad y=1.2336,
\qquad x-y=0.0008.
$$

In four-significant-digit arithmetic, both inputs round to $1.234$. Their computed difference is therefore zero.

| Quantity | Exact value | Stored value | Absolute error |
| --- | --- | --- | --- |
| $x$ | $1.2344$ | $1.234$ | $0.0004$ |
| $y$ | $1.2336$ | $1.234$ | $0.0004$ |
| $x-y$ | $0.0008$ | $0$ | $0.0008$ |

The errors in $x$ and $y$ are small compared with numbers of size 1. They are not small compared with the desired difference, $0.0008$. The relative error in that difference is 100%.

We can reproduce this with Python's `Decimal` arithmetic. Creating a `Decimal` from a string preserves the stated decimal value; unary `+` then rounds it to the chosen precision. Each `localcontext` confines the precision setting to that example.

```{code-cell} ipython3
from decimal import Decimal, localcontext

with localcontext() as ctx:
    ctx.prec = 4
    x = +Decimal('1.2344')
    y = +Decimal('1.2336')
    print("Stored inputs:", x, y)
    print("Computed difference:", x - y)

print("Exact difference:", Decimal('1.2344') - Decimal('1.2336'))
```

### The Subtraction Itself Can Be Exact

In this example, the subtraction $1.234-1.234=0$ is exact. The error arose when the inputs were rounded. Subtracting them exposed that error because their difference was so small.

By contrast, if the inputs are exactly $1.234$ and $1.233$, then their difference $0.001$ is exact in four-digit decimal arithmetic. Cancellation alone does not make an answer inaccurate. What matters is whether the quantities being subtracted already contain errors that are significant compared with their difference. This distinction between harmless and harmful cancellation is discussed in [Goldberg's treatment of floating-point arithmetic](https://docs.oracle.com/cd/E19957-01/806-3568/ncg_goldberg.html).

Let $x$ and $y$ denote the exact values we want to subtract, and let $\widehat{x}$ and $\widehat{y}$ denote their stored approximations. The hats mark quantities that may already contain rounding errors. The answer we want is $x-y$, but the computer subtracts $\widehat{x}$ and $\widehat{y}$ instead.

For now, assume that this subtraction introduces no further rounding error, as in the example above. The computed answer is then exactly $\widehat{x}-\widehat{y}$. Its difference from the desired answer is

$$
(\widehat{x}-\widehat{y})-(x-y)
=(\widehat{x}-x)-(\widehat{y}-y).
$$

The left-hand side is the error in the answer: the computed difference minus the exact difference we wanted. On the right, $\widehat{x}-x$ and $\widehat{y}-y$ are the errors already present in the two inputs. Thus, even an exact subtraction can give an inaccurate answer because it subtracts approximate inputs.

In our four-digit example, the two input errors are $-0.0004$ and $+0.0004$. Their difference is $-0.0008$, which is precisely the error in returning 0 instead of $0.0008$.

The important point is that the input errors need not become smaller when the leading digits cancel. If $x$ and $y$ have magnitude about $M$ and errors of size about $uM$, their difference $s=x-y$ can have relative error of size

$$
\boxed{\frac{uM}{|s|}.}
$$

This is a scale estimate, not a prediction of the exact error. The input errors might be smaller, or might cancel too. But it identifies when a calculation can lose accuracy.

For example, suppose two intermediate quantities of size $10^8$ have been computed in double precision and their difference is about 1. Rounding errors at the intermediate scale can be about $10^{-8}$. The difference may therefore have only about eight accurate decimal digits, even though the intermediate quantities have about sixteen.

### How to Examine a Calculation

Before introducing individual rounding-error variables:

1. **Estimate the sizes of the intermediate quantities.** A rounded quantity of size $M$ can carry an absolute error of order $uM$.
2. **Look for cancellation.** Compare the size of an addition or subtraction result with the quantities being combined. Check intermediate differences as well as the final answer.
3. **Compare the error scale with the desired result.** If the result is comparable to, or smaller than, the errors that can reach it, few or no correct digits may remain.

For a short calculation with accurate inputs and no cancellation, basic arithmetic generally preserves small relative errors, provided overflow and underflow are avoided. Multiplication and division by themselves do not have the cancellation mechanism: they scale the operands and their errors together. For example, $10^{-20}\times10^{-20}=10^{-40}$ is small but need not be inaccurate.

## A Small Term Lost Before Cancellation

In four-digit decimal arithmetic,

$$
\operatorname{fl}(10000+1)=10000.
$$

The absolute error is 1, but the relative error in this sum is only about $10^{-4}$. Now subtract 10000:

$$
\operatorname{fl}\bigl(\operatorname{fl}(10000+1)-10000\bigr)=0,
$$

although the exact answer is 1. The rounding error that was small relative to 10001 is as large as the final answer.

```{code-cell} ipython3
with localcontext() as ctx:
    ctx.prec = 4
    a = Decimal('10000')
    b = Decimal('1')
    c = Decimal('-10000')
    print("a + b:", a + b)
    print("(a + b) + c:", (a + b) + c)
    print("(a + c) + b:", (a + c) + b)
```

The two groupings give 0 and 1. In the second, the exactly represented values $10000$ and $-10000$ cancel before 1 is added. Floating-point addition is therefore not associative: changing the order can change the answer.

## Equivalent Formulas Can Give Different Accuracy

Consider

$$
P(x)=(x-1)^2=x^2-2x+1
$$

at $x=1+10^{-13}$. The exact answer is $10^{-26}$.

**Compare the scales first.** The expanded formula combines quantities of size 1 to produce a result of size $10^{-26}$. At sixteen-digit precision, errors near $10^{-16}$ in those quantities can overwhelm the answer. The factored formula first computes $x-1$, then squares that small difference.

Using sixteen significant decimal digits, the expanded calculation proceeds as follows:

$$
\begin{aligned}
x^2&=1+2\times10^{-13}+10^{-26},\\
\operatorname{fl}(x^2)&=1+2\times10^{-13},\\
\operatorname{fl}\bigl(\operatorname{fl}(x^2)-2x\bigr)&=-1,\\
\operatorname{fl}(-1+1)&=0.
\end{aligned}
$$

The small term is lost when $x^2$ is rounded. The final addition cancels the remaining terms and exposes that loss.

```{code-cell} ipython3
with localcontext() as ctx:
    ctx.prec = 16
    x = Decimal('1.0000000000001')
    expanded = x*x - 2*x + 1
    factored = (x - 1)**2
    print("Expanded form:", expanded)
    print("Factored form:", factored)
```

The factored form gives $10^{-26}$ exactly in this example. Here $x$ and 1 are exactly represented, so $x-1=10^{-13}$ is accurate despite the cancellation. Squaring it then preserves the answer. If $x$ itself were uncertain, the factored formula could not recover information already lost from the input.

The same issue becomes more pronounced for $(x-1)^{10}$. At $x=1.001$, the exact answer is $10^{-30}$, while the expanded polynomial contains terms of size up to a few hundred:

```{code-cell} ipython3
from math import comb

with localcontext() as ctx:
    ctx.prec = 8
    x = Decimal('1.001')
    expanded = Decimal(0)
    for k in range(11):
        expanded += Decimal((-1)**(10-k) * comb(10, k)) * x**k
    print("Expanded form:", expanded)
    print("Factored form:", (x - 1)**10)
```

These examples suggest a practical remedy: **rearrange a formula so that a small answer is computed directly, without first forming and rounding much larger quantities that must cancel.**

## Summation: A Bound That Captures the Main Idea

For a sum $S=\sum_{i=1}^n x_i$, the relevant comparison is between $|S|$ and $\sum_{i=1}^n|x_i|$. The latter measures the total size of the terms before cancellation.

````{prf:theorem} Error Bound for Summation
:label: thm:summation_error

Let $\widehat{S}$ be the result of adding $n\ge2$ floating-point numbers sequentially, starting with $x_1$. Assume rounding to nearest, no overflow or underflow, and $(n-1)u<1$. Then

$$
|\widehat{S}-S|
\le \gamma_{n-1}\sum_{i=1}^n|x_i|,
\qquad
\gamma_{n-1}=\frac{(n-1)u}{1-(n-1)u}.
$$

Here $S$ is the exact sum of the stored inputs. If $(n-1)u\ll1$, then $\gamma_{n-1}\approx(n-1)u$.
````

For $S\ne0$, dividing by $|S|$ gives

$$
\frac{|\widehat{S}-S|}{|S|}
\le \gamma_{n-1}
\frac{\sum_{i=1}^n|x_i|}{\left|\sum_{i=1}^n x_i\right|}.
$$

The two factors describe different effects. The factor $\gamma_{n-1}$ accounts for the number of additions. The ratio compares the sizes of the terms with the final answer:

- If all terms have the same sign, the ratio is 1. There is no cancellation, and the relative error is bounded by $\gamma_{n-1}$.
- If positive and negative terms nearly cancel, the ratio can be very large. Small rounding errors can then be large relative to the sum.

This is an upper bound; cancellation signals a risk, not a guarantee of a large error. It also shows why very long sums deserve care even without cancellation: errors can accumulate over many additions.

The same comparison applies to an inner product $x^Ty=\sum_i x_i y_i$: compare $\sum_i|x_i y_i|$ with $|x^Ty|$. Rounding in the products adds another source of error, but the question about scale is the same.

## Example: The Taylor Series for $e^{-20}$

Consider evaluating

$$
e^{-20}=\sum_{k=0}^{\infty}\frac{(-20)^k}{k!}.
$$

The answer is approximately $2.06\times10^{-9}$. Yet the largest terms, at $k=19$ and $k=20$, have magnitude about $4.31\times10^7$. We are combining numbers separated from the answer by about sixteen orders of magnitude.

With ten significant decimal digits, rounding a term of size $10^7$ can introduce an error of order $10^{-2}$. That is millions of times larger than the answer. We can recognize the danger without expanding a single rounding-error expression.

The sum of the magnitudes makes the comparison precise:

$$
\frac{\displaystyle\sum_{k=0}^{\infty}20^k/k!}{e^{-20}}
=\frac{e^{20}}{e^{-20}}
=e^{40}\approx2.35\times10^{17}.
$$

A better approach is to compute $e^{20}$ using its positive-term series and then take the reciprocal. There is no cancellation in that sum, and a small relative error in $e^{20}$ produces a comparably small relative error in its reciprocal.

```{code-cell} ipython3
def taylor_exp(x, degree=100):
    """Sum through x**degree / degree! at the current precision."""
    term = Decimal(1)
    total = term
    for k in range(1, degree + 1):
        term = term * x / Decimal(k)
        total += term
    return total

with localcontext() as ctx:
    ctx.prec = 10
    direct = taylor_exp(Decimal('-20'))
    reciprocal = Decimal(1) / taylor_exp(Decimal('20'))

with localcontext() as ctx:
    ctx.prec = 50
    reference = Decimal('-20').exp()
    print(f"Reference:         {reference:.9e}")
    print(f"Alternating sum:   {direct:.9e}")
    print(f"Reciprocal method: {reciprocal:.9e}")
    print(f"Relative error, alternating sum:   {abs((direct-reference)/reference):.2e}")
    print(f"Relative error, reciprocal method: {abs((reciprocal-reference)/reference):.2e}")
```

Both series are summed through $k=100$. In exact arithmetic, each omitted tail has magnitude less than $4\times10^{-29}$, so truncation does not explain the failure of the alternating sum. The problem is rounding at the scale of its large terms, followed by cancellation. In practice, use the exponential function provided by a numerical library; these series illustrate why the choice of formula matters.

## What to Look For

**Ask how large the quantities being combined are, how accurately they are known, and how small their sum or difference is.** If rounding errors at the intermediate scale are comparable to the desired answer, rewrite the calculation or use more precision. A large intermediate value alone is not evidence of inaccuracy; the danger is that its error survives after the large quantities cancel.

This comparison is the starting point for analyzing rounding errors in matrix algorithms. We will use it again when studying [LU factorization with pivoting](lu_pivoting.md).
