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

A floating-point number is the computer's way of representing a real number using a binary version of scientific notation. It balances the need to store numbers of vast a range (from the infinitesimally small to the astronomically large) with a fixed, finite number of bits.

This system is composed of three parts: a sign, a significand (the significant digits), and an exponent.

## The Big Idea: Scientific Notation for Computers 🖥️

You're likely familiar with scientific notation in base 10. For instance, we can write the Avogadro constant as approximately $6.022 \times 10^{23}$. This representation has three key pieces:

* **Sign:** Positive (+)
* **Significand:** 6.022 (the core digits)
* **Exponent:** 23 (how far to shift the decimal point)

Floating-point numbers apply this exact same concept, but in **base 2 (binary)**, the native language of computers. By storing numbers in this format, the binary point can "float" to whatever position is needed, controlled by the exponent.

## The Three Parts of a Floating-Point Number

A binary floating-point number is formally defined as:

$$\pm (1 + \sum_{i=1}^{p-1} d_i 2^{-i} )\; 2^e$$

Let's break this down into its three components.

1.  **The Sign (S):** This is the simplest part, represented by a single bit. It's 0 for positive (+) and 1 for negative (-).

2.  **The Exponent (E):** This component, $e$, determines the number's magnitude (its size). A large positive exponent means a large number, and a large negative exponent means a very small number close to zero.

3.  **The Significand (or Mantissa):** This part contains the number's significant digits, which determines its **precision**. It's represented by the term $(1 + \sum_{i=1}^{p-1} d_i 2^{-i} )$.
    * **The Hidden Bit:** In binary scientific notation, any non-zero number can be written with a leading "1" before the binary point (e.g., $1.01101 \times 2^5$). Since this leading digit is *always* 1, there's no need to waste a bit storing it. This is known as the "hidden bit." The `1 + ...` in the formula reflects this.
    * **The Fraction:** The summation term, $\sum d_i 2^{-i}$, represents the fractional part of the significand, which is what's actually stored in memory. The $d_i$ are the bits (0 or 1).


### An Example Decoded: 3.140625

Let's trace how the number `3.140625` is stored.

**1. Sign:** The number is positive, so the sign bit is **0 (+)**.

**2. Find the Exponent:** We locate the number between powers of 2.

$$2^1 \le 3.140625 < 2^2$$

Since it falls in the range of $2^1$, the exponent $e$ is **1**.

**3. Find the Significand:** To find the significand, we scale the number down by its exponent power:

$$\frac{3.140625}{2^1} = 1.5703125$$

So, the significand is **1.5703125**.

**4. Decompose the Significand:** Now, we break the significand into the `1 + fraction` format required by the formula:

$$1.5703125 = 1 + 0.5703125$$

The fractional part, `0.5703125`, is what we need to represent as a sum of negative powers of 2.

$$0.5703125 = 0.5 + 0.0625 + 0.0078125 = \frac{1}{2} + \frac{1}{16} + \frac{1}{128} = 2^{-1} + 2^{-4} + 2^{-7}$$

This tells us that the bits for the fractional part are $d_1=1, d_4=1, d_7=1$, and the others in between are 0.

Putting it all together, the number `3.140625` is represented by these three components.

## The Inevitable Trade-Off: Range vs. Precision

This system is an elegant compromise. With a fixed number of bits (e.g., 32 or 64), you must decide how many to allocate to the exponent versus the significand.
* More bits for the **exponent** gives you a larger **range** of numbers.
* More bits for the **significand** gives you higher **precision**.

You can't have unlimited of both. This limitation is why computers must round numbers, leading to the small but potentially consequential errors that are central to the study of numerical stability.

## The Fundamental Rule: Round to Nearest 🎯

Floating-point arithmetic is the set of rules computers use to perform calculations on real numbers. Because computers have finite memory, they can't store numbers with infinite precision. Consequently, every operation can introduce a tiny **rounding error**. While these errors are small, they cause floating-point arithmetic to behave differently from the perfect arithmetic you learn in math class.

The core principle of floating-point arithmetic is that after every operation (like addition or multiplication), if the true mathematical result isn't a number that can be stored exactly, the computer must **round it to the nearest representable floating-point number**.

This can be expressed with a simple formula:

$$fl(a \text{ op } b) = (a \text{ op } b)(1+\epsilon)$$

Here:

* `a op b` is the exact mathematical result.
* `fl(...)` is the final result the computer stores.
* $\epsilon$ (`epsilon`) is the small **relative error** introduced by rounding.

The maximum possible size of this error for any single operation is called the **unit roundoff**, denoted by $u$. For standard 64-bit "double precision" numbers, this value is incredibly small, around $10^{-16}$, which means you can expect about 15-16 correct decimal digits for any single calculation.

## The Surprising Consequences of Rounding

These small, seemingly harmless rounding errors lead to several counter-intuitive results.

### 1\. Equality Is Not What It Seems (`0.1 + 0.2 ≠ 0.3`)

One of the most famous quirks is that in nearly all programming languages, `0.1 + 0.2` is not equal to `0.3`. This is because numbers like 0.1 and 0.2, which are simple in base 10, are infinitely repeating fractions in base 2, much like 1/3 is repeating in base 10 (0.333...).

* The computer stores the closest binary approximation of 0.1 and 0.2.
* It adds these two approximations together.
* The result is a number that is extremely close to, but not *exactly* the same as, the computer's approximation for 0.3.

Let's see this in action with Python.

```{code-cell} ipython3
a = 0.1
b = 0.2
c = 0.3

# The sum looks correct when printed with default precision
sum_ab = a + b
print(f"The sum of {a} and {b} is: {sum_ab}")
print(f"Is the sum equal to {c}? {sum_ab == c}")

# To see the true stored values, let's print with more digits
print("\n--- Print more digits ---")
print(f"Value of sum_ab: {sum_ab:.20f}")
print(f"Value of c:      {c:.20f}")
```

The output clearly shows that the stored value for `0.1 + 0.2` is a number slightly larger than 0.3, while the stored value for `0.3` is slightly smaller.

This is why you should **never use a direct equality check `==` with floating-point numbers.** Instead, you should check if they are close enough: `abs(a - b) < tolerance`.

### 2\. The Rules of Arithmetic Are Bent

While floating-point addition is **commutative** (`a + b` is the same as `b + a`), it is not always **associative**. This means the order of operations can change the final answer.

$$(a + b) + c \neq a + (b + c)$$

This happens because the rounding error at each step depends on the magnitude of the numbers involved. Adding two small numbers together first and then adding a large one can produce a different result than adding a small number to a large one, where its precision might be lost.

### 3\. Numbers Are Not Evenly Spaced

Floating-point numbers are not distributed uniformly on the real number line. They are densest around zero and become progressively sparser as their magnitude increases.

Think of it like a ruler where the tick marks get farther apart the further you move from zero. The gap between 1 and the next representable number is tiny (this gap defines the unit roundoff $u$). However, the gap between 1,000,000,000 and the next representable number is much larger. This means you have far more precision for small numbers than for very large ones.

## Numerical Explorations: Understanding Roundoff Errors

The behavior of floating-point arithmetic is fundamental to everything we do in this course. An algorithm that is perfectly sound in exact arithmetic can produce meaningless results if it's unstable in the face of the small, inevitable errors introduced by computation.

To explore this, let's use Python's `Decimal` library. This is a powerful tool because it allows us to simulate a computer with a precision that we can define, letting us see these errors in a controlled and exaggerated way.

### Modeling Finite Precision

Let's begin by configuring our environment to operate with a very low precision of just four significant digits. This will make the effects of roundoff error obvious.

```{code-cell} ipython3
from decimal import Decimal, getcontext
# Set the precision to 4 significant digits
getcontext().prec = 4
```

Every number and every result of an operation will now be rounded to fit within these four digits.

The first type of error is **representation error**. Most real numbers cannot be stored perfectly. For example, even with our low precision, let's see what happens to $\pi$.

```{code-cell} ipython3
# The true value of pi cannot be stored, so it's rounded.
pi_approx = +Decimal(3.141592653589793)
# + needed to apply the context precision
print(f"Approximation of pi: {pi_approx}")
```

The value is immediately rounded to `3.142`. This initial error is the starting point for all subsequent numerical inaccuracies.

### The Dangers of Addition and Subtraction

Addition and subtraction are the primary sources of significant roundoff error, but they manifest in different ways depending on the numbers involved.

#### **Case 1: Adding Numbers of Drastically Different Magnitudes**

When you add a very large number to a very small one, the smaller number's information can be completely lost. This is called **swamping**.

Consider adding a large number to our `pi_approx`.

```{code-cell} ipython3
large_number = Decimal('1234')
small_number = Decimal('3.142') # Our 4-digit pi

# The true sum is 1237.142
result = large_number + small_number

print(f"'{large_number}' + '{small_number}' = '{result}'")
```

To perform the addition, the computer must align the decimal points. The sum `1237.142` requires seven significant digits to store. Our machine, with a precision of four, must round it, resulting in `1237`. The entire contribution of $\pi$ has vanished.

#### **Case 2: Subtracting Nearly Equal Numbers**

The most insidious error is **catastrophic cancellation**, which occurs when subtracting two numbers that are very close to each other. The operation itself may be exact, but the result loses a massive amount of relative precision.

Let's define two numbers that differ only in their fourth significant digit.

```{code-cell} ipython3
x = Decimal('1.234')
y = Decimal('1.233')

# The subtraction itself is simple and the result is exact.
difference = x - y

print(f"'{x}' - '{y}' = '{difference}'")
print(f"Exponent of the result's leading digit: {difference.adjusted()}")
```

Observe what has happened. The `.adjusted()` method returns the exponent of the number's most significant digit, which tells us its magnitude. For our inputs `x` and `y`, this value would be `0` (since their leading digit is in the $10^0$ place). However, for the `difference`, the exponent is `-3`. This dramatic drop in magnitude is the signature of catastrophic cancellation.

We started with two numbers, `x` and `y`, each containing four significant digits of information. The leading digits `1.23` cancelled each other out, leaving a result of `0.001`.

While this result is stored exactly as `1.000E-3` in our 4-digit system, we have gone from four figures of precision to just **one**. The trailing zeros in `1.000` are meaningless; they are artifacts of the arithmetic, not real information from our original numbers. We have lost 75% of our information in a single operation. This is why it is called "catastrophic."

### The Relative Safety of Multiplication

In contrast, multiplication is a far more benign operation, especially when dealing with numbers of different magnitudes. Where addition would "swamp" a small number, multiplication preserves its contribution.

Let's use the same large and small numbers from our addition example to see the difference.

```{code-cell} ipython3
large_number = Decimal('1234')
small_number = Decimal('3.142') # Our 4-digit pi

# The true product is 3877.228
product = large_number * small_number

print(f"'{large_number}' * '{small_number}' = '{product}'")
```

The exact answer, `3877.228`, was rounded to `3877`. Notice the crucial difference: the information from the `small_number` was not lost. Its value was essential in determining the final result. We began with two numbers, each with four significant digits, and our result still has four meaningful significant digits. A small amount of information was lost to rounding, but there was no catastrophic loss of precision.

Understanding these different error modes is critical. A stable algorithm must be structured to avoid subtractions of nearly equal quantities whenever possible, as this is the primary mechanism by which numerical error grows uncontrollably.


## The Accumulation of Error in Summation 📈

Let's build on our understanding of roundoff error. A single rounding error is typically harmless, but the central problem in numerical computation is how these small errors **accumulate** over the course of millions or billions of operations. In some cases, this accumulation can lead to a complete loss of accuracy in the final result.

A classic scenario where this occurs is the summation of a long sequence of numbers.

Consider the computation of a sum, $S_n = \sum_{i=1}^n x_i$. In exact arithmetic, the order of summation does not matter. In floating-point arithmetic, it does. The standard algorithm computes a sequence of partial sums:

$\tilde{S}_1 = x_1$

$\tilde{S}_2 = fl(\tilde{S}_1 + x_2)$

...

$\tilde{S}_k = fl(\tilde{S}_{k-1} + x_k)$

Each step introduces a small roundoff error. A rigorous analysis shows that the error in the final computed sum, $\Delta S_n = |S_n - \tilde{S}_n|$, is bounded by approximately:

$$|\Delta S_n| \le n u \sum_{i=1}^n |x_i| + O(u^2)$$

Here, $u$ is the unit roundoff of the machine. Let us dissect this formula. It tells us that the absolute error of the computed sum depends on:

1.  **$n$**: The number of terms. The error grows, at worst, linearly with the length of the sum.
2.  **$u$**: The machine precision. A more precise machine reduces the error.
3.  **$\sum_{i=1}^n |x_i|$**: The sum of the *magnitudes* of the terms.

This last point is the most critical. The error bound does not depend on the magnitude of the final sum, $|S_n|$, but on the sum of the absolute values of the terms that created it. This implies that if we compute a sum where the intermediate partial sums are large but the final answer is small (due to cancellation), the **relative error** can be enormous.

### Pathological Example: The Taylor Series for $e^{-x}$

A perfect illustration of this principle is the evaluation of the Taylor series for $e^{-x}$ for a large, positive $x$. The series is an alternating sum:

$$e^{-x} = \sum_{k=0}^{\infty} \frac{(-x)^k}{k!} = 1 - x + \frac{x^2}{2!} - \frac{x^3}{3!} + \dots$$

Let's attempt to calculate $e^{-20}$. The true answer is very small, approximately $2.06 \times 10^{-9}$. However, the individual terms in the series first become enormous before decaying. For instance, the term for $k=19$ is $\frac{20^{19}}{19!} \approx 4.3 \times 10^{8}$.

We are adding and subtracting gigantic numbers to produce a tiny one. This is the exact scenario where our error bound predicts a catastrophic loss of precision.

Let's observe this using a simulated machine with 10-digit precision.

```{code-cell} ipython3
import math
# Set precision to 10 significant digits
getcontext().prec = 10

def exp_taylor_unstable(x):
    """Computes e^(-x) using the standard Taylor series."""
    x = Decimal(x)
    term = Decimal(1)
    current_sum = Decimal(1)
    k = 1
    # Sum until the terms are too small to matter
    while abs(term) > 1e-15:
        term = term * (-x / Decimal(k))
        current_sum += term
        k += 1
        # Stop if terms get huge to prevent infinite loop
        if k > 200: break
    return current_sum

# The value to test
x_val = 20

# Calculate with our unstable algorithm
computed_val = exp_taylor_unstable(x_val)

# For comparison, calculate a "true" value with high precision
getcontext().prec = 50
true_val = Decimal(-x_val).exp()

print(f"Precision: {10} digits")
print(f"Computed e^({-x_val}):           {computed_val}")
print(f"Extended precision e^({-x_val}): {true_val:.10e}")
print(f"Math library e^({-x_val}):       {math.exp(-x_val):.10e}")
```

The result is completely wrong. This is a direct consequence of the catastrophic cancellation among the large intermediate terms, where the accumulated roundoff error is larger than the final true sum.

### A More Subtle Case: Evaluating Polynomials

This same issue arises when evaluating a polynomial in its standard form, $P(x) = \sum a_i x^i$, particularly near a root.

Consider the polynomial $P(x) = (x-1)^{10}$. Clearly, $P(1)=0$. Let's evaluate it at $x=1.001$.

The numerically stable way is to compute $(1.001-1)^{10} = (0.001)^{10} = 10^{-30}$.

The unstable way is to first expand the polynomial via the binomial theorem and then sum the resulting terms. The expanded form is $P(x) = x^{10} - 10x^9 + 45x^8 - \dots + 1$. When $x$ is close to 1, this involves summing large, alternating terms.

Let's observe the failure with 8-digit precision.

```{code-cell} ipython3
from math import comb
# Set precision to 8 significant digits
getcontext().prec = 8

# The value to test, close to the root at x=1
x = Decimal('1.001')
n = 10

# Method 1: The stable calculation
stable_result = (x - 1)**n

# Method 2: The unstable summation of the expanded form
unstable_result = Decimal(0)
for k in range(n + 1):
    # Binomial expansion: C(n,k) * x^k * (-1)^(n-k)
    term = Decimal(comb(n, k)) * (x**k) * ((-1)**(n - k))
    unstable_result += term

print(f"Precision: {8} digits")
print(f"Stable calculation (x-1)^10:   {stable_result:e}")
print(f"Unstable sum of terms:        {unstable_result:e}")
```

The result from the expanded sum is incorrect. Again, this is a classic case of catastrophic cancellation. The lesson here is profound: two formulas that are mathematically identical can have vastly different numerical properties. A key task for a numerical analyst is to identify and implement the stable formulation.

### Even Simpler Case: A Quadratic Polynomial

The issue is not limited to high-degree polynomials. It can appear even in a simple quadratic. Consider the two equivalent forms of the polynomial $P(x) = (x-1)^2 = x^2 - 2x + 1$.

Let's evaluate this for a value of $x$ extremely close to 1, using a precision of 16 digits, which is typical for standard 64-bit "double" floating-point numbers. The true answer should be $((1+10^{-13})-1)^2 = (10^{-13})^2 = 10^{-26}$.

```{code-cell} ipython3
# Set precision to 16 significant digits (like double precision)
getcontext().prec = 16

# Define x very close to 1
x = Decimal('1') + Decimal('1e-13')

# Method 1: The stable, factored form
stable_result = (x - 1)**2

# Method 2: The unstable, expanded form
unstable_result = x**2 - 2*x + 1

print(f"x = {x}")
print("-" * 40)
print(f"True Answer:                  1.000000000000000e-26")
print(f"Stable calculation (x-1)^2:   {stable_result:e}")
print(f"Unstable sum x^2 - 2x + 1:    {unstable_result:e}")
```

The stable form gives the correct answer, but the expanded form is completely wrong! It's equal to `0` in our 16-digit precision system.

The failure happens because the calculation of `x**2 - 2*x` involves subtracting two numbers that differ by exactly 1:

*   `x**2` is `1 + 2e-13`.
*   `2*x` is `2 + 2e-13`.

When `x**2` is computed, the term $(10^{-13})^2 = 10^{-26}$ is too small to be represented in the 16 available digits of `1.000...` and is lost to rounding. The subsequent subtraction `(x**2 - 2*x)` then suffers from catastrophic cancellation, destroying the accuracy of the result. This simple example powerfully demonstrates that the mathematical form of an expression is critically important in numerical computation.