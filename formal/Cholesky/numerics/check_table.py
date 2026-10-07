"""Numerical cross-check of the table "Numerical Examples" in content/cholesky.md.

This is not part of the Lean proof. It reruns the routine `cholesky_in_place` of the
notes on the counterexamples and checks:

* the failing step and computed pivot of each row of the table;
* that the trailing block after the 7k leading steps equals the matrix S-hat
  (this is proved in Lean as `CE.Par.trailing_eq_Shat`);
* that accumulating each inner product from zero and then subtracting it succeeds.

The 9216 x 9216 binary64 row is checked by factoring S-hat (2048 x 2048) only. Lean
proves that the leading 7k steps produce S-hat exactly, see `CE.Par.trailing_eq_Shat`,
and `run_natAdd` shows that the remaining steps act on S-hat alone.
"""
import math
from fractions import Fraction

import numpy as np


def cholesky_in_place(A):
    n = A.shape[0]
    for k in range(n):
        if not np.isfinite(A[k, k]) or A[k, k] <= 0:
            return k + 1, A[k, k]
        A[k, k] = np.sqrt(A[k, k])
        A[k + 1:, k] /= A[k, k]
        for j in range(k + 1, n):
            A[j:, j] -= A[j:, k] * A[j, k]
    return None, None


def cholesky_accumulate(A):
    """Accumulate each inner product from zero, then subtract it."""
    n = A.shape[0]
    L = np.zeros_like(A)
    for j in range(n):
        s = np.zeros(n - j, dtype=A.dtype)
        for k in range(j):
            s += L[j:, k] * L[j, k]
        d = A[j, j] - s[0]
        if not d > 0:
            return j + 1
        L[j, j] = np.sqrt(d)
        L[j + 1:, j] = (A[j + 1:, j] - s[1:]) / L[j, j]
    return None


def breakdown_matrix(q, dtype):
    p = np.finfo(dtype).nmant + 1
    u = 2.0**-p
    k = 4**q
    H = np.ones((1, 1), dtype=dtype)
    while H.shape[0] < k:
        H = np.block([[H, H], [H, -H]])
    Q = H / dtype(2**q)
    h = 2 * u / 2**q
    eta = 6 * k * h
    delta = k * eta
    t = math.ldexp(3.0, -(p + q + 2) // 2)
    assert Fraction(t)**2 == Fraction(9, 8) * Fraction(h)
    I = np.eye(k, dtype=dtype)
    T = np.vstack((np.full((6 * k, 2 * k), t, dtype=dtype),
                   dtype(1.5) * np.hstack((I, Q))))
    G = dtype(2.25) * np.block([[I, Q], [Q.T, I]])
    C = G + dtype(eta) * np.ones((2 * k, 2 * k), dtype=dtype)
    C += dtype(delta) * np.eye(2 * k, dtype=dtype)
    return np.block([[dtype(4) * np.eye(7 * k, dtype=dtype), dtype(2) * T],
                     [dtype(2) * T.T, C]])


def shat(q, dtype):
    p = np.finfo(dtype).nmant + 1
    u = 2.0**-p
    k = 4**q
    eta = 12 * 2**q * u
    a = (k + 9 / 8) * eta
    J = np.ones((k, k), dtype=dtype)
    B = dtype(a) * np.eye(k, dtype=dtype) - dtype(eta / 8) * J
    return np.block([[B, dtype(-eta) * J], [dtype(-eta) * J, B]])


def leading_steps(A, m):
    n = A.shape[0]
    for k in range(m):
        A[k, k] = np.sqrt(A[k, k])
        A[k + 1:, k] /= A[k, k]
        for j in range(k + 1, n):
            A[j:, j] -= A[j:, k] * A[j, k]
    return A


if __name__ == "__main__":
    for q, dtype, step, pivot in [(3, np.float64, 565, -1.66255e-12),
                                  (4, np.float32, 2252, -1.14055e-3)]:
        A = breakdown_matrix(q, dtype)
        k = 4**q
        W = leading_steps(A.copy(), 7 * k)
        same = np.array_equal(np.tril(W[7 * k:, 7 * k:]), np.tril(shat(q, dtype)))
        idx, piv = cholesky_in_place(A.copy())
        acc = cholesky_accumulate(A)
        print(f"n={9 * k} {dtype.__name__}: trailing block == S-hat: {same}; "
              f"fails at step {idx} (table {step}), pivot {piv:.5e} (table {pivot:.5e}); "
              f"accumulate-then-subtract: {'succeeds' if acc is None else 'fails'}")
    q, dtype = 5, np.float64
    k = 4**q
    idx, piv = cholesky_in_place(shat(q, dtype))
    print(f"n={9 * k} {dtype.__name__}: S-hat fails at its step {idx}, i.e. step {7 * k + idx} "
          f"(table 9002), pivot {piv:.5e} (table -5.04263e-11)")
