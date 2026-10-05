# Cholesky breakdown experiments

These runs test the matrix family in the optional section of `content/cholesky.md`.
The dense matrix was allocated and touched in full, including diagonal padding.
The C++ factorization skips zero contributions, keeps each effective update in its original order, and disables fused multiply-add.

All nine runs failed. Both small cases matched the notes’ NumPy construction and factorization bit for bit.
Every run matched the predicted stored Schur complement in every lower-triangular entry.
The rational value classes of the input entries were checked for exact representability.
The bounds below come from the proved formula; they are not computed condition numbers.

| Arithmetic | Size | Base size | Matrix GiB | Failing step | Pivot | Factor seconds |
|---|---:|---:|---:|---:|---:|---:|
| Binary64 | 576 | 576 | 0.0025 | 565 | -1.66255476e-12 | 0.001 |
| Binary32 | 2304 | 2304 | 0.0198 | 2252 | -1.14055234e-03 | 0.029 |
| Binary64 | 9216 | 9216 | 0.6328 | 9002 | -5.04262592e-11 | 5.379 |
| Binary64 | 12288 | 9216 | 1.1250 | 9002 | -5.04262592e-11 | 5.723 |
| Binary64 | 14336 | 9216 | 1.5312 | 9002 | -5.04262592e-11 | 5.707 |
| Binary64 | 16384 | 9216 | 2.0000 | 9002 | -5.04262592e-11 | 6.310 |
| Binary64 | 17408 | 9216 | 2.2578 | 9002 | -5.04262592e-11 | 6.205 |
| Binary64 | 18432 | 9216 | 2.5312 | 9002 | -5.04262592e-11 | 6.678 |
| Binary64 | 19456 | 9216 | 2.8203 | 9002 | -5.04262592e-11 | 6.726 |

For the 576, 2304, and 9216 base matrices, the respective condition-number upper bounds are
3.531 × 10¹³, 8239, and 5.517 × 10¹¹, rounded upward.
All larger tested matrices are the 9216 matrix padded with a 4I block. Their condition number and failing pivot are unchanged.

The machine has 32 GiB of RAM. The largest completed matrix, of size 19456, occupied 2.8203125 GiB.
The initial sweep imposed a 3 GiB allocation cap. Separate attempts to raise that cap and allocate sizes 20480 and 23040 were rejected by the available-memory check.
The 20480 attempt needed 3.125 GiB for the matrix alone, while only 4.78 GiB was estimated available. The check reserves 2 GiB for other work plus 128 MiB of process overhead.
No allocation failure was deliberately triggered. The largest completed size is a limit of these runs under current memory availability, not the maximum possible on an otherwise idle machine.

The next admissible unpadded binary64 matrix has size 147456 and requires 162 GiB in full dense storage. It was not allocated.
Padding does not test the next construction parameter or improve the condition number.

Timing covers the compiled factorization and the Schur-complement comparison, excluding allocation and NumPy reference validation.
Peak process RSS in `results.json` is cumulative within the initial sweep, so it is not a per-case matrix-storage measurement.
Swap counters are system-wide; changes cannot be attributed solely to this experiment. No swap-out increments were observed during completed cases.

## Reproduce

On macOS, from the repository root, with NumPy, psutil, and clang++ available:

```sh
.venv/bin/python experiments/cholesky/run.py --output /tmp/cholesky-results.json
```

To run one larger case in a fresh process while retaining the available-memory guard:

```sh
.venv/bin/python experiments/cholesky/run.py --case 64 5 20480 --max-gib 24 --reserve-gib 2 --output /tmp/cholesky-20480.json
```

Each candidate is checked against both the allocation cap and estimated available RAM before allocation.
The output records skipped cases as well as completed ones. Increasing the allocation cap does not disable the reserve check.
