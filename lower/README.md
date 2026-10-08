# Lower-bound certificate

This directory contains the numerical certificate for
`K_G > 1773/1000`.
The executable mathematical sources and finite witness inputs are unchanged;
`check-tail.py` is the additional rational rotation-weight check supplied in
that archive. Source hashes and packaging changes are recorded in
`SOURCE.json` and `MANIFEST.sha256`.

## Run

Install the dependencies following the [repository README](../README.md).
From the repository root, run:

```sh
python verify.py --lower --workers 4
```

The root command checks file hashes, uses a fresh working copy, checks the
70 exact rational rotation-weight inequalities, and invokes this directory's
unchanged `verify.py`. The verifier builds the arithmetic backends, regenerates
all 35 universal fibre bounds and 1,820 non-fibre matrix checks, runs the
negative controls, and performs the final exact comparison. The certified
norm target is `1000/1773 - 1/1000000000`.

The lower calculation requires Python without `-O`, NumPy, Cython, setuptools,
a C++17 compiler with signed 128-bit integers, GMP development headers, and
MPFR. The pinned Python versions are in [requirements.txt](requirements.txt).
It builds its own MPFR interface named `flint`; an installed Python-FLINT
package does not replace this backend. See [ARITHMETIC.md](ARITHMETIC.md) for
the arithmetic and acceptance conditions.

## Files

| Path | Purpose |
| --- | --- |
| `verify.py` | Complete lower-bound verification driver. |
| `check-tail.py` | Standard-library check of 70 rational rotation-weight inequalities; the root driver supplies `data/lower_rows.json`. |
| `backend/` | Directed MPFR endpoint arithmetic and its Cython interface. |
| `tools/analytic.py` | Regeneration of interpolation, Gaussian-tail and rounding bounds. |
| `tools/exact_graph.cpp`, `tools/graph.py` | Exact integer graph calculations. |
| `tools/replay.py` | Verification of every finite proof tree and its analytic input. |
| `tools/assemble_final.py` | Exact assembly and final strict comparison. |
| `tools/test_*.py`, `backend/test_backend.py` | Arithmetic tests and rejection of invalid evidence. |
| `tools/frozen_row_limit.py` | Supplementary limitation of this fixed row family; not an upper bound on `K_G`. |
| `baseline/` | Frozen rational row data and original non-fibre checker; paths are retained because the verifier uses them. |
| `inputs/` | Fixed-point and analytic inputs regenerated and checked during verification. |
| `certificates/` | The 35 finite proof trees. |
| `recorded_replay/` | Historical output for comparison, never accepted in place of a fresh calculation. |

## Historical outputs

Files in `recorded_replay/` were included in the source archive. Their status
fields describe that recorded run, not a new verification of this repository.
The historical report retains its original filenames and path references.
Fresh outputs are written under the run directory chosen by the root driver.
The historical `FROZEN_ROW_LIMIT.json` concerns only refinement of this fixed
family of auxiliary inequalities and does not bound the Grothendieck constant
from above.
