# Upper certificate

This directory contains the complete certificate for

$$
K_G^{\mathbb R}<1779893/1000000=1.779893.
$$

From the repository root, after installing the dependencies as described in the
[root README](../README.md), run:

```sh
python verify.py --upper
```

The root command runs the complete certificate and the additional scalar check
identified in the paper. It writes fresh outputs under `runs/`. For the preserved
standalone certificate driver alone, use:

```sh
python upper/verify.py --output runs/upper
```

The standalone driver creates a new `run-*` directory and accepts only after all
12 generation stages, coverage checks, and the strict endpoint comparison pass.
A successful run ends with:

```text
FRESH_COMPLETE_UPPER_CERTIFICATE_VERIFIED=1
CERTIFIED_KG_UPPER=1779893/1000000
RESULT=.../VERIFIED.json
```

Interrupted or failed runs leave logs without a `VERIFIED.json` acceptance file.
Run ordinary Python without `-O`; assertions are part of the verification.

## Files

| Path | Purpose |
| --- | --- |
| `verify.py` | Checks hashes, copies the 18 frozen source/input files to a fresh working tree, regenerates the recipe and all certificate stages, and checks the final result. |
| `source/certification/reproduce.py` | Executes the 12 mathematical generation stages at degree 681. |
| `source/certification/upper_head_681/` | Certified coefficient matrix, contraction, and explicit error budget. |
| `source/certification/weighted_staircase/` | Rounding recipe, staircase head, and nonadditive tails. |
| `source/certification/degree681/` | Linear-tail bounds. |
| `source/certification/upper_tail/` | Continuation, factorized energies, and heavy additive tail. |
| `source/certification/upper_head/arc/` | Complex-arc interval routine used by the tail bounds. |
| `source/certification/compose_head.py` | Assembly of the certified head. |
| `source/certification/certify_rare_additive_tails.py` | Rare additive tails. |
| `source/certification/certify_upper.py` | Exact final assembly and endpoint comparison. |
| `source/upper/refined_explicit_mixture.json` | Exact rounding mixture. |
| `source/competitors/staircase_3_5.json` | Exact staircase input. |
| `requirements.txt` | Exact public dependency versions. |
| `MANIFEST.sha256` | Hashes of all distributed files in this directory except the manifest itself. |

The mathematical source, both inputs, and `verify.py` are unchanged from the
packaged complete upper certificate. Internal source paths are preserved. No
stored success flag or generated numerical array is needed for verification.
See [PAPER_MAP.md](../PAPER_MAP.md) for the correspondence with the current paper.

## Arithmetic and resources

Use Python 3.12 with NumPy 2.3.5 and python-flint 0.9.0. The driver enforces the
exact dependency versions, IEEE-754 binary64, and round-to-nearest when the
platform exposes its rounding mode. The computation combines directed Arb
intervals with explicit bounds for binary64 operations. It uses the usual trust
in these libraries, Python, and the hardware; it is not a proof-assistant
formalization. After dependency installation the computation is offline.

Allow approximately 12 GB of available memory; 16 GB or more is recommended.
Runtime depends on the machine. The original package records a **historical**
268-second full replay on 8 September 2026 with macOS ARM64, Python 3.12.14,
NumPy 2.3.5, and python-flint 0.9.0. That report is provenance, not a substitute
for a fresh run of this repository. Generated work files can be large and should
remain outside `source/`.
