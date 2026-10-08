# Grothendieck constant: numerical certificates

Companion code for *Primal/dual method for the Grothendieck constant*, by Steven Heilman, Chris Jones and Giulio Malavolta.  
https://arxiv.org/abs/2610.10477

This repository contains independently runnable numerical certificates for

$$
K_G \ge \frac{1773000000000}{999999998227} > 1.773,
\qquad
K_G < \frac{1779893}{1000000}=1.779893.
$$

Start with [the paper-to-code map](PAPER_MAP.md) for the correspondence between mathematical statements, source routines, arithmetic precision and acceptance conditions.

## Install

Use Python 3.12 or a compatible newer version. From the repository root, create a virtual environment and install the pinned dependencies:

```sh
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

The lower verifier additionally requires a C++17 compiler with signed 128-bit integers, GMP development headers, and the GMP/MPFR libraries available to the compiler and linker. Install these using your system's package manager before running the lower bound. Its MPFR interface is built during verification; the upper bound uses the installed Python-FLINT package. See the [lower](lower/README.md) and [upper](upper/README.md) documentation for their arithmetic details.

Run ordinary Python **without `-O`**: assertions are part of verification. The calculation runs offline after dependencies have been installed.

## Verify

A quick check of hashes, exact input identities, interval coverage, normalization and all 70 rational rotation-tail inequalities requires only the Python standard library:

```sh
python verify.py --data-only
```

It prints `DATA_LINKS_CHECKED=1` and `NO_ANALYTIC_NUMERICAL_REPLAY=1`. This mode does not certify either numerical bound.

To regenerate a complete certificate, run one of:

```sh
python verify.py --lower --workers 4
python verify.py --upper
python verify.py --all --workers 4
```

With no mode specified, the driver runs both sides. `--workers` controls the lower verifier's parallel row computations. The lower command includes the exact tail check; the upper command includes the additional 256-bit Fourier-tail scalar inequality printed in the paper.

Each invocation creates a new directory under `runs/` and prints its location. Use `--output /path/to/results` to choose another parent. Source and input directories are not used for generated output.

The upper calculation historically required about 12 GB of available memory; 16 GB or more is recommended. Runtime depends on the machine. Lower-row computations also use additional memory as worker count increases. The [validation record](validation/README.md) states exactly which checks were freshly rerun while preparing this repository; its [machine-readable summary](validation/summary.json) is evidence of those runs, not a substitute for your own execution.

The wrapper's fast regression checks can also be run with the standard library:

```sh
python -m unittest discover -s tests -v
```

These test relocation, missing or altered inputs, unexpected source files, disabled assertions, output locations and 70 invalid tail-weight examples. They do not replay the numerical proofs.

## What counts as acceptance

Only a completed numerical run produces `runs/complete-*/ACCEPTED.json`. Success ends with:

```text
NUMERICAL_CERTIFICATES_ACCEPTED=lower,upper
ACCEPTANCE_FILE=.../ACCEPTED.json
```

For a single-side invocation, the first line lists only that side. Check `certificates_replayed` in the JSON: both bounds have been verified in that run only when it contains both `lower` and `upper` and `both_certificates_and_endpoint_checks_passed` is `true`. Failed or interrupted runs leave diagnostic logs without a root acceptance file.

The fresh run contains:

| File or directory | Meaning |
| --- | --- |
| `DATA_CHECKS.json`, `TAIL_CHECK.json` | Exact data and rotation-tail checks. |
| `lower_results/VERIFICATION_REPORT.json` | All 35 lower rows, regenerated analytic bounds, finite proof checks and exact assembly. |
| `upper_results/run-*/VERIFIED.json` | Fresh completion of all 12 upper stages and the strict endpoint comparison. |
| `UPPER_PAPER_CHECK.json` | The additional 256-bit Fourier-tail scalar check, when the upper side is requested. |
| `lower.log`, `upper.log` | Driver output for the requested sides; further logs reside in their result directories. |
| `ACCEPTED.json` | Final root acceptance, including the manuscript's printed numerical comparisons and source-manifest identity. |

The lower proof uses outward-rounded 384-bit MPFR intervals and checked integer/rational calculations. The upper proof uses directed Arb intervals at the precisions specified in the paper, together with explicit error bounds for binary64 operations. The programs verify the finite numerical hypotheses of the paper's analytic reductions; they are not a proof-assistant formalization.

## Repository guide

| Path | Contents |
| --- | --- |
| [PAPER_MAP.md](PAPER_MAP.md) | Stable manuscript labels mapped to source, inputs, reports and precise inequalities. |
| [verify.py](verify.py), [data_link.py](data_link.py) | Unified orchestration, file-identity checks and exact linkage to the manuscript data. |
| [data/](data/) | Readable exact mathematical definitions for both bounds. |
| [lower/](lower/) | Lower numerical source, arithmetic backend, frozen inputs and 35 finite proof trees. |
| [upper/](upper/) | Upper numerical source and exact rounding inputs. |
| [rounding/](rounding/) | Description of the rounding functions output by the primal/dual search. |
| [MANIFEST.sha256](MANIFEST.sha256) | Frozen file hashes for code, mathematical inputs and component package files used by the root driver; root documentation and validation records are outside this manifest. |
| [PROVENANCE.json](PROVENANCE.json) | Source-archive identities, manuscript snapshot and packaging changes. |
| [validation/](validation/) | Records of preparation-time checks, with their scope stated explicitly. |
| `runs/` | Locally generated output, excluded from Git. |

The files in `lower/certificates/` are finite witnesses: the verifier checks their branches, graph bounds and coverage. They are proof inputs. Files in `lower/recorded_replay/` are historical output, retained for comparison; stored success flags never replace recomputation. Fresh results belong under `runs/`.
