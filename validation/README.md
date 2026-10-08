# Validation of the GitHub package

On 6 October 2026, the code and exact inputs in this folder were copied to an
independent directory and both numerical certificates were regenerated using
the root driver. The lower and upper commands were run separately; a combined
`--all` invocation was not used. Both full commands exited successfully and
wrote their own root `ACCEPTED.json`.

The frozen source/input manifest has SHA-256
`75a1912de8dd33ef2e14a16f14fd8b53a8a27e52b52f08a1754f4b0e0bf688a4`.

| Check | Fresh result |
| --- | --- |
| Lower | All 35 rows; 1820 matrix checks; 1450 integer flows; 32 original negative controls rejected. |
| Added lower tail check | All 70 exact rational inequalities passed. |
| Upper | All 12 stages and all 32 components; the five printed comparisons and the explicit 256-bit Fourier-tail inequality passed. |
| Relocation and wrapper | Seven regression tests passed, including rejection of all 70 corrupted leading tail weights. |
| Original numerical source | 125 lower files and 20 upper runtime files retain their source-edition bytes. |

The verified endpoints are
`1773000000000/999999998227 <= K_G < 1779893/1000000`.

The lower invocation took 170.7 seconds and the upper
invocation took 308.6 seconds on the recorded machine.
The invocations overlapped, so these are observed elapsed times rather than
isolated performance benchmarks. Python packages and native arithmetic-library
versions are recorded in [summary.json](summary.json).

- [fresh-lower/](fresh-lower/) contains the new root acceptance, lower report,
  environment, exact tail result and driver output.
- [fresh-upper/](fresh-upper/) contains the new root acceptance, upper report,
  pipeline result, final interval certificate and added scalar check.
- [packaging-tests.log](packaging-tests.log) records the fast regression tests.
- [source-edition/](source-edition/) preserves the historical reports shipped
  with `lower-bound-shorter-cert.zip`, clearly separate from these new runs.

Absolute temporary paths in the unmodified reports identify their actual run
locations. They are provenance only; the runnable code does not depend on those
locations. Large regenerated arrays, object files, native libraries and virtual
environments are not distributed. A new run creates its own evidence under
`runs/`.

These checks establish successful fresh execution of the stated numerical
pipelines and manuscript comparisons. They do not replace the analytic proofs
in the paper and are not a proof-assistant formalization. The root manifest
covers verifier source, input and component package files; root documentation
and this validation directory are outside that manifest.
