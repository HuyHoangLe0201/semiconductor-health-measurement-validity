# Semiconductor-health measurement validity

Reproducibility materials for **Identifiability and uncertainty of
semiconductor-health measurements in multilevel inverter drives**.

This study examines four distinct requirements for an indirect health
measurement: structural support, reference- and nuisance-dependent precision,
model/calibration adequacy, and detection capability.

## Evidence and scope

The repository contains exact model fixtures, Gaussian Monte Carlo evidence,
independent electrical and electrothermal simulations, circuit-replay diagnostics,
and compiled offline closed-loop software-in-the-loop (SIL) records.
**No author-run physical converter bench or hardware-in-the-loop measurements
are included.** Physical accuracy, measured thermal calibration and real-time
feasibility remain unestablished.

The independent-record study retains 3,200 primary closed-loop SIL records and
200 explicitly labelled common-noise comparator records. For each ADC setting,
200 null calibration pairs precede 400 held-out evaluation triplets. A supported
50 mV contrast is detected in 400/400 injected comparisons at both resolutions;
the supported 0.2 mOhm resistance contrast is detected in 19.25% (16-bit) and
6.75% (12-bit). These results describe the declared simulated acquisition context.
Null calibration is per-target and marginal under exchangeable simulated records;
it is not a guarantee for repeated online monitoring or physical calibration.

## Added operating evidence

Release v1.1.0 adds 144 controller-generated records, complete-operator support
and conditioning diagnostics, and a local observation audit of a published
estimator. See [OPERATING_EVIDENCE.md](OPERATING_EVIDENCE.md) for scope and commands.
The audit is not a reproduction of the published HIL implementation.

## Quick verification

Use Python with the dependencies in `requirements.txt`:

```shell
python -m pip install -r requirements.txt
python scripts/verify_upload_bundle.py
python scripts/check_independent_sil.py
```

The first check verifies all release payload hashes. The second recomputes the
reported decision rates from the 1,200 trial groups, checks disjoint primary noise
streams and frozen calibration partitions, and verifies retained native parity.
It writes `validation/independent_sil_audit.json`; verify checksums first.
These checks do not rerun the expensive Monte Carlo experiment.

## Rebuilding and rerunning

See `INDEPENDENT_SIL_REPRODUCE.txt` for commands, seeds and restart semantics.
Native rebuilding uses Windows x64 and MSVC BuildTools, with the compiler root
configured through `-VsRoot` in the PowerShell build scripts. Retained DLLs are
host-specific evidence artifacts, not a portable runtime or timing qualification.
The full Monte Carlo run is a substantial offline computation. Preserve released
outputs and use a separate working copy for regeneration.

## Contents

| Directory | Contents |
| --- | --- |
| `data/` | Operators, summary results, simulated traces and independent SIL trials |
| `scripts/` | Model, plant, simulation, reporting and verification code |
| `figures/` | Data-derived plots and detailed monochrome LaTeX mathematical plates |
| `validation/` | Numerical audits, circuit diagnostics, capture fixtures and native artifacts |
| `external/` | Third-party provenance and retrieval information |

All independent-record trial summaries are retained. Full arrays are retained
for parity and representative trials; remaining traces can be regenerated from
the documented seeds and source. The manuscript and author declarations are not
distributed in this repository.

## Third-party materials

NASA raw MAT files, manufacturer PDFs and the bundled ngspice distribution are
not redistributed. Retrieve NASA `Device2__1.mat` through `Device5__1.mat` using
the links in `external/nasa/provenance.json`. Current SI calibration and upstream
byte matching remain unverified. Datasheet anchors are documented in
`external/datasheets/Infineon_IKx40N65H5.json`; ngspice retrieval is described by
`scripts/fetch_ngspice.ps1`. Full-window ngspice refinements remain unresolved.
External materials retain their respective terms.

`redistribution_audit.json` records exclusions and sanitization inherited from
the audited export. Local filesystem identity was removed; numerical arrays and
scientific formulas were preserved.

## Version and citation

The current reproducibility snapshot is release **v1.1.0**; v1.0.0 remains available. Cite its release URL
and the corresponding full commit SHA; see `CITATION.md`. No dataset DOI is
assigned. Publication of this repository does not establish manuscript acceptance.

## Reuse terms and assistance

No additional reuse license has been assigned; see `REUSE_TERMS.txt`. Public
availability does not grant an MIT, Apache, Creative Commons or other license.
Third-party inputs retain their own terms.

OpenAI Codex assisted preparation of code and visualizations. Numerical outputs
were produced by executing the documented methods on retained inputs. Human
authors retain responsibility for verification of methods, sources and results.
