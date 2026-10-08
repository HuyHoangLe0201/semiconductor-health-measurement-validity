# Added operating evidence (v1.1.0)

The predeclared ensemble contains 144 noiseless closed-loop electrothermal
records and four observation prefixes. It varies load resistance, back-EMF,
current amplitude, reference phase, case temperature and steady/envelope demand.
The complete 50-column operator is analyzed without prescribing a 36-column
quotient. Unsupported targets receive no finite reference standard error.
Reference covariance is a separate 2 mA Gaussian model, not an empirical
feedback/acquisition likelihood. Support fractions are finite-ensemble summaries.

At 4096 periods the voltage target is supported in 119/144 records; the resistance
target in 144/144. Global rank ranges from 24 to 36. Original source units are
unchanged; configurable derivatives retain exact substitution provenance. Gates,
currents and junction hashes are retained for every case; four cases retain arrays.

```shell
python scripts/build_trajectory_native.py
python scripts/trajectory_campaign.py
python scripts/report_trajectory_campaign.py
python scripts/published_estimator_audit.py
```

Native rebuilding requires Windows x64 and MSVC BuildTools. Rebuilding/rerunning
changes generated outputs; use a separate checkout after verifying release hashes.

The additional published-method audit uses equations (1)--(10) of Ou et al.,
APEC 2025, DOI 10.1109/APEC48143.2025.10977401. Accepted manuscript:
https://vbn.aau.dk/ws/portalfiles/portal/779663726/2024360959.pdf
Only our local first-order averaged operator and perturbations are included.
This does not reproduce the PI/PWM controller, mission sampling or HIL results.
The restricted nominated-resistance model is supported; allowing changing
threshold/filter/common-drop terms makes the scalar reduction insufficient.
Stationary nuisance cancels in matched records, while changing nuisance does not.

No author-run converter/HIL measurements, new physical calibration or real-time
qualification are added. Previous scientific data remain unchanged. Release
v1.0.0 and its commit remain available; no dataset DOI or reuse license is assigned.
