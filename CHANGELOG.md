# Changelog

## 1.1.0 — 2026-11-18

- Package the completed 2026-10-30 manuscript calculation sources, including the current geometric chemical porosity solver and the seven accompanying cycle/first-event cases.
- Provide the portable `run_thmc` entry point.
- Use exact Dirichlet boundary values by default; expose finite Robin exchange as an explicit option.
- Retain the fixed reference fluid storage coefficient and the current material parameters and constitutive laws.
- Correct software citation metadata to Yuhan Wang, ORCID 0000-0003-0490-5977, in both GitHub and Zenodo metadata.
- Add the MIT license file and consistent version metadata. Validate boundaries and all runtime entry points locally.
- Distribute only runtime code, parameter files, and necessary documentation; omit plotting, postprocessing, outputs, and old source archives.

## 1.0.0 — 2026-08-22

Initial public release.
