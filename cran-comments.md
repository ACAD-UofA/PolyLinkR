## Submission notes (template)

Edit this file before a CRAN submission and keep a copy for the submission email.

## R CMD check environments

- Local: R release and R devel with `devtools::check()` (or `R CMD check`) after `renv::restore()`.
- CI: GitHub Actions `R-CMD-check` workflow (Ubuntu, R release and R devel) with `r-lib/actions/setup-renv`.

## CRAN policies and extras

- Ran through the repository **CRAN extra checks** checklist (see `.agents/skills/cran-extrachecks/SKILL.md`): `DESCRIPTION` metadata, `NEWS.md`, exported documentation (`@return`, examples), and URL hygiene (`urlchecker::url_check()` when preparing a release).
- **rhub** / win-builder: run before first release or substantive compiled-dependency changes, per CRAN guidance.

## Downstream dependencies

- None on CRAN yet (first release pending).

## Special notes

- Optional: summarise any `--as-cran` NOTEs you accept and why.
