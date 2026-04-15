# polylinkR (development version)

## polylinkR 0.0.0.9000 (in progress)

### Package layout

- Modularised the `new_code/` sources into topic files under `R/` with roxygen-managed `NAMESPACE` and `man/`.
- Added **testthat** edition 3 tests for argument validation, S3 helpers, and selected numerical internals.
- Added **renv** lockfile scaffolding, GitHub Actions (`R-CMD-check` with **renv** restore on Ubuntu for R release and R devel), and **pkgdown** site configuration with a getting-started vignette.
- Documented contribution, release, and agent expectations in `CONTRIBUTING.md` and `AGENTS.md`, and added a pull request template aligned with CRAN-oriented R packages.

### User-visible behaviour

- User-facing prose follows **Australian English** spelling where practical; function names remain stable.
- **Active development** disclaimers appear in `DESCRIPTION`, the pkgdown home blurb, and the README.
