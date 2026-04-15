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

### Documentation and tests (follow-up)

- Added **`inst/extdata/tiny_polylinkR/`** for a fast `plR_read()` smoke test and documented it in the inputs vignette.
- New tests: `.path_check`, `.par_params`, `.get_quant` (when **tdigest** is available), and `plR_read()` on the tiny fixture.
- Additional vignettes: input formats, parallel execution and seeds, and relationship to legacy PolyLink / PolyLinkR.
- pkgdown reference index now groups **S3** print/summary helpers.
