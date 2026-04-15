# Agent instructions (PolyLinkR)

## Read first

1. [`CONTRIBUTING.md`](CONTRIBUTING.md) — branches (`dev` as integration), **renv** + **pak**, Australian English, algorithm scope, and **release / versioning** policy.
2. [`.github/PULL_REQUEST_TEMPLATE.md`](.github/PULL_REQUEST_TEMPLATE.md) — required PR sections and the AI-contributor checklist.
3. [`.hippo/index.json`](.hippo/index.json) — durable coordinator notes and policies.

## Algorithm source of truth

`new_code/` contains the original, working, algorithmically-sound R scripts (`main_plR.R`, `class_plR.R`, `internal_plR.R`) that serve as the authoritative reference for the package's statistical machinery. The modular files in `R/` are derived from these. When the package and `new_code/` disagree on behaviour, the package is the source of truth (but treat algorithmic correctness as frozen per the constraints below).

## Public workflow entrypoints

The package exposes four public functions following a sequential pipeline: `plR_read()` → `plR_permute()` → `plR_rescale()` → `plR_prune()`. Start with `?plR_read` for required table formats and column conventions.

## Verifiable commands

```bash
# Install dependencies (pak-enabled; also set RENV_CONFIG_PAK_ENABLED=true in CI)
renv::restore()

# Document (roxygen2) → NAMESPACE + man/
devtools::document()

# Full check (run before opening a PR)
devtools::check()

# Run tests only
devtools::test()

# Local pkgdown preview
pkgdown::build_site()
```

CI runs two jobs: **R release + renv** and **R devel + DESCRIPTION deps** (see `.github/workflows/R-CMD-check.yaml`). Report actual output — do not fabricate `R CMD check` results.

## Branch and PR conventions

- `dev` is the integration branch. Branch `feature/*` from `dev`, open PRs into `dev`.
- `master` holds production history; `release/*` only when cutting a release.
- Use `gh pr create --base dev`.

## Algorithm constraints

The permutation null, set-score definitions, tail/GPD logic, pruning semantics, rescaling, and extreme-value machinery are **frozen**. Do not change them except to fix a demonstrable bug, paired with a regression test and explicit maintainer sign-off.

## Language

Australian English in user-facing prose (roxygen, messages, README, vignettes). Keep function and argument names stable unless this is an agreed breaking release.

## pkgdown and CRAN

- After changing exports or user-facing documentation, run `devtools::document()` and consider `pkgdown::build_site()`.
- Public site: **https://acad-uofa.github.io/PolyLinkR/** — keep `DESCRIPTION` and `_pkgdown.yml` `url` fields aligned.
- Before claiming CRAN readiness, walk **`.agents/skills/cran-extrachecks/SKILL.md`** and update `cran-comments.md`.

## Secrets

Do not commit tokens, passwords, personal paths, or large generated artefacts.

## Versioning

Bump `Version` in `DESCRIPTION` only when cutting a deliberate release. Day-to-day merges to `dev` do not require a version change. Follow the [R Packages lifecycle](https://r-pkgs.org/lifecycle.html) (patch / minor / major).
