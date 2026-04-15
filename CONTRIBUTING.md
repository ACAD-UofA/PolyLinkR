# Contributing to polylinkR

Thank you for helping improve **polylinkR**. This document summarises how we work day to day and how we cut releases.

## Branches and pull requests

This repository follows a **git-flow–style** workflow:

- **`dev`** is the integration branch (equivalent to git-flow `develop`). Open pull requests into `dev` unless a maintainer directs otherwise.
- **`master`** holds production-ready history; releases merge from `release/*` or `hotfix/*` per team practice.
- Use **`feature/*`** branches for new work. Prefer the GitHub CLI (`gh pr create`) for pull requests.

See also `AGENTS.md` and `.github/PULL_REQUEST_TEMPLATE.md`.

## Development environment

- **R**: install current **R release** and **R devel** locally. The [**rig**](https://github.com/r-lib/rig) tool is convenient, for example `rig add release` and `rig add devel`, then invoke `R-*` binaries when running checks.
- **Dependencies**: this project uses [**renv**](https://rstudio.github.io/renv/). After cloning, open the project in R and run `renv::restore()`. Pak is enabled in `.Rprofile` via `options(renv.config.pak.enabled = TRUE)`; CI sets `RENV_CONFIG_PAK_ENABLED=true` for the same behaviour.
- **Checks**: run `devtools::document()` (or `roxygen2::roxygenise()`) after changing roxygen comments, then `devtools::check()` on at least R release before opening a PR.

## Documentation and language

- Use **Australian English** in user-facing strings: roxygen, vignettes, README, messages, and GitHub templates.
- Keep **function and argument names** stable unless you are intentionally shipping a breaking change with `NEWS.md` coverage.

## Algorithms and scope

- Do **not** change the permutation null, set-score definitions, tail/GPD logic, or pruning semantics unless you are fixing a demonstrable bug. Pair any such fix with a **regression test** that locks intended behaviour, and call it out explicitly in the PR.

## Releases, versioning, and NEWS

We follow the lifecycle guidance in [**R Packages** (2e), Chapter 21](https://r-pkgs.org/lifecycle.html):

- **Commits** are for developers; a **released version** is the user-facing signal. Bump `Version` in `DESCRIPTION` only when cutting a deliberate release—not on every merged PR.
- Use **patch** for internal fixes and non-API changes, **minor** for backward-compatible features, and **major** when breaking behaviour or renaming stable API.
- When you *do* cut a release: update `NEWS.md`, tag the repository, run full `R CMD check` on release and devel, and apply the **CRAN extra checks** checklist (see `.agents/skills/cran-extrachecks/SKILL.md` and `cran-comments.md`).
- Before a first CRAN submission or a major compiled dependency change, run **rhub** or **win-builder** as appropriate.

## pkgdown

- Local preview: `pkgdown::build_site()`.
- The published manual lives at **[https://acad-uofa.github.io/PolyLinkR/](https://acad-uofa.github.io/PolyLinkR/)** (GitHub Pages from `gh-pages`, built by the **`pkgdown` GitHub Actions workflow**). The same URL is listed in `DESCRIPTION` and `_pkgdown.yml`.

## Record keeping

Material decisions (dependency policy, CI matrices, release bar) should be noted under **`.hippo/`** so future contributors (human or automated) can rely on a single durable log.
