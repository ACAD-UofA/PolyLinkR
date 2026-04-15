# Agent instructions (Cursor, Copilot, and other assistants)

This file is the **authoritative behaviour contract** for automated contributors. Humans should read it too; it mirrors the pull request template and contribution guide.

## Read first

1. [`CONTRIBUTING.md`](CONTRIBUTING.md) — branches (`dev` as integration), **renv** + **rig**, Australian English, algorithm scope, and **release / versioning** policy (R Packages lifecycle).
2. [`.github/PULL_REQUEST_TEMPLATE.md`](.github/PULL_REQUEST_TEMPLATE.md) — required PR sections and the **AI contributors / “If you are an AI agent”** checklist (interim standard until the organisation’s verbatim block is pasted there).
3. [`.hippo/index.json`](.hippo/index.json) — durable coordinator notes and policies (for example **renv + pak** in CI).

## Operating rules

- **Honesty about verification**: never invent `R CMD check` output, URLs, or CI logs. Run checks or cite real workflow runs.
- **Algorithms**: treat the permutation / linkage-preserving null, rescaling, pruning, and extreme-value machinery as frozen unless a maintainer explicitly approves a change backed by regression tests.
- **Language**: Australian English in user-facing prose; keep stable identifiers unless a breaking release is planned.
- **Scope**: avoid drive-by refactors and unrelated files. Match existing style and dependencies.
- **Secrets**: do not commit tokens, passwords, or personal paths.

## pkgdown and CRAN

- After changing exports or user-facing documentation, regenerate `NAMESPACE` / `man/` and consider `pkgdown::build_site()` for layout issues.
- The public site is **[https://acad-uofa.github.io/PolyLinkR/](https://acad-uofa.github.io/PolyLinkR/)**; keep `DESCRIPTION` and `_pkgdown.yml` `url` fields aligned with whatever GitHub Pages serves.
- Before claiming CRAN readiness, walk the **CRAN extra checks** skill at `.agents/skills/cran-extrachecks/SKILL.md` and update `cran-comments.md` when submitting.
