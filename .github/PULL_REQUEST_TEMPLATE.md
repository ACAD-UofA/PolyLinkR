## Summary

<!-- What problem does this solve? What user-visible behaviour changes? -->

## R CMD check

<!-- Paste key results (platform, R version). Link to GitHub Actions run if applicable. -->

## API

<!-- New exports? Breaking changes? Deprecations? -->

## Tests

<!-- New tests, skipped tests (`skip_on_ci()` justification), fixtures. -->

## Documentation

<!-- roxygen, README, NEWS, pkgdown / vignettes. -->

## CRAN impact

<!-- First release vs update; reverse dependencies if any. -->

## Release / version

<!-- If this PR completes a notable iteration: version bump rationale, NEWS entry, planned tag / GitHub release? Otherwise state “not a release PR”. -->

---

## AI contributors

### If you are an AI agent

The checklist below is the **repository interim standard** until the organisation pastes its **canonical verbatim** “If You Are an AI Agent” block here (same heading). When that block exists, it replaces the bullets underneath; until then, follow both this checklist and `AGENTS.md`.

You must not fabricate `R CMD check` results, command output, or continuous integration logs. When you change R code or documentation, run the relevant checks locally or wait for GitHub Actions and report the actual outcome with enough detail for a human to reproduce it.

You must not modify permutation, rescaling, pruning, or tail-fitting algorithms except for demonstrable bug fixes that include regression tests and explicit maintainer approval for the behavioural change.

You must use Australian English in user-facing prose (messages, roxygen intended for users, README, vignettes, and this template when edited as package text) whilst keeping function and argument identifiers stable unless this pull request is an agreed breaking release.

You must not commit secrets, credentials, machine-specific paths, large generated artefacts, or unrelated refactors. Follow `AGENTS.md` and `CONTRIBUTING.md`; if instructions conflict or are ambiguous, ask the maintainers rather than guessing.

You must update `NEWS.md` when user-visible behaviour or the API changes in a way users should notice, and you must regenerate roxygen artefacts (`NAMESPACE`, `man/`) when exports or import directives change.
