# Copilot Cloud Agent Instructions for `dominoSignal`

## Repository at a glance
- This is an **R package** (Bioconductor-oriented) for ligand-receptor-transcription factor signaling analysis.
- Core code is in `/R`.
- Tests are in `/tests/testthat`.
- Generated docs are in `/man` and `NAMESPACE` (via roxygen2).
- Package metadata is in `DESCRIPTION` and release notes in `NEWS.md`.
- CI workflows are in `/.github/workflows`.

## High-value workflow for agents
1. Read `DESCRIPTION`, `README.md`, and relevant files in `/R` + matching tests in `/tests/testthat`.
2. Keep changes minimal and scoped to the requested behavior.
3. If you change exported functions, signatures, or roxygen blocks in `/R`, regenerate docs:
   - `Rscript -e 'roxygen2::roxygenise()'`
4. Add/update tests in `/tests/testthat` when logic changes.
5. Re-run build/check/tests before finalizing.

## Validation commands
Run from repository root:
- Build package: `R CMD build .`
- Run checks: `R CMD check --no-manual dominoSignal_*.tar.gz`
- Run tests: `Rscript -e 'testthat::test_dir("tests/testthat")'`
- Optional full package test entrypoint: `Rscript -e 'testthat::test_check("dominoSignal")'`

## Documentation/site commands
- Regenerate Rd + NAMESPACE after roxygen changes:
  - `Rscript -e 'roxygen2::roxygenise()'`
- Build pkgdown site (as in CI):
  - `Rscript -e 'pkgdown::build_site_github_pages(new_process = FALSE, install = TRUE)'`

## Coding/testing patterns to follow
- Input validation is centralized in helper utilities such as `check_arg`; preserve existing error message style because tests assert exact text.
- Keep S4 class/object behavior consistent (`domino`, `linkage_summary` in `R/class_definitions.R`).
- Prefer extending existing utility and processing functions over introducing new patterns.

## Errors encountered in this onboarding run + workarounds
- Error: `R: command not found` when attempting `R --version` and `R CMD build .` in the current sandbox.
- Workaround used/documented:
  1. Use the repository `Dockerfile` to get a working R environment (`rocker/tidyverse:4` base).
  2. Install Bioconductor dependencies listed in Dockerfile (`biomaRt`, `ComplexHeatmap`, `S4Arrays`, `SingleCellExperiment`, `SummarizedExperiment`).
  3. Install package dependencies with `devtools::install_deps(".", dependencies = TRUE)`.
  4. Install the package with `devtools::install(".", dependencies = TRUE)`.
  5. Re-run build/check/test commands inside that prepared environment.

## CI references
- `r-build-check.yml`: PR checks for `master`/`dev` (reusable workflow).
- `document.yaml`: auto-runs roxygen update when `R/**` changes.
- `pkgdown.yaml`: builds/deploys documentation site.
