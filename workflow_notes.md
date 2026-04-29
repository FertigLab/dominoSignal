# Developer Guidelines: GitHub and Bioconductor Synchronization

## Purpose

This document describes how to keep the GitHub repository and
Bioconductor package repository in sync while preserving:

- a clear development path for new features
- a safe path for stable-release bug fixes
- Bioconductor-compliant versioning
- automatic documentation deployment

This repository is the active collaboration mirror. The Bioconductor Git
server is the canonical distribution source for package builds installed
through BiocManager.

## Repository and Branch Mapping

Use the following branch mapping:

- Bioconductor devel branch: bioc-upstream/devel
- GitHub development mirror: origin/master
- Bioconductor stable release branch: bioc-upstream/RELEASE_X_YY
- GitHub stable mirror branch: origin/bioconductor/RELEASE_X_YY

Supporting branches used in this repository:

- dev: integration branch before promotion to master
- feature branches: short-lived branches created from dev
- archive/\*: retained historical branches
- gh-pages: site deployment branch maintained by workflow

## One-Time Bioconductor Remote Setup

If you do not already have the Bioconductor remote configured:

``` bash
git remote add bioc-upstream git@git.bioconductor.org:packages/dominoSignal.git
```

Notes:

- Bioconductor SSH credentials are required.
- The remote name bioc-upstream is used in this document but not
  required.

## Versioning Rules (Bioconductor)

Bioconductor uses x.y.z versioning with the following conventions:

- odd y: devel series
- even y: stable release series
- z: patch increment within the same series

How to implement version numbers:

- Increment z for any changes pushed to the stable or development
  releases (no changes necessary for the dev branch or feature branches
  as they are not user facing).
- Keep DESCRIPTION and NEWS.md consistent for each pushed version.
- For major transitions, Bioconductor convention uses x.99.z in devel
  before rolling to (x+1).0.z for the next stable release. Use this when
  breaking changes are released.

## Existing GitHub Actions in This Repository

Current workflows relevant to this runbook:

- r-build-check workflow: runs on pull requests to master and dev
- pkgdown workflow: runs on push and pull request for master, plus
  manual dispatch

## Standard Feature Development Workflow

Use this for new functionality and non-urgent maintenance.

1.  Create a feature branch from dev.
2.  Implement code, tests, and documentation updates.
3.  Add changelog notes under an Unreleased section while feature
    branches are in progress.
4.  Open a pull request into dev.
5.  Resolve review comments and ensure checks pass.
6.  When ready to promote, open pull request dev -\> master.
7.  Before merging to master, finalize version bump in DESCRIPTION and
    consolidate NEWS entries for that version.
8.  Merge to master.
9.  Push master updates to Bioconductor devel:

``` bash
git fetch bioc-upstream
git checkout master
git merge --ff-only bioc-upstream/devel
git push bioc-upstream master:devel
```

## Stable Release Bug-Fix Workflow

Use this only for fixes appropriate for stable Bioconductor release
branches.

1.  Check out the stable mirror branch on GitHub:

``` bash
git checkout bioconductor/RELEASE_X_YY
git pull --ff-only origin bioconductor/RELEASE_X_YY
```

2.  Create a bug-fix branch from that stable branch.
3.  Implement minimal bug fix and tests.
4.  Bump patch version (z) in DESCRIPTION.
5.  Add NEWS entry for the patched stable version.
6.  Open pull request into bioconductor/RELEASE_X_YY.
7.  After merge, push branch to Bioconductor stable release branch:

``` bash
git checkout bioconductor/RELEASE_X_YY
git push bioc-upstream bioconductor/RELEASE_X_YY:RELEASE_X_YY
```

8.  Apply the same fix to devel by opening a second pull request to
    master, then push to bioc-upstream/devel after merge.

## Documentation Site Workflow

This repository uses pkgdown with explicit development mode in
\_pkgdown.yml.

Guidelines:

- master should represent the devel package line and keep pkgdown
  development mode set for devel docs behavior.
- release branches should set pkgdown development mode to release when
  building stable site content.
- pkgdown deploy is automated to gh-pages on non-PR runs and can also be
  triggered manually.

Recommended release-branch docs update sequence:

1.  Switch to release branch and verify \_pkgdown.yml mode.
2.  Trigger pkgdown workflow manually if needed.
3.  Confirm expected pages are present on gh-pages deployment.

## Bioconductor Biannual Release Procedure (April/October)

Right before Bioconductor creates the new release branch (after current
release changes are disabled):

1.  Create a pre-release branch of the form `bioconductor/prep_X_YY`
    (this is a promotion branch that will not be merged to master):

``` bash
git checkout -b bioconductor/prep_X_YY bioc-upstream/devel
```

2.  Update versioning for NEWS.md, CITATION, README.md, index.Rmd, etc,
    to match new stable release version.
3.  Update DESCRIPTION by bumping patch number but leave as devel
    versioning (odd y).
4.  Push to Bioconductor devel so that upon release, when DESCRIPTION
    version is bumped, documentation matches new stable version.

``` bash
git push bioc-upstream bioconductor/prep_X_YY:devel
```

Once Bioconductor creates a new release branch:

1.  Fetch Bioconductor remotes:

``` bash
git fetch bioc-upstream
```

2.  Create local GitHub mirror branch for new release and push:

``` bash
git checkout -b bioconductor/RELEASE_X_NEW bioc-upstream/RELEASE_X_NEW
git push -u origin bioconductor/RELEASE_X_NEW
```

3.  On that new release branch, ensure pkgdown mode is set to release
    and trigger site deployment.
4.  Sync master with Bioconductor devel using fast-forward only:

``` bash
git fetch bioc-upstream
git checkout master
git merge --ff-only bioc-upstream/devel
git push origin master
```

5.  Verify DESCRIPTION version progression is correct (release branch
    even y, devel branch odd y).

## Verification Checklist for Any Push to Bioconductor

Before pushing:

- R CMD check passes locally and in CI
- tests pass
- DESCRIPTION version updated appropriately
- NEWS updated for the corresponding version

After pushing:

- Bioconductor build report shows successful build (check [the
  Bioconductor Build Report](https://bioconductor.org/checkResults/))
- GitHub mirror branch and Bioconductor branch remain aligned
- pkgdown deployment status is green when docs changed

## R Versioning

R releases an update each April. Bioconductor developers should use the
R-devel version from October through April to ensure the Bioconductor
release aligns with the new R version and then use the R-release version
from April through October.
