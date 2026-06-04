# Master Release Pull Request

## Summary

Describe what is being promoted to `master`/`main` and why.

## Master Release Checklist (Required)

- [ ] [release-master] DESCRIPTION version updated per Bioconductor conventions
- [ ] [release-master] NEWS.md updated for the same version and formatted appropriately
- [ ] [release-master] inst/CITATION software entry is current (authors/year/version)
- [ ] [release-master] README citation/version text is current
- [ ] [release-master] index.Rmd reviewed/updated/knit again if any metadata/package information changed
- [ ] [release-master] Documentation has been updated and includes example usage
- [ ] [release-master] _pkgdown.yml file has been updated to include all functions and vignettes; check that the development mode is set appropriately (release vs devel)
- [ ] [release-master] `R CMD check`, `BiocCheck`, and tests all pass

## Bioconductor Sync Plan

- [ ] [release-master] post-merge push plan to `bioc-upstream/devel` confirmed

## Notes

Add any reviewer guidance or release caveats.
