# The domino linkage summary class

The linkage summary class contains linkages established in multiple
domino objects through gene regulatory network inference and reference
to receptor- ligand databases. A data frame summarizing meta features
that describe the
[`domino()`](https://FertigLab.github.io/dominoSignal/dev/reference/domino-class.md)
objects compared in the linkage summary facilitates comparisons of
established linkages and differential signaling interactions across
categorical sample covariates.

## Value

an instance of class `linkage_summary`

## Slots

- `subject_names`:

  unique names for each domino result included in the summary

- `subject_meta`:

  data.frame with each row describing one subject and columns describing
  features of the subjects by which to draw comparisons of signaling
  networks. The first column holds the subject names, in the same order
  as `subject_names`.

- `subject_linkages`:

  nested list of linkages inferred for each subject, named by and in the
  same order as `subject_names`. Lists are stored in a hierarchical
  structure of subject-cluster-linkages, including transcription factors
  (tfs), links between transcription factors and receptors (tfs_rec),
  active receptors (rec), links between receptors and ligands(rec_lig),
  and incoming ligands (incoming_lig)
