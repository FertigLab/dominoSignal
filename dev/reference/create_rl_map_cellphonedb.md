# Create a receptor-ligand map from a CellPhoneDB signaling database

Generates a data frame of receptor-ligand interactions from a
CellPhoneDB database annotating the genes encoding the interacting
ligands and receptors to be queried in transcriptomic data.

## Usage

``` r
create_rl_map_cellphonedb(
  genes,
  proteins,
  interactions,
  complexes = NULL,
  database_name = "CellPhoneDB",
  gene_conv = NULL,
  gene_conv_host = "https://www.ensembl.org",
  alternate_convert = FALSE,
  alternate_convert_table = NULL
)
```

## Arguments

- genes:

  data frame or file path to table of gene names in uniprot,
  hgnc_symbol, or ensembl format in CellPhoneDB database format

- proteins:

  data frame or file path to table of protein features in CellPhoneDB
  format

- interactions:

  data frame or file path to table of protein-protein interactions in
  CellPhoneDB format

- complexes:

  optional: data frame or file path to table of protein complexes in
  CellPhoneDB format

- database_name:

  name of the database being used, stored in output

- gene_conv:

  character vector of length 2 formatted as (from, to) or (source,
  target) if gene conversion to orthologs is desired; options are
  ENSMUSG, ENSG, MGI, or HGNC (for from) and MGI or HGNC (for to).
  Identifiers in the `gene_name` column of `genes` (gene symbols, such
  as HGNC for CellPhoneDB human databases). Interactions with any gene
  lacking an ortholog are skipped.

- gene_conv_host:

  host for conversion; default ensembl, could also use mirrors if
  desired

- alternate_convert:

  boolean if you would like to use a non-ensembl method of conversion
  (must supply table; not recommended, use only if ensembl is down)

- alternate_convert_table:

  a data frame with column names corresponding to gene symbol types
  (mm.ens, hs.ens, mgi, hgnc) and rows corresponding to the gene symbols
  themselves for use with alternate_convert = TRUE

## Value

Data frame where each row describes a possible receptor-ligand
interaction

## See also

[`create_domino()`](https://FertigLab.github.io/dominoSignal/dev/reference/create_domino.md)
where output of this function can be used

## Examples

``` r
data(CellPhoneDB)
rl_map_tiny <- create_rl_map_cellphonedb(genes = CellPhoneDB$genes_tiny,
 proteins = CellPhoneDB$proteins_tiny,
 interactions = CellPhoneDB$interactions_tiny,
 complexes = CellPhoneDB$complexes_tiny)

if (FALSE) { # \dontrun{
rl_map_tiny_conv <- create_rl_map_cellphonedb(genes = CellPhoneDB$genes_tiny,
  proteins = CellPhoneDB$proteins_tiny,
  interactions = CellPhoneDB$interactions_tiny,
  gene_conv = c("HGNC", "MGI"),
  gene_conv_host = "https://beta.ensembl.org")
} # }

# Using alternate conversion table instead of biomaRt
# (this simplified table is for illustration only, this is
# NOT the correct way to convert between HGNC and MGI symbols)
hs_genes <- unique(CellPhoneDB$genes_tiny$gene_name)
ortho_table <- data.frame(
  hgnc = hs_genes,
  mgi = paste0(substr(hs_genes, 1, 1), tolower(substr(hs_genes, 2, nchar(hs_genes)))))

rl_map_tiny_alt <- create_rl_map_cellphonedb(genes = CellPhoneDB$genes_tiny,
  proteins = CellPhoneDB$proteins_tiny,
  interactions = CellPhoneDB$interactions_tiny,
  complexes = CellPhoneDB$complexes_tiny, gene_conv = c("HGNC", "MGI"),
  alternate_convert = TRUE, alternate_convert_table = ortho_table)
```
