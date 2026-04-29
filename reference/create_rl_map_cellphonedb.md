# Create a receptor - ligand map from a CellPhoneDB signaling database

Generates a data frame of ligand-receptor interactions from a
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
  ENSMUSG, ENSG, MGI, or HGNC

- gene_conv_host:

  host for conversion; default ensembl, could also use mirrors if
  desired

- alternate_convert:

  boolean if you would like to use a non-ensembl method of conversion
  (must supply table; not recommended, use only if ensembl is down)

- alternate_convert_table:

  supplied table for non-ensembl method of conversion

## Value

Data frame where each row describes a possible receptor-ligand
interaction

## Examples

``` r
data(CellPhoneDB)
rl_map_tiny <- create_rl_map_cellphonedb(genes = CellPhoneDB$genes_tiny,
 proteins = CellPhoneDB$proteins_tiny,
 interactions = CellPhoneDB$interactions_tiny,
 complexes =CellPhoneDB$complexes_tiny)
```
