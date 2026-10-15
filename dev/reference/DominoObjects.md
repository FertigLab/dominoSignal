# Example domino objects

A list of two domino objects built from the example data subsets
([PBMC](https://FertigLab.github.io/dominoSignal/dev/reference/PBMC.md),
[SCENIC](https://FertigLab.github.io/dominoSignal/dev/reference/SCENIC.md),
and
[CellPhoneDB](https://FertigLab.github.io/dominoSignal/dev/reference/CellPhoneDB.md)).

- dom_tiny:

  A
  [`domino()`](https://FertigLab.github.io/dominoSignal/dev/reference/domino-class.md)
  object created using
  [`create_domino()`](https://FertigLab.github.io/dominoSignal/dev/reference/create_domino.md).

- built_dom_tiny:

  The same object after running
  [`build_domino()`](https://FertigLab.github.io/dominoSignal/dev/reference/build_domino.md).

## Usage

``` r
data("DominoObjects")
```

## Format

An object of class `list` of length 2.

## Source

Generated with `data-raw/tiny_data.R`.
