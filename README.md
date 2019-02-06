# `biomartExons`

#### (**R** **P**ackage **T**emplate)
[![DOI](https://zenodo.org/badge/149768973.svg)](https://zenodo.org/badge/latestdoi/149768973)

&nbsp;

###### [Emily Ayala], &lt;emily.ayala@mail.utoronto.ca&gt;


----

<!-- TOCbelow -->
1. About this package:<br/>
&nbsp;&nbsp;&nbsp;&nbsp;1.1. What it is ...<br/>
&nbsp;&nbsp;&nbsp;&nbsp;1.2. How it works ...<br/>
&nbsp;&nbsp;&nbsp;&nbsp;1.3. How to use the package ...<br/>
4. Notes<br/>
5. Acknowledgements<br/>
<!-- TOCabove -->

----


# 1 About this package:

## 1.1 What it is ...

This package takes an input of HNGC symbols in a list and outputs a dataframe containing the start and end co-ordinates of the selected isoform for each gene.

## 1.2 How it works ...

Uisng the ensp2hngc.Rdata produced by Professor Steipe in BCB420.2019.STRING I can retrieve the ENSP ids for the HGNC symbols given.

Using biomart you can then get the ENSG id's as well as the transcript ids.

The isoforms were selected by choosing the longest principal1 isoform as given by APPRIS.

The exons were then retrived using biomaRt and the start and end co-ordinates given.

Return the start and end coordinates coresponding to each HNGC symbol.


## 1.3 How to use the package ...

Source the script file
```R

source('inst/scripts/getStartAndEndCoordinates.R')

```

wheere v is a vector of HNGC symbols

```R
getStartAndEndCOordinates(v)
```
----

# 4 Notes

----

&nbsp;

# 5 Acknowledgements

Thanks to Prof. Steipe, who provided the basis for this project and the rpt base package.

&nbsp;

&nbsp;

<!-- END -->
