## StAMPP: Statistical Analysis of Mixed Ploidy Populations <img align="right" src="inst/StAMPP.svg" height="100">

<!-- [Luke W Pembleton](https://website.com) personal website to be added -->

---

StAMPP is an [R](https://www.r-project.org) package for the Statistical Analysis of Mixed Ploidy Populations. It allows users to calculate pairwise Nei's Genetic Distances (Nei 1972), pairwise Fixation Indexes (Fst) (Weir & Cockerham 1984) and also Genomic Relationship matrixes following Yang et al. (2010) in mixed and single ploidy populations. Bootstrapping across loci is implemented during Fst calculation to generate confidence intervals and p-values around pairwise Fst values. StAMPP utilises SNP genotype data of any ploidy level (with the ability to handle missing data) and is coded to utilise multithreading where available to allow efficient analysis of large datasets. StAMPP is able to handle genotype data from genlight objects allowing integration with other packages such adegenet. 

---

### Installation

You can install StAMPP from [CRAN](https://cran.r-project.org/package=StAMPP):

```r
install.packages("StAMPP")
```

---

### Vignette

A vignette describing the use of the package is available from within
R (and [also here](https://cran.r-project.org/package=StAMPP/StAMPP.pdf)). Load the package
and then use the `vignette` function.

```r
library(StAMPP)
vignette("StAMPP", package="StAMPP")
```

---

### Citation

LW Pembleton, NOI Cogan & JW Forster, 2013, Molecular Ecology Resources, 13(5), 946-952. [doi:10.1111/1755-0998.12129](https://doi.org/10.1111/1755-0998.12129).

Thank you in advance.
