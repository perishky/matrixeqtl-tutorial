# eQTL and meQTL practical in R

Nov 2024

Created by Matthew Suderman
[matthew.suderman@bristol.ac.uk](mailto:matthew.suderman@bristol.ac.uk)

## Learning objectives

1. Prepare the data for analysis#

2. Conduct an eQTL analysis using the 
   [MatrixEQTL R package](https://cran.r-project.org/web/packages/MatrixEQTL/index.html)
   
3. Produce plots to visually inspect findings

4. (Optional) Conduct an meQTL analysis

5. (Optional) Identify trios of associated genes, CpG sites and SNPs

## Generate HTML version

In R:

```
library(rmarkdown)
rmarkdown::render("otheromics.rmd")
```

To include answers, rerun `render` after changing `show.answers` to 
`TRUE` in [otheromics.rmd](otheromics.rmd). 

## Student files

* [otheromics.rmd](otheromics.rmd)
* [answers.r](answers.r)
* [plink](plink) 
* [expected-outputs](expected-outputs)/*.*

## Dataset

The dataset was prepared by running [scripts/create-dataset.r](scripts/create-dataset.r). 


