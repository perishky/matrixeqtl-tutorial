# eQTL and meQTL practical in R

Nov 2021

Created by Matthew Suderman
[matthew.suderman@bristol.ac.uk](mailto:matthew.suderman@bristol.ac.uk)

* * * 

## Objectives

1. Prepare the data for analysis.<br>
2. Conduct an eQTL analysis.<br>
3. Produce plots to visually inspect findings.<br>
4. (Optional) Conduct an meQTL analysis.<br>
5. (Optional) Identify trios of associated genes, CpG sites and SNPs.

* * * 

## Preliminaries

Login to BlueCrystal using PuTTY. 

Clone a copy of this tutorial in your home directory.
```
git clone perishky/matrixeqtl-tutorial
```

Run the following command to request an interactive session on a compute node:
```
salloc --nodes=1 --ntasks=1 --mem=3G  --time=02:00:00
srun --job-name="example" --pty bash -i
```

Change to cloned repository for this practical:
```
cd ~/matrixeqtl-tutorial
```

It is good practice to save outputs to another directory
so that it is easy to distinguish output files
from other types of files such as 
data, script or documentation files.
We create the output directory as follows:
```
mkdir output
```

---

## Data

In this practical we will run eQTL and meQTL analyses
on a small publicly available dataset:

> Wagner, et al.
> [The relationship between DNA methylation, genetic and expression inter-individual variation in untransformed human fibroblasts.](http://www.ncbi.nlm.nih.gov/pubmed/24555846)
> Genome Biol. 2014 Feb 20;15(2):R37.
> [GEO link](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE53261)

It includes:
- genotype (Illumina Human1M-Duov3 DNA Analysis BeadChip),
- DNA methylation (Illumina HumanMethylation450 BeadChip), and
- gene expression (Illumina HumanRef-8 v3.0 Expression Beadchip) profiles
from 57 human skin fibroblast samples obtained from Coriell and McGill Cellbank.

---

## Analysis steps

1. Preprocessing.<br>
   a. Population stratification.<br>
   b. Gene expression outliers.<br>
   c. Matching samples between datasets.<br>
   d. Non-genetic variation removal.
2. Software and data.<br>
   a. R package MatrixEQTL.<br>
   b. Genetic data.<br>
   c. Population stratification covariates.<br>
   d. Expression data.<br>
   e. Gene and SNP locations.
3. eQTL analysis.
4. Aftermath.<br>
   a. Number of associations.<br>
   b. Top associations.<br>
   c. Association plotting.
5. (Optional) meQTL analysis.
6. (Optional) Associated trios.

---

## 1a. Preprocessing: Population stratification

This is detected by visualizing
genetic distances between all pairs of individuals
using multi-dimensional scaling (MDS).

We will use PLINK for these computations.
Please ensure that PLINK is installed before continuing.
```
module add apps/plink/2.00
```

Calculate genetic distances:
```
plink --bfile data/snp --genome --out output/plink
```

The output `output/plink.genome` file in the current directory
gives all pairwise distances (see the `DST` column).

You can type the following command to see the first few lines
in the file:
```
head output/plink.genome
```

---

## 1a. Preprocessing: Population stratification in 2-dimensions

We represent these distances in two dimensions as follows:
```
plink --bfile data/snp --read-genome output/plink.genome --cluster --mds-plot 2 --out output/plink
```

This generates a file `output/plink.mds` that looks something like this:
```
 FID       IID    SOL           C1           C2 
   1   GM02704      0  -0.00860421   -0.0114082 
   2   GM02706      0   -0.0145491    0.0005202 
   3   GM01650      0  -0.00431194  -0.00171059 
   4   GM01653      0  -0.00306369  -0.00175129 
   5   GM02640      0   -0.0035343  -0.00811608 
   6   GM02641      0  -0.00876442  -0.00616735 
```

FID - Family ID;<br>
IID - Individual ID (or sample name);<br>
SOL - Assigned solution code;<br>
C1 - Position on first dimension;<br>
C2 - Position on second dimension.

---

## 1a. Preprocessing: Population stratification visually

Start `R`. First you'll need to load the module:
```
module add languages/r/4.1.0
```

To start R, simply type `R`.
```
R
```

We'll send outputs from R to the output directory we set up earlier.
For convenience, we?ll save the location as a variable.
```r
out.dir <- "output"
```

We now plot the coordinates in `output/plink.mds` in `R` as follows:

```r
snp.mds <- read.table(file.path(out.dir, "plink.mds"), header=T, sep="", stringsAsFactors=F)
pdf(file.path(out.dir, "mds.pdf"))
plot(snp.mds[,"C1"], snp.mds[,"C2"], main="SNP MDS plot", xlab="first", ylab="second", pch=19)
dev.off()
```

> *Question*: Can you identify the four potential outlier samples in the plot?

Answer:

```r
snp.mds[which(snp.mds[,"C1"] > 0.06 | snp.mds[,"C2"] > 0.06),]
```

```
   FID    IID SOL         C1        C2
39  39 WG2120   0 -0.0336041 0.0728445
40  40 WG2121   0 -0.0323253 0.0719037
43  43 WG1977   0  0.0725762 0.0344586
44  44 WG1978   0  0.0713593 0.0332588
```

Given how far these samples differ from the rest,
we might actually consider removing them from the analysis.
In this case, we'll attempt to adjust for these differences.

---

## 1b. Preprocessing: Gene expression outliers

We could use a similar approach to identify gene expression outliers.

```r
rna <- read.csv("data/rna-data.csv", row.names=1)
rna.dist <- 1-abs(cor(rna))
rna.mds <- cmdscale(rna.dist,eig=TRUE, k=2)
pdf(file.path(out.dir, "mds-rna.pdf"))
plot(rna.mds$points[,1], rna.mds$points[,2],
     main="RNA MDS plot", xlab="first", ylab="second", pch=19)
dev.off()
```

---

## 1c. Preprocessing: Matching samples between datasets

All too often samples are mislabelled in the lab.
Fortunately, it is possible to computationally detect mismatches between
genetic and gene expression data using MixupMapper:

http://genenetwork.nl/wordpress/mixupmapper/

In this case, we will skip this step, assuming that the study authors
haven't made such errors!

It is possible to do the same for genetic and DNA methylation data,
particularly for data generated using the Illumina HumanMethylation450 BeadChip.
It actually measures a selection of SNPs in addition to the
methylation levels.

---

## 1d. Preprocessing: Non-genetic variation removal

Gene expression measurements can be influenced
by factors other than genetic variation.
To maximize power to detect associations between genes and genetic variants,
we can remove all non-genetic variation from the gene expression data.
To save time, we will not do that here. In case you're curious, here are the basic steps:

1. Calculate the principal components of the gene expression data.
2. Test associations between principal components and genetic variants.
3. Remove any components with strong genetic associations.
4. Remove variation of the remaining principal components from the gene expression data.

---

## 2a. Software and data: R package MatrixEQTL

If `MatrixEQTL` has not already been installed, we install it as follows:

```r
BiocManager::install("MatrixEQTL")
# Respond by typing 'yes' to any questions.
```


We then load the library in R as follows:

```r
library(MatrixEQTL)
```

---

## 2b. Software and data: Genetic data

MatrixEQTL operates on `SlicedData` objects.
These objects provide access to large datasets
without using huge amounts of computer memory.

The genetic data is contained in the file `data/snp-data.csv`.
We load it as follows:

```r
snp.data <- SlicedData$new()            ## create a new SlicedData object
snp.data$fileDelimiter <- ","           ## columns in the file are separated by commas
snp.data$fileOmitCharacters <- "NA"     ## missing values are represented by "NA"
snp.data$fileSkipRows <- 1              ## the first column provides row names
snp.data$fileSkipColumns <- 1           ## the first row provides column names
snp.data$fileSliceSize <- 2000          ## only load the data 2000 rows at a time to save memory
snp.data$LoadFile("data/snp-data.csv") 
```
Loading the data takes about half a minute.

---

## 2c. Software and data: Population stratification covariates

When we test associations between genetic variants,
we will need to adjust for population stratification.
This is encoded by the multi-dimensional scaling coordinates
that we created earlier.

As required by `MatrixEQTL`,
we load them as `SlicedData` objects:

```r
covariates.data <- SlicedData$new()
covariates.data$initialize(t(snp.mds[,c("C1","C2")]))
```

---

## 2d. Software and data: Expression data

Next we convert the gene expression data to `SlicedData` objects
and verify that samples match the genotype data samples.

```r
rna.data <- SlicedData$new()            ## create a new SlicedData object
rna.data$fileDelimiter <- ","           ## columns in the file are separated by commas
rna.data$fileOmitCharacters <- "NA"     ## missing values are represented by "NA"
rna.data$fileSkipRows <- 1              ## the first column provides row names
rna.data$fileSkipColumns <- 1           ## the first row provides column names
rna.data$fileSliceSize <- 2000          ## only load the data 2000 rows at a time to save memory
rna.data$LoadFile("data/rna-data.csv") 
```

Don't forget to make sure that the gene expression and genetic
datasets have the samples in the same order!

```r
stopifnot(all(colnames(snp.data) == colnames(rna.data)))
```

---

## 2e. Software and data: Gene and SNP locations

> *Question*: If we tested associations between all SNP-gene pairs,
> how many tests would we have to perform?

Answer:

```r
nrow(snp.data) * nrow(rna.data)
```

```
## [1] 869593928
```

That's a lot of tests. To save time, we will only consider
*cis* pairs (i.e. SNPs that are with 1 million bases of the gene).
`MatrixEQTL` can do this for us if we provide
SNP and gene locations in the genome:

```r
snp.loc <- read.csv("data/snp-features.csv", row.names=1)
rna.loc <- read.csv("data/rna-features.csv", row.names=1)
```

---

## 3. eQTL analysis: p-value thresholds

> *Question*: `MatrixEQTL` requires a p-value threshold for significance.
> To adjust for the number of tests, we could divide 0.05 by the number of
> tests to be performed. Because we only test *cis* pairs, the exact
> number of difficult to estimate.
> Can you think of a simple lower bound and upper bound on the number of tests?

Answer:

```r
lower.bound <- nrow(snp.data)
upper.bound <- nrow(snp.data) * nrow(rna.data)
```

We'll use the lower bound to estimate a p-value threshold
so that we don't miss any associations
with statistically significant p-values.

```r
threshold <- 0.05/lower.bound
```

---


## 3. eQTL analysis: finally!

Finally we are ready to run the eQTL analysis.


```r
eqtl <- Matrix_eQTL_main(snps = snp.data,
                         gene = rna.data,
                         cvrt = covariates.data, ## population stratification
                         pvOutputThreshold = 0, ## consider only cis pairs
                         pvOutputThreshold.cis = threshold, 
                         output_file_name.cis = file.path(out.dir, "eqtl-cis.txt"),
                         snpspos = snp.loc, 
                         genepos = rna.loc[,c("geneid","chr","left","right")],
                         cisDist = 1e6, ## define cis as < 1Mb
                         useModel = modelLINEAR, ## test using linear models 
                         verbose = TRUE,
                         pvalue.hist = TRUE,
                         min.pv.by.genesnp = FALSE,
                         noFDRsaveMemory = FALSE)
```

This takes about 30 seconds.

---

## 4a. Aftermath: Number of associations

Here is the exact number of tests that were performed:

```r
eqtl$cis$ntests
```

```
## [1] 4863929
```

> *Question*: Above we guessed 1 test per SNP.
> How many tests were performed on average for each SNP?
> What p-value threshold should we apply given this number?

Answer:

```r
eqtl$cis$ntests/nrow(snp.data)     ## number of tests per SNP
```

```
## [1] 27.28987
```

```r
threshold <- 0.05/eqtl$cis$ntests  ## correct p-value threshold
```

At this threshold (p < 1.03e-08), we identify hundreds of
associations.

```r
eqtls <- eqtl$cis$eqtls
eqtls <- eqtls[which(eqtls$pvalue < threshold),]
nrow(eqtls)
```

```
## [1] 468
```

> *Question*: How many eQTLs did we identify?

Answer:

```r
length(unique(eqtls$snps))
```

```
## [1] 403
```

---

## 4b. Aftermath: Top associations

We will plot the strongest association:

```r
idx <- which.min(eqtls$pvalue)
eqtls[idx,]
```

```
##        snps         gene statistic       pvalue          FDR       beta
## 1 rs1025996 ILMN_2100085 -21.33979 1.096262e-27 8.063919e-22 -0.4323303
```

We save the corresponding gene and SNP:

```r
top.gene <- as.character(eqtls$gene[idx])
top.snp <- as.character(eqtls$snps[idx])
top.gene
```

```
## [1] "ILMN_2100085"
```

```r
top.snp
```

```
## [1] "rs1025996"
```

---

## 4c. Aftermath: Association plotting (data extraction)

To plot the association, we will need to extract the
SNP and gene expression values for the pair.
This can be a bit painful using `SlicedData` objects
because they don't support direct access to individual
rows or columns.
Here we create our own function to do this.
Don't worry about how it works, just copy and paste into R.

```r
sliced.data.row <- function(object, row) {
    if (is.factor(row))
        row <- as.character(row)
    if (is.character(row)) {
        row <- which(rownames(object) == row)
        if (length(row) == 0)
            stop("Invalid row")
    }
    slice.size <- nrow(object[[1]])
    slice.idx <- floor(row/slice.size) + 1
    slice <- object[[slice.idx]]
    row <- row - slice.size * (slice.idx - 1)
    ret <- slice[row,]
    names(ret) <- colnames(object)
    ret             
}
```

---

## 4c. Aftermath: Association plotting

We can use this function to obtain the genotypes for any given SNP.
For example, below we extract data for the top SNP-gene pair.


```r
genotypes <- sliced.data.row(snp.data, top.snp)
genotypes <- factor(genotypes, levels=0:2, labels=c("AA","AB","BB"))
expression.levels <- sliced.data.row(rna.data, top.gene)
```

We plot their association as follows:


```r
pdf(file.path(out.dir, "top-snp-gene.pdf"))
boxplot(expression.levels ~ genotypes,
        main=paste(top.snp,top.gene),
        ylim=range(expression.levels),
        outline=F)
stripchart(expression.levels ~ genotypes,
           method="jitter", add=TRUE, vertical=TRUE,
           col=c("blue","orange","black"), pch=19)
dev.off()
```


---

## 5. (Optional) meQTL analysis


```r
cpg.data <- SlicedData$new()
cpg.data$fileDelimiter <- ","           ## columns are separated by commas
cpg.data$fileOmitCharacters <- "NA"     ## missing values are represented by "NA"
cpg.data$fileSkipRows <- 1              ## the first column provides row names
cpg.data$fileSkipColumns <- 1           ## the first row provides column names
cpg.data$fileSliceSize <- 2000          ## only load the data 2000 rows at a time
cpg.data$LoadFile("data/cpg-data.csv") 

cpg.loc <- read.csv("data/cpg-features.csv", row.names=1)

threshold <- 0.05/nrow(snp.data)

meqtl <- Matrix_eQTL_main(snps = snp.data,
                          gene = cpg.data,
                          cvrt = covariates.data,
                          pvOutputThreshold = 0,  ## ignore trans associations
                          pvOutputThreshold.cis = threshold,
                          output_file_name.cis = file.path(out.dir, "meqtl-cis.txt"),
                          snpspos = snp.loc, 
                          genepos = cpg.loc[,c("geneid","chr","left","right")],
                          cisDist = 1e6,                  ## define cis as < 1Mb
                          useModel = modelLINEAR,
                          verbose = TRUE,
                          pvalue.hist = TRUE,
                          min.pv.by.genesnp = FALSE,
                          noFDRsaveMemory = FALSE)

threshold <- 0.05/meqtl$cis$ntests
meqtls <- meqtl$cis$eqtls
meqtls <- meqtls[which(meqtls$pvalue < threshold),]
```

---

## 6. (Optional) Associated trios

Such trios are of interest because they could indicate cases where
a SNP modifies methylation levels which in turn modify gene expression levels.

Construct the trios:

```r
trios <- merge(data.frame(snp=eqtls$snps, gene=eqtls$gene, eqtl.statistic=eqtls$statistic),
               data.frame(snp=meqtls$snps, cpg=meqtls$gene, meqtl.statistic=meqtls$statistic))
```

There are 168 trios.

Calculate p-values for each gene/CpG site association:

```r
trios$p.value <- sapply(1:nrow(trios), function(i) {
    expr <- sliced.data.row(rna.data, trios$gene[i])
    meth <- sliced.data.row(cpg.data, trios$cpg[i])
    cor.test(meth, expr)$p.value
})
```

Identify the strongest gene/CpG site association.

```r
idx <- which.min(trios$p.value)
trios[idx,]
```

```
##          snp         gene eqtl.statistic        cpg meqtl.statistic
## 5 rs10220917 ILMN_1656045       9.664163 cg25118879       -13.50116
##        p.value
## 5 2.513491e-15
```

