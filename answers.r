## ---- question1 -------------------------

snp.mds[which(snp.mds[,"C1"] > 0.1),]

## ---- question2 -------------------------

nrow(snp.data) * nrow(rna.data)

## ---- question3 -------------------------

lower.bound <- nrow(snp.data)
upper.bound <- nrow(snp.data) * nrow(rna.data)

## ---- question4 -------------------------

eqtl$cis$ntests/nrow(snp.data)     ## number of tests per SNP
threshold <- 0.05/eqtl$cis$ntests  ## correct p-value threshold

## ---- question5 -------------------------

length(unique(eqtls$snps))

