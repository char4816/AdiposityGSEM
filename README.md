# AdiposityGSEM

Title: Modeling the genomic architecture of adiposity and anthropometrics across the lifespan
https://doi.org/10.1101/2024.08.14.24312003

`GenomicSEM_efa_cfa_multivariateGWAS.sh`

A vignette-style script including code used to:
1. wrangle the input GWAS summary statistics
2. perform the Genomic SEM exploratory factor analysis (EFA)
3. perform the Genomic SEM confirmatory factor analysis (CFA)
4. perform the multivariate GWAS for the 4 factors

`DrugGeneNetworks.R`

An R script that was used to:
1. construct the drug-gene networks
2. make the plots to visualize the network

`PheWAS_Results.R`

An R script that was used to:
1. plot the phewas associations

Implementation of Genomic SEM EFA and CFA was based on the software wiki:
[https://github.com/GenomicSEM/GenomicSEM/wiki](https://github.com/GenomicSEM/GenomicSEM/wiki)

Implementation of LDpred2 to derive PRS weights was based on the software vignette:
[https://privefl.github.io/bigsnpr/articles/LDpred2.html](https://privefl.github.io/bigsnpr/articles/LDpred2.html)

Implementation of DEPICT to analyze the multivariate GWAS was based on the software example:
[https://github.com/perslab/depict](https://github.com/perslab/depict)

Implementation of FOCUS to perform TWAS analysis was based on the software wiki:
[https://github.com/mancusolab/ma-focus](https://github.com/mancusolab/ma-focus/wiki)


