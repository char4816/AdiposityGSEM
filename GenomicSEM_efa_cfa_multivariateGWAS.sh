
cat << 'EOF' >gather_sumstats.R
#load in packages
library(data.table)
library(tidyverse)
library(GenomicSEM)
library(bigsnpr)

########## CHILDHOOD OBESITY

# http://egg-consortium.org/childhood-obesity-2019.html
# wget http://egg-consortium.org/Childhood_Obesity_2019/CHILDHOOD_OBESITY.TRANS_ANCESTRAL.RESULTS.txt.gz

childhoodObesity <- fread("./CHILDHOOD_OBESITY.TRANS_ANCESTRAL.RESULTS.txt", header=T, stringsAsFactor=F)
summary(as.numeric(childhoodObesity$AFR_N)) # median: 4592
summary(as.numeric(childhoodObesity$ASN_N)) # median: 578
summary(as.numeric(childhoodObesity$EUR_N)) # median: 21211
summary(as.numeric(childhoodObesity$AMR_N)) # median: 1979
childhoodObesity2 <- childhoodObesity[, c("CHR", "POS", "EA", "OA", "EUR_BETA", "EUR_P", "EUR_SE", "EUR_N", "EUR_FRQ")]
colnames(childhoodObesity2) <- c("chr", "pos", "a0", "a1", "beta", "p", "se", "N", "EUR_AF")
# a0 is being encoded as the effect allele for the matching to rsid
childhoodObesity2$a0 <- toupper(childhoodObesity2$a0)
childhoodObesity2$a1 <- toupper(childhoodObesity2$a1)
childhoodObesity2$chr <- as.numeric(childhoodObesity2$chr)
childhoodObesity2$pos <- as.numeric(childhoodObesity2$pos)
childhoodObesity2 <- subset(childhoodObesity2, childhoodObesity2$beta != "-")
dim(childhoodObesity2) #from 15,504,218 locations to the 8,411,444 locations with non NA effect sizers
# childhoodObesity2$beta[childhoodObesity2$beta == "-"] <- NA
childhoodObesity2$beta <- as.numeric(as.character(childhoodObesity2$beta))
summary(childhoodObesity2$beta)

# we need to get the rsids for the EGG Childhood Obesity summary stats:
# CHILDHOOD_OBESITY.TRANS_ANCESTRAL.RESULTS.txt

rsidmap <- fread("./rsid_position_allele_map.txt", header = F, stringsAsFactor=F)
# ./rsid_position_allele_map.txt
# 1	10019	rs775809821	TA	T
# 1	10039	rs978760828	A	C
# 1	10043	rs1008829651	T	A
# 1	10051	rs1052373574	A	G
head(rsidmap)
colnames(rsidmap) <- c("chr", "pos", "id", "a0", "a1")
rsidmap2 <- rsidmap[, c("id", "chr", "pos", "a0", "a1")]
str(rsidmap2)
str(childhoodObesity2)
rsidmap2$chr <- as.numeric(rsidmap2$chr)
# rsidmap2$pos <- as.numeric(rsidmap2$pos)
sumstats_matched <- bigsnpr::snp_match(childhoodObesity2, rsidmap2, match.min.prop = 0.5)
head(sumstats_matched)
sumstats_matched$`_NUM_ID_.ss` <- NULL
sumstats_matched$`_NUM_ID_` <- NULL
# let's remove the allele frequency column since we don't need it
# the allele frequency would be confusing in this new dataframe since the effect allele
# has been flipped for some of the SNPs.
sumstats_matched$`EUR_AF` <- NULL
head(sumstats_matched)

# after matching, a0 is still the effect allele, 
# but the betas were just multiplied by -1 if the alleles needed to be flipped

# # # We want to format the summary stats so that there are 5 columns:
# The rsid of the SNP.
# An A1 allele column, with A1 indicating the effect allele.
# An A2 allele column, with A2 indicating the non-effect allele.
# Either a logistic or continuous regression effect.
# The p-value associated with this effect.
final_ss <- sumstats_matched[, c("id", "a0", "a1", "beta", "p", "chr", "pos", "se", "N")]
colnames(final_ss) <- c("SNP", "A1", "A2", "effect", "P", "CHR", "POS", "SE", "N")

final_ss <- subset(final_ss, SE != "Inf")
print(paste0(round(nrow(final_ss)/nrow(childhoodObesity2)*100,0),"% SNPs retained: Nsnp = ",nrow(final_ss)))
# "77% SNPs retained: Nsnp = 6,461,721"
head(final_ss)
write.table(final_ss,"ChildhoodObesity_SEM.txt", sep = "\t", col.names=T, row.names=F, quote=F)
# these effect sizes for childhood obesity are from a LOGISTIC REGRESSION framework
# https://doi.org/10.1093/hmg/ddz161

########## WHR adj bmi
# https://portals.broadinstitute.org/collaboration/giant/index.php/GIANT_consortium_data_files#2016_GIANT_Body_Shape_Meta-analysis
# https://zenodo.org/record/1251813#.ZAe7PezMK3J
# whradjbmi.giant-ukbb.meta-analysis.combined.23May2018.txt
# Meta-analysis of waist-to-hip ratio adjusted for body mass index (whradjbmi) 
# in UK Biobank and GIANT data. Combined set of samples, max N = 694,649.
# wget https://zenodo.org/record/1251813/files/whradjbmi.giant-ukbb.meta-analysis.combined.23May2018.txt.gz?download=1
# wget https://zenodo.org/record/1251813/files/whradjbmi.giant-ukbb.meta-analysis.females.23May2018.txt.gz?download=1
# wget https://zenodo.org/record/1251813/files/whradjbmi.giant-ukbb.meta-analysis.males.23May2018.txt.gz?download=1

whr <- fread("./whradjbmi.giant-ukbb.meta-analysis.combined.23May2018.txt", header=T, stringsAsFactor=F)
# whr <- fread("./whradjbmi.giant-ukbb.meta-analysis.females.23May2018.txt", header=T, stringsAsFactor=F)
# whr <- fread("./whradjbmi.giant-ukbb.meta-analysis.males.23May2018.txt", header=T, stringsAsFactor=F)

final_ss <- whr[, c("SNP", "Tested_Allele", "Other_Allele", "BETA", "P", "CHR", "POS", "SE", "N")]
colnames(final_ss) <- c("SNP", "A1", "A2", "effect", "P", "CHR", "POS", "SE", "N")
final_ss <- final_ss %>% separate(SNP,into=c("SNP","del1","del2"), convert=TRUE, sep=":")
final_ss$del1 <- NULL
final_ss$del2 <- NULL
# sort( sapply(ls(),function(x){object.size(get(x))})) 
str(final_ss)
summary(final_ss$effect)
summary(final_ss$P)
final_ss <- subset(final_ss, SE != "Inf")
print(paste0(round(nrow(final_ss)/nrow(whr)*100,0),"% SNPs retained: Nsnp = ",nrow(final_ss)))
#  "100% SNPs retained: Nsnp = 27,375,636"
head(final_ss)
write.table(final_ss,"WHRadjBMI_SEM.txt", sep = "\t", col.names=T, row.names=F, quote=F)
# write.table(final_ss,"WHRadjBMI_females_SEM.txt", sep = "\t", col.names=T, row.names=F, quote=F)
# write.table(final_ss,"WHRadjBMI_males_SEM.txt", sep = "\t", col.names=T, row.names=F, quote=F)


########## WC adj bmi
# https://portals.broadinstitute.org/collaboration/giant/index.php/GIANT_consortium_data_files#GWAS_Anthropometric_2015_Waist_Summary_Statistics
# wget https://portals.broadinstitute.org/collaboration/giant/images/7/73/GIANT_2015_WCadjBMI_COMBINED_EUR.txt.gz
# wget https://portals.broadinstitute.org/collaboration/giant/images/9/92/GIANT_2015_WCadjBMI_FEMALES_EUR.txt.gz
# wget https://portals.broadinstitute.org/collaboration/giant/images/9/93/GIANT_2015_WCadjBMI_MALES_EUR.txt.gz

wc <- fread("./GIANT_2015_WCadjBMI_COMBINED_EUR.txt", header=F, stringsAsFactor=F, skip = 8)
# wc <- fread("./GIANT_2015_WCadjBMI_FEMALES_EUR.txt", header=F, stringsAsFactor=F, skip = 8)
# wc <- fread("./GIANT_2015_WCadjBMI_MALES_EUR.txt", header=F, stringsAsFactor=F, skip = 8)

colnames(wc) <- c("MarkerName", "Chr", "Pos", "Allele1", "Allele2", "FreqAllele1HapMapCEU", "b", "se", "p", "N")
head(wc)
# Allele1: The first allele (hg19 + strand). Where the regression coefficients (betas) are provided, 
# the first allele is the effect allele. Where betas are not provided (typically the 2010 data), 
# the first allele is the trait-increasing allele. 
final_ss <- wc[, c("MarkerName", "Allele1", "Allele2", "b", "p", "se", "N")]
# An A1 allele column, with A1 indicating the effect allele.
# An A2 allele column, with A2 indicating the non-effect allele.
colnames(final_ss) <- c("SNP", "A1", "A2", "effect", "P", "SE", "N")
str(final_ss)
summary(final_ss$effect)
summary(final_ss$P)
final_ss <- subset(final_ss, SE != "Inf")
print(paste0(round(nrow(final_ss)/nrow(wc)*100,0),"% SNPs retained: Nsnp = ",nrow(final_ss)))
#  "100% SNPs retained: Nsnp = 2,546,073"
head(final_ss)
write.table(final_ss,"WCadjBMI_SEM.txt", sep = "\t", col.names=T, row.names=F, quote=F)
# write.table(final_ss,"WCadjBMI_females_SEM.txt", sep = "\t", col.names=T, row.names=F, quote=F)
# write.table(final_ss,"WCadjBMI_males_SEM.txt", sep = "\t", col.names=T, row.names=F, quote=F)

########## HIP
# https://portals.broadinstitute.org/collaboration/giant/index.php/GIANT_consortium_data_files#GWAS_Anthropometric_2015_Waist_Summary_Statistics
# wget https://portals.broadinstitute.org/collaboration/giant/images/5/52/GIANT_2015_HIPadjBMI_COMBINED_EUR.txt.gz
# wget https://portals.broadinstitute.org/collaboration/giant/images/a/a4/GIANT_2015_HIPadjBMI_FEMALES_EUR.txt.gz
# wget https://portals.broadinstitute.org/collaboration/giant/images/5/5f/GIANT_2015_HIPadjBMI_MALES_EUR.txt.gz

# hip <- fread("./GIANT_2015_HIPadjBMI_COMBINED_EUR.txt", header=F, stringsAsFactor=F, skip = 8)
# hip <- fread("./GIANT_2015_HIPadjBMI_FEMALES_EUR.txt", header=T, stringsAsFactor=F)
hip <- fread("./GIANT_2015_HIPadjBMI_MALES_EUR.txt", header=T, stringsAsFactor=F)

colnames(hip) <- c("MarkerName", "Chr", "Pos", "Allele1", "Allele2", "FreqAllele1HapMapCEU", "b", "se", "p", "N")
head(hip)
# Allele1: The first allele (hg19 + strand). Where the regression coefficients (betas) are provided, 
# the first allele is the effect allele. Where betas are not provided (typically the 2010 data), 
# the first allele is the trait-increasing allele. 
final_ss <- hip[, c("MarkerName", "Allele1", "Allele2", "b", "p", "se", "N")]
# An A1 allele column, with A1 indicating the effect allele.
# An A2 allele column, with A2 indicating the non-effect allele.
colnames(final_ss) <- c("SNP", "A1", "A2", "effect", "P", "SE", "N")
str(final_ss)
summary(final_ss$effect)
summary(final_ss$P)
final_ss <- subset(final_ss, SE != "Inf")
print(paste0(round(nrow(final_ss)/nrow(hip)*100,0),"% SNPs retained: Nsnp = ",nrow(final_ss)))
#  "100% SNPs retained: Nsnp = 2,540,925"
head(final_ss)
# write.table(final_ss,"HIPadjBMI_SEM.txt", sep = "\t", col.names=T, row.names=F, quote=F)
# write.table(final_ss,"HIPadjBMI_females_SEM.txt", sep = "\t", col.names=T, row.names=F, quote=F)
write.table(final_ss,"HIPadjBMI_males_SEM.txt", sep = "\t", col.names=T, row.names=F, quote=F)

########## BMI
# https://portals.broadinstitute.org/collaboration/giant/index.php/GIANT_consortium_data_files#GWAS_Anthropometric_2015_Waist_Summary_Statistics
# https://zenodo.org/record/1251813#.ZAe-8uzMK3J
# bmi.giant-ukbb.meta-analysis.combined.23May2018.txt
# Meta-analysis of body mass index (bmi) in UK Biobank and GIANT data. 
# Combined set of samples, max N = 806,834.
# wget https://zenodo.org/record/1251813/files/bmi.giant-ukbb.meta-analysis.combined.23May2018.txt.gz?download=1
# Meta-analysis of bmi in UK Biobank and GIANT data. Female samples only, max N = 434,794.
# wget https://zenodo.org/record/1251813/files/bmi.giant-ukbb.meta-analysis.females.23May2018.txt.gz?download=1
# Meta-analysis of bmi in UK Biobank and GIANT data. Male samples only, max N = 374,756.
# wget https://zenodo.org/record/1251813/files/bmi.giant-ukbb.meta-analysis.males.23May2018.txt.gz?download=1

bmi <- fread("./bmi.giant-ukbb.meta-analysis.combined.23May2018.txt", header=T, stringsAsFactor=F)
# bmi <- fread("./bmi.giant-ukbb.meta-analysis.females.23May2018.txt", header=T, stringsAsFactor=F)
# bmi <- fread("./bmi.giant-ukbb.meta-analysis.males.23May2018.txt", header=T, stringsAsFactor=F)

head(bmi)
final_ss <- bmi[, c("SNP", "Tested_Allele", "Other_Allele", "BETA", "P", "CHR", "POS", "SE", "N")]
# An A1 allele column, with A1 indicating the effect allele.
# An A2 allele column, with A2 indicating the non-effect allele.
colnames(final_ss) <- c("SNP", "A1", "A2", "effect", "P", "CHR", "POS", "SE", "N")
final_ss <- final_ss %>% separate(SNP,into=c("SNP","del1","del2"), convert=TRUE, sep=":")
final_ss$del1 <- NULL
final_ss$del2 <- NULL
str(final_ss)
summary(final_ss$effect)
summary(final_ss$P)
final_ss <- subset(final_ss, SE != "Inf")
print(paste0(round(nrow(final_ss)/nrow(bmi)*100,0),"% SNPs retained: Nsnp = ",nrow(final_ss)))
#  "100% SNPs retained: Nsnp = 27,381,302"
head(final_ss)
write.table(final_ss,"BMI_SEM.txt", sep = "\t", col.names=T, row.names=F, quote=F)
# write.table(final_ss,"BMI_females_SEM.txt", sep = "\t", col.names=T, row.names=F, quote=F)
# write.table(final_ss,"BMI_males_SEM.txt", sep = "\t", col.names=T, row.names=F, quote=F)

########## HEIGHT
# https://www.joelhirschhornlab.org/giant-consortium-results
# wget https://504394d8-624a-4827-9f25-95a83cd9675a.filesusr.com/archives/1d0101_394c2d3120ba4d0a9f6326ed56ff8854.gz?dn=GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_EUR.gz

# Yengo et al., 2022 https://doi.org/10.1038/s41586-022-05275-y
# Ancestry-group-specific GWAS meta-analyses of 173 studies of EUR, 
# 56 studies of EAS, 29 studies of AFR, 11 studies of HIS and 12 studies of SAS. 

# height <- fread("./GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_EUR.txt", header=T, stringsAsFactor=F)
# head(height)
# # I am assuming EFFECT_ALLELE is the effect allele per the readme
# final_ss <- height[, c("RSID", "EFFECT_ALLELE", "OTHER_ALLELE", "BETA", "P", "SE", "N")]
# colnames(final_ss) <- c("SNP", "A1", "A2", "effect", "P", "SE", "N")
# summary(final_ss$N)
# # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# #  344 1540077 1587709 1452363 1596610 1597374
# dim(final_ss)
# # 1,373,020 is the number of SNPs
# str(final_ss)
# summary(final_ss$effect)
# final_ss$P <- as.numeric(final_ss$P)
# summary(final_ss$P)
# final_ss <- subset(final_ss, SE != "Inf")
# print(paste0(round(nrow(final_ss)/nrow(height)*100,0),"% SNPs retained: Nsnp = ",nrow(final_ss)))
# #  "100% SNPs retained: Nsnp = 1,372,608"
# head(final_ss)
# write.table(final_ss,"HEIGHT_SEM.txt", sep = "\t", col.names=T, row.names=F, quote=F)

# other height GWAS worth considering to boost SNP count
# https://portals.broadinstitute.org/collaboration/giant/index.php/GIANT_consortium_data_files#BMI_and_Height_GIANT_and_UK_BioBank_Meta-analysis_Summary_Statistics
# wget https://portals.broadinstitute.org/collaboration/giant/images/c/c8/Meta-analysis_Locke_et_al%2BUKBiobank_2018_UPDATED.txt.gz
# wget https://portals.broadinstitute.org/collaboration/giant/images/6/63/Meta-analysis_Wood_et_al%2BUKBiobank_2018.txt.gz
# cite: Yengo L, Sidorenko J, Kemper KE, Zheng Z, Wood AR, Weedon MN, Frayling TM, Hirschhorn J, Yang J, Visscher PM, GIANT Consortium. (2018). Meta-analysis of genome-wide association studies for height and body mass index in ~700,000 individuals of European ancestry. Biorxiv.

# Wood et al., 2014 https://doi.org/10.1038/ng.3097
# Defining the role of common variation in the genomic and biological architecture of adult human height
# combined the height summary association statistics from 79 GWAS in a meta-analysis of 253,288 individuals
height <- fread("./Meta-analysis_Wood_et_al+UKBiobank_2018.txt", header=T, stringsAsFactor=F)
dim(height)
# 2,334,001 is the number of SNPs
summary(height$N)
 #   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 # 485983  695748  704526  695648  708089  709706
head(height)
# I am assuming Tested_Allele is the effect allele per the readme
final_ss <- height[, c("SNP", "Tested_Allele", "Other_Allele", "BETA", "P", "SE", "N")]
colnames(final_ss) <- c("SNP", "A1", "A2", "effect", "P", "SE", "N")
summary(final_ss$N)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  344 1540077 1587709 1452363 1596610 1597374
dim(final_ss)
# 1,373,020 is the number of SNPs
str(final_ss)
summary(final_ss$effect)
final_ss$P <- as.numeric(final_ss$P)
summary(final_ss$P)
final_ss <- subset(final_ss, SE != "Inf")
print(paste0(round(nrow(final_ss)/nrow(height)*100,0),"% SNPs retained: Nsnp = ",nrow(final_ss)))
#  "100% SNPs retained: Nsnp = 1,372,608"
head(final_ss)
write.table(final_ss,"HEIGHT_SEM.txt", sep = "\t", col.names=T, row.names=F, quote=F)


# Locke et al., 2015 https://doi.org/10.1038/nature14177
# Genetic studies of body mass index yield new insights for obesity biology
# In stage 1 we performed meta-analysis of 80 GWAS (n = 234,069); 
# and stage 2 incorporated data from 34 additional studies (n = 88,137) genotyped using Metabochip
# height <- fread("./Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt", header=T, stringsAsFactor=F)
# dim(height)
# # 2,336,269 is the number of SNPs
# summary(height$N)
# # 688,566 is the median N
#  #   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  # 485648  679885  688566  683365  691547  795640


########## INFANT HEAD CIRCUMFERENCE
# http://egg-consortium.org/Head-circumference-2022.html
# wget http://egg-consortium.org/HC2/HC_INFANT.filteredDf0.txt.gz

ihc <- fread("./HC_INFANT.filteredDf0.txt", header=T, stringsAsFactor=F)
head(ihc)
# Allele1: Effect allele
# Allele2: Other allele
final_ss <- ihc[, c("ID", "Allele1", "Allele2", "Effect", "P.value", "StdErr", "TotalSampleSize")]
# An A1 allele column, with A1 indicating the effect allele.
# An A2 allele column, with A2 indicating the non-effect allele.
colnames(final_ss) <- c("SNP", "A1", "A2", "effect", "P", "SE", "N")
str(final_ss)
summary(final_ss$effect)
summary(final_ss$P)
final_ss$A1 <- toupper(final_ss$A1)
final_ss$A2 <- toupper(final_ss$A2)
final_ss <- subset(final_ss, SE != "Inf")
final_ss <- subset(final_ss, SNP != ".")
print(paste0(round(nrow(final_ss)/nrow(ihc)*100,0),"% SNPs retained: Nsnp = ",nrow(final_ss)))
#  "88% SNPs retained: Nsnp = 10,143,709
head(final_ss)
write.table(final_ss,"IHC_SEM.txt", sep = "\t", col.names=T, row.names=F, quote=F)

########## BIRTH LENGTH
# http://egg-consortium.org/birth-length.html
# wget http://egg-consortium.org/Birth_Length/EGG-GWAS-BL.txt.gz

bl <- fread("./EGG-GWAS-BL.txt", header=T, stringsAsFactor=F)
head(bl)
# Allele1: Effect allele
# Allele2: Other allele
final_ss <- bl[, c("RSID", "EA", "NEA", "BETA", "P", "SE", "N")]
# An A1 allele column, with A1 indicating the effect allele.
# An A2 allele column, with A2 indicating the non-effect allele.
colnames(final_ss) <- c("SNP", "A1", "A2", "effect", "P", "SE", "N")
str(final_ss)
summary(final_ss$effect)
summary(final_ss$P)
final_ss$A1 <- toupper(final_ss$A1)
final_ss$A2 <- toupper(final_ss$A2)
final_ss <- subset(final_ss, SE != "Inf")
final_ss <- subset(final_ss, SNP != ".")
print(paste0(round(nrow(final_ss)/nrow(bl)*100,0),"% SNPs retained: Nsnp = ",nrow(final_ss)))
#  "100% SNPs retained: Nsnp = 2,201,971
head(final_ss)
write.table(final_ss,"BL_SEM.txt", sep = "\t", col.names=T, row.names=F, quote=F)

########## BIRTH WEIGHT
# http://egg-consortium.org/birth-weight-2019.html
# European-only meta-analysis of own birth weight in up to 298,142 individuals;
# wget http://egg-consortium.org/BW5/Fetal_BW_European_meta.NG2019.txt.gz

bw <- fread("./Fetal_BW_European_meta.NG2019.txt", header=T, stringsAsFactor=F)
head(bw)
# Allele1: Effect allele
# Allele2: Other allele
final_ss <- bw[, c("rsid", "ea", "nea", "beta", "p", "se", "n")]
# An A1 allele column, with A1 indicating the effect allele.
# An A2 allele column, with A2 indicating the non-effect allele.
colnames(final_ss) <- c("SNP", "A1", "A2", "effect", "P", "SE", "N")
str(final_ss)
summary(final_ss$effect)
summary(final_ss$P)
final_ss$A1 <- toupper(final_ss$A1)
final_ss$A2 <- toupper(final_ss$A2)
final_ss <- subset(final_ss, SE != "Inf")
final_ss <- subset(final_ss, SNP != ".")
print(paste0(round(nrow(final_ss)/nrow(bw)*100,0),"% SNPs retained: Nsnp = ",nrow(final_ss)))
#  "100% SNPs retained: Nsnp = 13,536,081
head(final_ss)
write.table(final_ss,"BW_SEM.txt", sep = "\t", col.names=T, row.names=F, quote=F)

########## LONGITUDINAL BMI
# https://www.sciencedirect.com/science/article/pii/S0002929722000490

# primary goal is to test 
# (1) the mean effect of genotype, beta_g, i.e., whether a genotype shifts the mean of the biomarker trajectory; 
# (2) the WS variance effect of genotype, tau_g, i.e., whether a genotype changes the within-subject variation of the biomarker trajectory around its mean;

# Table 2 statistics for the BMI GWAS
# BMI mean (SD): 28.3 (5.7)
# sample size: 144,414
# female %: 54.9
# age (SD): 57.3 (9.9)

# https://ucla.app.box.com/v/trajgwassummary
# hit the download button at https://ucla.app.box.com/v/trajgwassummary/file/875801611119

# chr: chromosome
# pos: position
# snpid: RefSNP ID of each SNP
# varid: Variant ID denoting chromosome , position, and reference allele, and alternative allele
# hwepval: p-value of Hardy-Weinberg equilibrium test. SNPs with hwepval less than 10^{-10} have been removed.
# maf: minor allele frequency. SNPs with maf less than 0.002 have been removed.
# infoscore: information score. SNPs with infoscore below 0.3 have been removed.
# betapval: p-value for beta_g adjusted by the saddlepoint approximation.
# betachisqpval: p-value of the ordinary chi-square test for beta_g
# betadir: direction of effect toward beta_g, i.e., sign of the score. 1 for positive and -1 for negative.
# taupval: p-value for tau_g adjusted by the saddlepoint approximation.
# tauchisqpval: p-value of the ordinary chi-square test for tau_g
# taudir: direction of effect toward tau_g, i.e., sign of the score. 1 for positive and -1 for negative.

library(data.table)
library(tidyverse)
bmitraj <- fread("./ukbbgen_bmi_clean.csv", header=T, stringsAsFactor=F, sep = ",")
bmitraj <- bmitraj %>% separate(varid, into=c("location", "Ref", "Alt"), convert=TRUE, sep="_")
tau_ukbbgen <- bmitraj[,c("snpid","Ref", "Alt","maf","infoscore","taupval","taudir")]
tau_ukbbgen <- tau_ukbbgen %>% rename("SNP" = "snpid","REF" = "Ref","ALT" = "Alt","MAF" = "maf","INFO" = "infoscore","PVAL" = "taupval","BETAdir" = "taudir")
# munge will do:
# Interpreting the REF column as the A1 column, with A1 indicating the effect allele.
# Interpreting the ALT column as the A2 column, with A2 indicating the non-effect allele.
n_bmi <- 144414
tau_ukbbgen$Z <- sign(tau_ukbbgen$BETAdir)*qnorm(1 - tau_ukbbgen$PVAL/2)
# https://www.biostars.org/p/319584/
# Beta = z / sqrt(2p(1 − p)(n + z^2))
tau_ukbbgen$BETA <- tau_ukbbgen$Z/sqrt(2*tau_ukbbgen$MAF*(1-tau_ukbbgen$MAF)*(n_bmi+tau_ukbbgen$Z^2))
# SE = 1 / sqrt(2p(1 − p)(n + z^2))
tau_ukbbgen$SE <- 1/sqrt(2*tau_ukbbgen$MAF*(1-tau_ukbbgen$MAF)*(n_bmi+tau_ukbbgen$Z^2))
# you get the same answer if you se = beta/Z
# tau_ukbbgen$SE <- tau_ukbbgen$BETA/tau_ukbbgen$Z
summary(tau_ukbbgen$Z)
summary(tau_ukbbgen$BETA)
summary(tau_ukbbgen$SE)

head(tau_ukbbgen)
tail(tau_ukbbgen)

write.table(tau_ukbbgen, file = "bmiTrajTau_SEM.txt", col.names=TRUE,quote = FALSE, row.names = FALSE)

## also output summmary stats for bmiTrajBeta
# (1) the mean effect of genotype, beta_g, i.e., whether a genotype shifts the mean of the biomarker trajectory; 

bmitraj <- fread("./ukbbgen_bmi_clean.csv", header=T, stringsAsFactor=F, sep = ",")
bmitraj <- bmitraj %>% separate(varid, into=c("location", "Ref", "Alt"), convert=TRUE, sep="_")
beta_ukbbgen <- bmitraj[,c("snpid","Ref", "Alt","maf","infoscore","betapval","betadir")]
beta_ukbbgen <- beta_ukbbgen %>% rename("SNP" = "snpid","REF" = "Ref","ALT" = "Alt","MAF" = "maf","INFO" = "infoscore","PVAL" = "betapval","BETAdir" = "betadir")
# munge will do:
# Interpreting the REF column as the A1 column, with A1 indicating the effect allele.
# Interpreting the ALT column as the A2 column, with A2 indicating the non-effect allele.
n_bmi <- 144414
beta_ukbbgen$Z <- sign(beta_ukbbgen$BETAdir)*qnorm(1 - beta_ukbbgen$PVAL/2)
# https://www.biostars.org/p/319584/
# Beta = z / sqrt(2p(1 − p)(n + z^2))
beta_ukbbgen$BETA <- beta_ukbbgen$Z/sqrt(2*beta_ukbbgen$MAF*(1-beta_ukbbgen$MAF)*(n_bmi+beta_ukbbgen$Z^2))
# SE = 1 / sqrt(2p(1 − p)(n + z^2))
beta_ukbbgen$SE <- 1/sqrt(2*beta_ukbbgen$MAF*(1-beta_ukbbgen$MAF)*(n_bmi+beta_ukbbgen$Z^2))
# you get the same answer if you se = beta/Z
# beta_ukbbgen$SE <- beta_ukbbgen$BETA/beta_ukbbgen$Z
summary(beta_ukbbgen$Z)
summary(beta_ukbbgen$BETA)
beta_ukbbgen2 <- beta_ukbbgen[which(is.na(beta_ukbbgen$BETA) == F),]
summary(beta_ukbbgen2$Z)
summary(beta_ukbbgen2$SE)

head(beta_ukbbgen2)
tail(beta_ukbbgen2)

write.table(beta_ukbbgen2, file = "bmiTrajBeta_SEM.txt", col.names=TRUE,quote = FALSE, row.names = FALSE)

########## BODY FAT DISTRIBUTION 
# segmental bio-electrical impedance analysis (sBIA)
# arm fat ratio (AFR): The proportions of body fat distributed to the arms
# leg fat ratio (LFR): The proportions of body fat distributed to the legs
# trunk fat ratio (TFR): The proportions of body fat distributed to the trunk
# calculated by dividing the fat mass per compartment with the total body fat mass for each participant
# https://www.nature.com/articles/s41467-018-08000-4#data-availability

# wget https://myfiles.uu.se/ssf/s/readFile/share/3993/1270878243748486898/publicLink/GWAS_summary_stats_ratios.zip
# had to unzip the files on a mac, rather than on RC

# "After filtering, 116,138 participants remained in the discovery cohort and 246,361 in the replication cohort"
# Full Cohort N = 116138 + 246361 = 362499
# But we see from the NMISS column that the sample size is closer to 191163 females + 163201 males = 354364 total

# CHR - Chromosome number
# SNP - RSID
# BP - Basepair
# A1 - Minor/effect allele
# TEST - Type of test (ADD = Additive)
# NMISS - Number of informative participants 
# BETA - Regression coefficient
# STAT - Coefficient t-statistic
# P - Asymptotic p-value for t-statistic

ref <- fread("reference.1000G.maf.0.005.txt")
# 9,667,224 variants are in ref
colnames(ref) <- paste0(colnames(ref), "_ref")

# # setDTthreads(16)

# impedance_gwas <- fread("impedance_gwas/Full_cohort/AFR_females_full.sumstats", header=T, stringsAsFactor=F)
# impedance_gwas <- fread("impedance_gwas/Full_cohort/AFR_males_full.sumstats", header=T, stringsAsFactor=F)
# impedance_gwas <- fread("impedance_gwas/Full_cohort/TFR_females_full.sumstats", header=T, stringsAsFactor=F)
# impedance_gwas <- fread("impedance_gwas/Full_cohort/TFR_males_full.sumstats", header=T, stringsAsFactor=F)
# impedance_gwas <- fread("impedance_gwas/Full_cohort/LFR_females_full.sumstats", header=T, stringsAsFactor=F)
impedance_gwas <- fread("impedance_gwas/Full_cohort/LFR_males_full.sumstats", header=T, stringsAsFactor=F)
head(impedance_gwas)
dim(impedance_gwas)

#subset impedance_gwas to reference file, so we have a smaller dataset (and things will go faster)
impedance_gwas_sub <- subset(impedance_gwas, impedance_gwas$SNP %in% ref$SNP)
dim(impedance_gwas_sub)
# 8,842,481 variants in impedance_gwas_sub

# get the sumstats ready for the merge
summary(as.factor(impedance_gwas_sub$TEST))
impedance_gwas_sub$TEST <- NULL
# se = beta/Z
impedance_gwas_sub$se <- impedance_gwas_sub$BETA/impedance_gwas_sub$STAT
impedance_gwas_sub$STAT <- NULL
impedance_gwas_sub$nonEffectAllele <- NA
impedance_gwas_sub <- impedance_gwas_sub[,c("SNP","CHR", "BP","A1","nonEffectAllele","BETA","P","se","NMISS")]
impedance_gwas_sub <- impedance_gwas_sub %>% rename("SNP" = "SNP", "chr" = "CHR", "pos" = "BP", "A1" = "A1", "nonEffectAllele" = "nonEffectAllele", "beta" = "BETA", "P" = "P",  "se" = "se", "N" = "NMISS")

# merge with the tgp reference file
m_impedance_gwas <- merge(impedance_gwas_sub, ref, by.x = "SNP", by.y = "SNP_ref")
dim(m_impedance_gwas)
head(m_impedance_gwas)
# if A1 and A1_ref match, set nonEffectAllele equal to A2_ref
m_impedance_gwas$nonEffectAllele[m_impedance_gwas$A1 == m_impedance_gwas$A1_ref] <- m_impedance_gwas$A2_ref[m_impedance_gwas$A1 == m_impedance_gwas$A1_ref]
m_impedance_gwas[sample(x = 1:nrow(m_impedance_gwas), size = 10, replace = F),c("SNP", "A1", "nonEffectAllele", "A1_ref", "A2_ref")]
# if A1 and A2_ref match, set nonEffectAllele equal to A1_ref
m_impedance_gwas$nonEffectAllele[m_impedance_gwas$A1 == m_impedance_gwas$A2_ref] <- m_impedance_gwas$A1_ref[m_impedance_gwas$A1 == m_impedance_gwas$A2_ref]
m_impedance_gwas[sample(x = 1:nrow(m_impedance_gwas), size = 10, replace = F),c("SNP", "A1", "nonEffectAllele", "A1_ref", "A2_ref")]
# remove the rows when there wasn't a match of the alleles
sum(is.na(m_impedance_gwas$nonEffectAllele))
m_impedance_gwas <- m_impedance_gwas[which(is.na(m_impedance_gwas$nonEffectAllele) == FALSE),]

# # # We want to format the summary stats:
# An A1 allele column, with A1 indicating the effect allele.
# An A2 allele column, with A2 indicating the non-effect allele.
final_ss <- m_impedance_gwas[, c("SNP", "A1", "nonEffectAllele", "beta", "P", "chr", "pos", "se", "N")]
colnames(final_ss) <- c("SNP", "A1", "A2", "effect", "P", "CHR", "POS", "SE", "N")

head(final_ss)
dim(final_ss)
summary(final_ss$N)
# females
 #   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 # 182115  190161  191163  190525  191590  191738 
# males
 #    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 # 155462  162345  163201  162656  163565  163692

# write.table(final_ss, "AFR_females_SEM.txt", sep = "\t", col.names=T, row.names=F, quote=F)
# write.table(final_ss, "AFR_males_SEM.txt", sep = "\t", col.names=T, row.names=F, quote=F)
# write.table(final_ss, "TFR_females_SEM.txt", sep = "\t", col.names=T, row.names=F, quote=F)
# write.table(final_ss, "TFR_males_SEM.txt", sep = "\t", col.names=T, row.names=F, quote=F)
# write.table(final_ss, "LFR_females_SEM.txt", sep = "\t", col.names=T, row.names=F, quote=F)
write.table(final_ss, "LFR_males_SEM.txt", sep = "\t", col.names=T, row.names=F, quote=F)


# awk -F'\t' '{print $20}' CAD_GWAS_primary_discovery_meta.tsv | sort -n | awk '{
#     numbers[NR] = $1
# }

# END {
#     count = NR
#     middle = int((count + 1) / 2)
#     if (count % 2 == 0) {
#         median = (numbers[middle] + numbers[middle + 1]) / 2
#     } else {
#         median = numbers[middle]
#     }
#     print "Median:", median
# }'
EOF

####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################

# need to do download the other files needed to run Genomic SEM
# https://utexas.app.box.com/s/vkd36n197m8klbaio3yzoxsee6sxo11v

cat << 'EOF' >munge_script.R
#load in packages
library(data.table)
library(tidyverse)
library(GenomicSEM)
library(bigsnpr)

# # N_CO
# Values aggregated from Table S1 from https://doi.org/10.1093/hmg/ddz161
nr.cases <- 9116
nr.controls <- 13292
N_CO <- nr.cases + nr.controls
# "Obesity is having a dramatic impact on modern societies, 
# leading to substantial health issues, with an overall 
# prevalence among children already >20% in many populations, including the USA"
CO_P1 <- 0.2
# CO_K1 <- nr.cases/(nr.cases+nr.controls) #0.406819
CO_K1 <- 0.5 #since it was a meta-analysis we are using the sum of effective sample sizes, and therefore put 0.5 here.

# Removed because h2 Z: 3.91 which is less than 4.
######## bmiTrajTau (longitudinal BMI Tau, WI variability) - continuous

awk -F'\t' '{print $9}' AFR_males_SEM.txt | sort -n | awk '{
    numbers[NR] = $1
}

END {
    count = NR
    middle = int((count + 1) / 2)
    if (count % 2 == 0) {
        median = (numbers[middle] + numbers[middle + 1]) / 2
    } else {
        median = numbers[middle]
    }
    print "Median N:", median
}'

# 1. BMI_females (body mass index females) - continuous - Median N : 262817
# 2. BMI_males (body mass index males) - continuous - Median N: 221863
# 3. HIPadjBMI_females (hip adjusted for BMI females) - continuous - Median N: 86748
# 4. HIPadjBMI_males (hip adjusted for BMI males) - continuous - Median N: 56944
# 5. WCadjBMI_females (waist circumference adjusted for BMI females) - continuous - Median N: 91298
# 6. WCadjBMI_males (waist circumference adjusted for BMI males) - continuous - Median N: 61075
# 7. WHRadjBMI_females (waist-to-hip ratio adjusted for BMI females) - continuous - Median N: 262759
# 8. WHRadjBMI_males (waist-to-hip ratio adjusted for BMI males) - continuous - Median N: 221804
# 9. BL (birth length) - continuous - Median N: 22183.9
# 10. BW (birth weight) - continuous - Median N: 281141
# 11. CO (childhood obesity) - binary - Median N: 21224
# 12. IHC (infant head circumference) - continuous - Median N: 20904
# 13. Height - continuous - Median N: 704526
# 14. bmiTrajBeta (longitudinal BMI Beta, shifts the mean of the biomarker trajectory) - continuous - N: 144414
# 15. TFR_females (bio-electrical impedance trunk fat ratio females) - continuous - Median N: 191132
# 16. TFR_males (bio-electrical impedance trunk fat ratio males) - continuous - Median N: 163158
# 17. AFR_females (bio-electrical impedance arm fat ratio females) - continuous - Median N: 191163
# 18. AFR_males (bio-electrical impedance arm fat ratio males) - continuous - Median N: 163201

N_sumstat_files <- 18

# # Step 1: Munge the summary statistics

# All have SNP-specific sum of effective sample sizes as a column in the summary stats 
NN <- rep(NA, N_sumstat_files)
NN[c(14)] <- 144414

munge(c("BMI_females_SEM.txt", "BMI_males_SEM.txt", "HIPadjBMI_females_SEM.txt", "HIPadjBMI_males_SEM.txt", 
	"WCadjBMI_females_SEM.txt", "WCadjBMI_males_SEM.txt", "WHRadjBMI_females_SEM.txt", "WHRadjBMI_males_SEM.txt", 
	"BL_SEM.txt", "BW_SEM.txt", "ChildhoodObesity_SEM.txt", "IHC_SEM.txt", "HEIGHT_SEM.txt", "bmiTrajBeta_SEM.txt", 
	"TFR_females_SEM.txt", "TFR_males_SEM.txt", "AFR_females_SEM.txt", "AFR_males_SEM.txt"),
        "reference.1000G.maf.0.005.txt", 
        trait.names = c("BMI_females", "BMI_males", "HIPadjBMI_females", "HIPadjBMI_males", 
	"WCadjBMI_females", "WCadjBMI_males", "WHRadjBMI_females", "WHRadjBMI_males", 
	"BL", "BW", "CO", "IHC", "Height", "bmiTrajBeta", 
	"TFR_females", "TFR_males", "AFR_females", "AFR_males"), 
        N = NN, info.filter = 0.9, maf.filter = 0.01)

# # Step 2: Run multivariable LDSC

#enter sample prevalence of .5 to reflect that all traits were munged using the sum of effective sample size
sample.prev <- rep(NA, N_sumstat_files); sample.prev[11] <- CO_K1

#vector of population prevalences
population.prev <- rep(NA, N_sumstat_files); population.prev[11] <- CO_P1

#the folder of LD scores
ld <- "eur_w_ld_chr/"

#the folder of LD weights [typically the same as folder of LD scores]
wld <- "eur_w_ld_chr/"

trait.names <- c("BMI_females", "BMI_males", "HIPadjBMI_females", "HIPadjBMI_males", 
	"WCadjBMI_females", "WCadjBMI_males", "WHRadjBMI_females", "WHRadjBMI_males", 
	"BL", "BW", "CO", "IHC", "Height", "bmiTrajBeta", 
	"TFR_females", "TFR_males", "AFR_females", "AFR_males")
LDSCoutput <- ldsc(traits = c("BMI_females.sumstats.gz", "BMI_males.sumstats.gz", "HIPadjBMI_females.sumstats.gz", "HIPadjBMI_males.sumstats.gz", 
	"WCadjBMI_females.sumstats.gz", "WCadjBMI_males.sumstats.gz", "WHRadjBMI_females.sumstats.gz", "WHRadjBMI_males.sumstats.gz", 
	"BL.sumstats.gz", "BW.sumstats.gz", "CO.sumstats.gz", "IHC.sumstats.gz", "Height.sumstats.gz", "bmiTrajBeta.sumstats.gz", 
	"TFR_females.sumstats.gz", "TFR_males.sumstats.gz", "AFR_females.sumstats.gz", "AFR_males.sumstats.gz"), 
	sample.prev = sample.prev, population.prev = population.prev, 
	ld = ld, wld = wld, trait.names = trait.names, stand = T)

save(LDSCoutput, file="LDSCoutput.Rdata")
# load("LDSCoutput.Rdata")

# Save the genetic covariance matrix
print(LDSCoutput$S)
genetic_covariance_matrix <- LDSCoutput$S
rownames(genetic_covariance_matrix) <- colnames(genetic_covariance_matrix)
write.table(genetic_covariance_matrix, "genetic_covariance_matrix.txt", row.names = T, col.names = T, sep = "\t", quote = F)

# Save the genetic correlation matrix
print(LDSCoutput$S_Stand)
genetic_correlation_matrix <- LDSCoutput$S_Stand
rownames(genetic_correlation_matrix) <- colnames(genetic_correlation_matrix)
write.table(genetic_correlation_matrix, "genetic_correlation_matrix.txt", row.names = T, col.names = T, sep = "\t", quote = F)


# # Step 3: Prepare the summary statistics for GWAS
files <- c("BMI_females_SEM.txt", "BMI_males_SEM.txt", "HIPadjBMI_females_SEM.txt", "HIPadjBMI_males_SEM.txt", 
	"WCadjBMI_females_SEM.txt", "WCadjBMI_males_SEM.txt", "WHRadjBMI_females_SEM.txt", "WHRadjBMI_males_SEM.txt", 
	"BL_SEM.txt", "BW_SEM.txt", "ChildhoodObesity_SEM.txt", "IHC_SEM.txt", "HEIGHT_SEM.txt", "bmiTrajBeta_SEM.txt", 
	"TFR_females_SEM.txt", "TFR_males_SEM.txt", "AFR_females_SEM.txt", "AFR_males_SEM.txt")
ref = "reference.1000G.maf.0.005.txt"
se.logit = rep(F, N_sumstat_files); se.logit[11] <- T
linprob = rep(F, N_sumstat_files)
info.filter=.6
maf.filter=0.01
# OLS=NULL
# N=NULL
# betas=NULL

setDTthreads(N_sumstat_files)
OVERWEIGHT_sumstats <- sumstats(files=files, ref=ref, trait.names=trait.names, 
	se.logit=se.logit, OLS=NULL, linprob=linprob, N=NULL, betas=NULL, 
	info.filter=info.filter, maf.filter=maf.filter, keep.indel=FALSE, parallel=F,
	cores=N_sumstat_files)

save(OVERWEIGHT_sumstats, file="OVERWEIGHT_sumstats.Rdata")
# load("OVERWEIGHT_sumstats.Rdata")
EOF

### ---------------------------------------------------------------- ###

cat << 'EOF' >run_munge_script.sh
#!/bin/bash
#SBATCH --qos=preemptable
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --mem=110000
#SBATCH --time=08:30:00
#SBATCH -J c_GSEM
#SBATCH -o run_munge_script.out


date
startTime="$(date +%s)"

module load anaconda
conda activate nmr_prs
date

Rscript munge_script.R

scontrol show job $SLURM_JOB_ID

date
endTime="$(date +%s)"
diff=$(echo "$endTime-$startTime" |bc)
diff=$(($diff + 0))

if [ $diff -lt "$((60))" ]
then
  printf %.3f "$((1000 * $diff/1))e-3"; echo " seconds"
elif [ $diff -lt "$((60 * 60))" ]
then
  printf %.3f "$((1000 * $diff/60))e-3"; echo " minutes"
elif [ $diff -lt "$((60 * 60 * 24))" ]
then
  printf %.3f "$((1000 * $diff/3600))e-3"; echo " hours"
else
  printf %.3f "$((1000 * $diff/86400))e-3"; echo " days"
fi

# * # * # SCRIPT END # * # * #
EOF
# sbatch run_munge_script.sh

grep "SNPs are left in the summary statistics" BMI_females_BMI_males_HIPadjBMI_females_HIPadjBMI_males_WCadjBMI_females_WCadjBMI_males_WHRadjBMI_fe_sumstats.log

# 7,475,483 SNPs are left in the summary statistics file BMI_females_SEM.txt after QC and merging with the reference file.
# 7,479,319 SNPs are left in the summary statistics file BMI_males_SEM.txt after QC and merging with the reference file.
# 2,359,655 SNPs are left in the summary statistics file HIPadjBMI_females_SEM.txt after QC and merging with the reference file.
# 2,123,803 SNPs are left in the summary statistics file HIPadjBMI_males_SEM.txt after QC and merging with the reference file.
# 2,362,480 SNPs are left in the summary statistics file WCadjBMI_females_SEM.txt after QC and merging with the reference file.
# 2,223,277 SNPs are left in the summary statistics file WCadjBMI_males_SEM.txt after QC and merging with the reference file.
# 7,468,408 SNPs are left in the summary statistics file WHRadjBMI_females_SEM.txt after QC and merging with the reference file.
# 7,472,397 SNPs are left in the summary statistics file WHRadjBMI_males_SEM.txt after QC and merging with the reference file.
# 2,137,955 SNPs are left in the summary statistics file BL_SEM.txt after QC and merging with the reference file.
# 7,943,952 SNPs are left in the summary statistics file BW_SEM.txt after QC and merging with the reference file.
# 6,086,078 SNPs are left in the summary statistics file ChildhoodObesity_SEM.txt after QC and merging with the reference file.
# 7,900,160 SNPs are left in the summary statistics file IHC_SEM.txt after QC and merging with the reference file.
# 2,259,169 SNPs are left in the summary statistics file HEIGHT_SEM.txt after QC and merging with the reference file.
# 5,702,132 SNPs are left in the summary statistics file bmiTrajBeta_SEM.txt after QC and merging with the reference file.
# 7,798,333 SNPs are left in the summary statistics file TFR_females_SEM.txt after QC and merging with the reference file.
# 7,800,162 SNPs are left in the summary statistics file TFR_males_SEM.txt after QC and merging with the reference file.
# 7,798,331 SNPs are left in the summary statistics file AFR_females_SEM.txt after QC and merging with the reference file.
# 7,800,162 SNPs are left in the summary statistics file AFR_males_SEM.txt after QC and merging with the reference file.

cat << 'EOF' >EFA_CFA.R

library(data.table)
library(R.utils)
require(Matrix)
require(stats)
require(dplyr)
require(stringr)
library(GenomicSEM)

print("loading LDSCoutput.Rdata")
load("LDSCoutput.Rdata")
print("loading OVERWEIGHT_sumstats.Rdata")
load("OVERWEIGHT_sumstats.Rdata")


N_sumstat_files <- 18

# Values aggregated from Table S1 from https://doi.org/10.1093/hmg/ddz161
nr.cases <- 9116
nr.controls <- 13292
N_CO <- nr.cases + nr.controls
# "Obesity is having a dramatic impact on modern societies, 
# leading to substantial health issues, with an overall 
# prevalence among children already >20% in many populations, including the USA"
CO_P1 <- 0.2
# CO_K1 <- nr.cases/(nr.cases+nr.controls) #0.406819
CO_K1 <- 0.5 #since it was a meta-analysis we are using the sum of effective sample sizes, and therefore put 0.5 here.

# All have SNP-specific sum of effective sample sizes as a column in the summary stats 
NN <- rep(NA, N_sumstat_files)
NN[c(14)] <- 144414

#enter sample prevalence of .5 to reflect that all traits were munged using the sum of effective sample size
sample.prev <- rep(NA, N_sumstat_files); sample.prev[11] <- CO_K1

#vector of population prevalences
population.prev <- rep(NA, N_sumstat_files); population.prev[11] <- CO_P1


ld <- "eur_w_ld_chr/"
wld <- "eur_w_ld_chr/"

trait.names <- c("BMI_females", "BMI_males", "HIPadjBMI_females", "HIPadjBMI_males", 
	"WCadjBMI_females", "WCadjBMI_males", "WHRadjBMI_females", "WHRadjBMI_males", 
	"BL", "BW", "CO", "IHC", "Height", "bmiTrajBeta", 
	"TFR_females", "TFR_males", "AFR_females", "AFR_males")

# Full ldsc
# sink("LDSCoutput.txt")
# LDSCoutput <- ldsc(traits = c("BMI_females.sumstats.gz", "BMI_males.sumstats.gz", "HIPadjBMI_females.sumstats.gz", "HIPadjBMI_males.sumstats.gz", 
# 	"WCadjBMI_females.sumstats.gz", "WCadjBMI_males.sumstats.gz", "WHRadjBMI_females.sumstats.gz", "WHRadjBMI_males.sumstats.gz", 
# 	"BL.sumstats.gz", "BW.sumstats.gz", "CO.sumstats.gz", "IHC.sumstats.gz", "Height.sumstats.gz", "bmiTrajBeta.sumstats.gz", 
# 	"TFR_females.sumstats.gz", "TFR_males.sumstats.gz", "AFR_females.sumstats.gz", "AFR_males.sumstats.gz"), 
# 	sample.prev = sample.prev, population.prev = population.prev, 
# 	ld = ld, wld = wld, trait.names = trait.names, stand = T)
# sink()
# save(LDSCoutput, file="LDSCoutput.Rdata")

# grep "h2 Z" run_munge_script.out
# 1.  h2 Z: 32.8 BMI_females (body mass index females) - continuous
# 2.  h2 Z: 31.6 BMI_males (body mass index males) - continuous
# 3.  h2 Z: 12.4 HIPadjBMI_females (hip adjusted for BMI females) - continuous
# 4.  h2 Z: 10.7 HIPadjBMI_males (hip adjusted for BMI males) - continuous
# 5.  h2 Z: 11.2 WCadjBMI_females (waist circumference adjusted for BMI females) - continuous
# 6.  h2 Z: 12.5 WCadjBMI_males (waist circumference adjusted for BMI males) - continuous
# 7.  h2 Z: 16.4 WHRadjBMI_females (waist-to-hip ratio adjusted for BMI females) - continuous
# 8.  h2 Z: 20.8 WHRadjBMI_males (waist-to-hip ratio adjusted for BMI males) - continuous
# 9.  h2 Z: 7.09 BL (birth length) - continuous
# 10. h2 Z: 18.3 BW (birth weight) - continuous
# 11. h2 Z: 10.2 CO (childhood obesity) - binary
# 12. h2 Z: 6.94 IHC (infant head circumference) - continuous
# 13. h2 Z: 23.5 Height - continuous
# 14. h2 Z: 22.7 bmiTrajBeta (longitudinal BMI Beta, shifts the mean of the biomarker trajectory) - continuous
# 15. h2 Z: 20.9 TFR_females (bio-electrical impedance trunk fat ratio females) - continuous
# 16. h2 Z: 18.7 TFR_males (bio-electrical impedance trunk fat ratio males) - continuous
# 17. h2 Z: 27.5 AFR_females (bio-electrical impedance arm fat ratio females) - continuous
# 18. h2 Z: 18.6 AFR_males (bio-electrical impedance arm fat ratio males) - continuous

######## Removed bmiTrajTau: h2 Z: 3.91

# Odd ldsc
anthro_ODD <- ldsc(traits = c("BMI_females.sumstats.gz", "BMI_males.sumstats.gz", "HIPadjBMI_females.sumstats.gz", "HIPadjBMI_males.sumstats.gz", 
	"WCadjBMI_females.sumstats.gz", "WCadjBMI_males.sumstats.gz", "WHRadjBMI_females.sumstats.gz", "WHRadjBMI_males.sumstats.gz", 
	"BL.sumstats.gz", "BW.sumstats.gz", "CO.sumstats.gz", "IHC.sumstats.gz", "Height.sumstats.gz", "bmiTrajBeta.sumstats.gz", 
	"TFR_females.sumstats.gz", "TFR_males.sumstats.gz", "AFR_females.sumstats.gz", "AFR_males.sumstats.gz"), 
	sample.prev = sample.prev, population.prev = population.prev, 
	ld = ld, wld = wld, trait.names = trait.names, stand = T, select="ODD")
save(anthro_ODD, file="anthro_ODD.Rdata")

# Even ldsc
anthro_EVEN <- ldsc(traits = c("BMI_females.sumstats.gz", "BMI_males.sumstats.gz", "HIPadjBMI_females.sumstats.gz", "HIPadjBMI_males.sumstats.gz", 
	"WCadjBMI_females.sumstats.gz", "WCadjBMI_males.sumstats.gz", "WHRadjBMI_females.sumstats.gz", "WHRadjBMI_males.sumstats.gz", 
	"BL.sumstats.gz", "BW.sumstats.gz", "CO.sumstats.gz", "IHC.sumstats.gz", "Height.sumstats.gz", "bmiTrajBeta.sumstats.gz", 
	"TFR_females.sumstats.gz", "TFR_males.sumstats.gz", "AFR_females.sumstats.gz", "AFR_males.sumstats.gz"), 
	sample.prev = sample.prev, population.prev = population.prev, 
	ld = ld, wld = wld, trait.names = trait.names, stand = T, select="EVEN")
save(anthro_EVEN, file="anthro_EVEN.Rdata")

load("LDSCoutput.Rdata")
load("anthro_ODD.Rdata")
load("anthro_EVEN.Rdata")
anthro_Full <- LDSCoutput

#the n_factors function I'm providing to run the kaiser, acceleration, and optimal coordinates tests
n_factors<-function(S_LD){
eig <- eigen(cov2cor(S_LD), only.values = TRUE)
eig <- eig$values; nk <- length(eig); k <- 1:nk
criteria <- mean(eig)

# Compute the Kaiser rule.
nkaiser <- sum(eig >= rep(criteria, nk))
# Compute the acceleration factor.
aparallel <- rep(criteria, length(eig))
pred.eig <- af <- rep(NA, nk)
for (j in 2:(nk - 1)) {
  if (eig[j - 1] >= aparallel[j - 1]) {
    af[j] <- (eig[j + 1] - 2 * eig[j]) + eig[j - 1]
  }
}
naf <- which(af == max(af, na.rm = TRUE), TRUE) - 1
# Compute the optimal coordinates.
proportion <- eig/sum(eig)
cumulative <- proportion
for (i in 2:nk) cumulative[i] = cumulative[i - 1] + proportion[i]
proportion[proportion < 0] <- 0
cond1 <- TRUE
cond2 <- TRUE
i <- 0
pred.eig <- af <- rep(NA, nk)
while ((cond1 == TRUE) && (cond2 == TRUE) && (i < nk)) {
  i <- i + 1
  ind <- k[c(i + 1, nk)]
  vp.p <- lm(eig[c(i + 1, nk)] ~ ind)
  vp.prec <- pred.eig[i] <- sum(c(1, i) * coef(vp.p))
  cond1 <- (eig[i] >= vp.prec)
  cond2 <- (eig[i] >= aparallel[i])
  noc <- i - 1
}

cat("Number of factors to retain according to the Kaiser rule:", nkaiser, "\n",
    "Number of factors to retain according to the acceleration factor:", naf, "\n",
    "Number of factors to retain according to the optimal coordinates:", noc, "\n")
}

#run the function in the odd chromosome matrix
n_factors(anthro_ODD$S)
#  Number of factors to retain according to the Kaiser rule: 4 
#  Number of factors to retain according to the acceleration factor: 4 
#  Number of factors to retain according to the optimal coordinates: 4 

#EFA in odd chromosomes
# EFA_2_pro <- factanal(factors = 4, covmat = anthro_ODD$S, rotation = "promax")

#print the eigen values
# eigen(anthro_ODD$S)$values
#we see that there is one small negative eigenvalue

#smooth the matrix
S_LD <- as.matrix((nearPD(anthro_ODD$S, corr = FALSE))$mat)
eigen(S_LD)$values

#try running EFA again
nfactors <- 4
# EFA_2_pro <- factanal(factors = 5, covmat = S_LD, rotation = "varimax")
# EFA_2_pro <- factanal(factors = nfactors, covmat = S_LD, rotation = "promax", lower = 0.009) 
EFA_2_pro <- factanal(factors = nfactors, covmat = S_LD, rotation = "promax", lower = 0.01)  

#save the EFa loadings for the N_sumstat_files variables and the nfactors factors
Loadings_2_pro <- (EFA_2_pro$loadings[1:N_sumstat_files,1:nfactors])
write.table(as.data.frame(Loadings_2_pro), "EFA_loadings.csv", col.names=T, row.names=T, quote=F, sep=",")

# ?usermodel
#use the write.model function to write out the model
Model_2_pro <- write.model(Loadings_2_pro, S_LD, 0.5, mustload=TRUE)
# Model_2_pro <- write.model(Loadings_2_pro, S_LD, 0.5, mustload=FALSE)
save(Model_2_pro, file="Model_2_pro.Rdata")
# load("Model_2_pro.Rdata")

Model_2_pro <- "
F1=~c_BL*BL + c_BW*BW + c_IHC*IHC
c_BL < 0.9999
F2=~c_WHRadjBMI_females*WHRadjBMI_females + c_WHRadjBMI_males*WHRadjBMI_males + c_WCadjBMI_females*WCadjBMI_females
c_WCadjBMI_females < 0.9999
F3=~c_WCadjBMI_males*WCadjBMI_males + c_Height*Height + c_HIPadjBMI_females*HIPadjBMI_females + c_HIPadjBMI_males*HIPadjBMI_males + c_TFR_females*TFR_females + c1_TFR_males*TFR_males
F4=~c_AFR_females*AFR_females + c_AFR_males*AFR_males + c_bmiTrajBeta*bmiTrajBeta + c_BMI_females*BMI_females + c_BMI_males*BMI_males + c_CO*CO + c2_TFR_males*TFR_males
WCadjBMI_females~~WCadjBMI_males
F1~~F2 \n F1~~F3 \n F1~~F4 \n F2~~F3 \n F2~~F4 \n F3~~F4
BMI_females ~~ ceio*BMI_females \n ceio > .0001 \n BMI_males ~~ acmr*BMI_males \n acmr > .0001 \n HIPadjBMI_females ~~ aekl*HIPadjBMI_females \n aekl > .0001 \n HIPadjBMI_males ~~ hiky*HIPadjBMI_males \n hiky > .0001 \n WCadjBMI_females ~~ adkq*WCadjBMI_females \n adkq > .0001 \n WCadjBMI_males ~~ ctuz*WCadjBMI_males \n ctuz > .0001 \n WHRadjBMI_females ~~ bgqz*WHRadjBMI_females \n bgqz > .0001 \n WHRadjBMI_males ~~ crsy*WHRadjBMI_males \n crsy > .0001 \n BL ~~ mnrt*BL \n mnrt > .0001 \n BW ~~ cefx*BW \n cefx > .0001 \n CO ~~ aciz*CO \n aciz > .0001 \n IHC ~~ bhsu*IHC \n bhsu > .0001 \n Height ~~ cdvw*Height \n cdvw > .0001 \n bmiTrajBeta ~~ fijx*bmiTrajBeta \n fijx > .0001 \n TFR_females ~~ cirw*TFR_females \n cirw > .0001 \n AFR_females ~~ fjps*AFR_females \n fjps > .0001 \n AFR_males ~~ cmuz*AFR_males \n cmuz > .0001
"

# BMI_females~~BMI_males 0.93, F5
# HIPadjBMI_females~~HIPadjBMI_males 0.89, F2
# WCadjBMI_females~~WCadjBMI_males 0.77, F3 F2
# WHRadjBMI_females~~WHRadjBMI_males 0.66, F3
# TFR_females~~TFR_males 0.43, F2 F4
# AFR_females~~AFR_males 0.42, F5
# LFR_females~~LFR_males 0.35, F2 F4

CFA2_EVEN <- usermodel(anthro_EVEN, model=Model_2_pro, std.lv=TRUE, imp_cov=TRUE)
CFA2_EVEN$modelfit
#       chisq  df p_chisq      AIC      CFI      SRMR
# df 21460.36 127       0 21548.36 0.892376 0.1382072

commonfactor_EVEN <- commonfactor(anthro_EVEN)
commonfactor_EVEN$modelfit
#       chisq  df p_chisq      AIC       CFI      SRMR
# df 38157.88 135       0 38229.88 0.8081796 0.2606608
CFA2_EVEN$modelfit

CFA2_Full <- usermodel(anthro_Full, model=Model_2_pro, std.lv=TRUE, imp_cov=TRUE)
CFA2_Full$modelfit
#       chisq  df p_chisq      AIC       CFI      SRMR
# df 50482.39 127       0 50570.39 0.8878071 0.1300755
sink("CFA2_Full.txt")
CFA2_Full$results
sink()
write.table(as.data.frame(CFA2_Full$results), "CFA2_Full.csv", col.names=T, row.names=F, quote=F, sep=",")
EOF

cat << 'EOF' >GenomicSEM_commonfactorGWAS_CHUNKS.sh
#!/bin/bash
#SBATCH --qos=normal
#SBATCH --partition=amilan
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=3700
#SBATCH --time=04:30:00
#SBATCH -J c_GSEM
#SBATCH --array=4-1000
#SBATCH -o ./GenomicSEM_commonfactorGWAS_CHUNKS_%A-%a.out


date
startTime="$(date +%s)"

arrayCount=1000

module load anaconda
conda activate nmr_prs
date

Rscript GenomicSEM_commonfactorGWAS_CHUNKS.R "${SLURM_ARRAY_TASK_ID}" "${SLURM_NTASKS}" "${arrayCount}"

scontrol show job $SLURM_JOB_ID

date
endTime="$(date +%s)"
diff=$(echo "$endTime-$startTime" |bc)
diff=$(($diff + 0))

if [ $diff -lt "$((60))" ]
then
  printf %.3f "$((1000 * $diff/1))e-3"; echo " seconds"
elif [ $diff -lt "$((60 * 60))" ]
then
  printf %.3f "$((1000 * $diff/60))e-3"; echo " minutes"
elif [ $diff -lt "$((60 * 60 * 24))" ]
then
  printf %.3f "$((1000 * $diff/3600))e-3"; echo " hours"
else
  printf %.3f "$((1000 * $diff/86400))e-3"; echo " days"
fi

# * # * # SCRIPT END # * # * #
EOF

cat << 'EOF' >GenomicSEM_commonfactorGWAS_CHUNKS.R
library(GenomicSEM)
library(data.table)
library(R.utils)
library(gdata)
`%ni%` <- Negate(`%in%`)
setwd("working/directory")

## Run the multivariate GWAS
print("loading LDSCoutput.Rdata")
load("LDSCoutput.Rdata")
print("loading OVERWEIGHT_sumstats.Rdata")
load("OVERWEIGHT_sumstats.Rdata")

args <- commandArgs(trailingOnly = TRUE)
job_id <- as.numeric(args[1])
ncores <- as.numeric(args[2])
arrayCount <- as.numeric(args[3])

# job_id = 1
# ncores = 24
# arrayCount = 4

############ standardize the ldsc matrix ############

# LDSCoutputStand <- LDSCoutput
# LDSCoutputStand$V <- LDSCoutputStand$V_Stand
# LDSCoutputStand$S <- LDSCoutputStand$S_Stand



############ set up snp chunk ############

print(paste0("Using ", ncores, " cores"))

print(paste0("there are ", nrow(OVERWEIGHT_sumstats), " rows in sumstats object"))
# nrow(OVERWEIGHT_sumstats) ##954,086

snpj <- 1:nrow(OVERWEIGHT_sumstats)
chunks <- split(snpj, ceiling(seq_along(snpj)/1000))
chunk_idx <- job_id
chunk <- chunks[[chunk_idx]]

# useful for troubleshooting with a smaller dataset
# chunk <- sample(snpj, size = 24*4, replace = F)
# chunk <- 1:100

# SPECIFY THE MODEL and prevent residual variance from being estimated below 0 for waist circumference
modelSNP <- paste0("
F1=~BL + BW + IHC
F2=~WHRadjBMI_females + WHRadjBMI_males + WCadjBMI_females
F3=~WCadjBMI_females + WCadjBMI_males + Height + HIPadjBMI_females + HIPadjBMI_males + TFR_females + TFR_males
F4=~AFR_females + AFR_males + bmiTrajBeta + BMI_females + BMI_males + CO + TFR_males
F1~~F2 \n F1~~F3 \n F1~~F4 \n F2~~F3 \n F2~~F4 \n F3~~F4
BMI_females ~~ ceio*BMI_females \n ceio > .0001 \n BMI_males ~~ acmr*BMI_males \n acmr > .0001 \n HIPadjBMI_females ~~ aekl*HIPadjBMI_females \n aekl > .0001 \n HIPadjBMI_males ~~ hiky*HIPadjBMI_males \n hiky > .0001 \n WCadjBMI_females ~~ adkq*WCadjBMI_females \n adkq > .0001 \n WCadjBMI_males ~~ ctuz*WCadjBMI_males \n ctuz > .0001 \n WHRadjBMI_females ~~ bgqz*WHRadjBMI_females \n bgqz > .0001 \n WHRadjBMI_males ~~ crsy*WHRadjBMI_males \n crsy > .0001 \n BL ~~ mnrt*BL \n mnrt > .0001 \n BW ~~ cefx*BW \n cefx > .0001 \n CO ~~ aciz*CO \n aciz > .0001 \n IHC ~~ bhsu*IHC \n bhsu > .0001 \n Height ~~ cdvw*Height \n cdvw > .0001 \n bmiTrajBeta ~~ fijx*bmiTrajBeta \n fijx > .0001 \n TFR_females ~~ cirw*TFR_females \n cirw > .0001 \n AFR_females ~~ fjps*AFR_females \n fjps > .0001 \n AFR_males ~~ cmuz*AFR_males \n cmuz > .0001
F1~SNP
F2~SNP
F3~SNP
F4~SNP
")
cat(modelSNP)

print("---------- FactorModel via userGWAS(...) ----------")
sub <- c("F1 ~ SNP", "F2 ~ SNP", "F3 ~ SNP", "F4 ~ SNP", "F5 ~ SNP")
FactorModel <- userGWAS(covstruc = LDSCoutput, SNPs = OVERWEIGHT_sumstats[chunk, ], 
                        estimation = "DWLS", model = modelSNP, cores = ncores, parallel = TRUE, 
                        std.lv = TRUE, printwarn = TRUE, 
                        sub=sub, 
                        GC="standard", fix_measurement=T,
                        MPI=FALSE, SNPSE = FALSE, smooth_check=T)
FactorModel_df <- FactorModel
F1df <- as.data.frame(FactorModel_df[[1]])
F2df <- as.data.frame(FactorModel_df[[2]])
F3df <- as.data.frame(FactorModel_df[[3]])
F4df <- as.data.frame(FactorModel_df[[4]])
F5df <- as.data.frame(FactorModel_df[[5]])

pdf(paste0("./Qsnp_Plots/Qsnp_Plots", chunk_idx, ".pdf"), width = 7, height = 5)

# ## F1 Qsnp --#####
print("---------- F1 Qsnp ----------")
model_F1Q <- paste0("
F1=~BL + BW + IHC
F2=~WHRadjBMI_females + WHRadjBMI_males + WCadjBMI_females
F3=~WCadjBMI_females + WCadjBMI_males + Height + HIPadjBMI_females + HIPadjBMI_males + TFR_females + TFR_males
F4=~AFR_females + AFR_males + bmiTrajBeta + BMI_females + BMI_males + CO + TFR_males
F1~~F2 \n F1~~F3 \n F1~~F4 \n F2~~F3 \n F2~~F4 \n F3~~F4
BMI_females ~~ ceio*BMI_females \n ceio > .0001 \n BMI_males ~~ acmr*BMI_males \n acmr > .0001 \n HIPadjBMI_females ~~ aekl*HIPadjBMI_females \n aekl > .0001 \n HIPadjBMI_males ~~ hiky*HIPadjBMI_males \n hiky > .0001 \n WCadjBMI_females ~~ adkq*WCadjBMI_females \n adkq > .0001 \n WCadjBMI_males ~~ ctuz*WCadjBMI_males \n ctuz > .0001 \n WHRadjBMI_females ~~ bgqz*WHRadjBMI_females \n bgqz > .0001 \n WHRadjBMI_males ~~ crsy*WHRadjBMI_males \n crsy > .0001 \n BL ~~ mnrt*BL \n mnrt > .0001 \n BW ~~ cefx*BW \n cefx > .0001 \n CO ~~ aciz*CO \n aciz > .0001 \n IHC ~~ bhsu*IHC \n bhsu > .0001 \n Height ~~ cdvw*Height \n cdvw > .0001 \n bmiTrajBeta ~~ fijx*bmiTrajBeta \n fijx > .0001 \n TFR_females ~~ cirw*TFR_females \n cirw > .0001 \n AFR_females ~~ fjps*AFR_females \n fjps > .0001 \n AFR_males ~~ cmuz*AFR_males \n cmuz > .0001
BL + BW + IHC ~ SNP
F2~SNP
F3~SNP
F4~SNP
")

cat(model_F1Q)

subQsnp <- c("F1~~F2", "F1~~F3", "F1~~F4", "F1~~F5", "F2~~F3", "F2~~F4", "F2~~F5", "F3~~F4", "F3~~F5", "F4~~F5")
FactorModel_F1Q <- userGWAS(covstruc = LDSCoutput, SNPs = OVERWEIGHT_sumstats[chunk, ],
                            estimation = "DWLS", model = model_F1Q, cores = ncores, parallel = TRUE,
                            sub=subQsnp, fix_measurement = T,
                            std.lv = TRUE, printwarn = TRUE, toler = FALSE, SNPSE = FALSE)
# dfc <- as.data.frame((FactorModel_F1Q))
# plot_SEM_loadings(whichValueToPlot = "est", df = dfb)
# plot_SEM_loadings(whichValueToPlot = "est", df = dfc)

Q_chisq.F1 <- FactorModel[[1]]$chisq - FactorModel_F1Q[[1]]$chisq
hist(FactorModel[[1]]$chisq, breaks = 30, main = "F1 Step 1 Model (Red Model) chisq", xlim = c(min(c(FactorModel[[1]]$chisq,FactorModel_F1Q[[1]]$chisq)),max(c(FactorModel[[1]]$chisq,FactorModel_F1Q[[1]]$chisq))))
hist(FactorModel_F1Q[[1]]$chisq, breaks = 30, main = "F1 Step 2 Indep. Pathways Model (Black Model) chisq", xlim = c(min(c(FactorModel[[1]]$chisq,FactorModel_F1Q[[1]]$chisq)),max(c(FactorModel[[1]]$chisq,FactorModel_F1Q[[1]]$chisq))))
df.F1 <- FactorModel[[1]]$chisq_df[1] - FactorModel_F1Q[[1]]$chisq_df[1]
hist(Q_chisq.F1, breaks = 30, main = paste0("F1 df = ",df.F1,"\nStep 1 Model (Red) chisq - Step 2 Indep. Pathways Model (Black) chisq"))
Q_chisq_pval.F1 <- pchisq(Q_chisq.F1, df.F1, lower.tail=FALSE)
hist(Q_chisq_pval.F1, breaks = 30, main = paste0("F1 df = ",df.F1,", p-value for difference in chisq"))
QsnpIndices.F1 <- which(Q_chisq_pval.F1 < 5E-08)
nonQsnpIndices.F1 <- which(Q_chisq_pval.F1 >= 5E-08)
print(paste0("F1 N Genome Wide Sig Qsnps: ", length(QsnpIndices.F1)))
F1df$CommonPath_chisq <- FactorModel[[1]]$chisq
F1df$IndepPath_chisq <- FactorModel_F1Q[[1]]$chisq
F1df$Q_chisq_pval <- Q_chisq_pval.F1

# ## F2 Qsnp --#####
print("---------- F2 Qsnp ----------")
model_F2Q <- paste0("
F1=~BL + BW + IHC
F2=~WHRadjBMI_females + WHRadjBMI_males + WCadjBMI_females
F3=~WCadjBMI_females + WCadjBMI_males + Height + HIPadjBMI_females + HIPadjBMI_males + TFR_females + TFR_males
F4=~AFR_females + AFR_males + bmiTrajBeta + BMI_females + BMI_males + CO + TFR_males
F1~~F2 \n F1~~F3 \n F1~~F4 \n F2~~F3 \n F2~~F4 \n F3~~F4
BMI_females ~~ ceio*BMI_females \n ceio > .0001 \n BMI_males ~~ acmr*BMI_males \n acmr > .0001 \n HIPadjBMI_females ~~ aekl*HIPadjBMI_females \n aekl > .0001 \n HIPadjBMI_males ~~ hiky*HIPadjBMI_males \n hiky > .0001 \n WCadjBMI_females ~~ adkq*WCadjBMI_females \n adkq > .0001 \n WCadjBMI_males ~~ ctuz*WCadjBMI_males \n ctuz > .0001 \n WHRadjBMI_females ~~ bgqz*WHRadjBMI_females \n bgqz > .0001 \n WHRadjBMI_males ~~ crsy*WHRadjBMI_males \n crsy > .0001 \n BL ~~ mnrt*BL \n mnrt > .0001 \n BW ~~ cefx*BW \n cefx > .0001 \n CO ~~ aciz*CO \n aciz > .0001 \n IHC ~~ bhsu*IHC \n bhsu > .0001 \n Height ~~ cdvw*Height \n cdvw > .0001 \n bmiTrajBeta ~~ fijx*bmiTrajBeta \n fijx > .0001 \n TFR_females ~~ cirw*TFR_females \n cirw > .0001 \n AFR_females ~~ fjps*AFR_females \n fjps > .0001 \n AFR_males ~~ cmuz*AFR_males \n cmuz > .0001
F1~SNP
WHRadjBMI_females + WHRadjBMI_males + WCadjBMI_females ~SNP
F3~SNP
F4~SNP
")

cat(model_F2Q)

subQsnp <- c("F1~~F2", "F1~~F3", "F1~~F4", "F1~~F5", "F2~~F3", "F2~~F4", "F2~~F5", "F3~~F4", "F3~~F5", "F4~~F5")
FactorModel_F2Q <- userGWAS(covstruc = LDSCoutput, SNPs = OVERWEIGHT_sumstats[chunk, ],
                            estimation = "DWLS", model = model_F2Q, cores = ncores, parallel = TRUE,
                            sub=subQsnp, fix_measurement = T,
                            std.lv = TRUE, printwarn = TRUE, toler = FALSE, SNPSE = FALSE)
# dfc <- as.data.frame((FactorModel_F2Q))
# plot_SEM_loadings(whichValueToPlot = "est", df = dfb)
# plot_SEM_loadings(whichValueToPlot = "est", df = dfc)

Q_chisq.F2 <- FactorModel[[1]]$chisq - FactorModel_F2Q[[1]]$chisq
hist(FactorModel[[1]]$chisq, breaks = 30, main = "F2 Step 1 Model (Red Model) chisq", xlim = c(min(c(FactorModel[[1]]$chisq,FactorModel_F2Q[[1]]$chisq)),max(c(FactorModel[[1]]$chisq,FactorModel_F2Q[[1]]$chisq))))
hist(FactorModel_F2Q[[1]]$chisq, breaks = 30, main = "F2 Step 2 Indep. Pathways Model (Black Model) chisq", xlim = c(min(c(FactorModel[[1]]$chisq,FactorModel_F2Q[[1]]$chisq)),max(c(FactorModel[[1]]$chisq,FactorModel_F2Q[[1]]$chisq))))
df.F2 <- FactorModel[[1]]$chisq_df[1] - FactorModel_F2Q[[1]]$chisq_df[1]
hist(Q_chisq.F2, breaks = 30, main = paste0("F2 df = ",df.F2,"\nStep 1 Model (Red) chisq - Step 2 Indep. Pathways Model (Black) chisq"))
Q_chisq_pval.F2 <- pchisq(Q_chisq.F2, df.F2, lower.tail=FALSE)
hist(Q_chisq_pval.F2, breaks = 30, main = paste0("F2 df = ",df.F2,", p-value for difference in chisq"))
QsnpIndices.F2 <- which(Q_chisq_pval.F2 < 5E-08)
nonQsnpIndices.F2 <- which(Q_chisq_pval.F2 >= 5E-08)
print(paste0("F2 N Genome Wide Sig Qsnps: ", length(QsnpIndices.F2)))
F2df$CommonPath_chisq <- FactorModel[[1]]$chisq
F2df$IndepPath_chisq <- FactorModel_F2Q[[1]]$chisq
F2df$Q_chisq_pval <- Q_chisq_pval.F2

# ## F3 Qsnp --#####
print("---------- F3 Qsnp ----------")
model_F3Q <- paste0("
F1=~BL + BW + IHC
F2=~WHRadjBMI_females + WHRadjBMI_males + WCadjBMI_females
F3=~WCadjBMI_females + WCadjBMI_males + Height + HIPadjBMI_females + HIPadjBMI_males + TFR_females + TFR_males
F4=~AFR_females + AFR_males + bmiTrajBeta + BMI_females + BMI_males + CO + TFR_males
F1~~F2 \n F1~~F3 \n F1~~F4 \n F2~~F3 \n F2~~F4 \n F3~~F4
BMI_females ~~ ceio*BMI_females \n ceio > .0001 \n BMI_males ~~ acmr*BMI_males \n acmr > .0001 \n HIPadjBMI_females ~~ aekl*HIPadjBMI_females \n aekl > .0001 \n HIPadjBMI_males ~~ hiky*HIPadjBMI_males \n hiky > .0001 \n WCadjBMI_females ~~ adkq*WCadjBMI_females \n adkq > .0001 \n WCadjBMI_males ~~ ctuz*WCadjBMI_males \n ctuz > .0001 \n WHRadjBMI_females ~~ bgqz*WHRadjBMI_females \n bgqz > .0001 \n WHRadjBMI_males ~~ crsy*WHRadjBMI_males \n crsy > .0001 \n BL ~~ mnrt*BL \n mnrt > .0001 \n BW ~~ cefx*BW \n cefx > .0001 \n CO ~~ aciz*CO \n aciz > .0001 \n IHC ~~ bhsu*IHC \n bhsu > .0001 \n Height ~~ cdvw*Height \n cdvw > .0001 \n bmiTrajBeta ~~ fijx*bmiTrajBeta \n fijx > .0001 \n TFR_females ~~ cirw*TFR_females \n cirw > .0001 \n AFR_females ~~ fjps*AFR_females \n fjps > .0001 \n AFR_males ~~ cmuz*AFR_males \n cmuz > .0001
F1~SNP
F2~SNP
WCadjBMI_males + Height + HIPadjBMI_females + HIPadjBMI_males + TFR_females + TFR_males~SNP
F4~SNP
")

cat(model_F3Q)

subQsnp <- c("F1~~F2", "F1~~F3", "F1~~F4", "F1~~F5", "F2~~F3", "F2~~F4", "F2~~F5", "F3~~F4", "F3~~F5", "F4~~F5")
FactorModel_F3Q <- userGWAS(covstruc = LDSCoutput, SNPs = OVERWEIGHT_sumstats[chunk, ],
                            estimation = "DWLS", model = model_F3Q, cores = ncores, parallel = TRUE,
                            sub=subQsnp, fix_measurement = T,
                            std.lv = TRUE, printwarn = TRUE, toler = FALSE, SNPSE = FALSE)
# dfc <- as.data.frame((FactorModel_F3Q))
# plot_SEM_loadings(whichValueToPlot = "est", df = dfb)
# plot_SEM_loadings(whichValueToPlot = "est", df = dfc)

Q_chisq.F3 <- FactorModel[[1]]$chisq - FactorModel_F3Q[[1]]$chisq
hist(FactorModel[[1]]$chisq, breaks = 30, main = "F3 Step 1 Model (Red Model) chisq", xlim = c(min(c(FactorModel[[1]]$chisq,FactorModel_F3Q[[1]]$chisq)),max(c(FactorModel[[1]]$chisq,FactorModel_F3Q[[1]]$chisq))))
hist(FactorModel_F3Q[[1]]$chisq, breaks = 30, main = "F3 Step 2 Indep. Pathways Model (Black Model) chisq", xlim = c(min(c(FactorModel[[1]]$chisq,FactorModel_F3Q[[1]]$chisq)),max(c(FactorModel[[1]]$chisq,FactorModel_F3Q[[1]]$chisq))))
df.F3 <- FactorModel[[1]]$chisq_df[1] - FactorModel_F3Q[[1]]$chisq_df[1]
hist(Q_chisq.F3, breaks = 30, main = paste0("F3 df = ",df.F3,"\nStep 1 Model (Red) chisq - Step 2 Indep. Pathways Model (Black) chisq"))
Q_chisq_pval.F3 <- pchisq(Q_chisq.F3, df.F3, lower.tail=FALSE)
hist(Q_chisq_pval.F3, breaks = 30, main = paste0("F3 df = ",df.F3,", p-value for Difference in chisq"))
QsnpIndices.F3 <- which(Q_chisq_pval.F3 < 5E-08)
nonQsnpIndices.F3 <- which(Q_chisq_pval.F3 >= 5E-08)
print(paste0("F3 N Genome Wide Sig Qsnps: ", length(QsnpIndices.F3)))
F3df$CommonPath_chisq <- FactorModel[[1]]$chisq
F3df$IndepPath_chisq <- FactorModel_F3Q[[1]]$chisq
F3df$Q_chisq_pval <- Q_chisq_pval.F3

# ## F4 Qsnp --#####
print("---------- F4 Qsnp ----------")
model_F4Q <- paste0("
F1=~BL + BW + IHC
F2=~WHRadjBMI_females + WHRadjBMI_males + WCadjBMI_females
F3=~WCadjBMI_females + WCadjBMI_males + Height + HIPadjBMI_females + HIPadjBMI_males + TFR_females + TFR_males
F4=~AFR_females + AFR_males + bmiTrajBeta + BMI_females + BMI_males + CO + TFR_males
F1~~F2 \n F1~~F3 \n F1~~F4 \n F2~~F3 \n F2~~F4 \n F3~~F4
BMI_females ~~ ceio*BMI_females \n ceio > .0001 \n BMI_males ~~ acmr*BMI_males \n acmr > .0001 \n HIPadjBMI_females ~~ aekl*HIPadjBMI_females \n aekl > .0001 \n HIPadjBMI_males ~~ hiky*HIPadjBMI_males \n hiky > .0001 \n WCadjBMI_females ~~ adkq*WCadjBMI_females \n adkq > .0001 \n WCadjBMI_males ~~ ctuz*WCadjBMI_males \n ctuz > .0001 \n WHRadjBMI_females ~~ bgqz*WHRadjBMI_females \n bgqz > .0001 \n WHRadjBMI_males ~~ crsy*WHRadjBMI_males \n crsy > .0001 \n BL ~~ mnrt*BL \n mnrt > .0001 \n BW ~~ cefx*BW \n cefx > .0001 \n CO ~~ aciz*CO \n aciz > .0001 \n IHC ~~ bhsu*IHC \n bhsu > .0001 \n Height ~~ cdvw*Height \n cdvw > .0001 \n bmiTrajBeta ~~ fijx*bmiTrajBeta \n fijx > .0001 \n TFR_females ~~ cirw*TFR_females \n cirw > .0001 \n AFR_females ~~ fjps*AFR_females \n fjps > .0001 \n AFR_males ~~ cmuz*AFR_males \n cmuz > .0001
F1~SNP
F2~SNP
F3~SNP
AFR_females + AFR_males + bmiTrajBeta + BMI_females + BMI_males + CO + TFR_males~SNP
")

cat(model_F4Q)

subQsnp <- c("F1~~F2", "F1~~F3", "F1~~F4", "F1~~F5", "F2~~F3", "F2~~F4", "F2~~F5", "F3~~F4", "F3~~F5", "F4~~F5")
FactorModel_F4Q <- userGWAS(covstruc = LDSCoutput, SNPs = OVERWEIGHT_sumstats[chunk, ],
                            estimation = "DWLS", model = model_F4Q, cores = ncores, parallel = TRUE,
                            sub=subQsnp, fix_measurement = T,
                            std.lv = TRUE, printwarn = TRUE, toler = FALSE, SNPSE = FALSE)
# dfc <- as.data.frame((FactorModel_F4Q))
# plot_SEM_loadings(whichValueToPlot = "est", df = dfb)
# plot_SEM_loadings(whichValueToPlot = "est", df = dfc)

Q_chisq.F4 <- FactorModel[[1]]$chisq - FactorModel_F4Q[[1]]$chisq
hist(FactorModel[[1]]$chisq, breaks = 30, main = "F4 Step 1 Model (Red Model) chisq", xlim = c(min(c(FactorModel[[1]]$chisq,FactorModel_F4Q[[1]]$chisq)),max(c(FactorModel[[1]]$chisq,FactorModel_F4Q[[1]]$chisq))))
hist(FactorModel_F4Q[[1]]$chisq, breaks = 30, main = "F4 Step 2 Indep. Pathways Model (Black Model) chisq", xlim = c(min(c(FactorModel[[1]]$chisq,FactorModel_F4Q[[1]]$chisq)),max(c(FactorModel[[1]]$chisq,FactorModel_F4Q[[1]]$chisq))))
df.F4 <- FactorModel[[1]]$chisq_df[1] - FactorModel_F4Q[[1]]$chisq_df[1]
hist(Q_chisq.F4, breaks = 30, main = paste0("F4 df = ",df.F4,"\nStep 1 Model (Red) chisq - Step 2 Indep. Pathways Model (Black) chisq"))
Q_chisq_pval.F4 <- pchisq(Q_chisq.F4, df.F4, lower.tail=FALSE)
hist(Q_chisq_pval.F4, breaks = 30, main = paste0("F4 df = ",df.F4,", p-value for difference in chisq"))
QsnpIndices.F4 <- which(Q_chisq_pval.F4 < 5E-08)
nonQsnpIndices.F4 <- which(Q_chisq_pval.F4 >= 5E-08)
print(paste0("F4 N Genome Wide Sig Qsnps: ", length(QsnpIndices.F4)))
F4df$CommonPath_chisq <- FactorModel[[1]]$chisq
F4df$IndepPath_chisq <- FactorModel_F4Q[[1]]$chisq
F4df$Q_chisq_pval <- Q_chisq_pval.F4

dev.off()

# ## saveRDS --#####
print("---------- saveRDS GWAS Factors ----------")
saveRDS(F1df, file = paste0("./output_4F/F1/OVERWEIGHT_F1_", chunk_idx, ".rds"))
saveRDS(F2df, file = paste0("./output_4F/F2/OVERWEIGHT_F2_", chunk_idx, ".rds"))
saveRDS(F3df, file = paste0("./output_4F/F3/OVERWEIGHT_F3_", chunk_idx, ".rds"))
saveRDS(F4df, file = paste0("./output_4F/F4/OVERWEIGHT_F4_", chunk_idx, ".rds"))

print("-----END-----")
EOF


cat << 'EOF' >combineOutputFiles.R
library(data.table)
library(plyr)
library(gtools)

args <- commandArgs(trailingOnly = TRUE)
whichFactor <- as.numeric(args[1])

# whichFactor <- 4

if(whichFactor == 1){
	Flist <- list.files("./output_4F/F1")
}else if(whichFactor == 2){
	Flist <- list.files("./output_4F/F2")
}else if(whichFactor == 3){
	Flist <- list.files("./output_4F/F3")
}else if(whichFactor == 4){
	Flist <- list.files("./output_4F/F4")
}

Flist <- mixedsort(sort(Flist))
# Flist <- Flist[1:500]

length(Flist)
head(Flist); tail(Flist)
c_Fdf_list <- list()
rowCount <- 0
told <- Sys.time(); tstart <- Sys.time()
tsum = 0
for(i in 1:length(Flist)){
	
	if(i %in% round(seq(from = 1, to = length(Flist), length.out = 20))){
		tnew <- Sys.time()
		print(paste0(round(100*i/length(Flist)),"%"))
		print(tnew - told)
		told <- Sys.time()
	}

	if(whichFactor == 1){
		Fdf <- readRDS(paste0("./output_4F/F1/OVERWEIGHT_F1_",i,".rds"))
	}else if(whichFactor == 2){
		Fdf <- readRDS(paste0("./output_4F/F2/OVERWEIGHT_F2_",i,".rds"))
	}else if(whichFactor == 3){
		Fdf <- readRDS(paste0("./output_4F/F3/OVERWEIGHT_F3_",i,".rds"))
	}else if(whichFactor == 4){
		Fdf <- readRDS(paste0("./output_4F/F4/OVERWEIGHT_F4_",i,".rds"))
	}

	c_Fdf_list[[i]] <- Fdf
	rowCount <- rowCount + nrow(Fdf)
	rm(Fdf)
}
c_Fdf <- rbind.fill(c_Fdf_list)
tsum = Sys.time() - tstart
print(tsum)

# this value should be 1000, since that was the chunk of SNPs in each array job
print(rowCount/length(Flist))

summary(c_Fdf$Z_smooth)
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.000e+00 6.470e-06 1.949e-05 4.704e-05 5.154e-05 9.723e-03 

summary(as.factor(c_Fdf$warning))
summary(as.factor(c_Fdf$error))
c_Fdf_qc1 <- subset(c_Fdf, error == 0)
c_Fdf_qc1$lhs <- NULL
c_Fdf_qc1$op <- NULL
c_Fdf_qc1$rhs <- NULL
c_Fdf_qc1$free <- NULL
c_Fdf_qc1$label <- NULL
c_Fdf_qc1$error <- NULL
c_Fdf_qc1$warning <- NULL

names(c_Fdf_qc1)

summary(c_Fdf_qc1$Q_chisq_pval)

# how many Qsnps are there
c_Fdf_qc1_Qsnp <- subset(c_Fdf_qc1, Q_chisq_pval < 5E-08)
nrow(c_Fdf_qc1_Qsnp)
nrow(c_Fdf_qc1)
paste0(prettyNum(nrow(c_Fdf_qc1_Qsnp), big.mark=",", scientific=FALSE), ", or ", round(100*nrow(c_Fdf_qc1_Qsnp)/nrow(c_Fdf_qc1),5), "% of SNPs are Q SNPs")

# how many significant SNPs are there
c_Fdf_qc1_GWsig <- subset(c_Fdf_qc1, Pval_Estimate < 5E-08)
nrow(c_Fdf_qc1_GWsig)
nrow(c_Fdf_qc1)
paste0(prettyNum(nrow(c_Fdf_qc1_GWsig), big.mark=",", scientific=FALSE), ", or ", round(100*nrow(c_Fdf_qc1_GWsig)/nrow(c_Fdf_qc1),5), "% of SNPs are Genome Wide Significant")

# make indicators for the significant SNPs
c_Fdf_qc1$isQsnp <- 0
c_Fdf_qc1$isQsnp[c_Fdf_qc1$Q_chisq_pval < 5E-08] <- 1
summary(as.factor(c_Fdf_qc1$isQsnp))

c_Fdf_qc1$isGWS <- 0
c_Fdf_qc1$isGWS[c_Fdf_qc1$Pval_Estimate < 5E-08] <- 1
summary(as.factor(c_Fdf_qc1$isGWS))

c_Fdf_qc1$isQsnp_and_isGWS <- 0
c_Fdf_qc1$isQsnp_and_isGWS[(c_Fdf_qc1$Q_chisq_pval < 5E-08 & c_Fdf_qc1$Pval_Estimate < 5E-08)] <- 1
summary(as.factor(c_Fdf_qc1$isQsnp_and_isGWS))
paste0(prettyNum(sum(c_Fdf_qc1$isQsnp_and_isGWS), big.mark=",", scientific=FALSE), ", or ", round(100*sum(c_Fdf_qc1$isQsnp_and_isGWS)/nrow(c_Fdf_qc1),5), "% of SNPs are both Qsnps and Genome Wide Significant")

# Write the files out in tab-delimited format
c_Fdf_qc1_noQsnp <- subset(c_Fdf_qc1, c_Fdf_qc1$isQsnp == 0)
Qsnp_rsids <- c_Fdf_qc1$SNP[c_Fdf_qc1$isQsnp == 1]
if(whichFactor == 1){
	write.table(c_Fdf_qc1, "./output_4F/F1_GSEM_Birth_Size_sumstats.txt", sep = "\t", quote = F, col.names = T, row.names = F)
	write.table(Qsnp_rsids, "./output_4F/F1_GSEM_Birth_Size_Qsnps.txt", sep = "\t", quote = F, col.names = F, row.names = F)
	# write.table(c_Fdf_qc1_noQsnp, "./output_4F/F1_GSEM_Birth_Size_noQsnp_sumstats.txt", sep = "\t", quote = F, col.names = T, row.names = F)
}else if(whichFactor == 2){
	write.table(c_Fdf_qc1, "./output_4F/F2_GSEM_Abdominal_Size_sumstats.txt", sep = "\t", quote = F, col.names = T, row.names = F)
	write.table(Qsnp_rsids, "./output_4F/F2_GSEM_Abdominal_Size_Qsnps.txt", sep = "\t", quote = F, col.names = F, row.names = F)
	# write.table(c_Fdf_qc1_noQsnp, "./output_4F/F2_GSEM_Abdominal_Size_noQsnp_sumstats.txt", sep = "\t", quote = F, col.names = T, row.names = F)
}else if(whichFactor == 3){
	write.table(c_Fdf_qc1, "./output_4F/F3_GSEM_Body_Size_and_Adipose_Distribution_sumstats.txt", sep = "\t", quote = F, col.names = T, row.names = F)
	write.table(Qsnp_rsids, "./output_4F/F3_GSEM_Body_Size_and_Adipose_Distribution_Qsnps.txt", sep = "\t", quote = F, col.names = F, row.names = F)
	# write.table(c_Fdf_qc1_noQsnp, "./output_4F/F3_GSEM_Body_Size_and_Adipose_Distribution_noQsnp_sumstats.txt", sep = "\t", quote = F, col.names = T, row.names = F)
}else if(whichFactor == 4){
	write.table(c_Fdf_qc1, "./output_4F/F4_GSEM_Adiposity_sumstats.txt", sep = "\t", quote = F, col.names = T, row.names = F)
	write.table(Qsnp_rsids, "./output_4F/F4_GSEM_Adiposity_Qsnps.txt", sep = "\t", quote = F, col.names = F, row.names = F)
	# write.table(c_Fdf_qc1_noQsnp, "./output_4F/F4_GSEM_Adiposity_noQsnp_sumstats.txt", sep = "\t", quote = F, col.names = T, row.names = F)
}
EOF

Rscript combineOutputFiles.R 1
Rscript combineOutputFiles.R 2
Rscript combineOutputFiles.R 3
Rscript combineOutputFiles.R 4


