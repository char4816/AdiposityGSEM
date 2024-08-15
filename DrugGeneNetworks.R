setwd("working/directory")

library(data.table)
library(tidyverse)
library(UpSetR)
library(ComplexHeatmap)
library(PheWAS)
library(igraph)
library(corrr)
set.seed(22)

range01 <- function(x){(x-min(x))/(max(x)-min(x))}
'%ni%' <- Negate('%in%')

st_Red <- "#e6194B"; st_Green <- "#3cb44b"; st_Yellow <- "#ffe119"; st_Blue <- "#4363d8"
st_Orange <- "#f58231"; st_Purple <- "#911eb4"; st_Cyan <- "#42d4f4"; st_Magenta <- "#f032e6"
st_Lime <- "#bfef45"; st_Pink <- "#fabed4"; st_Teal <- "#469990"; st_Lavender <- "#dcbeff"
st_Brown <- "#9A6324"; st_Beige <- "#fffac8"; st_Maroon <- "#800000"; st_Mint <- "#aaffc3"
st_Olive <- "#808000"; st_Apricot <- "#ffd8b1"; st_Navy <- "#000075"; st_brick <- "firebrick"
st_Chartreuse <- "chartreuse1"; st_Grey <- "mistyrose4"; st_White <- "grey85"; st_Black <- "grey22"
st_Cyan2 <- "cyan2"; st_Pink2 <- "hotpink"; st_Green2 <- "springgreen2"; st_Blue2 <- "skyblue1"
st_colors <- c(st_Red, st_Green, st_Yellow, st_Blue,
               st_Orange, st_Purple, st_Cyan, st_Magenta,
               st_Lime, st_Pink, st_Teal, st_Lavender,
               st_Brown, st_Beige, st_Maroon, st_Mint,
               st_Olive, st_Apricot, st_Navy, st_brick,
               st_Chartreuse, st_Grey, st_White, st_Black,
               st_Cyan2, st_Pink2, st_Green2, st_Blue2)


# The Drug Repurposing Hub: a next-generation drug library and information resource
# "Please cite usage as: Corsello SM, et al. Nature Medicine. 2017 Apr 7;23(4):405-408. doi: 10.1038/nm.4306."	
# !Contact	repurposing@broadinstitute.org	
# https://repo-hub.broadinstitute.org/repurposing#home
# https://repo-hub.broadinstitute.org/repurposing#download-data
# Drug information
# Latest version: 3/24/2020
# wget https://s3.amazonaws.com/data.clue.io/repurposing/downloads/repurposing_drugs_20200324.txt
# wget https://s3.amazonaws.com/data.clue.io/repurposing/downloads/repurposing_samples_20200324.txt

# Read in the DRH ############

DRH_drugs <- fread("./DRH/repurposing_drugs_20200324.txt", header = T, skip =9)
DRH_samples <- fread("./DRH/repurposing_samples_20200324.txt", header = T, skip =9)
head(DRH_drugs)
DRH_drugs_launched <- subset(DRH_drugs, clinical_phase == "Launched")
DRH_drugs_launched <- unique(DRH_drugs_launched$pert_iname)
DRH_drugs_launched <- toupper(DRH_drugs_launched)
length(DRH_drugs_launched)
# 2427 unique drugs that "Launched"

# Drug-Gene Interaction Database (DGIdb) ############
# https://www.dgidb.org/downloads
# 2023-Dec
# https://www.dgidb.org/data/2023-Dec/interactions.tsv
# https://www.dgidb.org/data/2023-Dec/genes.tsv
# https://www.dgidb.org/data/2023-Dec/drugs.tsv
# https://www.dgidb.org/data/2023-Dec/categories.tsv
# DGIdb v5
# https://github.com/dgidb/dgidb-v5

DGIdb_interactions <- fread("./DGIdb/interactions.tsv", header = T)
DGIdb_genes <- fread("./DGIdb/genes.tsv", header = T)
DGIdb_drugs <- fread("./DGIdb/drugs.tsv", header = T)
DGIdb_categories <- fread("./DGIdb/categories.tsv", header = T)

DGIdb_drugs_launched <- subset(DGIdb_drugs, approved == "TRUE")
DGIdb_drugs_launched <- unique(DGIdb_drugs_launched$drug_claim_name)
DGIdb_drugs_launched <- toupper(DGIdb_drugs_launched)
length(DGIdb_drugs_launched)

# Ensembl to HGNC ############
# MAPPING ENSEMBL GENE IDs (that are missing)
# https://www.genenames.org/download/custom/
# Approved Symbol
# The official gene symbol that has been approved by the HGNC and is publicly available. 
# Symbols are approved based on specific HGNC nomenclature guidelines. In the HTML results page this ID links to the HGNC Symbol Report for that gene.
# Ensembl Gene ID
# This column contains a manually curated Ensembl Gene ID.
HGNC_Ensembl <- fread("./GeneMaps/HGNC_Ensembl_Map.txt", header = T)
colnames(HGNC_Ensembl) <- c("HGNC_map", "Ensembl_map")

mapHGNC_to_Ensembl <- function(df = df, map = HGNC_Ensembl, replace = "Gene", mergex = "EnsemblGene", mergey = "Ensembl_map"){
  mdf <- left_join(df, map, by = setNames(nm=mergex, mergey))
  mdf[[replace]][which(mdf[[replace]] == "-")] <- mdf$HGNC_map[which(mdf[[replace]] == "-")]
  mdf[[replace]][is.na(mdf[[replace]])] <- mdf[[mergex]][is.na(mdf[[replace]])]
  mdf$HGNC_map <- NULL
  return(mdf)
}

# The 4 factor's gene lists ############
F1_SNP_to_Gene <- fread("./DEPICT/F1_GSEM_Birth_Size/F1_GSEM_Birth_Size_geneprioritization.txt", header = T)
F1_SNP_to_Genem <- F1_SNP_to_Gene[,c("Locus", "Ensembl gene ID", "Gene symbol", "Nominal P value")]; F1_SNP_to_Genem$Locus <- as.character(F1_SNP_to_Genem$Locus); F1_SNP_to_Genem$`Gene symbol` <- as.character(F1_SNP_to_Genem$`Gene symbol`); F1_SNP_to_Genem$`Nominal P value` <- as.numeric(F1_SNP_to_Genem$`Nominal P value`); colnames(F1_SNP_to_Genem) <- c("SNP", "EnsemblGene", "Gene", "Gpval")
F1_SNP_to_Genem <- mapHGNC_to_Ensembl(df = F1_SNP_to_Genem, map = HGNC_Ensembl, replace = "Gene", mergex = "EnsemblGene", mergey = "Ensembl_map")
F1_SNP_PrioritizedGene <- fread("./DEPICT/F1_GSEM_Birth_Size/F1_GSEM_Birth_Size_PrioritizedSNP2Gene.txt", header = T)
F1_FUMA_Genes <- fread("./DEPICT/F1_GSEM_Birth_Size/FUMA_gene2func127955/geneIDs.txt", header = T)

F2_SNP_to_Gene <- fread("./DEPICT/F2_GSEM_Abdominal_Size/F2_GSEM_Abdominal_Size_geneprioritization.txt", header = T)
F2_SNP_to_Genem <- F2_SNP_to_Gene[,c("Locus", "Ensembl gene ID", "Gene symbol", "Nominal P value")]; F2_SNP_to_Genem$Locus <- as.character(F2_SNP_to_Genem$Locus); F2_SNP_to_Genem$`Gene symbol` <- as.character(F2_SNP_to_Genem$`Gene symbol`); F2_SNP_to_Genem$`Nominal P value` <- as.numeric(F2_SNP_to_Genem$`Nominal P value`); colnames(F2_SNP_to_Genem) <- c("SNP", "EnsemblGene", "Gene", "Gpval")
F2_SNP_to_Genem <- mapHGNC_to_Ensembl(df = F2_SNP_to_Genem, map = HGNC_Ensembl, replace = "Gene", mergex = "EnsemblGene", mergey = "Ensembl_map")
F2_SNP_PrioritizedGene  <- fread("./DEPICT/F2_GSEM_Abdominal_Size/F2_GSEM_Abdominal_Size_PrioritizedSNP2Gene.txt", header = T)
F2_FUMA_Genes <- fread("./F2_GSEM_Abdominal_Size/FUMA_gene2func127956/geneIDs.txt", header = T)

F3_SNP_to_Gene <- fread("./DEPICT/F3_GSEM_Body_Size_and_Adipose_Distribution/F3_GSEM_Body_Size_and_Adipose_Distribution_geneprioritization.txt", header = T)
F3_SNP_to_Genem <- F3_SNP_to_Gene[,c("Locus", "Ensembl gene ID", "Gene symbol", "Nominal P value")]; F3_SNP_to_Genem$Locus <- as.character(F3_SNP_to_Genem$Locus); F3_SNP_to_Genem$`Gene symbol` <- as.character(F3_SNP_to_Genem$`Gene symbol`); F3_SNP_to_Genem$`Nominal P value` <- as.numeric(F3_SNP_to_Genem$`Nominal P value`); colnames(F3_SNP_to_Genem) <- c("SNP", "EnsemblGene", "Gene", "Gpval")
F3_SNP_to_Genem <- mapHGNC_to_Ensembl(df = F3_SNP_to_Genem, map = HGNC_Ensembl, replace = "Gene", mergex = "EnsemblGene", mergey = "Ensembl_map")
F3_SNP_PrioritizedGene  <- fread("./DEPICT/F3_GSEM_Body_Size_and_Adipose_Distribution/F3_GSEM_Body_Size_and_Adipose_Distribution_PrioritizedSNP2Gene.txt", header = T)
F3_FUMA_Genes <- fread("./DEPICT/F3_GSEM_Body_Size_and_Adipose_Distribution/FUMA_gene2func128062/geneIDs.txt", header = T)

F4_SNP_to_Gene <- fread("./DEPICT/F4_GSEM_Adiposity/F4_GSEM_Adiposity_geneprioritization.txt", header = T)
F4_SNP_to_Genem <- F4_SNP_to_Gene[,c("Locus", "Ensembl gene ID", "Gene symbol", "Nominal P value")]; F4_SNP_to_Genem$Locus <- as.character(F4_SNP_to_Genem$Locus); F4_SNP_to_Genem$`Gene symbol` <- as.character(F4_SNP_to_Genem$`Gene symbol`); F4_SNP_to_Genem$`Nominal P value` <- as.numeric(F4_SNP_to_Genem$`Nominal P value`); colnames(F4_SNP_to_Genem) <- c("SNP", "EnsemblGene", "Gene", "Gpval")
F4_SNP_to_Genem <- mapHGNC_to_Ensembl(df = F4_SNP_to_Genem, map = HGNC_Ensembl, replace = "Gene", mergex = "EnsemblGene", mergey = "Ensembl_map")
F4_SNP_PrioritizedGene <- fread("./DEPICT/F4_GSEM_Adiposity/F4_GSEM_Adiposity_PrioritizedSNP2Gene.txt", header = T)
F4_FUMA_Genes <- fread("./DEPICT/F4_GSEM_Adiposity/FUMA_gene2func127958/geneIDs.txt", header = T)

# QC for DRH and DGIdb ############
# These databases include information on the MOA for each drug (eg, antagonist). This allowed for matching drugs to gene products that are likely to have therapeutic, as opposed to adverse, effects based on whether upward or downward expression was associated with psychiatric risk.

# In line with prior quality control of the DGIdb resource, only drug-gene pairs with interaction scores greater than 0.5 were retained.
# The interaction score reflects the product of supporting publications and the relative drug and gene specificity, such that higher values indicate greater support of a drug-gene interaction.
# DGIdb_interactions_05 <- subset(DGIdb_interactions, interaction_score > 0.5)
DGIdb_interactions$interaction_score <- as.numeric(DGIdb_interactions$interaction_score)
DGIdb_interactions_05 <- DGIdb_interactions[which(DGIdb_interactions$interaction_score > 0.5)]
dim(DGIdb_interactions_05)
DGIdb_interactions_05 <- DGIdb_interactions_05[,c("gene_name", "gene_claim_name", "drug_name", "drug_claim_name", "interaction_type", "interaction_score")]

F4_merge <- left_join(DGIdb_interactions_05, F4_SNP_PrioritizedGene, by = c("gene_name" = "F4_Gene"))

# Reformat the DRH database, split the pipe separated gene lists into columns
# This is important to do because it makes our search match exactly the gene name
DRH_drugs_targets <- as.data.frame(DRH_drugs$target)
colnames(DRH_drugs_targets) <- c("target")
# DRH_drugs_targets_split <- within(DRH_drugs_targets, FOO<-data.frame(do.call('rbind', strsplit(as.character(FOO), '|', fixed=TRUE))))
DRH_drugs_targets_split <- DRH_drugs_targets %>% separate(target, into = paste0("target_", 1:144), sep = "\\|", extra = "merge", fill = "right")
summary(is.na(DRH_drugs_targets_split))
# looks like the longest number of splits is 144
DRH_drugs_targets_split[which(is.na(DRH_drugs_targets_split$target_144) == F),]
DRH_drugs_targets[which(is.na(DRH_drugs_targets_split$target_144) == F),]

# Bring in the FOCUS results ##############

# tissueMap <- fread("tissueMap.txt", header = T)
# tissueMap <- as.data.frame(rbind(as.matrix(tissueMap),cbind(
#   c("stomach"),
#   c("stomach")
# )))
# colnames(tissueMap) <- c("tissue", "general_tissue")
# write.table(tissueMap, "tissueMap.txt", sep = "\t", quote = F, col.names = T, row.names = F)
tissueMap <- fread("tissueMap.txt", header=T)

# low PIP: >0.10
# moderate PIP: >0.25 
# high PIP: >0.75
FOCUS_PIP_THRESH <- 0.1

F1_FOCUS <- fread(".//ma-focus/F1/F1_Significant_NonNull_Blocks_CS.txt", header = T)
hist(F1_FOCUS$pips_pop1, breaks = 50)
F1_FOCUS$ens_gene_id <- gsub("\\..*","", F1_FOCUS$ens_gene_id)
unique(F1_FOCUS$tissue)
nrow(F1_FOCUS) #1584
length(unique(F1_FOCUS$block)) #69
F1_FOCUS <- subset(F1_FOCUS, pips_pop1 > FOCUS_PIP_THRESH)
F1_FOCUS <- left_join(F1_FOCUS, tissueMap, by = "tissue")
nrow(F1_FOCUS) #215
length(unique(F1_FOCUS$block)) #69
F1_FOCUS$general_tissue <- as.factor(F1_FOCUS$general_tissue)

length(table(F1_FOCUS$general_tissue)[order(table(F1_FOCUS$general_tissue))])
table(F1_FOCUS$general_tissue)[order(table(F1_FOCUS$general_tissue))]

focus_pip_F1 <- ggplot(F1_FOCUS, aes_string(y = "general_tissue", x = "pips_pop1", fill = "general_tissue")) + 
  # geom_histogram(alpha = 0.5, position = position_dodge(width = 0.1), binwidth = .5) + 
  geom_boxplot(show.legend = F, outlier.size = 0.6) + geom_point(size = .6, show.legend = F) + 
  # scale_fill_manual(values = st_colors) +
  scale_x_continuous(breaks=seq(0,1,by=.1)) +
  ggtitle(paste0("F1: ", length(unique(F1_FOCUS$block)), " blocks with non-null 90% CS, ", length(unique(F1_FOCUS$ens_gene_id)), " genes with PIP > ", FOCUS_PIP_THRESH)) +
  ylab("General Tissue") +
  xlab("PIP") +
  theme_bw() 

png("./F1_FOCUS_Tissues_PIPs_Above_Thresh_Distribution.png", width = 8, height = 5, unit = "in", res = 300)
focus_pip_F1
dev.off()

F2_FOCUS <- fread(".//ma-focus/F2/F2_Significant_NonNull_Blocks_CS.txt", header = T)
hist(F2_FOCUS$pips_pop1, breaks = 50)
F2_FOCUS$ens_gene_id <- gsub("\\..*","", F2_FOCUS$ens_gene_id)
unique(F2_FOCUS$tissue)
nrow(F2_FOCUS) #5242
length(unique(F2_FOCUS$block)) #243
F2_FOCUS <- subset(F2_FOCUS, pips_pop1 > FOCUS_PIP_THRESH)
F2_FOCUS <- left_join(F2_FOCUS, tissueMap, by = "tissue")
nrow(F2_FOCUS) #862
length(unique(F2_FOCUS$block)) #243
F2_FOCUS$general_tissue <- as.factor(F2_FOCUS$general_tissue)

length(table(F2_FOCUS$general_tissue)[order(table(F2_FOCUS$general_tissue))])
table(F2_FOCUS$general_tissue)[order(table(F2_FOCUS$general_tissue))]

focus_pip_F2 <- ggplot(F2_FOCUS, aes_string(y = "general_tissue", x = "pips_pop1", fill = "general_tissue")) + 
  # geom_histogram(alpha = 0.5, position = position_dodge(width = 0.1), binwidth = .5) + 
  geom_boxplot(show.legend = F, outlier.size = 0.6) + geom_point(size = .6, show.legend = F) + 
  # scale_fill_manual(values = st_colors) +
  scale_x_continuous(breaks=seq(0,1,by=.1)) +
  ggtitle(paste0("F2: ", length(unique(F2_FOCUS$block)), " blocks with non-null 90% CS, ", length(unique(F2_FOCUS$ens_gene_id)), " genes with PIP > ", FOCUS_PIP_THRESH)) +
  ylab("General Tissue") +
  xlab("PIP") +
  theme_bw() 

png("./F2_FOCUS_Tissues_PIPs_Above_Thresh_Distribution.png", width = 8, height = 5, unit = "in", res = 300)
focus_pip_F2
dev.off()

F3_FOCUS <- fread(".//ma-focus/F3/F3_Significant_NonNull_Blocks_CS.txt", header = T)
hist(F3_FOCUS$pips_pop1, breaks = 50)
F3_FOCUS$ens_gene_id <- gsub("\\..*","", F3_FOCUS$ens_gene_id)
unique(F3_FOCUS$tissue)
nrow(F3_FOCUS) #8011
length(unique(F3_FOCUS$block)) #690
F3_FOCUS <- subset(F3_FOCUS, pips_pop1 > FOCUS_PIP_THRESH)
F3_FOCUS <- left_join(F3_FOCUS, tissueMap, by = "tissue")
nrow(F3_FOCUS) #2944
length(unique(F3_FOCUS$block)) #689
F3_FOCUS$general_tissue <- as.factor(F3_FOCUS$general_tissue)

F4_FOCUS <- fread(".//ma-focus/F4/F4_Significant_NonNull_Blocks_CS.txt", header = T)
hist(F4_FOCUS$pips_pop1, breaks = 50)
F4_FOCUS$ens_gene_id <- gsub("\\..*","", F4_FOCUS$ens_gene_id)
unique(F4_FOCUS$tissue)
nrow(F4_FOCUS) #1014
length(unique(F4_FOCUS$block)) #335
F4_FOCUS <- subset(F4_FOCUS, pips_pop1 > FOCUS_PIP_THRESH)
F4_FOCUS <- left_join(F4_FOCUS, tissueMap, by = "tissue")
nrow(F4_FOCUS) #864
length(unique(F4_FOCUS$block)) #335
F4_FOCUS$general_tissue <- as.factor(F4_FOCUS$general_tissue)

length(table(F4_FOCUS$general_tissue)[order(table(F4_FOCUS$general_tissue))])
table(F4_FOCUS$general_tissue)[order(table(F4_FOCUS$general_tissue))]

focus_pip_F4 <- ggplot(F4_FOCUS, aes_string(y = "general_tissue", x = "pips_pop1", fill = "general_tissue")) + 
  # geom_histogram(alpha = 0.5, position = position_dodge(width = 0.1), binwidth = .5) + 
  geom_boxplot(show.legend = F, outlier.size = 0.6) + geom_point(size = .6, show.legend = F) + 
  # scale_fill_manual(values = st_colors) +
  scale_x_continuous(breaks=seq(0,1,by=.1)) +
  ggtitle(paste0("F4: ", length(unique(F4_FOCUS$block)), " blocks with non-null 90% CS, ", length(unique(F4_FOCUS$ens_gene_id)), " genes with PIP > ", FOCUS_PIP_THRESH)) +
  ylab("General Tissue") +
  xlab("PIP") +
  theme_bw() 

png("./F4_FOCUS_Tissues_PIPs_Above_Thresh_Distribution.png", width = 8, height = 5, unit = "in", res = 300)
focus_pip_F4
dev.off()

BMI_FOCUS <- fread(".//ma-focus/BMI/BMI_Significant_NonNull_Blocks_CS.txt", header = T)
hist(BMI_FOCUS$pips_pop1, breaks = 50)
BMI_FOCUS$ens_gene_id <- gsub("\\..*","", BMI_FOCUS$ens_gene_id)
unique(BMI_FOCUS$tissue)
nrow(BMI_FOCUS) #1085
length(unique(BMI_FOCUS$block)) #376
BMI_FOCUS <- subset(BMI_FOCUS, pips_pop1 > FOCUS_PIP_THRESH)
BMI_FOCUS <- left_join(BMI_FOCUS, tissueMap, by = "tissue")
nrow(BMI_FOCUS) #940
length(unique(BMI_FOCUS$block)) #376
BMI_FOCUS$general_tissue <- as.factor(BMI_FOCUS$general_tissue)

length(table(BMI_FOCUS$general_tissue)[order(table(BMI_FOCUS$general_tissue))])
table(BMI_FOCUS$general_tissue)[order(table(BMI_FOCUS$general_tissue))]

focus_pip_BMI <- ggplot(BMI_FOCUS, aes_string(y = "general_tissue", x = "pips_pop1", fill = "general_tissue")) + 
  # geom_histogram(alpha = 0.5, position = position_dodge(width = 0.1), binwidth = .5) + 
  geom_boxplot(show.legend = F, outlier.size = 0.6) + geom_point(size = .6, show.legend = F) + 
  # scale_fill_manual(values = st_colors) +
  scale_x_continuous(breaks=seq(0,1,by=.1)) +
  ggtitle(paste0("BMI: ", length(unique(BMI_FOCUS$block)), " blocks with non-null 90% CS, ", length(unique(BMI_FOCUS$ens_gene_id)), " genes with PIP > ", FOCUS_PIP_THRESH)) +
  ylab("General Tissue") +
  xlab("PIP") +
  theme_bw() 

png("./BMI_FOCUS_Tissues_PIPs_Above_Thresh_Distribution.png", width = 8, height = 5, unit = "in", res = 300)
focus_pip_BMI
dev.off()

# stop()
# popout genes and upset plot ############

upsetList <- list(
  # DEPICT possible genes that are worth considering at each independent GWAS locus
  "DEPICT_S2G" = F1_SNP_to_Gene$`Ensembl gene ID`[which(F1_SNP_to_Gene$`False discovery rate` %in% c("<=0.01",  "<0.05"))],
  # significant genes from FOCUS non-Null credible sets
  "TWAS_FOCUS" = F1_FOCUS$ens_gene_id
)
F1_upsetMatrix <- list_to_matrix(upsetList)
# Plot
F1upset <- upset(fromList(upsetList), 
                 nintersects = 40, nsets = 7, order.by = "degree", decreasing = T, 
                 mb.ratio = c(0.6, 0.4),number.angles = 0, text.scale = 1.1, 
                 point.size = 2.8, line.size = 1
)
F1upset

upsetList <- list(
  # DEPICT possible genes that are worth considering at each independent GWAS locus
  "DEPICT_S2G" = F2_SNP_to_Gene$`Ensembl gene ID`[which(F2_SNP_to_Gene$`False discovery rate` %in% c("<=0.01",  "<0.05"))],
  # significant genes from FOCUS non-Null credible sets
  "TWAS_FOCUS" = F2_FOCUS$ens_gene_id
)
F2_upsetMatrix <- list_to_matrix(upsetList)
# Plot
F2upset <- upset(fromList(upsetList), 
                 nintersects = 40, nsets = 6, order.by = "degree", decreasing = T, 
                 mb.ratio = c(0.6, 0.4),number.angles = 0, text.scale = 1.1, 
                 point.size = 2.8, line.size = 1
)
F2upset

upsetList <- list(
  # DEPICT possible genes that are worth considering at each independent GWAS locus
  "DEPICT_S2G" = F3_SNP_to_Gene$`Ensembl gene ID`[which(F3_SNP_to_Gene$`False discovery rate` %in% c("<=0.01",  "<0.05"))],

  # significant genes from FOCUS non-Null credible sets
  "TWAS_FOCUS" = F3_FOCUS$ens_gene_id
)
F3_upsetMatrix <- list_to_matrix(upsetList)
# Plot
F3upset <- upset(fromList(upsetList), 
                 nintersects = 40, nsets = 6, order.by = "degree", decreasing = T, 
                 mb.ratio = c(0.6, 0.4),number.angles = 0, text.scale = 1.1, 
                 point.size = 2.8, line.size = 1
)
F3upset

upsetList <- list(
  # DEPICT possible genes that are worth considering at each independent GWAS locus
  "DEPICT_S2G" = F4_SNP_to_Gene$`Ensembl gene ID`[which(F4_SNP_to_Gene$`False discovery rate` %in% c("<=0.01",  "<0.05"))],
  # significant genes from FOCUS non-Null credible sets
  "TWAS_FOCUS" = F4_FOCUS$ens_gene_id
)
F4_upsetMatrix <- list_to_matrix(upsetList)
# Plot
F4upset <- upset(fromList(upsetList), 
                 nintersects = 40, nsets = 6, order.by = "degree", decreasing = T, 
                 mb.ratio = c(0.6, 0.4),number.angles = 0, text.scale = 1.1, 
                 point.size = 2.8, line.size = 1
)
F4upset

# Read in the novel genes to F4, that are not in BMI GWAS
NovelGenes <- fread("./DEPICT/Genes_in_F4_not_in_BMI.txt", header = T)
BMI_SNP_to_Gene <- fread("./DEPICT/BMI/BMI_geneprioritization.txt", header = T)
BMI_GeneHits <- BMI_SNP_to_Gene$`Ensembl gene ID`[which(BMI_SNP_to_Gene$`False discovery rate` %in% c("<=0.01",  "<0.05"))]

upsetList <- list(
  "F4_TWAS_FOCUS" = F4_FOCUS$ens_gene_id,
  "BMI_TWAS_FOCUS" = BMI_FOCUS$ens_gene_id
)
F4_BMI_FOCUS_upsetMatrix <- list_to_matrix(upsetList)
F4_BMI_FOCUS_upset <- upset(fromList(upsetList), 
                     nintersects = 40, nsets = 16, order.by = "degree", decreasing = T, 
                     mb.ratio = c(0.6, 0.4),number.angles = 0, text.scale = 1.1, 
                     point.size = 2.8, line.size = 1
)
png("./F4_BMI_FOCUS_Upset.png", width = 6, height = 4, unit = "in", res = 300)
F4_BMI_FOCUS_upset
dev.off()

upsetList <- list(
  "F4_GWAS_DEPICT" = F4_SNP_to_Gene$`Ensembl gene ID`[which(F4_SNP_to_Gene$`False discovery rate` %in% c("<=0.01",  "<0.05"))],
  "BMI_GWAS_DEPICT" = BMI_GeneHits
)
F4_BMI_DEPICT_upsetMatrix <- list_to_matrix(upsetList)
F4_BMI_DEPICT_upset <- upset(fromList(upsetList), 
                            nintersects = 40, nsets = 16, order.by = "degree", decreasing = T, 
                            mb.ratio = c(0.6, 0.4),number.angles = 0, text.scale = 1.1, 
                            point.size = 2.8, line.size = 1
)
png("./F4_BMI_DEPICT_Upset.png", width = 6, height = 4, unit = "in", res = 300)
F4_BMI_DEPICT_upset
dev.off()

upsetList <- list(
  "F4_GWAS_DEPICT" = F4_SNP_to_Gene$`Ensembl gene ID`[which(F4_SNP_to_Gene$`False discovery rate` %in% c("<=0.01",  "<0.05"))],
  "F4_TWAS_FOCUS" = F4_FOCUS$ens_gene_id,
  "BMI_TWAS_FOCUS" = BMI_FOCUS$ens_gene_id,
  "BMI_GWAS_DEPICT" = BMI_GeneHits
   )
F4_BMI_upsetMatrix <- list_to_matrix(upsetList)
# Plot
F4_BMIupset <- upset(fromList(upsetList), 
      nintersects = 40, nsets = 16, order.by = "degree", decreasing = T, 
      mb.ratio = c(0.6, 0.4),number.angles = 0, text.scale = 1.1, 
      point.size = 2.8, line.size = 1
)
# stop("bingo")
png("./F4_BMI_DEPCIT_FOCUS_Upset.png", width = 6, height = 4, unit = "in", res = 300)
F4_BMIupset
dev.off()

upsetList <- list(
  "F4_GWAS_DEPICT" = F4_SNP_to_Gene$`Ensembl gene ID`[which(F4_SNP_to_Gene$`False discovery rate` %in% c("<=0.01",  "<0.05"))],
  "F4_TWAS_FOCUS" = F4_FOCUS$ens_gene_id,
  "F3_GWAS_DEPICT" = F3_SNP_to_Gene$`Ensembl gene ID`[which(F3_SNP_to_Gene$`False discovery rate` %in% c("<=0.01",  "<0.05"))],
  "F3_TWAS_FOCUS" = F3_FOCUS$ens_gene_id,
  "F2_GWAS_DEPICT" = F2_SNP_to_Gene$`Ensembl gene ID`[which(F2_SNP_to_Gene$`False discovery rate` %in% c("<=0.01",  "<0.05"))],
  "F2_TWAS_FOCUS" = F2_FOCUS$ens_gene_id,
  "F1_GWAS_DEPICT" = F1_SNP_to_Gene$`Ensembl gene ID`[which(F1_SNP_to_Gene$`False discovery rate` %in% c("<=0.01",  "<0.05"))],
  "F1_TWAS_FOCUS" = F1_FOCUS$ens_gene_id
  # "BMI_GWAS_DEPICT" = BMI_SNP_to_Gene$`Ensembl gene ID`[which(BMI_SNP_to_Gene$`False discovery rate` %in% c("<=0.01",  "<0.05"))],
  # "BMI_TWAS_FOCUS" = BMI_FOCUS$ens_gene_id
)
lapply(upsetList, length)

# F4_GWAS_DEPICT: 437
write.table(F4_SNP_to_Gene[which(F4_SNP_to_Gene$`False discovery rate` %in% c("<=0.01",  "<0.05")),c("Ensembl gene ID","Gene symbol")], "F4_DEPICT_Identified_Genes.txt", sep = "\t", quote = F, col.names = T, row.names = F)
# F4_TWAS_FOCUS: 864
write.table(F4_FOCUS[,c("ens_gene_id","mol_name")], "F4_FOCUS_Identified_Genes.txt", sep = "\t", quote = F, col.names = T, row.names = F)
# F3_GWAS_DEPICT: 1864
write.table(F3_SNP_to_Gene[which(F3_SNP_to_Gene$`False discovery rate` %in% c("<=0.01",  "<0.05")),c("Ensembl gene ID","Gene symbol")], "F3_DEPICT_Identified_Genes.txt", sep = "\t", quote = F, col.names = T, row.names = F)
# F3_TWAS_FOCUS: 2944
write.table(F3_FOCUS[,c("ens_gene_id","mol_name")], "F3_FOCUS_Identified_Genes.txt", sep = "\t", quote = F, col.names = T, row.names = F)
# F2_GWAS_DEPICT: 319
write.table(F2_SNP_to_Gene[which(F2_SNP_to_Gene$`False discovery rate` %in% c("<=0.01",  "<0.05")),c("Ensembl gene ID","Gene symbol")], "F2_DEPICT_Identified_Genes.txt", sep = "\t", quote = F, col.names = T, row.names = F)
# F2_TWAS_FOCUS: 862
write.table(F2_FOCUS[,c("ens_gene_id","mol_name")], "F2_FOCUS_Identified_Genes.txt", sep = "\t", quote = F, col.names = T, row.names = F)
# F1_GWAS_DEPICT: 24
write.table(F1_SNP_to_Gene[which(F1_SNP_to_Gene$`False discovery rate` %in% c("<=0.01",  "<0.05")),c("Ensembl gene ID","Gene symbol")], "F1_DEPICT_Identified_Genes.txt", sep = "\t", quote = F, col.names = T, row.names = F)
# F1_TWAS_FOCUS: 215
write.table(F1_FOCUS[,c("ens_gene_id","mol_name")], "F1_FOCUS_Identified_Genes.txt", sep = "\t", quote = F, col.names = T, row.names = F)
# BMI_GWAS_DEPICT: 1109
write.table(BMI_SNP_to_Gene[which(BMI_SNP_to_Gene$`False discovery rate` %in% c("<=0.01",  "<0.05")),c("Ensembl gene ID","Gene symbol")], "BMI_DEPICT_Identified_Genes.txt", sep = "\t", quote = F, col.names = T, row.names = F)
# BMI_TWAS_FOCUS: 940
write.table(BMI_FOCUS[,c("ens_gene_id","mol_name")], "BMI_FOCUS_Identified_Genes.txt", sep = "\t", quote = F, col.names = T, row.names = F)

FAll_upsetMatrix <- list_to_matrix(upsetList)
# Plot
FAllupset <- upset(fromList(upsetList), 
                 nintersects = 40, nsets = 10, order.by = "freq", decreasing = T, 
                 mb.ratio = c(0.6, 0.4),number.angles = 0, text.scale = 1.1, 
                 point.size = 2.8, line.size = 1
)
png("./F1_F2_F3_F4_DEPCIT_FOCUS_Upset.png", width = 14, height = 6, unit = "in", res = 300)
FAllupset
dev.off()


upsetList <- list(
  "F4_GWAS_DEPICT" = F4_SNP_to_Gene$`Ensembl gene ID`[which(F4_SNP_to_Gene$`False discovery rate` %in% c("<=0.01",  "<0.05"))],
  "F3_GWAS_DEPICT" = F3_SNP_to_Gene$`Ensembl gene ID`[which(F3_SNP_to_Gene$`False discovery rate` %in% c("<=0.01",  "<0.05"))],
  "F2_GWAS_DEPICT" = F2_SNP_to_Gene$`Ensembl gene ID`[which(F2_SNP_to_Gene$`False discovery rate` %in% c("<=0.01",  "<0.05"))],
  "F1_GWAS_DEPICT" = F1_SNP_to_Gene$`Ensembl gene ID`[which(F1_SNP_to_Gene$`False discovery rate` %in% c("<=0.01",  "<0.05"))]
)
FAll_upsetMatrix <- list_to_matrix(upsetList)
# Plot
FAllupset <- upset(fromList(upsetList), 
                   nintersects = 40, nsets = 10, order.by = "freq", decreasing = T, 
                   mb.ratio = c(0.6, 0.4),number.angles = 0, text.scale = 1.1, 
                   point.size = 2.8, line.size = 1
)
png("./F1_F2_F3_F4_DEPCIT_Upset.png", width = 14, height = 6, unit = "in", res = 300)
FAllupset
dev.off()

upsetList <- list(
  "F4_TWAS_FOCUS" = F4_FOCUS$ens_gene_id,
  "F3_TWAS_FOCUS" = F3_FOCUS$ens_gene_id,
  "F2_TWAS_FOCUS" = F2_FOCUS$ens_gene_id,
  "F1_TWAS_FOCUS" = F1_FOCUS$ens_gene_id
)
FAll_upsetMatrix <- list_to_matrix(upsetList)
# Plot
FAllupset <- upset(fromList(upsetList), 
                   nintersects = 40, nsets = 10, order.by = "freq", decreasing = T, 
                   mb.ratio = c(0.6, 0.4),number.angles = 0, text.scale = 1.1, 
                   point.size = 2.8, line.size = 1
)
png("./F1_F2_F3_F4_FOCUS_Upset.png", width = 14, height = 6, unit = "in", res = 300)
FAllupset
dev.off()

# stop("pause")

# Search the DRH and the DGIdb ############
search_Drug_Databases <- function(geneVector = c(NA, NA), 
                                  funDRH_drugs_targets_split = DRH_drugs_targets_split,
                                  funDRH_drugs = DRH_drugs,
                                  funDGIdb_interactions_05 = DGIdb_interactions_05
){
  funDRH_drugs_targets_split$rownum <- 1:nrow(funDRH_drugs_targets_split)
  funDRH_drugs_targets_split <- mutate_all(funDRH_drugs_targets_split, .funs=toupper)
  funDGIdb_interactions_05$rownum <- 1:nrow(funDGIdb_interactions_05)
  funDGIdb_interactions_05$gene_name <- toupper(funDGIdb_interactions_05$gene_name)
  if(sum(is.na(geneVector)) > 0){
    stop("You need to provide a non-NA vector of HGNC gene names in the geneVector argument")
  }else{
    
    DRH_output_list <- list()
    DGIdb_output_list <- list()
    
    for(g in 1:length(geneVector)){
      curgene <- geneVector[g]
      print(paste0("#--",g,"--",curgene,"--###################################"))
      DRH_rows <- c()
      DGIdb_rows <- c()
      
      # Identify rows of the DRH that match our gene
      print("----------rowsDRH----------------------------------------------")
      # https://stackoverflow.com/questions/28432119/using-grep-to-search-along-row-in-r
      DRH_rows <- which(rowSums(t(apply(funDRH_drugs_targets_split, 1, function(u) grepl(toupper(curgene),u)))) > 0)
      # for(drhg in 1:nrow(funDRH_drugs_targets_split)){
      #   if(toupper(curgene) %in% toupper(funDRH_drugs_targets_split[drhg,])){
      #     DRH_rows <- c(DRH_rows, drhg)
      #   }
      # }
      
      print("--------------rowsDGIdb----------------------------------------")
      # Identify rows of the DGIdb that match our gene
      funDGIdb_interactions_05_sub <- subset(funDGIdb_interactions_05, gene_name == toupper(curgene) )
      DGIdb_rows <- funDGIdb_interactions_05_sub$rownum
      # for(dgidbg in 1:nrow(funDGIdb_interactions_05)){
      #   if(toupper(curgene) %in% toupper(funDGIdb_interactions_05[dgidbg,c("gene_name")])){
      #     DGIdb_rows <- c(DGIdb_rows, dgidbg)
      #   }
      # }
      
      if(length(DRH_rows) > 0){
        print("----------DRH----------------------------------------------")
        DRH_output_list[[curgene]] <- funDRH_drugs[DRH_rows,]
      }else{
        # DRH_output_list[[curgene]] <- NA
      }
      
      if(length(DGIdb_rows) > 0){
        print("--------------DGIdb----------------------------------------")
        DGIdb_output_list[[curgene]] <- funDGIdb_interactions_05[DGIdb_rows,]
      }else{
        # DGIdb_output_list[[curgene]] <- NA
      }
    }
    returnList <- list("DRH_output_list" = DRH_output_list, "DGIdb_output_list" = DGIdb_output_list)
    return(returnList)
  }
}

print("F1 Search Drug Databases")
F1_upsetMatrix <- as.data.frame(F1_upsetMatrix)
F1_crossRef <- as.data.frame(rownames(F1_upsetMatrix)[which(F1_upsetMatrix$DEPICT_S2G == 1 | F1_upsetMatrix$TWAS_FOCUS == 1)])
colnames(F1_crossRef) <- c("Ensembl_map")
nrow(F1_crossRef)
F1_crossRef <- merge(F1_crossRef, HGNC_Ensembl, by = "Ensembl_map")
nrow(F1_crossRef)
F1_drugMatch <- search_Drug_Databases(geneVector = F1_crossRef$HGNC_map)
F1_drugMatch_DRH_output_list <- F1_drugMatch$DRH_output_list
F1_drugMatch_DRH <- bind_rows(F1_drugMatch_DRH_output_list, .id = "column_label")
F1_drugMatch_DGIdb_output_list <- F1_drugMatch$DGIdb_output_list
F1_drugMatch_DGIdb <- bind_rows(F1_drugMatch_DGIdb_output_list, .id = "column_label")
# incororate the upset plot, showing DEPICT and FOCUS hits
m_F1_upsetMatrix <- F1_upsetMatrix
m_F1_upsetMatrix$DEPICT_and_FOCUS <- 0; m_F1_upsetMatrix$DEPICT_and_FOCUS[which(m_F1_upsetMatrix$DEPICT_S2G == 1 & m_F1_upsetMatrix$TWAS_FOCUS == 1)] <- 1
m_F1_upsetMatrix$Ensembl_map <- rownames(m_F1_upsetMatrix)
m_F1_upsetMatrix <- left_join(m_F1_upsetMatrix, HGNC_Ensembl, by = "Ensembl_map")
f_F1_drugMatch_DRH <- left_join(F1_drugMatch_DRH, m_F1_upsetMatrix, by = c("column_label" = "HGNC_map"))
f_F1_drugMatch_DGIdb <- left_join(F1_drugMatch_DGIdb, m_F1_upsetMatrix, by = c("column_label" = "HGNC_map"))
f_F1_drugMatch_DRH <- f_F1_drugMatch_DRH %>% arrange(desc(TWAS_FOCUS)) %>% arrange(desc(DEPICT_S2G)) %>% arrange(desc(DEPICT_and_FOCUS))
f_F1_drugMatch_DGIdb <- f_F1_drugMatch_DGIdb %>% arrange(desc(TWAS_FOCUS)) %>% arrange(desc(DEPICT_S2G)) %>% arrange(desc(DEPICT_and_FOCUS))
# 
print("F2 Search Drug Databases")
F2_upsetMatrix <- as.data.frame(F2_upsetMatrix)
F2_crossRef <- as.data.frame(rownames(F2_upsetMatrix)[which(F2_upsetMatrix$DEPICT_S2G == 1 | F2_upsetMatrix$TWAS_FOCUS == 1)])
colnames(F2_crossRef) <- c("Ensembl_map")
nrow(F2_crossRef)
F2_crossRef <- merge(F2_crossRef, HGNC_Ensembl, by = "Ensembl_map")
nrow(F2_crossRef)
F2_drugMatch <- search_Drug_Databases(geneVector = F2_crossRef$HGNC_map)
F2_drugMatch_DRH_output_list <- F2_drugMatch$DRH_output_list
F2_drugMatch_DRH <- bind_rows(F2_drugMatch_DRH_output_list, .id = "column_label")
F2_drugMatch_DGIdb_output_list <- F2_drugMatch$DGIdb_output_list
F2_drugMatch_DGIdb <- bind_rows(F2_drugMatch_DGIdb_output_list, .id = "column_label")
# incororate the upset plot, showing DEPICT and FOCUS hits
m_F2_upsetMatrix <- F2_upsetMatrix
m_F2_upsetMatrix$DEPICT_and_FOCUS <- 0; m_F2_upsetMatrix$DEPICT_and_FOCUS[which(m_F2_upsetMatrix$DEPICT_S2G == 1 & m_F2_upsetMatrix$TWAS_FOCUS == 1)] <- 1
m_F2_upsetMatrix$Ensembl_map <- rownames(m_F2_upsetMatrix)
m_F2_upsetMatrix <- left_join(m_F2_upsetMatrix, HGNC_Ensembl, by = "Ensembl_map")
f_F2_drugMatch_DRH <- left_join(F2_drugMatch_DRH, m_F2_upsetMatrix, by = c("column_label" = "HGNC_map"))
f_F2_drugMatch_DGIdb <- left_join(F2_drugMatch_DGIdb, m_F2_upsetMatrix, by = c("column_label" = "HGNC_map"))
f_F2_drugMatch_DRH <- f_F2_drugMatch_DRH %>% arrange(desc(TWAS_FOCUS)) %>% arrange(desc(DEPICT_S2G)) %>% arrange(desc(DEPICT_and_FOCUS))
f_F2_drugMatch_DGIdb <- f_F2_drugMatch_DGIdb %>% arrange(desc(TWAS_FOCUS)) %>% arrange(desc(DEPICT_S2G)) %>% arrange(desc(DEPICT_and_FOCUS))
# 
print("F3 Search Drug Databases")
F3_upsetMatrix <- as.data.frame(F3_upsetMatrix)
F3_crossRef <- as.data.frame(rownames(F3_upsetMatrix)[which(F3_upsetMatrix$DEPICT_S2G == 1 | F3_upsetMatrix$TWAS_FOCUS == 1)])
colnames(F3_crossRef) <- c("Ensembl_map")
nrow(F3_crossRef)
F3_crossRef <- merge(F3_crossRef, HGNC_Ensembl, by = "Ensembl_map")
nrow(F3_crossRef)
F3_drugMatch <- search_Drug_Databases(geneVector = F3_crossRef$HGNC_map)
F3_drugMatch_DRH_output_list <- F3_drugMatch$DRH_output_list
F3_drugMatch_DRH <- bind_rows(F3_drugMatch_DRH_output_list, .id = "column_label")
F3_drugMatch_DGIdb_output_list <- F3_drugMatch$DGIdb_output_list
F3_drugMatch_DGIdb <- bind_rows(F3_drugMatch_DGIdb_output_list, .id = "column_label")
# incororate the upset plot, showing DEPICT and FOCUS hits
m_F3_upsetMatrix <- F3_upsetMatrix
m_F3_upsetMatrix$DEPICT_and_FOCUS <- 0; m_F3_upsetMatrix$DEPICT_and_FOCUS[which(m_F3_upsetMatrix$DEPICT_S2G == 1 & m_F3_upsetMatrix$TWAS_FOCUS == 1)] <- 1
m_F3_upsetMatrix$Ensembl_map <- rownames(m_F3_upsetMatrix)
m_F3_upsetMatrix <- left_join(m_F3_upsetMatrix, HGNC_Ensembl, by = "Ensembl_map")
f_F3_drugMatch_DRH <- left_join(F3_drugMatch_DRH, m_F3_upsetMatrix, by = c("column_label" = "HGNC_map"))
f_F3_drugMatch_DGIdb <- left_join(F3_drugMatch_DGIdb, m_F3_upsetMatrix, by = c("column_label" = "HGNC_map"))
f_F3_drugMatch_DRH <- f_F3_drugMatch_DRH %>% arrange(desc(TWAS_FOCUS)) %>% arrange(desc(DEPICT_S2G)) %>% arrange(desc(DEPICT_and_FOCUS))
f_F3_drugMatch_DGIdb <- f_F3_drugMatch_DGIdb %>% arrange(desc(TWAS_FOCUS)) %>% arrange(desc(DEPICT_S2G)) %>% arrange(desc(DEPICT_and_FOCUS))
# 
# 
print("F4 Search Drug Databases")
F4_upsetMatrix <- as.data.frame(F4_upsetMatrix)
F4_crossRef <- as.data.frame(rownames(F4_upsetMatrix)[which(F4_upsetMatrix$DEPICT_S2G == 1 | F4_upsetMatrix$TWAS_FOCUS == 1)])
colnames(F4_crossRef) <- c("Ensembl_map")
nrow(F4_crossRef)
F4_crossRef <- merge(F4_crossRef, HGNC_Ensembl, by = "Ensembl_map")
nrow(F4_crossRef) #1086
F4_drugMatch <- search_Drug_Databases(geneVector = F4_crossRef$HGNC_map)
F4_drugMatch_DRH_output_list <- F4_drugMatch$DRH_output_list
F4_drugMatch_DRH <- bind_rows(F4_drugMatch_DRH_output_list, .id = "column_label")
F4_drugMatch_DGIdb_output_list <- F4_drugMatch$DGIdb_output_list
F4_drugMatch_DGIdb <- bind_rows(F4_drugMatch_DGIdb_output_list, .id = "column_label")
# incororate the upset plot, showing DEPICT and FOCUS hits
m_F4_upsetMatrix <- F4_upsetMatrix
m_F4_upsetMatrix$DEPICT_and_FOCUS <- 0; m_F4_upsetMatrix$DEPICT_and_FOCUS[which(m_F4_upsetMatrix$DEPICT_S2G == 1 & m_F4_upsetMatrix$TWAS_FOCUS == 1)] <- 1
m_F4_upsetMatrix$Ensembl_map <- rownames(m_F4_upsetMatrix)
m_F4_upsetMatrix <- left_join(m_F4_upsetMatrix, HGNC_Ensembl, by = "Ensembl_map")
f_F4_drugMatch_DRH <- left_join(F4_drugMatch_DRH, m_F4_upsetMatrix, by = c("column_label" = "HGNC_map"))
f_F4_drugMatch_DGIdb <- left_join(F4_drugMatch_DGIdb, m_F4_upsetMatrix, by = c("column_label" = "HGNC_map"))
f_F4_drugMatch_DRH <- f_F4_drugMatch_DRH %>% arrange(desc(TWAS_FOCUS)) %>% arrange(desc(DEPICT_S2G)) %>% arrange(desc(DEPICT_and_FOCUS))
f_F4_drugMatch_DGIdb <- f_F4_drugMatch_DGIdb %>% arrange(desc(TWAS_FOCUS)) %>% arrange(desc(DEPICT_S2G)) %>% arrange(desc(DEPICT_and_FOCUS))
# 
write.table(f_F1_drugMatch_DRH, "F1_drugMatch_DEPICT_FOCUS_DRH.txt", sep = "\t", quote = F, col.names = T, row.names = F)
write.table(f_F1_drugMatch_DGIdb, "F1_drugMatch_DEPICT_FOCUS_DGIdb.txt", sep = "\t", quote = F, col.names = T, row.names = F)
write.table(f_F2_drugMatch_DRH, "F2_drugMatch_DEPICT_FOCUS_DRH.txt", sep = "\t", quote = F, col.names = T, row.names = F)
write.table(f_F2_drugMatch_DGIdb, "F2_drugMatch_DEPICT_FOCUS_DGIdb.txt", sep = "\t", quote = F, col.names = T, row.names = F)
write.table(f_F3_drugMatch_DRH, "F3_drugMatch_DEPICT_FOCUS_DRH.txt", sep = "\t", quote = F, col.names = T, row.names = F)
write.table(f_F3_drugMatch_DGIdb, "F3_drugMatch_DEPICT_FOCUS_DGIdb.txt", sep = "\t", quote = F, col.names = T, row.names = F)
write.table(f_F4_drugMatch_DRH, "F4_drugMatch_DEPICT_FOCUS_DRH.txt", sep = "\t", quote = F, col.names = T, row.names = F)
write.table(f_F4_drugMatch_DGIdb, "F4_drugMatch_DEPICT_FOCUS_DGIdb.txt", sep = "\t", quote = F, col.names = T, row.names = F)

print("NovelGenes Search Drug Databases")
NG_crossRef <- as.data.frame(NovelGenes$EnsemblGene)
colnames(NG_crossRef) <- c("Ensembl_map")
nrow(NG_crossRef)
NG_crossRef <- merge(NG_crossRef, HGNC_Ensembl, by = "Ensembl_map")
nrow(NG_crossRef) #44
NG_drugMatch <- search_Drug_Databases(geneVector = NG_crossRef$HGNC_map)
NG_drugMatch_DRH_output_list <- NG_drugMatch$DRH_output_list
NG_drugMatch_DRH <- bind_rows(NG_drugMatch_DRH_output_list, .id = "column_label")
NG_drugMatch_DGIdb_output_list <- NG_drugMatch$DGIdb_output_list
NG_drugMatch_DGIdb <- bind_rows(NG_drugMatch_DGIdb_output_list, .id = "column_label")

# # incororate the upset plot, showing DEPICT and FOCUS hits
m_F4_upsetMatrix <- as.data.frame(F4_upsetMatrix)
m_F4_upsetMatrix$DEPICT_and_FOCUS <- 0; m_F4_upsetMatrix$DEPICT_and_FOCUS[which(m_F4_upsetMatrix$DEPICT_S2G == 1 & m_F4_upsetMatrix$TWAS_FOCUS == 1)] <- 1
m_F4_upsetMatrix$Ensembl_map <- rownames(m_F4_upsetMatrix)
m_F4_upsetMatrix <- left_join(m_F4_upsetMatrix, HGNC_Ensembl, by = "Ensembl_map")
f_NG_drugMatch_DRH <- left_join(NG_drugMatch_DRH, m_F4_upsetMatrix, by = c("column_label" = "HGNC_map"))
f_NG_drugMatch_DGIdb <- left_join(NG_drugMatch_DGIdb, m_F4_upsetMatrix, by = c("column_label" = "HGNC_map"))
f_NG_drugMatch_DRH <- f_NG_drugMatch_DRH %>% arrange(desc(TWAS_FOCUS)) %>% arrange(desc(DEPICT_S2G)) %>% arrange(desc(DEPICT_and_FOCUS))
f_NG_drugMatch_DGIdb <- f_NG_drugMatch_DGIdb %>% arrange(desc(TWAS_FOCUS)) %>% arrange(desc(DEPICT_S2G)) %>% arrange(desc(DEPICT_and_FOCUS))
write.table(f_NG_drugMatch_DRH, "The_45_nonBMI_Genes_DEPICT_FOCUS_DRH.txt", sep = "\t", quote = F, col.names = T, row.names = F)
write.table(f_NG_drugMatch_DGIdb, "The_45_nonBMI_Genes_DEPICT_FOCUS_DGIdb.txt", sep = "\t", quote = F, col.names = T, row.names = F)

f_F1_drugMatch_DRH <- fread("F1_drugMatch_DEPICT_FOCUS_DRH.txt", header = T)
f_F1_drugMatch_DGIdb <- fread("F1_drugMatch_DEPICT_FOCUS_DGIdb.txt", header = T)
f_F2_drugMatch_DRH <- fread("F2_drugMatch_DEPICT_FOCUS_DRH.txt", header = T)
f_F2_drugMatch_DGIdb <- fread("F2_drugMatch_DEPICT_FOCUS_DGIdb.txt", header = T)
f_F3_drugMatch_DRH <- fread("F3_drugMatch_DEPICT_FOCUS_DRH.txt", header = T)
f_F3_drugMatch_DGIdb <- fread("F3_drugMatch_DEPICT_FOCUS_DGIdb.txt", header = T)
f_F4_drugMatch_DRH <- fread("F4_drugMatch_DEPICT_FOCUS_DRH.txt", header = T)
f_F4_drugMatch_DGIdb <- fread("F4_drugMatch_DEPICT_FOCUS_DGIdb.txt", header = T)

# Plot the gene-drug network ############

# make different network plots for DEPICT identified genes, FOCUS identified genes, and both.
DEPICT_ONLY = F
FOCUS_ONLY = F

if(DEPICT_ONLY == T){
  d1 <- f_F4_drugMatch_DRH[which(f_F4_drugMatch_DRH$DEPICT_S2G==1), c("column_label", "pert_iname")]
  d1$weight <- 0.01
}else if(FOCUS_ONLY == T){
  d1 <- f_F4_drugMatch_DRH[which(f_F4_drugMatch_DRH$TWAS_FOCUS==1), c("column_label", "pert_iname")]
  d1$weight <- 0.01
}else{
  d1 <- f_F4_drugMatch_DRH[, c("column_label", "pert_iname")]
  d1$weight <- 0.01
}
colnames(d1) <- c("from", "to", "weight")
d1 <- mutate_all(d1, .funs=toupper)
d1$weight <- as.numeric(as.character(d1$weight))
d1$to <- gsub(" ", "-", d1$to)
d1 <- d1 %>% distinct()

hist(range01(log(f_F4_drugMatch_DGIdb$interaction_score)), breaks = 50)
hist(range01(log(f_F4_drugMatch_DGIdb$interaction_score[which(f_F4_drugMatch_DGIdb$DEPICT_S2G==1)])), breaks = 50)
hist(range01(log(f_F4_drugMatch_DGIdb$interaction_score[which(f_F4_drugMatch_DGIdb$TWAS_FOCUS==1)])), breaks = 50)
if(DEPICT_ONLY == T){
  d2 <- f_F4_drugMatch_DGIdb[which(f_F4_drugMatch_DGIdb$DEPICT_S2G==1),c("column_label", "drug_claim_name")]
  d2$weight <- 0.5
}else if(FOCUS_ONLY == T){
  d2 <- f_F4_drugMatch_DGIdb[which(f_F4_drugMatch_DGIdb$TWAS_FOCUS==1),c("column_label", "drug_claim_name")]
  d2$weight <- 0.5
}else{
  d2 <- f_F4_drugMatch_DGIdb[,c("column_label", "drug_claim_name")]
  d2$weight <- 0.5
}

colnames(d2) <- c("from", "to", "weight")
d2 <- mutate_all(d2, .funs=toupper)
d2$weight <- as.numeric(as.character(d2$weight))
d2$to <- gsub(" ", "-", d2$to)
d2 <- d2 %>% distinct()

# the drug-gene pairs that were in both DRH and DGIdb
overlap_DRH_and_DGIdb <- inner_join(d1[,c("from", "to")], d2[,c("from", "to")], by = c("from" = "from", "to" = "to"))
rowsToKeep <- c()
for(i in 1:nrow(d1)){
  curRow <- as.data.frame(d1[i,])
  isInBoth <- inner_join(curRow[,c("from", "to")], overlap_DRH_and_DGIdb[,c("from", "to")], by = c("from" = "from", "to" = "to"))
  if(nrow(isInBoth) > 0){
    print(curRow)
  }else{
    # let's use the DGIdb instead of the DRH, when the drug-gene pair is in both
    rowsToKeep <- c(rowsToKeep, i)
  }
}
d1 <- d1[rowsToKeep,]
for(i in 1:nrow(d2)){
  curRow <- as.data.frame(d2[i,])
  isInBoth <- inner_join(curRow[,c("from", "to")], overlap_DRH_and_DGIdb[,c("from", "to")], by = c("from" = "from", "to" = "to"))
  if(nrow(isInBoth) > 0){
    d2$weight[i] <- .99
    print(curRow)
  }
}

length(unique(d1$from))
summary(as.factor(d1$weight))
summary(as.factor(d2$weight))

edgeList <- rbind(d1, d2)
dim(edgeList)
edgeList$duped <- edgeList[,c("from", "to")] %>% duplicated()
sum(edgeList$duped)
edgeList$weight <- as.numeric(as.character(edgeList$weight))

### filter to only approved/launched drugs
edgeList_sub <- subset(edgeList, edgeList$to %in% c(DGIdb_drugs_launched, DRH_drugs_launched))
dim(edgeList); length(unique(edgeList$to)); length(unique(edgeList$from))
dim(edgeList_sub); length(unique(edgeList_sub$to)); length(unique(edgeList_sub$from))
edgeList <- edgeList_sub

# Highlight Drug-Vertices with ADE-Weight-Changing label (OnSIDES)
adverse_reactions <- fread("./releases/v2.0.0/20231113/adverse_reactions.csv.gz", header = T)
length(unique(adverse_reactions$pt_meddra_term))
# 4302 adverse reactions
length(unique(adverse_reactions$ingredients_names))
# 2020 unique drugs
hist(adverse_reactions$percent_labels)
adverse_reactions_thresh <- subset(adverse_reactions, percent_labels > 0.75)
# 127,135 to 94,879

unique_pt_meddra_term <- as.data.frame(unique(adverse_reactions$pt_meddra_term))
colnames(unique_pt_meddra_term) <- c("pt_meddra_term")
dim(unique_pt_meddra_term)
write.table(unique_pt_meddra_term, "OnSIDES_unique_pt_meddra_term.txt", sep = "\t", quote = F, col.names = T, row.names = F)
# View(unique_pt_meddra_term)

weight_related_ADEs <- c(
  "Obesity",
  "Central obesity",
  "Weight increased",
  "Weight decreased",
  "Weight fluctuation",
  "Abnormal loss of weight",
  "Abnormal weight gain",
  "Weight loss poor",
  "Decreased appetite",
  "Increased appetite",
  "Appetite disorder",
  "Hunger",
  "Early satiety",
  "Binge eating",
  "Sleep-related eating disorder",
  "Eating disorder"
)
weightADE_adverse_reactions_thresh <- subset(adverse_reactions_thresh, pt_meddra_term %in% weight_related_ADEs)
table(weightADE_adverse_reactions_thresh$pt_meddra_term)[order(table(weightADE_adverse_reactions_thresh$pt_meddra_term))]
grepVec <- toupper(weightADE_adverse_reactions_thresh$ingredients_names)
length(grepVec)
grepVec <- gsub(" ", "-", grepVec)
# 513

edgeList$WeightADE <- F
counter = 0
for(i in 1:length(edgeList$to)){
  curDrug <- edgeList$to[i]
  grepRows <- grep(pattern = paste0("^",curDrug,"$"), x = grepVec)
  if(length(grepRows) > 0){
    counter <- counter + 1
    edgeList$WeightADE[i] <- T
    # stop("got one!")
  }
}
counter

length(unique(edgeList$to)) #approved drugs
length(unique(edgeList$to[which(edgeList$WeightADE == T)])) #drugs without wADE
# summary(edgeList$WeightADE)

nrow(edgeList) #drug-gene pairs
summary(as.factor(edgeList$weight))
length(unique(edgeList$from)) #genes
length(unique(edgeList$to)) #drugs

### DEPICT
# “Of the 248 approved drugs with drug-gene interactions to DEPICT-identified genes,
# 82 had prior descriptions of weight-related effects in the OnSIDES database. 
# The bipartite network for F4 DEPICT-identified genes included 314 drug-gene pairs (179 identified in the DRH, 
# 101 identified in the DGIdb, and 34 identified by both), consisting of 75 genes and 248 approved drugs.

### FOCUS
# “Of the 360 approved drugs with drug-gene interactions to FOCUS-identified genes,
# 107 had prior descriptions of weight-related effects in the OnSIDES database.
# The bipartite network for F4 FOCUS-identified genes included 442 drug-gene pairs (288 identified in the DRH,
# 95 identified in the DGIdb, and 59 identified by both), consisting of 84 genes and 360 approved drugs.

### DEPICT & FOCUS
# “Of the 529 approved drugs with drug-gene interactions to DEPICT- and FOCUS-identified genes,
# 148 had prior descriptions of weight-related effects in the OnSIDES database. 
# The bipartite network for F4 included 733 drug-gene pairs (451 identified in the DRH, 
# 192 identified in the DGIdb, and 90 identified by both), consisting of 151 genes and 529 approved drugs.

WeightADE_subset <- subset(edgeList, WeightADE == T)
WeightADE_subset <- WeightADE_subset[,c("from", "to")]

if(DEPICT_ONLY == T){
  write.table(WeightADE_subset, "WeightADE_subset_DrugeGenePairs_DEPICT.txt", sep = "\t", quote = F, col.names = F, row.names = F)
}else if(FOCUS_ONLY == T){
  write.table(WeightADE_subset, "WeightADE_subset_DrugeGenePairs_FOCUS.txt", sep = "\t", quote = F, col.names = F, row.names = F)
}else{
  write.table(WeightADE_subset, "WeightADE_subset_DrugeGenePairs.txt", sep = "\t", quote = F, col.names = F, row.names = F)
}

#### for graphical illustration cushion the weights from the endpoints
hist(edgeList$weight, breaks = 50, xlim = c(0,1))
edgeList$weight <- range01(edgeList$weight)
edgeList$weight[which(edgeList$weight < 0.05)] <- .05
edgeList$weight[which(edgeList$weight > 0.95)] <- .95
hist(edgeList$weight, breaks = 50, xlim = c(0,1))

MEDI_C <- fread("./MEDI-C/MEDI_Combined.csv", header = T, sep = ",")
MEDI_C <- subset(MEDI_C, VOCABULARY == "ICD10CM")
# The MEDI high precision subset (MEDI-HPS) includes indications found within either RxNorm or at least two of the three other resources.
MEDI_C_HPS <- subset(MEDI_C, HIGH_PRECISION_SUBSET == 1)
# 38,378 high precision drug-indication pairs 
length(unique(MEDI_C_HPS$DRUG_DESCRIPTION))
rxnorm_product_to_ingredient <- fread("//releases/v2.0.0/20231113/rxnorm_product_to_ingredient.csv.gz", header = T)
ingredient_rx_cui_name_MAP <- rxnorm_product_to_ingredient[,c("ingredient_rx_cui","ingredient_name")]
ingredient_rx_cui_name_MAP <- mutate_all(ingredient_rx_cui_name_MAP, .funs=toupper)
ingredient_rx_cui_name_MAP$ingredient_name <- gsub(" ", "-", ingredient_rx_cui_name_MAP$ingredient_name)
ingredient_rx_cui_name_MAP <- ingredient_rx_cui_name_MAP %>% distinct()
MEDI_C_HPS$RXCUI_IN <- as.character(MEDI_C_HPS$RXCUI_IN)
ingredient_rx_cui_name_MAP$ingredient_rx_cui <- as.character(ingredient_rx_cui_name_MAP$ingredient_rx_cui)
ingred_MEDI_C_HPS <- inner_join(MEDI_C_HPS, ingredient_rx_cui_name_MAP, by = c("RXCUI_IN" = "ingredient_rx_cui"))
dim(MEDI_C_HPS) #20149
dim(ingred_MEDI_C_HPS) #19440
sub_ingred_MEDI_C_HPS <- subset(ingred_MEDI_C_HPS, ingredient_name %in% edgeList$to)
dim(sub_ingred_MEDI_C_HPS) #3948
sub_ingred_MEDI_C_HPS$CODE_sep <- sub_ingred_MEDI_C_HPS$CODE
sub_ingred_MEDI_C_HPS$ICD_Range_Desc <- NA
sub_ingred_MEDI_C_HPS <- sub_ingred_MEDI_C_HPS %>% separate(CODE_sep,into=c("ICDRange","ICDCode"), convert=TRUE, sep="\\.")
### use the phewas package to map icd codes to phecodes, and then use those phecode categories to classify indications
sub_ingred_MEDI_C_HPS$id <- 1:nrow(sub_ingred_MEDI_C_HPS)
icd_codes <- data.frame(id=1:nrow(sub_ingred_MEDI_C_HPS), vocabulary_id=rep("ICD10CM", nrow(sub_ingred_MEDI_C_HPS)), code=as.character(sub_ingred_MEDI_C_HPS$CODE))
phecodes <- mapCodesToPhecodes(icd_codes, make.distinct = T, rollup.map = NULL)
dim(sub_ingred_MEDI_C_HPS)
sub_ingred_MEDI_C_HPS <- inner_join(x = sub_ingred_MEDI_C_HPS, y = phecodes, by = c("id" = "id"))
dim(sub_ingred_MEDI_C_HPS)
phecodesInfo <- addPhecodeInfo(as.character(sub_ingred_MEDI_C_HPS$phecode), groupcolors = T)
dim(phecodesInfo)
phecodesInfo <- phecodesInfo %>% distinct()
dim(phecodesInfo)
dim(sub_ingred_MEDI_C_HPS)
sub_ingred_MEDI_C_HPS <- left_join(x = sub_ingred_MEDI_C_HPS, y = phecodesInfo, by = c("phecode" = "phenotype"))
dim(sub_ingred_MEDI_C_HPS)
cat(unique(sub_ingred_MEDI_C_HPS$group), sep = "\n")

# Coninue plotting ############
v1 <- unique(edgeList$from)
v2 <- unique(edgeList$to[edgeList$WeightADE == F])
v3 <- unique(edgeList$to[edgeList$WeightADE == T])
g <- graph.empty()
g <- add.vertices(g, nv=length(v1), attr=list(name=v1, indication=rep("NA",length(v1)), type=rep("Gene",length(v1))))
g <- add.vertices(g, nv=length(v2), attr=list(name=v2, indication=rep("NA",length(v2)), type=rep("Drug",length(v2))))
g <- add.vertices(g, nv=length(v3), attr=list(name=v3, indication=rep("NA",length(v3)), type=rep("DrugADE",length(v3))))

# print the drugs that aren't in the MEDI-C database
v2[which(v2 %ni% sub_ingred_MEDI_C_HPS$ingredient_name)]
v3[which(v3 %ni% sub_ingred_MEDI_C_HPS$ingredient_name)]

for(i in 1:length(V(g))){
  curDrug <- V(g)$name[i]
  cur_type <- vertex_attr(graph = g, name = "type", index = which(V(g)$name == curDrug))
  if(cur_type == "Drug" | cur_type == "DrugADE"){
    cur_MEDI_C_sub <- subset(sub_ingred_MEDI_C_HPS, ingredient_name == curDrug)
    if(nrow(cur_MEDI_C_sub) > 0){
      cur_top_indication <- names(table(cur_MEDI_C_sub$group)[order(table(cur_MEDI_C_sub$group),decreasing = T)])[1]
      vertex_attr(graph = g, name = "indication", index = which(V(g)$name == curDrug)) <- cur_top_indication
    }
  }
}
vertex_attr(graph = g, name = "indication", index = which(V(g)$name == "INSULIN")) <- "endocrine/metabolic"

V(g)$indication
table(as.factor(V(g)$indication))

edgeListVec <- as.vector(t(as.matrix(data.frame(S1=edgeList$from, S2=edgeList$to))))
g <- add.edges(g, edgeListVec, weight = edgeList$weight)
is.bipartite(g)
c_scale <- colorRamp(c("#B2A7F4", "#F5C249", "#ff1700"))

E(g)$color = apply(c_scale(E(g)$weight), 1, function(x) rgb(x[1]/255,x[2]/255,x[3]/255) )
E(g)$arrow.mode="-" 
E(g)$curved=F 
V(g)$label.cex <- igraph::degree(g) * 0.2
V(g)$size <- igraph::degree(g)
# split the graph into high-degree nodes and low-degree nodes for visualization
gHD <- g
gLD <- g
gHD <- delete.vertices(gHD, V(gHD)[igraph::degree(gHD)<=1])
gHD <- delete.vertices(gHD, V(gHD)[igraph::degree(gHD)<1])
# gLD <- delete.vertices(gLD, V(gLD)[igraph::degree(gLD)>4])
length(V(g))
length(V(gHD))
# length(V(gLD))
hist(igraph::degree(g), breaks = 50); table(igraph::degree(g))
# TRUE IS NUMBER OF GENES, FALSE IS NUMBER OF DRUGS
summary(V(g)$type)
hist(igraph::degree(gHD), breaks = 50); table(igraph::degree(gHD))
summary(V(gHD)$type)
# hist(degree(gLD), breaks = 50); table(degree(gLD))
length(E(g))
length(E(gHD))

colcols <- c("grey60", "ivory2", "hotpink", "tomato3", "lightseagreen", "#ffd8b1", "#ffe119", "slateblue1", 
             "#f58231", "#fabed4", "#e6194B", "steelblue3", "springgreen2", "forestgreen", 
             "darkgoldenrod1", "#f032e6", "darkolivegreen2", "#42d4f4", "#9A6324")

collabs <- c("gene",
             "indication missing",
             "neoplasms",
             "neurological",
             "hematopoietic",
             "endocrine/metabolic",
             "respiratory",
             "digestive",
             "dermatologic",
             "genitourinary",
             "pregnancy complications",
             "injuries & poisonings",
             "circulatory system",
             "infectious diseases",
             "sense organs",
             "musculoskeletal",
             "symptoms",
             "mental disorders",
             "congenital anomalies")

g_fill_vec <- V(g)$indication
g_fill_vec[which(V(g)$type == "Gene")] <- "grey70"
g_fill_vec <- gsub("NA", "ivory2", g_fill_vec)
g_fill_vec <- gsub("neoplasms", "hotpink", g_fill_vec)
g_fill_vec <- gsub("neurological", "tomato3", g_fill_vec)
g_fill_vec <- gsub("hematopoietic", "lightseagreen", g_fill_vec)
g_fill_vec <- gsub("endocrine/metabolic", "#ffd8b1", g_fill_vec)
g_fill_vec <- gsub("respiratory", "#ffe119", g_fill_vec)
g_fill_vec <- gsub("digestive", "slateblue1", g_fill_vec)
g_fill_vec <- gsub("dermatologic", "#f58231", g_fill_vec)
g_fill_vec <- gsub("genitourinary", "#fabed4", g_fill_vec)
g_fill_vec <- gsub("pregnancy complications", "#e6194B", g_fill_vec)
g_fill_vec <- gsub("injuries & poisonings", "steelblue3", g_fill_vec)
g_fill_vec <- gsub("circulatory system", "springgreen2", g_fill_vec)
g_fill_vec <- gsub("infectious diseases", "forestgreen", g_fill_vec)
g_fill_vec <- gsub("sense organs", "darkgoldenrod1", g_fill_vec)
g_fill_vec <- gsub("musculoskeletal", "#f032e6", g_fill_vec)
g_fill_vec <- gsub("symptoms", "darkolivegreen2", g_fill_vec)
g_fill_vec <- gsub("mental disorders", "skyblue", g_fill_vec)
g_fill_vec <- gsub("congenital anomalies", "#9A6324", g_fill_vec)

gHD_color_vec <- V(gHD)$type
gHD_color_vec <- gsub("Gene", "grey70", gHD_color_vec)
gHD_color_vec <- gsub("DrugADE", "cyan2", gHD_color_vec)
gHD_color_vec <- gsub("Drug", "chartreuse1", gHD_color_vec)

g_color_vec <- V(g)$type
g_color_vec <- gsub("Gene", "grey70", g_color_vec)
g_color_vec <- gsub("DrugADE", "cyan2", g_color_vec)
g_color_vec <- gsub("Drug", "chartreuse1", g_color_vec)

g_outline_vec <- V(g)$type
g_outline_vec <- gsub("Gene", "white", g_outline_vec)
g_outline_vec <- gsub("DrugADE", "black", g_outline_vec)
g_outline_vec <- gsub("Drug", "white", g_outline_vec)

collabs2 <- c("Gene", "Drug", "Drug with wADE")
colcols2 <- c("grey70", "chartreuse1", "cyan2")

ghostWeight <- log(1.0010)
groupingWeightClose <- log(1.004)
groupingWeightFar <- log(1.009)

# create a layout for the network that brings drugs with similar indications closer together spatially
set.seed(12)
G_Grouped <- g
# all other edges get really small weight
G_Grouped = add_edges(G_Grouped, combn(V(g), 2), attr=list(weight=ghostWeight))
# real edges get weight of 1
E(G_Grouped)$weight[which(E(g)$weight == 0.05)] <- log(10)
E(G_Grouped)$weight[which(E(g)$weight == 0.50)] <- log(20)
E(G_Grouped)$weight[which(E(g)$weight == 0.95)] <- log(30)
## Add small weight to edges in the same group
summary(as.factor(V(G_Grouped)$indication))[order(summary(as.factor(V(g)$indication)))]
gcomps <- components(g)$membership
for(i in unique(V(G_Grouped)$indication)) {
  # if(i %ni% c("neoplasms", "circulatory system", "mental disorders", "respiratory", "neurological", "genitourinary", "dermatologic", "endocrine/metabolic")){
  print(i)
  GroupV <- which(V(G_Grouped)$indication == i)
  print(length(GroupV))
  if(length(GroupV) > 1 & i %ni% c("NA")){
    vertexPairs <- combn(GroupV, 2)
    connected_weights <- rep(F, ncol(vertexPairs))
    for(j in 1:ncol(vertexPairs)){
      if(gcomps[[V(g)[vertexPairs[1,j]]]] == gcomps[[V(g)[vertexPairs[2,j]]]]){
        connected_weights[j] <- T
      }else{
      }
    }
    G_Grouped <- add_edges(G_Grouped, vertexPairs[,which(connected_weights == T)], attr=list(weight=groupingWeightClose))
    G_Grouped <- add_edges(G_Grouped, vertexPairs[,which(connected_weights == F)], attr=list(weight=groupingWeightFar))
  }
} 
hist(E(G_Grouped)$weight, breaks = 50)
# hist(E(G_Grouped)$weight[which(E(G_Grouped)$weight >= groupingWeightClose & E(G_Grouped)$weight <= groupingWeightFar)], breaks = 50)
hist(E(G_Grouped)$weight[which(E(G_Grouped)$weight >= groupingWeightClose)], breaks = 50)
# LO = layout_with_fr(G_Grouped)
# LO = layout_with_kk(g)
# LO = layout_with_drl(G_Grouped, weights = E(G_Grouped)$weight, options=list(edge.cut = 0, simmer.attraction = 0))
LO = layout_with_fr(G_Grouped)

if(DEPICT_ONLY == T){
  set.seed(12)
  pdf("./Drug_Gene_Bipartite_Network_DEPICT.pdf", width = 12, height = 12)
  par(mar=c(0,0.25,0,2)) #mai gives margin of bottom left top right in inches
  plot(g, layout=LO, vertex.label.cex=.32, vertex.size = 2.7, vertex.frame.color = g_outline_vec, vertex.frame.width = 1, vertex.color=g_color_vec, vertex.label.color = "grey20", asp=0, margin=0)
  legend("topright", collabs2, pch=21, col="white", pt.bg=colcols2, pt.cex=1.15, cex=1, bty="n", ncol=1)
  dev.off()
  set.seed(12)
  pdf("./Drug_Gene_Bipartite_Network_DEPICT_Indications.pdf", width = 12, height = 12)
  par(mar=c(0,0.25,0,2)) #mai gives margin of bottom left top right in inches
  plot(g, layout=LO, vertex.label.cex=.32, vertex.size = 2.7, vertex.frame.color = g_outline_vec, vertex.frame.width = 1, vertex.color=g_fill_vec, vertex.label.color = "grey20", asp=0, margin=0)
  legend("topright", collabs, pch=21, col="white", pt.bg=colcols, pt.cex=0.9, cex=.6, bty="n", ncol=1)
  dev.off()
}else if(FOCUS_ONLY == T){
  set.seed(12)
  pdf("./Drug_Gene_Bipartite_Network_FOCUS.pdf", width = 14, height = 14)
  par(mar=c(0,0.25,0,2)) #mai gives margin of bottom left top right in inches
  plot(g, layout=LO, vertex.label.cex=.32, vertex.size = 2.4, vertex.frame.color = g_outline_vec, vertex.frame.width = 1, vertex.color=g_color_vec, vertex.label.color = "grey20", asp=0, margin=0)
  legend("topright", collabs2, pch=21, col="white", pt.bg=colcols2, pt.cex=1.5, cex=1, bty="n", ncol=1)
  dev.off()
  set.seed(12)
  pdf("./Drug_Gene_Bipartite_Network_FOCUS_Indications.pdf", width = 14, height = 14)
  par(mar=c(0,0.25,0,2)) #mai gives margin of bottom left top right in inches
  plot(g, layout=LO, vertex.label.cex=.32, vertex.size = 2.4, vertex.frame.color = g_outline_vec, vertex.frame.width = 1, vertex.color=g_fill_vec, vertex.label.color = "grey20", asp=0, margin=0)
  legend("topright", collabs, pch=21, col="white", pt.bg=colcols, pt.cex=0.9, cex=.6, bty="n", ncol=1)
  dev.off()
}else{
  set.seed(12)
  # LO = layout_with_fr(g, weights = rep(1,length(E(g)$weight)))
  # LO = layout_with_fr(g)
  pdf("./Drug_Gene_Bipartite_Network.pdf", width = 14, height = 14)
  par(mar=c(0,0.25,0,2)) #mai gives margin of bottom left top right in inches
  plot(g, layout=LO, vertex.label.cex=.3, vertex.size = 2.2, vertex.frame.color = g_outline_vec, vertex.frame.width = 1, vertex.color=g_color_vec, vertex.label.color = "grey20", asp=0, margin=0)
  legend("topright", collabs2, pch=21, col="white", pt.bg=colcols2, pt.cex=1.5, cex=1, bty="n", ncol=1)
  dev.off()
  set.seed(12)
  pdf("./Drug_Gene_Bipartite_Network_Indications.pdf", width = 14, height = 14)
  par(mar=c(0,0.25,0,2)) #mai gives margin of bottom left top right in inches
  plot(g, layout=LO, vertex.label.cex=.3, vertex.size = 2.2, vertex.frame.color = g_outline_vec, vertex.frame.width = 1, vertex.color=g_fill_vec, vertex.label.color = "grey20", asp=0, margin=0)
  legend("topright", collabs, pch=21, col="white", pt.bg=colcols, pt.cex=0.9, cex=.6, bty="n", ncol=1)
  dev.off()
}

View(edgeList)

