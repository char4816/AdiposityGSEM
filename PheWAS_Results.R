setwd("working/directory")

library(PheWAS)
library(data.table)
library(RColorBrewer)
library(wpa)
'%ni%' <- Negate('%in%')
options(ggrepel.max.overlaps = Inf)

# PheWAS is done by SPAtest, and corrects for age, sex, batch, and top 10 PCs. 
# Any phenotype with case # <10 is censored and not run.
# It is done in unrelated EUR in freeze2. 
# The association effect size reflects that of 1 s.d. of the score.

# *_EUR_phewas_qPRS_bPRS.txt
# “phecode”, “OR_qPRS”, “pval_qPRS”. 
# “bPRS” means stratifying a score by top 10% vs. the rest

# *_phewas.svg
# A quick pheWAS plot for visualization

# *_prs.log
# This records the information of applying your weight files to the freeze2 genotype. It has the information of how many SNPs were included vs. tossed out if they were not found etc.

# In the plot I used FDR<10%, because a lot of phecode are related to each other 
# (like “breast cancer” vs. “female breast cancer”).
# or replot anything in the way you like, mine was just to get a quick visualization

fnames <- c("F1_GSEM_Birth_Size_sumstats_noQsnp_hg38phewasFormat_WEIGHTED_EUR_phewas_qPRS_bPRS.txt",
            "F2_GSEM_Abdominal_Size_sumstats_noQsnp_hg38phewasFormat_WEIGHTED_EUR_phewas_qPRS_bPRS.txt",
            "F3_GSEM_Body_Size_and_Adipose_Distribution_sumstats_noQsnp_hg38phewasFormat_WEIGHTED_EUR_phewas_qPRS_bPRS.txt",
            "F4_GSEM_Adiposity_sumstats_noQsnp_hg38phewasFormat_WEIGHTED_EUR_phewas_qPRS_bPRS.txt",
            "F4_GSEM_Adiposity_sumstats_noQsnp_hg38phewasFormat_WEIGHTED_bmi_adjusted_EUR_phewas_qPRS_bPRS.txt"
)

for(ff in 1:length(fnames)){
  fname <- fnames[ff]
  print(fname)
  # fname <- "F1_GSEM_Birth_Size_sumstats_noQsnp_hg38phewasFormat_WEIGHTED_EUR_phewas_qPRS_bPRS.txt"
  # fname <- "F2_GSEM_Abdominal_Size_sumstats_noQsnp_hg38phewasFormat_WEIGHTED_EUR_phewas_qPRS_bPRS.txt"
  # fname <- "F3_GSEM_Body_Size_and_Adipose_Distribution_sumstats_noQsnp_hg38phewasFormat_WEIGHTED_EUR_phewas_qPRS_bPRS.txt"
  # fname <- "F4_GSEM_Adiposity_sumstats_noQsnp_hg38phewasFormat_WEIGHTED_EUR_phewas_qPRS_bPRS.txt"
  # fname <- "F4_GSEM_Adiposity_sumstats_noQsnp_hg38phewasFormat_WEIGHTED_bmi_adjusted_EUR_phewas_qPRS_bPRS.txt"
  df <- fread(fname, header = T)
  df <- subset(df, is.na(df$pval_qPRS) == F)
  
  df$phecode <- as.character(df$phecode)
  head(df)
  dim(df); length(unique(df$phecode))
  
  dfInfo <- addPhecodeInfo(df$phecode)
  mdf <- merge(df, dfInfo, by.x = "phecode", by.y = "phenotype")
  
  mdf$phenotype = mdf$phecode
  mdf$p = mdf$pval_qPRS
  mdf$OR = mdf$OR_qPRS
  mdf <- as.data.frame(mdf)
  colnames(mdf)
  
  ncoul <- length(unique(mdf$group))
  if(ncoul <= 8){
    coul <- brewer.pal(ncoul, "Accent")
  }else if(ncoul <= 20){
    coul <- c(brewer.pal(8, "Accent"), brewer.pal(ncoul - 8, "Paired"))
  }else{
    stop("phecode groups > 20")
  }
  set.seed(12)
  coul <- sample(coul, size = ncoul, replace = F)
  coul[which(coul == "#FFFF99")] <- "gold"
  coul[which(coul == "#FDBF6F")] <- "cyan"
  coul[which(coul == "#FDC086")] <- "darkred"
  coul[which(coul == "#386CB0")] <- "darkblue"
  
  mdf$groupnum <- as.numeric(as.factor(mdf$group))
  mdf$color <- coul[mdf$groupnum]
  mdf$anno <- F
  pthresh <- 0.10/nrow(mdf)
  # If there's more than 50 annotations, prioritize the top 5 per group
  nprioritize <- 6
  if(length(which(mdf$p < pthresh)) > 50){
    for(g in unique(mdf$group)){
      indices <- which(mdf$group == g & mdf$p < pthresh)
      ovec <- order(mdf$pval_qPRS[indices], decreasing = F)
      # plot(mdf$pval_qPRS[indices[ovec]])
      if(length(ovec) > nprioritize){
        # mdf$anno[indices[ovec]][1:nprioritize] <- mdf$description[indices[ovec]][1:nprioritize] 
        mdf$anno[indices[ovec]][1:nprioritize] <- T
      }else{
        # mdf$anno[indices[ovec]] <- mdf$description[indices[ovec]]
        mdf$anno[indices[ovec]] <- T
      }
    }
  }else{
    mdf$anno[which(mdf$p < pthresh)] <- T
  }
  
  plotTitle <- us_to_space(unlist(strsplit(unlist(strsplit(fname, split = "\\."))[1], split = "_sumstats_"))[1])
  if(ff == 5){
    plotTitle <- paste0(plotTitle, " BMI-Adjusted")
  }
  
  mdf$description[which(mdf$phecode == 965.1)] <- "Opiates and related narcotics causing adverse effects"
  mdf$description[which(mdf$phecode == 457)] <- "Encounter for long-term (current) use of anticoagulants"
  
  
  writeoutdf <- mdf[order(mdf$p, decreasing = F),c("phecode", "description", "group", "OR_qPRS", "SPA_SEbeta_qPRS", "pval_qPRS", "ncase", "nctrl", "notes_qPRS")]     
  write.table(writeoutdf, paste0(gsub(".txt","",fname),"forSupp.txt"), col.names = T, row.names = F, quote = F, sep = "\t")
  
  phewas_plot <- NULL
  phewas_plot <- phenotypeManhattan(mdf, 
                                    OR.direction = T, 
                                    title=plotTitle, 
                                    annotate.size=3, 
                                    annotate.angle = 0,
                                    use.color = T,
                                    suggestive.line = NA,
                                    significant.line = NA,
                                    # annotate.level = 0.05/nrow(mdf),
                                    annotate.list = mdf$phenotype[mdf$anno],
                                    annotate.phenotype.description = mdf[,c("phenotype", "description")]
  )
  phewas_plot <- phewas_plot + theme(text=element_text(family="Helvetica")) +
    theme(axis.text.x=element_text(angle = 0, hjust = 0.5)) +
    theme(axis.text.y=element_text(family="Helvetica")) +
    geom_hline(yintercept = -log10(pthresh), color = "red", linetype = "dashed", alpha = 0.5) +
    # scale_y_continuous(limits = c(0, ceiling(max(-log10(mdf$p))/10)*10), breaks = round(seq(0, ceiling(max(-log10(mdf$p))/10)*10, length.out=7))) +
    scale_y_continuous() +
    ylab("-log10(p)") +
    coord_flip()
  # phewas_plot
  
  closeAllConnections()
  png(paste0(unlist(strsplit(fname, split = "\\."))[1],".png"), width = 10, height = 16, unit = "in", res = 300)
  print(phewas_plot)
  dev.off()
  
}


gsub("[^[:alnum:]]+","_",unlist(strsplit(unlist(strsplit(fname, split = "\\."))[1], split = "_sumstats_"))[1])
sub(" ", "_", unlist(strsplit(unlist(strsplit(fname, split = "\\."))[1], split = "_sumstats_"))[1])




fname <- "F4_GSEM_Adiposity_sumstats_noQsnp_hg38phewasFormat_WEIGHTED_EUR_phewas_qPRS_bPRS.txt"
df <- fread(fname, header = T)
df <- subset(df, is.na(df$pval_qPRS) == F)
df$phecode <- as.character(df$phecode)
head(df)
dim(df); length(unique(df$phecode))
dfInfo <- addPhecodeInfo(df$phecode)
mdf <- merge(df, dfInfo, by.x = "phecode", by.y = "phenotype")
mdf$phenotype = mdf$phecode
mdf$p = mdf$pval_qPRS
mdf$OR = mdf$OR_qPRS
mdf <- as.data.frame(mdf)
colnames(mdf)

fname2 <- "F4_GSEM_Adiposity_sumstats_noQsnp_hg38phewasFormat_WEIGHTED_bmi_adjusted_EUR_phewas_qPRS_bPRS.txt"
df2 <- fread(fname2, header = T)
df2 <- subset(df2, is.na(df2$pval_qPRS) == F)
df2$phecode <- as.character(df2$phecode)
head(df2)
dim(df2); length(unique(df2$phecode))
df2Info <- addPhecodeInfo(df2$phecode)
mdf2 <- merge(df2, df2Info, by.x = "phecode", by.y = "phenotype")
mdf2$phenotype = mdf2$phecode
mdf2$p = mdf2$pval_qPRS
mdf2$OR = mdf2$OR_qPRS
mdf2 <- as.data.frame(mdf2)
colnames(mdf2)

merge_pre_post_adjBMI <- merge(mdf[,c("description","group","phenotype","OR","p")], mdf2[,c("description","group","phenotype","OR","p")], by = c("description","group","phenotype"))

merge_pre_post_adjBMI$beta.x <- log(merge_pre_post_adjBMI$OR.x)
merge_pre_post_adjBMI$beta.y <- log(merge_pre_post_adjBMI$OR.y)

p1 <- ggplot(merge_pre_post_adjBMI, aes(x=OR.x, y=OR.y)) + 
  geom_point(color = "#00AFBB", size = 0.75) + 
  stat_smooth(method = "lm", formula = y ~ x, geom = "smooth", color = "darkblue") +
  theme_classic() +
  geom_abline (slope=1, linetype = "dashed", color="black") +
  xlab("OR without adjusting for BMI") +
  ylab("OR with adjusting for BMI")

png("pheWAS_adjBMI_compare_OR.png", width = 3.5, height = 3.5, unit = "in", res = 300)
print(p1)
dev.off()

p2 <- ggplot(merge_pre_post_adjBMI, aes(x=beta.x, y=beta.y)) + 
  geom_point(color = "#00AFBB", size = 0.75) + 
  stat_smooth(method = "lm", formula = y ~ x, geom = "smooth", color = "darkblue") +
  theme_classic() +
  geom_abline (slope=1, linetype = "dashed", color="black") +
  xlab("Beta without adjusting for BMI") +
  ylab("Beta with adjusting for BMI")

png("pheWAS_adjBMI_compare_Beta.png", width = 3.5, height = 3.5, unit = "in", res = 300)
print(p2)
dev.off()

merge_pre_post_adjBMI$neglog10p.x <- -log10(merge_pre_post_adjBMI$p.x)
merge_pre_post_adjBMI$neglog10p.y <- -log10(merge_pre_post_adjBMI$p.y)



p3 <- ggplot(merge_pre_post_adjBMI, aes(x=neglog10p.x, y=neglog10p.y)) + 
  geom_point(color = "salmon", size = 0.75) + 
  stat_smooth(method = "lm", formula = y ~ x, geom = "smooth", color = "orange") +
  theme_classic() +
  geom_abline (slope=1, linetype = "dashed", color="black") +
  xlab("-log10(p) without adjusting for BMI") +
  ylab("-log10(p) with adjusting for BMI")

png("pheWAS_adjBMI_compare_p_all.png", width = 3.5, height = 3.5, unit = "in", res = 300)
print(p3)
dev.off()

p4 <- ggplot(merge_pre_post_adjBMI, aes(x=neglog10p.x, y=neglog10p.y)) + 
  geom_point(color = "salmon", size = 0.75) + 
  stat_smooth(method = "lm", formula = y ~ x, geom = "smooth", color = "orange") +
  theme_classic() +
  geom_abline (slope=1, linetype = "dashed", color="black") +
  xlab("-log10(p) without adjusting for BMI") +
  ylab("-log10(p) wtih adjusting for BMI") +
  xlim(c(0,30)) + ylim(c(0,30))

png("pheWAS_adjBMI_compare_p_crop.png", width = 3.5, height = 3.5, unit = "in", res = 300)
print(p4)
dev.off()




