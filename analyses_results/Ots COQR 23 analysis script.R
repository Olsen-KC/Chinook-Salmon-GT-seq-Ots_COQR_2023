library(magrittr)
library(tidyverse)
library(adegenet)
library(ggplot2)
library(cowplot)
library(readxl)
library(Rmisc)
library(dplyr)






##load QC filtered genos

load("C:/Users/olsenk2/OneDrive - Oregon State University/Desktop/OtsAC22_23COQR/Ots COQR 22_23 genos_2.0.R")

#look at relatedness in filtered dataset for identifying duplicates##
##format vcf for coancestry##


replacedGT<-as.data.frame(genos_2.0[, c(9:(ncol(genos_2.0)-1))])

delim_gtypes <- replacedGT[, ] %>%
  mutate_all(funs(case_when(. == "AA" ~ "A/A",
                            . == "AT" ~ "A/T", 
                            . == "AC" ~ "A/C", 
                            . == "AG" ~ "A/G", 
                            . == "TT" ~ "T/T", 
                            . == "TA" ~ "T/A", 
                            . == "TC" ~ "T/C", 
                            . == "TG" ~ "T/G",
                            . == "CC" ~ "C/C", 
                            . == "CA" ~ "C/A", 
                            . == "CT" ~ "C/T", 
                            . == "CG" ~ "C/G", 
                            . == "GG" ~ "G/G", 
                            . == "GA" ~ "G/A", 
                            . == "GT" ~ "G/T", 
                            . == "GC" ~ "G/C", 
                            . == "--" ~ "-/-", 
                            . == "A-" ~ "A/-", 
                            . == "T-" ~ "T/-", 
                            . == "C-" ~ "C/-", 
                            . == "G-" ~ "G/-", 
                            . == "-A" ~ "-/A",
                            . == "-T" ~ "-/T", 
                            . == "-C" ~ "-/C", 
                            . == "-G" ~ "-/G", 
                            . == "00" ~ "0/0", 
                            . == "0" ~ "0/0"
  )))



library(splitstackshape)


vrtnames<-c(colnames(delim_gtypes))

##split genotypes into 2 alleles

dfspltgeno<-as.data.frame(cSplit(delim_gtypes, splitCols = vrtnames, sep = "/", type.convert = FALSE))


coancestrylbls<-data.frame(c(genos_2.0$sample_simple))

colnames(coancestrylbls)<-c("Ind")

coancestrylbls %<>%
  mutate(Ind = substr(Ind, 4, 16))


rownames(dfspltgeno)<-as.character(coancestrylbls$Ind)

##replace nucleotides with numeric values for coancestry## 
##A = 1; T = 2; C = 3; G = 4; - = 5 ###

dfspltgeno %<>%
  mutate_all(funs(case_when(. == "A" ~ 1,
                            . == "T" ~ 2,
                            . == "C" ~ 3,
                            . == "G" ~ 4, 
                            . == "-" ~ 5,
                            . == "0" ~ 0, 
  )))



write.table(dfspltgeno, file = "C:/Users/olsenk2/OneDrive - Oregon State University/Desktop/OtsAC22_23COQR/OtsAC2223COQR coancestry input.txt", quote = FALSE, sep = "\t", col.names = FALSE, row.names = TRUE)


##read in relatedness estimates##


rel<-read.table("C:/ZSL/Coancestry/OtsAC2223COQR/RelatednessEstimates.txt", sep = ",")

##grab wang estimator##

rel<-rel[, c(1:3, 6)]

colnames(rel)<-c("dyad", "Ind1", "Ind2", "wang")




##read in COQR 2023 meta##

meta<-read.csv("C:/Users/olsenk2/OneDrive - Oregon State University/Desktop/OtsAC22_23COQR/Ots23COQR_meta.csv", header = TRUE)

##subset dataset to just 2023 samples
coquille_23<-genos_2.0[genos_2.0$sample_simple %in% meta$IndividualName, ]

##subset meta to QC filtered samples#

meta<-meta[meta$IndividualName %in% coquille_23$sample_simple, ]

##put dataset and meta in same order##

meta<-meta[order(match(meta$IndividualName, coquille_23$sample_simple)), ]

identical(coquille_23$sample_simple, meta$IndividualName)

##read marker information for Ots GT-seq panel

ots_panel_info<-read_excel("C:/Users/olsenk2/OneDrive - Oregon State University/Desktop/OtsAC22_23COQR/Ots GT-seq panel/Ots GT-seq panel info.xlsx", sheet = 2)

##subset panel info to those markers in the SFGL panel

ots_SFGL_panel_info<-ots_panel_info[ots_panel_info$Assay %in% colnames(genos_2.0), ]


##PCA##

##look at just the 31 Ots run timing markers that made it through filtering 
##2 of the run timing markers (Ots28_11143508 & Ots37124-12281207) were removed during filtering
##1 run timing marker (Ots37124-12279142) was not in the panel
##1 run timing maker (Ots37124-12270118) may be problematic drop it


##subset to run timing markers#
rt_markers<-subset(ots_SFGL_panel_info, ots_SFGL_panel_info$PresumedType == "runtiming")

run_timing<-coquille_23[, c(1:8, which(colnames(coquille_23) %in% rt_markers$Assay))]

#drop Ots37124-12270118

run_timing %<>%
  select(!'Ots37124-12270118')

##put markers in chromosomal order##

run_timing<-run_timing[, c(1:12, 36:37, 13:16, 38, 17:23, 39, 24:35)]


##write run timing markers to csv for polarization

write.csv(run_timing, quote = FALSE, row.names = FALSE, file = "C:/Users/olsenk2/OneDrive - Oregon State University/Desktop/OtsAC22_23COQR/OtsAC23COQR_runtiming_genotypes.csv")


rt_matrix <- as.matrix(run_timing[,c(9:(ncol(run_timing)))]) 
row.names(rt_matrix) <- run_timing$sample_simple
genind_run_timing <- df2genind(rt_matrix, sep ="", ploidy=2,NA.char = "00")

##In the run timing markers, calculate proportional heterozygosity 

hets<-as.data.frame(apply(rt_matrix, MARGIN = 1, function(x){ sum( x == "AT" | x == "AC" | x == "AG" | x == "TA" | x == "TC" | x == "TG" | x == "CA" | x == "CT" | x == "CG" | x == "GA" | x == "GT" | x == "GC") } ))

colnames(hets)<-"het loci"

hets$prop_het<-hets$`het loci` / 31



##Individual level PCA###

X <- scaleGen(genind_run_timing,  NA.method="mean")

#then run pca, keep all PCs
pca1 <- dudi.pca(X, scale = FALSE, scannf = FALSE, nf = 300)

snp_pcs <- pca1$li

##calculate percentage of total variation explained by each eigenvalue#

eig.perc <- 100*pca1$eig/sum(pca1$eig)

#create column of sample IDs from rownames##
snp_pcs %<>%
  rownames_to_column(var="sample") 

##check both in same order##

identical(snp_pcs$sample, meta$IndividualName)

##add meta data back in##

snp_pcs<-cbind(snp_pcs, meta)

snp_pcs %<>%
  relocate(sample, IndividualName, WILDorHAT., Marks, Gender)

##add proportion heterozygous data

identical(snp_pcs$sample, rownames(hets))

snp_pcs<-cbind(snp_pcs, hets)

##set order of NOR, HOR, and unknown (UK) origins


snp_pcs$WILDorHAT.<-factor(snp_pcs$WILDorHAT., levels = c('W', 'H', 'UK'))


##plot PCA with HOR / NOR comparison##
jitter <- position_jitter(width = 0.5, height = 0.5)
rt_PCA<-ggplot(snp_pcs)+
  geom_point(aes(Axis1, Axis2, color = WILDorHAT.), position = jitter, alpha = 0.5, size = 5.0) +
  scale_color_manual(values = c("#404080", "#69b3a2", "orange"),  
                     name = "Origin:", labels = c("NOR", "HOR", "UK")) +
  xlab("PC1 4.7%")+
  ylab("PC2 1.4%")+
  theme_classic() +
  theme(text = element_text(family = "sans"), 
        legend.title = element_text(size = 22), 
        legend.text = element_text(size = 22),
        legend.position = 'top',
        legend.justification = 'left',
        axis.title = element_text(color = "black", size = 22),
        axis.ticks = element_line(size = 2, color = "black"), 
        axis.ticks.length = unit(0.35, "cm"), 
        axis.line = element_line(size = 1.75, color = "black"), 
        axis.text = element_text(color = "black", size = 20.0))

##plot PCA with heterozygosity 

het_PCA<-ggplot(snp_pcs)+
  geom_point(aes(Axis1, Axis2, color = prop_het), position = jitter, alpha = 0.5, size = 5.0) +
  scale_color_gradient(name = "Proportion of heterozygous adult migration timing markers       ", limits = c(0, 1), breaks = c(0, 0.5, 1), low = "darkblue", high = "red") +
  xlab("PC1 4.7%")+
  ylab("PC2 1.4%")+
  theme_classic() +
  theme(text = element_text(family = "sans"), 
        legend.title = element_text(size = 22), 
        legend.text = element_text(size = 18, family = "sans"),
        legend.position = 'top',
        legend.justification = 'left',
        axis.title = element_text(color = "black", size = 22),
        axis.ticks = element_line(size = 2, color = "black"), 
        axis.ticks.length = unit(0.35, "cm"), 
        axis.line = element_line(size = 1.75, color = "black"), 
        axis.text = element_text(color = "black", size = 20.0))

##use package cowplot to put figures on single page##
together<-plot_grid(rt_PCA, het_PCA, ncol = 1, labels = c("a)", "b)"), label_size = 18)


##PCA sex
jitter <- position_jitter(width = 0.4, height = 0.4)
ggplot(snp_pcs)+
  geom_point(aes(Axis1, Axis2, color = Gender), position = jitter, alpha = 0.5, size = 5.0) +
  scale_color_manual(values = c("#69b3a2", "#404080"), na.value = "orange", 
                     name = "Sex:", labels = c("Female", "Male")) +
  xlab("")+
  ylab("")+
  theme_classic() +
  theme(text = element_text(family = "sans"), 
        legend.title = element_text(size = 22), 
        legend.text = element_text(size = 22),
        legend.position = 'top',
        legend.justification = 'left',
        axis.ticks = element_line(size = 2, color = "black"), 
        axis.ticks.length = unit(0.35, "cm"), 
        axis.line = element_line(size = 1.75, color = "black"), 
        axis.text = element_text(color = "black", size = 20.0))




##look at 239 markers identified as Neutral by CRITFC##

nm<-which(ots_SFGL_panel_info$PresumedType == "Neutral")

neutral_marker_list<-ots_SFGL_panel_info[nm, ]

neutral_matrix<-as.matrix(coquille_23[, colnames(coquille_23) %in% neutral_marker_list$Assay])

rownames(neutral_matrix)<-coquille_23$sample_simple

neutral_genind <- df2genind(neutral_matrix, sep ="", ploidy=2,NA.char = "00")



##Individual level PCA###

X <- scaleGen(neutral_genind,  NA.method="mean")

#then run pca, keep all PCs
pca1 <- dudi.pca(X, scale = FALSE, scannf = FALSE, nf = 300)

snp_pcs <- pca1$li

##calculate percentage of total variation explained by each eigenvalue#

eig.perc <- 100*pca1$eig/sum(pca1$eig)

#create column of sample IDs from rownames##
snp_pcs %<>%
  rownames_to_column(var="sample") 

##check both in same order##

identical(snp_pcs$sample, meta$IndividualName)

##add meta data back in##

snp_pcs<-cbind(snp_pcs, meta)

snp_pcs %<>%
  relocate(sample, IndividualName, WILDorHAT., Marks, Gender)



##plot PCA with HOR / NOR comparison##
neutral_PCA<-ggplot(snp_pcs)+
  geom_point(aes(Axis1, Axis2, color = WILDorHAT.), alpha = 0.5, size = 5.0) +
  scale_color_manual(values = c("#69b3a2", "#404080"), na.value = "orange", 
                     name = "Origin:", labels = c("HOR", "NOR")) +
  xlab("PC1 4.0%")+
  ylab("PC2 3.6%")+
  theme_classic() +
  theme(text = element_text(family = "sans"), 
        legend.title = element_text(size = 22), 
        legend.text = element_text(size = 22),
        legend.position = 'top',
        legend.justification = 'left',
        axis.title = element_text(color = "black", size = 22),
        axis.ticks = element_line(size = 2, color = "black"), 
        axis.ticks.length = unit(0.35, "cm"), 
        axis.line = element_line(size = 1.75, color = "black"), 
        axis.text = element_text(color = "black", size = 20.0))

##use package cowplot to put figures on single page##
together<-plot_grid(rt_PCA, non_rt_PCA, ncol = 1, labels = c("a)", "b)"), label_size = 18)

##PCA sex 
jitter <- position_jitter(width = 0.4, height = 0.4)
ggplot(snp_pcs)+
  geom_point(aes(Axis1, Axis2, color = Gender), position = jitter, alpha = 0.5, size = 5.0) +
  scale_color_manual(values = c("#69b3a2", "#404080"), na.value = "orange", 
                     name = "Sex:", labels = c("Female", "Male")) +
  xlab("")+
  ylab("")+
  theme_classic() +
  theme(text = element_text(family = "sans"), 
        legend.title = element_text(size = 22), 
        legend.text = element_text(size = 22),
        legend.position = 'top',
        legend.justification = 'left',
        axis.ticks = element_line(size = 2, color = "black"), 
        axis.ticks.length = unit(0.35, "cm"), 
        axis.line = element_line(size = 1.75, color = "black"), 
        axis.text = element_text(color = "black", size = 20.0))



##look at all non run timing/sex markers##

non_run_timing<-coquille_23[, c(which(!colnames(coquille_23) %in% rt_markers$Assay))]

non_run_timing_matrix<-as.matrix(non_run_timing[, c(9:(ncol(non_run_timing)-1))])

rownames(non_run_timing_matrix)<-coquille_23$sample_simple

non_run_timing_genind <- df2genind(non_run_timing_matrix, sep ="", ploidy=2,NA.char = "00")


##Individual level PCA###

X <- scaleGen(non_run_timing_genind,  NA.method="mean")

#then run pca, keep all PCs
pca1 <- dudi.pca(X, scale = FALSE, scannf = FALSE, nf = 300)

snp_pcs <- pca1$li

##calculate percentage of total variation explained by each eigenvalue#

eig.perc <- 100*pca1$eig/sum(pca1$eig)

#create column of sample IDs from rownames##
snp_pcs %<>%
  rownames_to_column(var="sample") 

##check both in same order##

identical(snp_pcs$sample, meta$IndividualName)

##add meta data back in##

snp_pcs<-cbind(snp_pcs, meta)

snp_pcs %<>%
  relocate(sample, IndividualName, WILDorHAT., Marks, Gender)

##set order of NOR, HOR, and unknown (UK) origins

snp_pcs$WILDorHAT.<-factor(snp_pcs$WILDorHAT., levels = c('W', 'H', 'UK'))

##correct phenotypic -genotypic sex issue in 3 samples

snp_pcs %<>%
  mutate(Gender = case_when(sample %in% c("OtsAC23COQR_0107", "OtsAC23COQR_0281", "OtsAC23COQR_0282") ~ 'F', 
                            TRUE ~ as.character(Gender)
  ))



##plot PCA with HOR / NOR comparison##
non_run_timing_PCA<-ggplot(snp_pcs)+
  geom_point(aes(Axis1, Axis2, color = WILDorHAT.), alpha = 0.5, size = 5.0) +
  scale_color_manual(values = c("#404080", "#69b3a2", "orange"), 
                     name = "Origin:", labels = c("NOR", "HOR", "UK")) +
  xlab("PC1 3.7%")+
  ylab("PC2 3.4%")+
  theme_classic() +
  theme(text = element_text(family = "sans"), 
        legend.title = element_text(size = 22), 
        legend.text = element_text(size = 22),
        legend.position = 'top',
        legend.justification = 'left',
        axis.title = element_text(color = "black", size = 22),
        axis.ticks = element_line(size = 2, color = "black"), 
        axis.ticks.length = unit(0.35, "cm"), 
        axis.line = element_line(size = 1.75, color = "black"), 
        axis.text = element_text(color = "black", size = 20.0))

##PCA sex 

non_run_timing_sex_PCA<-ggplot(snp_pcs)+
  geom_point(aes(Axis1, Axis2, color = Gender), position = jitter, alpha = 0.5, size = 5.0) +
  scale_color_manual(values = c("#69b3a2", "#404080"), na.value = "orange", 
                     name = "Sex:", labels = c("Female", "Male")) +
  xlab("PC1 3.7%")+
  ylab("PC2 3.4%")+
  theme_classic() +
  theme(text = element_text(family = "sans"), 
        legend.title = element_text(size = 22), 
        legend.text = element_text(size = 22),
        legend.position = 'top',
        legend.justification = 'left',
        axis.title = element_text(color = "black", size = 22),
        axis.ticks = element_line(size = 2, color = "black"), 
        axis.ticks.length = unit(0.35, "cm"), 
        axis.line = element_line(size = 1.75, color = "black"), 
        axis.text = element_text(color = "black", size = 20.0))

##use package cowplot to put figures on single page##
together<-plot_grid(non_run_timing_PCA, non_run_timing_sex_PCA, ncol = 1, labels = c("a)", "b)"), label_size = 18)

##quantify FST across run timing and non-run timing markers##

library(pegas)

##with run timing markers##

##add NOR HOW classifications as population identifiers to the genind object

nm_df<-as.data.frame(rownames(genind_run_timing$tab))

colnames(nm_df)<-"Sample"

identical(nm_df$Sample, meta$IndividualName)

nm_df<-cbind(nm_df, meta$WILDorHAT.)


##add HOR NOR classification as population identifier to genind object of run timing markers
genind_run_timing$pop<-as.factor(nm_df$`meta$WILDorHAT.`)

run_timing_loci<-genind2loci(genind_run_timing) 

##estimate FST between NOR and HOR at run timing markers
differentiation_run_timing<-as.data.frame(Fst(run_timing_loci))

##make density plot of FST among run timing markers
diff_plot_run_timing<-ggplot(differentiation_run_timing, aes(x = Fst)) +
  geom_density(alpha = 0.4, fill = "grey", size = 1.5) +
  scale_y_continuous(expand = c(0.005, 0.0, 0.0, 0))+
  scale_x_continuous(breaks = c(0, 0.05, 0.1, 0.15, 0.2 ), limits = c(-0.02, 0.2))+
  xlab(expression(F["ST"]))+
  ylab("Density")+
  theme_classic()+
  theme(text = element_text(family = "sans"), 
        axis.title = element_text(color = "black", size = 22),
        axis.ticks = element_line(size = 1.5, color = "black"), 
        axis.ticks.length.x = unit(0.35, "cm"), 
        axis.ticks.length.y = unit(0, "cm"),
        axis.line = element_line(size = 1.25, color = "black"), 
        axis.text.x = element_text(color = "black", size = 16.0), 
        axis.text.y = element_blank())

##add HOR NOR classification as population identifier to genind object of non run timing markers
non_run_timing_genind$pop<-as.factor(nm_df$`meta$WILDorHAT.`)

non_run_timing_loci<-genind2loci(non_run_timing_genind) 

##estimate FST between NOR and HOR at run timing markers
differentiation_non_run_timing<-as.data.frame(Fst(non_run_timing_loci))

##make density plot of FST among run timing markers
diff_plot_non_run_timing<-ggplot(differentiation_non_run_timing, aes(x = Fst)) +
  geom_density(alpha = 0.4, fill = "grey", size = 1.5) +
  scale_y_continuous(expand = c(0.005, 0.0, 0.0, 0))+
  scale_x_continuous(breaks = c(0, 0.05, 0.1, 0.15, 0.2 ), limits = c(-0.02, 0.2))+
  xlab(expression(F["ST"]))+
  ylab("Density")+
  theme_classic()+
  theme(text = element_text(family = "sans"), 
        axis.title = element_text(color = "black", size = 22),
        axis.ticks = element_line(size = 1.5, color = "black"), 
        axis.ticks.length.x = unit(0.35, "cm"), 
        axis.ticks.length.y = unit(0, "cm"),
        axis.line = element_line(size = 1.25, color = "black"), 
        axis.text.x = element_text(color = "black", size = 16.0), 
        axis.text.y = element_blank())

##use package cowplot to put figures on single page##
together<-plot_grid(diff_plot_non_run_timing, diff_plot_run_timing, ncol = 1, labels = c("a)", "b)"), label_size = 18)


##calculate mean FST and 95% CI across run timing markers##

mean(differentiation_run_timing$Fst)

FST_mean_run_timing_ci<-summarySE(differentiation_run_timing, measurevar = "Fst")

mean(differentiation_non_run_timing$Fst)

FST_mean_non_run_timing_ci<-summarySE(differentiation_non_run_timing, measurevar = "Fst")


#calculate expected heterozygosities for HOR and NOR at run timing markers##

hwe.pop.runtiming <- seppop(genind_run_timing)

NOR_loci_run_timing<-genind2loci(hwe.pop.runtiming$W)

HOR_loci_run_timing<-genind2loci(hwe.pop.runtiming$H)

NOR_het_run_timing<-heterozygosity(NOR_loci_run_timing, variance = TRUE)

HOR_het_run_timing<-heterozygosity(HOR_loci_run_timing, variance = TRUE)

NOR_het_run_timing_ci<-summarySE(NOR_het_run_timing, measurevar = "Hs")
HOR_het_run_timing_ci<-summarySE(HOR_het_run_timing, measurevar = "Hs")

#calculate expected heterozygosities for HOR and NOR at non run timing markers##

hwe.pop.non_runtiming <- seppop(non_run_timing_genind)

NOR_loci_non_run_timing<-genind2loci(hwe.pop.non_runtiming$W)

HOR_loci_non_run_timing<-genind2loci(hwe.pop.non_runtiming$H)

NOR_het_non_run_timing<-heterozygosity(NOR_loci_non_run_timing, variance = TRUE)

HOR_het_non_run_timing<-heterozygosity(HOR_loci_non_run_timing, variance = TRUE)

NOR_het_non_run_timing_ci<-summarySE(NOR_het_non_run_timing, measurevar = "Hs")
HOR_het_non_run_timing_ci<-summarySE(HOR_het_non_run_timing, measurevar = "Hs")




##Sample same number of samples in both groups and calculate expected heterozygosity##

##run timing markers##

run_timing<-coquille_23[, c(1:8, which(colnames(coquille_23) %in% rt_markers$Assay))]

run_timing %<>%
  filter(!sample_simple == "OtsAC23COQR_0286")

#drop Ots37124-12270118

run_timing %<>%
  select(!'Ots37124-12270118')

##put markers in chromosomal order##

run_timing<-run_timing[, c(1:12, 36:37, 13:16, 38, 17:23, 39, 24:35)]

rt_matrix <- as.matrix(run_timing[,c(9:(ncol(run_timing)))]) 
row.names(rt_matrix) <- run_timing$sample_simple
genind_run_timing <- df2genind(rt_matrix, sep ="", ploidy=2,NA.char = "00")

##add NOR HOR classisfications as population identifiers to the genind object

nm_df<-as.data.frame(rownames(genind_run_timing$tab))

colnames(nm_df)<-"Sample"


##remove individual with unknown origin from meta

meta %<>%
  filter(!IndividualName == "OtsAC23COQR_0286")


identical(nm_df$Sample, meta$IndividualName)

nm_df<-cbind(nm_df, meta$WILDorHAT.)

##add HOR NOR classification as population identifier to genind object of run timing markers
genind_run_timing$pop<-as.factor(nm_df$`meta$WILDorHAT.`)

hwe.pop.runtiming <- seppop(genind_run_timing)

##sample 20 individuals with HOR and NOR origin
library(pegas)

NOR<-rep(0, times = 1000)

HOR<-rep(0, times = 1000)

for (i in 1:length(NOR)){
mySamp <- lapply(hwe.pop.runtiming, function(x) x[sample(1:nrow(x$tab), 20)])
NOR_loci_run_timing<-genind2loci(mySamp$W)
HOR_loci_run_timing<-genind2loci(mySamp$H)
NOR_het_run_timing<-heterozygosity(NOR_loci_run_timing, variance = FALSE)
HOR_het_run_timing<-heterozygosity(HOR_loci_run_timing, variance = FALSE)
NOR[i]<-mean(NOR_het_run_timing[, 1])
HOR[i]<-mean(HOR_het_run_timing[, 1])
}


##non run timing markers##
non_run_timing<-coquille_23[, c(which(!colnames(coquille_23) %in% rt_markers$Assay))]

non_run_timing %<>%
  filter(!sample_simple == "OtsAC23COQR_0286")


non_run_timing_matrix<-as.matrix(non_run_timing[, c(9:(ncol(non_run_timing)-1))])
rownames(non_run_timing_matrix)<-non_run_timing$sample_simple
non_run_timing_genind <- df2genind(non_run_timing_matrix, sep ="", ploidy=2,NA.char = "00")


##add NOR HOR classisfications as population identifiers to the genind object

nm_df<-as.data.frame(rownames(non_run_timing_genind$tab))

colnames(nm_df)<-"Sample"


##remove individual with unknown origin from meta

meta %<>%
  filter(!IndividualName == "OtsAC23COQR_0286")


identical(nm_df$Sample, meta$IndividualName)

nm_df<-cbind(nm_df, meta$WILDorHAT.)

##add HOR NOR classification as population identifier to genind object of non run timing markers
non_run_timing_genind$pop<-as.factor(nm_df$`meta$WILDorHAT.`)

hwe.pop.nonruntiming <- seppop(non_run_timing_genind)

NOR<-rep(0, times = 1000)

HOR<-rep(0, times = 1000)

for (i in 1:length(NOR)){
  mySamp <- lapply(hwe.pop.nonruntiming, function(x) x[sample(1:nrow(x$tab), 20)])
  NOR_loci<-genind2loci(mySamp$W)
  HOR_loci<-genind2loci(mySamp$H)
  NOR_het<-heterozygosity(NOR_loci, variance = FALSE)
  HOR_het<-heterozygosity(HOR_loci, variance = FALSE)
  NOR[i]<-mean(NOR_het[, 1])
  HOR[i]<-mean(HOR_het[, 1])
}

##look at linkage

library(poppr)

non_run_timing_poppr<-genind2genalex(non_run_timing_genind, filename = "OtsAC23_nonruntiming_genalex.csv")

non_run_timing_poppr<-read.genalex("C:/Users/olsenk2/OneDrive - Oregon State University/Desktop/OtsAC22_23COQR/OtsAC23_nonruntiming_genalex.csv")

NOR<-popsub(non_run_timing_poppr, sublist = "W")

res<-pair.ia(NOR)

HOR<-popsub(non_run_timing_poppr, sublist = "H")

res_HOR<-pair.ia(HOR)

ia(NOR)

resample.ia(NOR, n = 20, reps = 99)

ia(HOR)

##calculate sequencing depth statistics##

##read in marker info

marker_info<-read_csv("C:/Users/olsenk2/OneDrive - Oregon State University/Desktop/OtsAC22_23COQR/marker_info.csv")

marker_info %<>%
  mutate(a1_count =  as.numeric(substr(a1_count, 3, nchar(a1_count)))) %>%
  mutate(a2_count =  as.numeric(substr(a2_count, 3, nchar(a2_count)))) %>%
  mutate(ind = str_remove(ind, "^\\./")) %>%
  mutate(ind = str_remove(ind, "\\.genos"))

##filter marker info file to samples that passed quality filtering in the 2023 Coquille samples and the 326 markers used in analyses (31 run timing ; 295 non run timing)##

marker_info %<>%
  filter(ind %in% coquille_23$Sample) %>%
  mutate(marker = gsub("\\.", "_", marker)) %>%
  filter(marker %in% colnames(rt_matrix) | marker %in% colnames(non_run_timing_matrix))


marker_info$total_depth<-marker_info$a1_count+marker_info$a2_count

depth_ci<-summarySE(marker_info, measurevar = "total_depth")

marker_info$ind<-factor(marker_info$ind, levels = c(unique(marker_info$ind)))

mean_depths<-marker_info %>%
  group_by(ind) %>%
  summarise(mean = mean(total_depth, na.rm = TRUE))

##use dartR to calculate FIS##

library(dartR)

non_run_timing_genlight<-gi2gl(non_run_timing_genind)

gl.report.heterozygosity(non_run_timing_genlight, 
                         method = "pop", 
)



##evaluate effective size##
##put data into Ne estimator format##

##start with non run timing markers##

replacedGT<-as.data.frame(non_run_timing[, c(9:(ncol(non_run_timing)-1))])


##replace nucleotides with numeric values## 
##A = 1; T = 2; C = 3; G = 4; - = 5 ###


replacedGT[, ] %<>%
  mutate_all(funs(case_when(. == "AA" ~ 11,
                            . == "AT" ~ 12, 
                            . == "AC" ~ 13, 
                            . == "AG" ~ 14, 
                            . == "TT" ~ 22, 
                            . == "TA" ~ 21, 
                            . == "TC" ~ 23, 
                            . == "TG" ~ 24,
                            . == "CC" ~ 33, 
                            . == "CA" ~ 31, 
                            . == "CT" ~ 32, 
                            . == "CG" ~ 34, 
                            . == "GG" ~ 44, 
                            . == "GA" ~ 41, 
                            . == "GT" ~ 42, 
                            . == "GC" ~ 43, 
                            . == "--" ~ 55, 
                            . == "A-" ~ 15, 
                            . == "T-" ~ 25, 
                            . == "C-" ~ 35, 
                            . == "G-" ~ 45, 
                            . == "-A" ~ 51,
                            . == "-T" ~ 52, 
                            . == "-C" ~ 53, 
                            . == "-G" ~ 54, 
                            . == "00" ~ 00, 
                            . == "0" ~ 00
  )))






##add sample name and origin
replacedGT$Sample<-non_run_timing$sample_simple

identical(replacedGT$Sample, meta$IndividualName)

replacedGT$origin<-meta$WILDorHAT.

replacedGT %<>%
  relocate(Sample, origin)

##write file with locus names separated by commas 

marker_names<-colnames(replacedGT)[-c(1, 2)]

write.table(marker_names, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",", file = "Nonruntiming marker names Ne estimator.txt")

##write file with genotypes 
##subset hatchery fish first##

HOR<-subset(replacedGT, replacedGT$origin == "H")

HOR %<>%
  mutate(Sample = paste(Sample, ",", sep = ""))


##drop origin column

HOR %<>%
  select(!origin)


write.table(HOR, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t", file = "Nonruntiming genotypes HOR Ne estimator.txt")


##subset natural origin#

NOR<-subset(replacedGT, replacedGT$origin == "W")

NOR %<>%
  mutate(Sample = paste(Sample, ",", sep = ""))


##drop origin column 

NOR %<>%
  select(!origin)

write.table(NOR, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t", file = "Nonruntiming genotypes NOR Ne estimator.txt")



##read in F estimates from coancestry##

##read in relatedness estimates##


inbreeding_F<-read.table("C:/ZSL/Coancestry/OtsAC2223COQR/InbreedingEstimates.txt")

colnames(inbreeding_F)<-c("Sample", "Ritland", "LynchRD")

##subset inbreeding values for 2023 samples##

inbreeding_F %<>%
  mutate(Sample = paste("Ots", Sample, sep = ""))

inbreeding_F_23<-inbreeding_F[inbreeding_F$Sample %in% non_run_timing$sample_simple, ]

identical(inbreeding_F_23$Sample, meta$IndividualName)

inbreeding_F_23$origin<-meta$WILDorHAT.

inbreeding_F_23 %<>%
  relocate(Sample, origin)

##calculate mean SD for both groups 

NOR_F<-subset(inbreeding_F_23, inbreeding_F_23$origin == "W")
HOR_F<-subset(inbreeding_F_23, inbreeding_F_23$origin == "H")

NOR_F_non_run_timing_ci<-summarySE(NOR_F, measurevar = "LynchRD")
HOR_F_non_run_timing_ci<-summarySE(HOR_F, measurevar = "LynchRD")







##evaluate relatedness within NOR and HOR fish with non run timing markers

##load QC filtered genos

load("C:/Users/olsenk2/OneDrive - Oregon State University/Desktop/OtsAC22_23COQR/Ots COQR 22_23 genos_2.0.R")

##read in COQR 2023 meta##

meta<-read.csv("C:/Users/olsenk2/OneDrive - Oregon State University/Desktop/OtsAC22_23COQR/Ots23COQR_meta.csv", header = TRUE)

##subset dataset to just 2023 samples
coquille_23<-genos_2.0[genos_2.0$sample_simple %in% meta$IndividualName, ]

##subset meta to QC filtered samples#

meta<-meta[meta$IndividualName %in% coquille_23$sample_simple, ]


##read marker information for Ots GT-seq panel

ots_panel_info<-read_excel("C:/Users/olsenk2/OneDrive - Oregon State University/Desktop/OtsAC22_23COQR/Ots GT-seq panel/Ots GT-seq panel info.xlsx", sheet = 2)

##subset panel info to those markers in the SFGL panel

ots_SFGL_panel_info<-ots_panel_info[ots_panel_info$Assay %in% colnames(genos_2.0), ]

##identify run timing markers#
rt_markers<-subset(ots_SFGL_panel_info, ots_SFGL_panel_info$PresumedType == "runtiming")

##look at all non run timing/sex markers##

non_run_timing<-coquille_23[, c(which(!colnames(coquille_23) %in% rt_markers$Assay))]



replacedGT<-as.data.frame(non_run_timing[, c(9:(ncol(non_run_timing)-1))])

delim_gtypes <- replacedGT[, ] %>%
  mutate_all(funs(case_when(. == "AA" ~ "A/A",
                            . == "AT" ~ "A/T", 
                            . == "AC" ~ "A/C", 
                            . == "AG" ~ "A/G", 
                            . == "TT" ~ "T/T", 
                            . == "TA" ~ "T/A", 
                            . == "TC" ~ "T/C", 
                            . == "TG" ~ "T/G",
                            . == "CC" ~ "C/C", 
                            . == "CA" ~ "C/A", 
                            . == "CT" ~ "C/T", 
                            . == "CG" ~ "C/G", 
                            . == "GG" ~ "G/G", 
                            . == "GA" ~ "G/A", 
                            . == "GT" ~ "G/T", 
                            . == "GC" ~ "G/C", 
                            . == "--" ~ "-/-", 
                            . == "A-" ~ "A/-", 
                            . == "T-" ~ "T/-", 
                            . == "C-" ~ "C/-", 
                            . == "G-" ~ "G/-", 
                            . == "-A" ~ "-/A",
                            . == "-T" ~ "-/T", 
                            . == "-C" ~ "-/C", 
                            . == "-G" ~ "-/G", 
                            . == "00" ~ "0/0", 
                            . == "0" ~ "0/0"
  )))



library(splitstackshape)


vrtnames<-c(colnames(delim_gtypes))

##split genotypes into 2 alleles

dfspltgeno<-as.data.frame(cSplit(delim_gtypes, splitCols = vrtnames, sep = "/", type.convert = FALSE))


coancestrylbls<-data.frame(c(non_run_timing$sample_simple))

colnames(coancestrylbls)<-c("Ind")

coancestrylbls %<>%
  mutate(Ind = substr(Ind, 4, 16))


rownames(dfspltgeno)<-as.character(coancestrylbls$Ind)

##replace nucleotides with numeric values for coancestry## 
##A = 1; T = 2; C = 3; G = 4; - = 5 ###

dfspltgeno %<>%
  mutate_all(funs(case_when(. == "A" ~ 1,
                            . == "T" ~ 2,
                            . == "C" ~ 3,
                            . == "G" ~ 4, 
                            . == "-" ~ 5,
                            . == "0" ~ 0, 
  )))



write.table(dfspltgeno, file = "C:/Users/olsenk2/OneDrive - Oregon State University/Desktop/OtsAC22_23COQR/OtsAC23COQR non run timing coancestry input.txt", quote = FALSE, sep = "\t", col.names = FALSE, row.names = TRUE)


##read in relatedness estimates##


rel<-read.table("C:/ZSL/Coancestry/Ots23COQR non run timing markers/RelatednessEstimates.txt", sep = ",")

##grab wang estimator##

r<-c(1:3, 6)

rel<-rel[, r]

colnames(rel)<-c("dyad", "Ind1", "Ind2", "wang")


##compare relatedness within NOR and HOR

##add "Ots" back into sample name, add year sampled for ind 1 and 2

rel %<>%
  mutate(Ind1 = paste("Ots", Ind1, sep = "")) %>%
  mutate(Ind2 = paste("Ots", Ind2, sep = "")) %>%
  mutate(Ind1_samp_year = substr(Ind1, 6, 7)) %>%
  mutate(Ind2_samp_year = substr(Ind2, 6, 7))


##read in COQR 2023 meta##

meta<-read.csv("C:/Users/olsenk2/OneDrive - Oregon State University/Desktop/OtsAC22_23COQR/Ots23COQR_meta.csv", header = TRUE)

##add NOR HOR origin to relatedness table 

HOR_meta<-subset(meta, meta$WILDorHAT. == "H")
NOR_meta<-subset(meta, meta$WILDorHAT. == "W")

##remove individual with unknown origin 

rel %<>%
  filter(!Ind1 == "OtsAC23COQR_0286") %>%
  filter(!Ind2 == "OtsAC23COQR_0286")


rel %<>%
  mutate(Ind1_origin = case_when(Ind1 %in% NOR_meta$IndividualName ~ 'NOR', 
                                 Ind1 %in% HOR_meta$IndividualName ~ 'HOR')) %>%
  mutate(Ind2_origin = case_when(Ind2 %in% NOR_meta$IndividualName ~ 'NOR', 
                                 Ind2 %in% HOR_meta$IndividualName ~ 'HOR')) %>%
  mutate(Origins = paste(Ind1_origin, Ind2_origin, sep = ""))


##subset relatedness table to only HOR-HOR and NOR-NOR estimates

rel_HOR_NOR<-subset(rel, rel$Origins == "HORHOR" | rel$Origins == "NORNOR")

##bootstrap estimates of relatedness in NOR and HOR and make density plot

##dataframe with just origin and relatedness##
df<-rel_HOR_NOR[, c(9, 4)]

##subset each origin##

HORHOR<-subset(df, df$Origins == "HORHOR")
NORNOR<-subset(df, df$Origins == "NORNOR")

##for HORHOR ##
resamp_HORHOR <- lapply(1:10000, function(i) sample(HORHOR$wang, size = 1, replace = TRUE))
mean_HORHOR<-as.data.frame(sapply(resamp_HORHOR, mean))

mean_HORHOR$origin<-rep("HOR-HOR", times = 10000)
colnames(mean_HORHOR)[1]<-"estimate"


##for NORNOR
resamp_NORNOR <- lapply(1:10000, function(i) sample(NORNOR$wang, size = 1, replace = TRUE))
mean_NORNOR<-as.data.frame(sapply(resamp_NORNOR, mean))

mean_NORNOR$origin<-rep("NOR-NOR", times = 10000)
colnames(mean_NORNOR)[1]<-"estimate"

booted<-rbind(mean_HORHOR, mean_NORNOR)

booted$origin<-factor(booted$origin, levels = c("NOR-NOR", "HOR-HOR"))

##density plot##
ggplot(booted, aes(x = estimate, fill = origin))+
  geom_density(alpha = 0.6)+
  scale_fill_manual(values=c("#404080", "#69b3a2"), 
                    name = "Origin:", labels = c('NOR-NOR', 'HOR-HOR')) +
  scale_x_continuous(expand = c(0.0, 0.05)) +
  scale_y_continuous(expand = c(0.0, 0.0, 0.15, 0), breaks = c(0, 250, 500, 750, 1000))+
  xlab("")+
  ylab("")+
  theme_classic()+
  theme(text = element_text(family = "sans"),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.position = c(0.9, 0.8),
        axis.ticks.length.y = unit(0, "cm"),
        axis.ticks.x = element_line(size = 1.5, color = "black"),
        axis.ticks.length.x = unit(0.35, "cm"),
        axis.text.y = element_blank(),
        axis.text.x = element_text(color = "black", size = 16.0),
        axis.line = element_line(size = 1.25, color = "black"))



##raw distribution 

rel_HOR_NOR$Origins<-factor(rel_HOR_NOR$Origins, levels = c("NORNOR", "HORHOR"))


##density plot##
ggplot(rel_HOR_NOR, aes(x = wang, fill = Origins))+
  geom_density(alpha = 0.6)+
  scale_fill_manual(values=c("#404080", "#69b3a2"), 
                    name = "Origin:", labels = c('NOR-NOR', 'HOR-HOR')) +
  scale_x_continuous(expand = c(0.0, 0.05)) +
  scale_y_continuous(expand = c(0.0, 0.0, 0.15, 0))+
  xlab("")+
  ylab("")+
  theme_classic()+
  theme(text = element_text(family = "sans"),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.position = c(0.9, 0.8),
        axis.ticks.length.y = unit(0, "cm"),
        axis.ticks.x = element_line(size = 1.5, color = "black"),
        axis.ticks.length.x = unit(0.35, "cm"),
        axis.text.y = element_blank(),
        axis.text.x = element_text(color = "black", size = 16.0),
        axis.line = element_line(size = 1.25, color = "black"))



perind<-rel_HOR_NOR %>%
  group_by(Ind1) %>%
  summarise(mean = mean(wang))
