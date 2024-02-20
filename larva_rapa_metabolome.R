
# analyze targeted metabolomics run on larva samples (some of Mitchell Lee's rapa flies).

suppressPackageStartupMessages({
library(data.table)
library(tidyverse)
library(stringr)
library(RColorBrewer)
library(lme4)
library(lmerTest)
library(caret)
library(ggpubr)
library(RColorBrewer)
library(pheatmap)
library(sjPlot)
library(sjmisc)
library(pls)
library(igraph)
library(dendextend)
library(ggrepel)
library(effects)
library(emmeans)
library(effectsize)
library(sjstats)
library(gridExtra)
})

rm(list=ls())
gc()
dev.off()

setwd("~/Google Drive/My Drive/Documents/rapa/metabolomics/larva_metabolomics")
dir()

# skip to line 331 for processed data

# Danijel ran the targeted LC-MS analysis in the Raftery Lab.  He sent us a complicated excel file.  5 sheets, (check them out!) most with non R-frendly formatting.  I manually formatted them in the working directory for import into R.

# the file that Danijel sent is called:
# 2022-03-17_Promislow-54_Data.xlsx
dir('data')

######################################################################
# Sample ID
######################################################################
samps <- read.csv('~/Google Drive/My Drive/Documents/rapa/metabolomics/larva_metabolomics/data/Sample_ID.csv', header=T, stringsAsFactors = F)
head(samps)
colnames(samps)[2] <- 'BCA.assay'

head(samps)
samps$sampleNumber <- sprintf("%02d", as.numeric(rownames(samps))) # adds a leading zero so order() will work
head(samps)

table(grepl('QC', samps$Sample.ID)) # 54 biological samples, 14 QC samples

biosamps <- samps[!grepl('QC', samps$Sample.ID), ]
head(biosamps)

# the sample information returned from the Rafery lab is partially ambiguous.  The sample ids returned are the dgrp line number, and the treatment (E = ethanol, R = rapa). There were however 2 collection days and most line - treatment combinations were collected on both days.  I (Ben H.) had added a red stripe or a blue stripe to the tube tops and asked that the 28 blue-stripe samples be extracted together (from Mitchell's sample box A, see image in Ben's Notebook VI pg 2), and the 26 red-stripe samples be prepped together (image Notebook VI pg 2).

# Hayley J Purcell and Danijel worked to resolve how the preps were done and Haley said 'Looks like the study was prepared in two batches, first one had 28 samples and the second had 26.', which corresponds to the instructions on the sample submission sheet to prep the blue-stripee and red-stripe samples each together in batches.  

biosamps$prep.batch <- as.factor(c(rep('A', 28), rep('B', 26))) # assuming preps were done and run in order
biosamps$sample_ID <- paste0('sample_', biosamps$Sample.ID, '_', biosamps$prep.batch)
biosamps$BCA.assay <- as.numeric(biosamps$BCA.assay)
head(biosamps)
biosamps <- biosamps %>% separate(Sample.ID, c("line", "treatment"))
head(biosamps)

# [NOTE] samples from E and R of dgrp line 57 were collected by Mitchell on BOTH collection dates, but the ones from Oct 15th 2021 (the second day of collection, which ended up with red-stripes) were removed and run by Danijel to test the LCMS ability on larval samples.

table(biosamps$line, biosamps$treatment)
## problem with the sample IDs: there are three 'R' (rapa) samples for line 318, and only one 'E' (ethanol) sample, when there should be 2 of each.# I don't think I will be able to tell which of the R samples is the E sample for line 318...

table(biosamps$line, biosamps$treatment, biosamps$prep.batch) # looks like one of the two '318' samples from block A is mislabeled.
biosamps$treatment[biosamps$line == '318' & biosamps$prep.batch=='A'] <- NA
table(biosamps$treatment)
table(is.na(biosamps$treatment))

# rapa phenotype for each line:
meta <- read.csv('~/Google Drive/My Drive/Documents/rapa/metabolomics/larva_metabolomics/Mitchell_sample_meta.data.csv', header=T, stringsAsFactors = F)
head(meta)
meta$TREAT <- as.factor(meta$treatment)
meta$treatment <- as.factor(ifelse(meta$TREAT=='EtOH', 'E', 'R'))
meta$prep.batch <- as.factor(ifelse(meta$collection.date == '10_14_2021', 'A', 'B'))
head(meta)

head(biosamps)
meta$sample_ID <- paste0('sample_', meta$line, ' ', meta$treatment, '_', meta$prep.batch)
head(meta)

table(meta$sample_ID %in% biosamps$sample_ID)
meta[!meta$sample_ID %in% biosamps$sample_ID, ] # remember, the true 318 R sample is not known 

biosamps <- merge(biosamps, meta, all.x=T) 
head(biosamps)
str(biosamps)
biosamps$phenotype <- as.factor(biosamps$phenotype)
biosamps$treatment <- as.factor(biosamps$treatment)

table(biosamps$line, biosamps$treatment) # sample summary
table(biosamps$line, biosamps$phenotype) 

# This sheet provides sample numbers that we have used for your samples (NWMRC ID) and the corresponding sample ID. The order of sample numbers shows the order in which the samples were analysed.														
## Quality Control (QC) sample codes:
# QC(I) = Instrument QC (Pooled Serum Samples): 
# QC(I)#1, QC(I)#2, etc. indicate the same instrument quality control sample analyzed repeatedly, every 10 samples. A pooled commercial human plasma sample was used as the instrument quality control to monitor the instrument stability during the analysis.														

# QC(S) Sample QC (Pooled Study Samples): 
# QC(S)#1, QC(S)#2, etc. indicate the same sample quality control sample analyzed repeatedly, every 10 samples. A pooled sample made using a small portion from each of your samples was used as the sample quality control.			
######################################################################


######################################################################
# Data Reproducibility
######################################################################
repro <- read.csv('data/LC_MS_reproducibility.csv', header=T, stringsAsFactors = T)
head(repro)

# This sheet provides metabolite data for the instrument quality control samples [QC(I)'s] and sample quality control samples [QC(S)'s] used in the analysis as indicated in Sheet 1.																					

# Coefficient of variation (CV) for each metabolite, average CV and median CV are also provided for both QC(I) and QC(S) samples to assess the instrument performance during the analysis.																											
# Any metabolite not detected is indicated as N/A.							
# QC(S) reproducibility data could be used to normalize MS signal for those metabolites with CV values > 10%. 	
######################################################################

######################################################################
# Relative Quant Data
######################################################################
dat <- read.csv('data/lcms.data.csv', header=T, stringsAsFactors = T)
head(dat)
colnames(dat) <- gsub('X', 'sample_', colnames(dat))
dat[1:5,1:8]

# This sheet provides the full list of metabolites targeted in the analysis and their relative concentration.  
# HMDB ID and KEGG ID are also provided for each metabolite in columns 2 and 3, respectively, to help with further analysis of the data such as pathway analysis, if needed.
# Sample IDs shown in the first row are the numbers we have given (NWMRC IDs). Please refer to Sheet 1 for the corresponding ID for your samples.
# When isotope labeled internal standards are used for absolute quantitation or additional quality check, they are listed towards the end of the sheet.
# Any metabolite not detected in your sample is indicated as N/A.
# Please note, the data is not normalized either with reference to the protein amount or quality control data [QC(I) or QC(S)].
######################################################################


######################################################################
# Metabolite Information
######################################################################
mz.info <- read.csv('data/metabolite.info.csv', header=T, stringsAsFactors = F)
head(mz.info)

# This sheet provides some important information on the metabolites targeted in the analysis.
# Information provided includes HMDB ID, KEGG ID, molecular formula, CAS number and typical pathway each metabolite represents.
######################################################################


######################################################################
# Analyze missingness
######################################################################
dat[1:4,1:8]
dat[dat=='N/A'] <- NA # replace 'N/A' with NA
hist(rowSums(is.na(dat[ ,4:ncol(dat)])))

# some metabolites have patterns of missingness that may relate to the biological variables
############################################################################################################
dat[1:4,1:8]
dim(dat)
hist(rowSums(is.na(dat[ ,4:ncol(dat)])))
dat[rowSums(is.na(dat[ ,4:ncol(dat)])) < 40 & rowSums(is.na(dat[ ,4:ncol(dat)])) > 10, ] # extract metabolites that are measured in at least 11 samples, but fewer than 40
hit.n.miss.mzs <- as.character(dat[rowSums(is.na(dat[ ,4:ncol(dat)])) < 40 & rowSums(is.na(dat[ ,4:ncol(dat)])) > 10, ][ ,1])
hit.n.miss <- t(dat[rowSums(is.na(dat[ ,4:ncol(dat)])) < 40 & rowSums(is.na(dat[ ,4:ncol(dat)])) > 10, ][ ,4:ncol(dat)])
colnames(hit.n.miss) <- hit.n.miss.mzs
head(hit.n.miss) 
hit.n.miss <- cbind(biosamps, hit.n.miss)
head(hit.n.miss)

k <- !is.na(hit.n.miss$cAMP)
table(k, hit.n.miss$line) # pattern of missingness assoc with genotype
fisher.test(print(table(k, hit.n.miss$treatment)))
fisher.test(print(table(k, hit.n.miss$phenotype)))
print(table(k, hit.n.miss$treatment, hit.n.miss$phenotype))
# cAMP shows a pattern of missingness assoc with genotype, but no other variables


k <- !is.na(hit.n.miss$UDP)
table(k, hit.n.miss$line) # pattern of missingness assoc with genotype
fisher.test(print(table(k, hit.n.miss$treatment)))
fisher.test(print(table(k, hit.n.miss$phenotype)))
print(table(k, hit.n.miss$treatment, hit.n.miss$phenotype))
# UDP shows a strong pattern of missingness assoc with sex, and maybe sex by treatment, but no other variables

k <- !is.na(hit.n.miss$NADH)
table(k, hit.n.miss$line) # pattern of missingness assoc with genotype
fisher.test(print(table(k, hit.n.miss$treatment)))
fisher.test(print(table(k, hit.n.miss$phenotype)))
print(table(k, hit.n.miss$treatment, hit.n.miss$phenotype))
# NADH a pattern of missingness assoc with genotype, but no other variables

# remove metabolties with missing values, based on the above analysis, this is unlikely to remove interesting metabolites
############################################################################################################
rm(hit.n.miss)
rm(hit.n.miss.mzs)
rm(k)

table(rowSums(is.na(dat[ ,4:ncol(dat)])) < 1)  # 154 metabolites and 4 spiked isotopes detected in all samples; 361 metabolites measured (last 32 rows are isotopes, of which 4 are in the 158 rows w/o missing data)

dat <- dat[rowSums(is.na(dat[ ,4:ncol(dat)])) < 1, ] # keep 154 complete metabolites and 4 isotopes
 

#####################################################################
# data normalization and potential confounder analysis
#####################################################################
dat[1:4,1:6]
compounds <- as.character(dat$COMPOUND)
grepl('C13', compounds) # [NOTE] the last few rows of dat are the spiked isotopic standards. 
dat[153:nrow(dat), 1:4] # a look at the bottom
isotopes <- compounds[grepl('C13', compounds)]

dat[1:4,1:6]
mat <- as.matrix(t(dat[ ,4:ncol(dat)]))
mat <- apply(mat, 2, as.numeric)
table (is.na(mat))
mat[1:4,1:4]
colnames(mat) <- compounds
head(biosamps)
biosamps <- biosamps[order(biosamps$sampleNumber), ]
rownames(mat) <- biosamps$sampleNumber
mat[1:4,1:4]

mzs <- colnames(mat)[!grepl('C1', colnames(mat))] # grab the metabolite names (but not the spiked isotopes too)

par(mfrow=c(1,2))
boxplot(mat, las=1)
boxplot(log(mat), las=1)
logmat <- log(mat)

#######################################################
# protein quantity (via BCA assay) a likely reflection of sample mass/#larva shows a non-linear relationship with many metabolites:
#######################################################

par(mfrow=c(2,2))
plot(BCA.assay ~ n_larvae, biosamps, pch=ifelse(biosamps$n_larvae>20, 19, 1))
abline(v=22.5, lty=2)
boxplot(t(logmat)[ ,order(biosamps$BCA.assay)], xlab='samples (ordered by protein quantity)', ylab='log(mzs)', las=1, main='protein quantity predicts metaboloite levels') 

table(biosamps$n_larvae > 20, biosamps$line, biosamps$treatment)
table(biosamps$n_larvae > 20, biosamps$line, biosamps$collection.date)

# for 13 of 14 lines, there is at least one collection date with >20 larva per sample per treatment.  

# I will apply a cutoff to remove samples with <=20 larva across both treatments for a collection date
lines <- colnames(table(biosamps$n_larvae > 20, biosamps$line, biosamps$collection.date))

keepersA <- lines[table(biosamps$n_larvae > 20, biosamps$line, biosamps$collection.date)[2, ,1] == 2] # lines to keep from the first collection date
keepersB <- lines[table(biosamps$n_larvae > 20, biosamps$line, biosamps$collection.date)[2, ,2] == 2] # lines to keep from the second collection date

keepersA <- biosamps[biosamps$line %in% keepersA & biosamps$prep.batch == 'A', ]
keepersB <- biosamps[biosamps$line %in% keepersB & biosamps$prep.batch == 'B', ]

keepers <- rbind(keepersA, keepersB)

table(keepers$n_larvae > 20, keepers$line, keepers$collection.date) 
# of the 52 biological samples, this cutoff kept 42 (removing 10)
keepers



par(mfrow=c(2,2))
plot(BCA.assay ~ n_larvae, biosamps, pch=ifelse(biosamps$n_larvae>20, 19, 1), las=1)
abline(v=22.5, lty=2)
boxplot(t(logmat)[ ,order(biosamps$BCA.assay)], xlab='samples (ordered by protein quantity)', ylab='log(mzs)', las=1, main='protein quantity predicts metaboloite levels') 
plot(BCA.assay ~ n_larvae, keepers, pch=19, main='after cutoff (larvae > 20 per treatment per batch)', las=1)


#######################################################
# apply the cutoff
cutlogmat <- logmat[keepers$sampleNumber, ]
boxplot(cutlogmat)
boxplot(t(scale(t(cutlogmat[ ,mzs])))) # don't include the isotopes in this normalization
boxplot(scale(t(cutlogmat[ ,mzs])))

scat <- t(scale(t(cutlogmat[ ,mzs])))

mzdat <- scat
biosamps <- keepers
#######################################################

# look at prep-batch effects and run-order effects (at the same time?  correcting for batch might remove run order effects, and run order effects should be removed on their own)
mzdat[1:4,1:4]
dim(mzdat)
run.order <- 1:nrow(mzdat)

dat <- cbind(biosamps, mzdat[biosamps$sampleNumber, ])
dat[1:4,1:16]

par(mfrow=c(4,6))
for(i in 1:length(mzs)) {
  plot(dat[ , mzs[i]], pch=16, ylab=mzs[i], col=dat$prep.batch) }

combat <- sva::ComBat(t(mzdat), batch=biosamps$prep.batch) 
combat <- t(combat)
colnames(combat) <- mzs

round(apply(mzdat, 1, mean), 3)
round(apply(combat, 1, mean), 3) 

# alternative batch correction, a simple lm of mz ~ batch
batchcor <- matrix(nr=nrow(mzdat), nc=ncol(mzdat))
for(i in 1:length(mzs)) {
  batchcor[ ,i] <- lm(dat[ , mzs[i]] ~ dat$prep.batch)$residuals }
colnames(batchcor) <- mzs

plot(dat[ , mzs[i]] ~ dat$prep.batch)

par(mfrow=c(2,2))
boxplot(mzdat, main='dat')
boxplot(combat, main='combat')
boxplot(t(mzdat), main='dat')
boxplot(t(combat), main='combat')

combat[1:4,1:4]
mzdat[1:4,1:4]

dat[ ,mzs] <- combat # replace mzdat in dat with batch cor

##############################################################################################################
# normalization finished
# save data
save(dat, mz.info, mzs, file='data/processed.data')
rm(list=ls())
##############################################################################################################

load('data/processed.data')

 
mzs[grepl('Ace', mzs)]
mzs[grepl('Acid', mzs)]
mzs[grepl('ate', mzs)]


write.table(mz.info, file='mz.info.csv', sep=',', row.names = F, quote=F)

dat[1:5,1:15]
mzdat <- as.matrix(dat[ ,mzs])
## check out pca:
levels(dat$treatment)
levels(dat$treatment) <- c('(-)', 'rapa')
dat$phenotype.treat <- as.factor(paste(dat$phenotype, dat$treatment))
levels(dat$phenotype.treat) 

table(dat$line, dat$phenotype)

colors <- c('lightblue', 'dodgerblue', 'pink' ,'red')
colors[dat$phenotype.treat]

pca <- prcomp(dat[ ,mzs], scale=T)
pairs(pca$x[ ,1:5], col=colors[dat$phenotype.treat], pch=19)

pcs <- pca$x[ ,1:6]
rownames(pcs) <- dat$sampleNumber
rownames(dat) <- dat$sampleNumber
pca <- as.data.frame(cbind(dat, pcs))

p <- prcomp(dat[ ,mzs], scale=T)
eigs <- p$sdev^2
propVar <- eigs / sum(eigs)
propVar

loadings <- prcomp(dat[ ,mzs], scale=T)$rotation
loadPC1 <- loadings[ ,1][rev(order(loadings[ ,1]))]
par(mfrow=c(1,2))
barplot(loadPC1, horiz=T, las=1, xlab='PC1 loading', col=ifelse(abs(loadPC1) > 0.1, 'red', 'lightgrey'), border=0,cex.names=0.5)

head(loadPC1) # top 6
tail(loadPC1) # bottom 6

ggplot(pca, aes(x=PC1, y=PC6, color = phenotype.treat)) +
  ylab(paste0('PC6 (', round(100*propVar[6],1), '%)')) +
  xlab(paste0('PC1 (', round(100*propVar[1],1), '%)')) +
  stat_ellipse(level=0.5, geom = "polygon", aes(fill = phenotype.treat), alpha = 0.25) +
  scale_fill_manual(values=c('lightblue', 'dodgerblue', 'pink' ,'red')) +
  geom_point() +
  scale_color_manual(values=c('lightblue', 'dodgerblue', 'pink' ,'red')) +
  theme_classic(base_size = 16 ) 


head(pca)

means <- pca %>% group_by(line, treatment, phenotype) %>% summarise_at(vars(c(PC1)), mean)   
head(means)

levels(means$treatment) <- c('control', 'rapamycin')
levels(means$phenotype) <- c('resistant', 'sensitive')
means$phenotype

ggplot(means, aes(y=PC1, x=treatment, color=treatment, group=line))+
  geom_line(color='grey40')+
  scale_color_manual(values=c('blue', 'red'))+
  geom_point(size=3)+
  xlab('')+
  ylab('metabolome PC1 (18.7% var)')+
  theme_bw(base_size = 16)+
  facet_wrap(~phenotype)

# pretty dramatic within-line effect of rapa in S lines, not R lines. Random effect (1|line) model will scream.


###########################################################
# clustering mz data:
###########################################################

par(mfrow=c(1,1))
par(mar=c(5,4,4,2))

# Create a vector giving a color for each extraction batch a sample belongs to
par(mar=c(5,8,4,2))

d <- dist(mzdat, method = "euclidean") # distance matrix
dend <- as.dendrogram(hclust(d, method="ward.D2"))
dend <- dend %>%  set("labels", dat$line[order.dendrogram(dend)]) # sets line names as tip lables
plot(dend, ylab='euclidian distance', las=1, main='clustering larva metabolome')
colored_bars(cbind(as.numeric(dat$treatment)+2, as.numeric(dat$phenotype)), dend, rowLabels = c('treatment', 'phenotype'))

par(mar=c(5,4,4,2))
tmp <- mzdat
rownames(tmp) <- dat$line
pheatmap(tmp, fontsize=5, scale='column')

# add group annotation
dat[1:5,1:14]
anno <- dat[ ,c('line', 'phenotype.treat')]
head(anno)
anno$line <- as.factor(anno$line)
str(anno)
rownames(anno) -> rownames(tmp) # trick to align annotation with rows

COLORS = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]

line = sample(COLORS, length(unique(anno$line)))
names(line) = levels(anno$line)
phenotype_treatment = brewer.pal(11, "Spectral")[c(7, 10, 5, 2)]
names(phenotype_treatment) <- levels(anno$phenotype.treat)
ann_colors = list(line = line, phenotype_treatment = phenotype_treatment)

pheatmap(t(tmp), annotation = anno, fontsize=7, annotation_colors = ann_colors, scale='row')

# here is clustering both by genotype and phenotype, little by treatment.  Since phenotype is confounded with genotype, I suspect that the clustering by P is partly an artifact of clustering by G  


rm(list=ls())
###########################################################
## univariate analysis
###########################################################
load('data/processed.data')

# get rapa effects for each phenotype group separately:
rapList <- list()
singular <- logical()

for(k in 1:2) {
  pheno <- levels(dat$phenotype)[k]
tmp <- dat[dat$phenotype == pheno,  ]
tmp$phenotype <- droplevels(tmp$phenotype)
tmprap <- matrix(nr=length(mzs), nc=2)

for(i in 1:length(mzs)) {
m <- lmer(tmp[ ,mzs[i]] ~ treatment + (1|line) , tmp)
singular[i] <- isSingular(m)
m <- summary(m)
m$coefficients
tmprap[i, ] <- m$coefficients[2, c(1, 5)] }

rownames(tmprap) = mzs
colnames(tmprap) <- c('beta', 'P')
tmprap <- as.data.frame(tmprap)
tmprap$FDR <- p.adjust(tmprap$P, 'BH')
tmprap$singular <- singular
rapList[[k]] <- tmprap }

head(tmprap)

names(rapList) <- levels(dat$phenotype)
lapply(rapList, function(x) table(x$FDR<=0.05)) # significant rapamycin effects per phenotype group
lapply(rapList, function(x) table(x$FDR<=0.05, singular)) # Q: any relationship between significance and singularity? A:No


# plot
library(ggrepel)

R_labes <- mzs
R_labes[!rapList[['R']]$FDR <=0.05] = NA
R_labes

ggplot(rapList[['R']], aes(-log(P, 10), beta, label = R_labes)) + 
  geom_point(color=ifelse(rapList[['R']]$FDR < 0.05, 4, 'grey66')) + 
  geom_text_repel(size=3) + 
  labs(title = "resistant lines") + 
  ylab('rapamycin effect (Beta)') + 
  theme_classic() +
  xlim(0, 6.5) +
  ylim(-0.38, 0.35) +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=12))

S_labes <- mzs
S_labes[!rapList[['S']]$FDR <=0.005] = NA
S_labes

S_labes[S_labes=='Leucine /D-Norleucine'] <- 'Leucine'
S_labes[S_labes=='iso-Leucine /allo-isoLeucine'] <- 'iso-Leucine'

ggplot(rapList[['S']], aes(beta, -log(P, 10), label = S_labes)) + 
  geom_point(color=ifelse(rapList[['S']]$FDR < 0.005, 2, 'grey66'), alpha=1) + 
  geom_text_repel(size=4, max.overlaps = 5) + 
  labs(title = "sensitive lines") + 
  ylab('rapamycin effect (Beta)') + 
  xlim(-0.6, 0.6)+
  theme_classic(base_size = 16) +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14))


mzs[rapList[['R']]$FDR<=0.05]
mzs[rapList[['S']]$FDR<=0.05]
intersect(mzs[rapList[['R']]$FDR<=0.05], mzs[rapList[['S']]$FDR<=0.05]) # all mzs with rapa effects in R and found among the rapa mzs in S

save(dat, rapList, file='univariate_analysis_results')
rm(list=ls())
####################################################################

# plot BetaRapa_R v. BetaRapa_S

load('univariate_analysis_results')

head(rapList[['S']])
tmp <- rapList[['S']]

rapList[['S']]$phenotype <- c('sensitive')
rapList[['R']]$phenotype <- c('resistant')
rapList[['S']]$metabolite <- rownames(rapList[['S']])
rapList[['R']]$metabolite <- rownames(rapList[['R']])
rapList[['S']]$significant <- rapList[['S']]$FDR <=0.05
rapList[['R']]$significant <- rapList[['R']]$FDR <=0.05

betaOrder <- rapList[['S']]$metabolite[order(rapList[['S']]$beta)]

head(rapList[['S']])

rapList <- lapply(rapList, function(x) {x[6] <- lapply(x[6], as.factor);x })
rapList[['R']]$metabolite <- factor(rapList[['R']]$metabolite, levels=betaOrder)
rapList[['S']]$metabolite <- factor(rapList[['S']]$metabolite, levels=betaOrder)
rapList[['S']]$metabolite # these should be in decreasing order of Beta-rapa among S lines

r <- do.call(rbind, rapList)

r$phenotype <- as.factor(r$phenotype)
r$phenotype <- factor(r$phenotype, levels=c('sensitive', 'resistant'))

head(r)
w <- r %>% pivot_wider(id_cols = -c(singular, P, significant), names_from= c(phenotype), values_from=c(beta, FDR))
head(w)

write.table(w, row.names = F, quote=F, file="/Users/ben/Library/CloudStorage/GoogleDrive-ben6@uw.edu/My Drive/Documents/rapa/paper/SUPPLEMENTARY TABLE 4 univariate metabolome.csv", sep=',')

w$effect <- paste(w$FDR_sensitive<=0.05, w$FDR_resistant<=0.05)
table(w$effect, w$beta_sensitive>0)

ggplot(w, aes(x=beta_resistant, y=beta_sensitive, label = metabolite, color = effect))+
  geom_vline(xintercept=0, color='grey40', alpha=0.5)+
  geom_hline(yintercept=0, color='grey40', alpha=0.5)+
  geom_point()+
  geom_text_repel(size=3, max.overlaps = 9) + 
  scale_colour_manual(name = 'effect', values = setNames(c('purple','red', 'blue', 'grey'), c("TRUE TRUE", "TRUE FALSE", "FALSE TRUE" ,"FALSE FALSE"))) +
  theme_bw(base_size = 16)+
  theme(legend.position="none") +
  ylab(expression("sensitive " (beta)))+
  xlab(expression("resistant " (beta)))+
  ggtitle('effects of treatment')

head(w)
cor.test(w$beta_resistant, w$beta_sensitive, method='spear') # sign of rapa effect is correlated among S and R


rm(list=ls())
###########################################################################################################
# Enrichment analysis: try FELLA
###########################################################################################################
load('data/processed.data')
load('univariate_analysis_results')
mzdat <- as.matrix(dat[ ,mzs])

suppressPackageStartupMessages({
  library(igraph)
  library(FELLA)
  library(org.Dm.eg.db)
  library(KEGGREST)
  library(magrittr)
  library(biomaRt)
  library(rvest)
  library(XML)
  library(RCy3)
})

set.seed(1)

# get KEGG id for metabolites with significant rapa effects within R and S lines
head(mz.info)

Rmzs <- mzs[rapList[['R']]$FDR <= 0.05]
Rkegg <- mz.info$KEGG.ID[mz.info$Current.MS.Compounds %in% Rmzs]
Smzs <- mzs[rapList[['S']]$FDR <= 0.05]
Skegg <- mz.info$KEGG.ID[mz.info$Current.MS.Compounds %in% Smzs]
keggList <- list(Rkegg, Skegg)
names(keggList) <- c('R_mzs', 'S_mzs')
keggList

keggList <- lapply(keggList, function(x) x[!grepl('N/A', x)]) # remove NA kegg IDs
keggList <- lapply(keggList, function(x) sapply(strsplit(x, "/"), '[[', 1)) # take first ID of multiplue KEGG IDs
keggList

# get compounds from the analysis
mz.info <- mz.info[mz.info$Current.MS.Compounds %in% colnames(mzdat), ]
all_compounds <- mz.info$KEGG.ID
all_compounds <- all_compounds[!grepl('N/A', all_compounds)]
all_compounds <- sapply(strsplit(all_compounds, "/"), '[[', 1)
all_compounds

tmpdir <- "/Users/ben/Library/CloudStorage/GoogleDrive-ben6@uw.edu/My Drive/Documents/rapa/metabolomics/larva_metabolomics/FELLA_database"

#############################################################################################################
# [IF YOU ALREADY RAN THE LINES HERE, SKIP TO THE LINE: fella.data <- loadKEGGdata ]
graph <- buildGraphFromKEGGREST(organism = "dme")  # Build network; don't bother filtering large pathways out, they don't enrich unless they survive permutation:
unlink(tmpdir, recursive = TRUE) # this line allows you to rewrite the: buildataFromGraph()
buildDataFromGraph(graph, databaseDir = tmpdir, internalDir = F, matrices = 'diffusion', normality = "diffusion", niter = 1000) # iterate over the graph 1000 times using: 'diffusion'
#############################################################################################################

fella.data <- loadKEGGdata(databaseDir = tmpdir, internalDir = F, loadMatrix = 'diffusion')
getInfo(fella.data) # prints the version KEGG used (Sergio says it always uses the most recent release)
fella.data

## enrichment analysis
enrichmentAnalysisList <- list()
enrichmentTableList <- list()

enrichmentAnalysisList <- lapply(keggList, function(x) enrich(compounds = x, compoundsBackground = all_compounds, data = fella.data, method = 'diffusion', approx = "simulation", niter=10000))
enrichmentAnalysisList

###
enrichmentTableList <- lapply(enrichmentAnalysisList, function(x) generateResultsTable(object=x, data=fella.data))
lapply(enrichmentTableList, head)

FDRList <- lapply(enrichmentTableList, function(x) p.adjust(x$p.score, 'BH'))
enrichmentTableList <- Map(cbind, enrichmentTableList, FDR = FDRList)
lapply(enrichmentTableList, head) 

KEGG.nameList <- lapply(enrichmentTableList, function(x) gsub(',', '_', x$KEGG.name)) # removing commas to ease data wrangling
enrichmentTableList <- lapply(enrichmentTableList, function(x) { x["KEGG.name"] <- NULL; x })
enrichmentTableList <- Map(cbind, enrichmentTableList, KEGG.name = KEGG.nameList)


save(enrichmentTableList, keggList, file='FELLA_results/enrichmentTableList')



rm(list=ls())
################################################
load('FELLA_results/enrichmentTableList')

# look at only pathways and modules
pathwaysList <- lapply(enrichmentTableList, subset, Entry.type %in% c('pathway', 'module'))
pathwaysList


# make a table
################################################################################################
dir('FELLA_results')

# [NOTE]: open enrichment tables in excel, find and replace ',' with '_', then save as .csv, then you can open with:
lapply(enrichmentTableList, head)

sigPs <- lapply(enrichmentTableList, function(x) signif(x[ ,3], 3)) # limit to sig figs:
sigFDRs <- lapply(enrichmentTableList, function(x) signif(x[ ,4], 3))
simpKEGGs <- lapply(enrichmentTableList, function(x) sapply(str_split(x$KEGG.name, "\\ -"), '[[', 1)) # simplify kegg names
head(enrichmentTableList[[1]])

enrichmentTableList <- lapply(enrichmentTableList, function(x) { x[c("p.score", "FDR", "KEGG.name")] <- NULL; x })
enrichmentTableList <- Map(cbind, enrichmentTableList, KEGG.name = simpKEGGs, p.score = sigPs, FDR = sigFDRs)
lapply(enrichmentTableList, head)

dir()
lapply(enrichmentTableList, function(x) DT::datatable(x)) # save an HTML table


library(gt)
library(gapminder)

for(i in 1:2) {
  enrichmentTableList[[i]] %>% arrange(p.score) %>% filter(Entry.type %in% c('pathway', 'module')) %>% gt(groupname_col = "Entry.type", rowname_col = "KEGG.name") %>% tab_header(title = "Pathway Enrichment") %>% cols_align(align = "center") %>% row_group_order(groups = c('pathway')) }
names(enrichmentTableList)
x <- bind_rows(enrichmentTableList, .id = "group")
x$phenotype <- ifelse(grepl('R', x$group), 'resistant', 'sensitive')

head(x)
getwd()
write.table(x, row.names=F, quote=F, file='/Users/ben/Library/CloudStorage/GoogleDrive-ben6@uw.edu/My Drive/Documents/rapa/paper/SUPPLEMENTARY TABLE 5.csv', sep=',')

x <- x[x$Entry.type %in% c('pathway', 'module'), ] # limit tables to only module and pathways

head(x)
x[ ,1:6] %>% arrange(p.score) %>% gt(groupname_col = "group", rowname_col = "KEGG.name") %>% tab_header(title = "Table 4. Pathway Enrichment by Rapamycin Effect") %>% cols_align(align = "center") %>% tab_style(style = list(cell_text(weight = "bold")), locations = cells_row_groups())
x


