############################################################
# QGG package for GBLUP and Genome Feature Models 
############################################################

if (!requireNamespace("BiocManager", quietly = TRUE)) 
  install.packages("BiocManager")

# Remove hashes and excecute to install packages.
# BiocManager::install("org.Dm.eg.db") 
# BiocManager::install("AnnotationDbi")
# BiocManager::install("GOSim")
# BiocManager::install("KEGGREST")

library(org.Dm.eg.db)            
library(AnnotationDbi)
library(car) 
library(qgg)
library(parallel)
library(qqman)
library(KEGGREST)
library(GOSim)
library(stringr)
library(msm)
library(rlang) 
library(tidyr)
library(DT)

rm(list=ls())

# refer to the tutorials here:
# http://psoerensen.github.io/qgg/articles/gsea.html

# for use in multi-core processing
n.cores <- detectCores(all.tests = F, logical = T) # detect the number of available CPU cores

setwd("/Users/ben/Library/CloudStorage/GoogleDrive-ben6@uw.edu/My Drive/Documents/rapa/GWAS/GBLUP")
dir()

pheno <- read.table("GWASinput_FINAL.txt", header=T) # should be 140 DGRP lines
head(pheno)
covars <- read.table("covars.txt", header = T)
head(covars)

################################################################################
# to reduce DGRP snps to those that pass MAF and genotype rate thresholds, etc:

############################################################
# load genotype data (from processing steps done in 'QGG_GFsets.R')
############################################################

## see note above, you need to complete the work in 'QGG_GFsets.R' prior to running the code below:

load('W.maf05.geno20')
dim(W)
W[1:4,1:4]

# snps already processed in PLINK, first limited by MAF <5%, geno >95%, then LD-pruned: --indep-pairwise 200 5 0.8 
LDsnps.maf5.geno5 <- unlist(read.table("/Users/ben/Library/CloudStorage/GoogleDrive-ben6@uw.edu/My Drive/Documents/rapa/GWAS/plink/rapaLD.maf5.geno5.prune.in")) 

GRM.5.5 <- grm(W=W[ ,colnames(W) %in% LDsnps.maf5.geno5]) # create GRM from pruned markers present in W
rm(LDsnps.maf5.geno5)

############################################################
# fit GBLUP
############################################################
dat <- merge(pheno, covars)
head(dat)
colnames(dat)[3] <- 'pheno'

# make a design matrix of the covariates (ignore genotype PCs, the grm will be used in qgg)
fm <- pheno ~  wolbachia + In.2L.t + In.2R.NS + In.3R.P + In.3R.K + In.3R.Mo
X <- model.matrix(fm, data=dat)

pheno <- dat$pheno
names(pheno) <- dat$FID

################################################
# fit GBLUPs to data:
################################################
fit <- greml(y=pheno, X=X, GRM=list(GRM.5.5), verbose=T) #greml fit required for later associations: this is where the covarates and realtedness are accounted for in the model.

GBLUPs <- fit$Py # these are the 'genetic effects', the residual phenotype after the greml fit of y ~ Xb + GRM

#########################################################
# a few methods to evaluate 'snp heritability':
#########################################################

#########################################################
# plotting
plot(GBLUPs, pheno, las=1, pch=16, ylab='rapamycin sensitivity', cex=0.8)

# the delta method, gives Hsnp for all data, along with estimate of error
theta <- fit$theta
covar <- fit$asd
rownames(covar) <- colnames(covar) <- names(theta)
heritability_summary <- deltaMethod(theta, "G1/(G1+E)", vcov=covar)
heritability_summary

plot(GBLUPs, pheno, las=1, pch=16, ylab='rapamycin sensitivity', cex=0.8)
legend('bottomright', legend=paste0('G1/(G1 + E)', ' =', round(heritability_summary$Estimate, 3), ' +/-', round(heritability_summary$SE, 2)), bty='n')
#########################################################

# the wide variance (se) of Hsnp may be due to small number of lines
# see how sensitive estiamte is to held-out lines:

allboots <- list()
subsample <- c(139, 135, 130, 120, 100)

for(k in 1:length(subsample)) {
  subsize <- subsample[k]

heritability_boot <- list()
for(i in 1:200) {
  tryCatch({
  index <- sample(1:nrow(X), subsize, replace=F)
boot <- greml(y=pheno[ index], X=X[ index, ], GRM=list(GRM.5.5[ index, index]), verbose=F) 
theta <- boot$theta
covar <- boot$asd
rownames(covar) <- colnames(covar) <- names(theta)
heritability_boot[[i]] <- deltaMethod(theta, "G1/(G1+E)", vcov=covar)  }, error=function(e){})  # empty error function added to pass non-convergent iterations
}
allboots[[k]] <- heritability_boot }

names(allboots) <- subsample


x <- lapply(allboots, function(x) do.call(rbind, x))
x
boots <- do.call(rbind, x)
head(boots)
boots$subsample <- as.factor(as.numeric(sapply(strsplit(rownames(boots), split='\\.'), '[[', 1)))
head(boots)
heritability_summary$subsample <- 140
boots <- rbind(heritability_summary, boots)

ggplot(boots, aes(y=Estimate, x=subsample))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width=0.1, size=1)+
  theme_bw(base_size = 16)+
  labs(y=expression(H["snp"]))+
  xlab('subsample (n lines)')

save(boots, file='heritability.bootstrap.results')

##################################################################
# permutation testing, permute phenotype across covaraites and GRM:
##################################################################
heritability_perms <- list()
p <- pheno

for(i in 1:1000) {
  tryCatch({
  index <- sample(1:140)
  perm <- greml(y=p[index], X=X, GRM=list(GRM.5.5), verbose=F, maxit = 500) 
  theta <- perm$theta
  covar <- perm$asd
  rownames(covar) <- colnames(covar) <- names(theta)
  heritability_perms[[i]] <- deltaMethod(theta, "G1/(G1+E)", vcov=covar) }, error=function(e){}) 
  }

perms <- do.call(rbind, heritability_perms)
table(rowSums(is.na(perms)) ==0) # ____ of the n perms converged
perms <-perms[rowSums(is.na(perms)) ==0, ]
head(perms)

hist(perms$Estimate, 20, border=0, col='grey40', las=1, xlab='gBLUP Hsnp', ylab='permutations', main='')
plot(Estimate ~ SE, perms, pch=20, xlim=c(0, 1))

table(heritability_summary$Estimate <= perms$Estimate)/1000 # empirical P


p1 <- ggplot(boots, aes(y=Estimate, x=subsample))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width=0.1, size=1)+
  theme_bw(base_size = 16)+
  labs(y=expression(H["snp"]))+
  xlab('subsample (n lines)')

p2 <- ggplot(perms) + 
  geom_histogram(binwidth=0.05, aes(x=Estimate), colour="Black")+
  theme_bw(base_size = 16)+
  labs(x=expression(H["snp"]))+
  ylab('frequency')

ggarrange(p1, p2)

save(perms, file='heritability.permutation.results')


#########################################################
# cross-validation method (based on scheme on ?greml(); help page for greml)

sets <- list(G1 = colnames(W)) # make a list of all snps 
GB <- lapply(sets, function(x) {grm(W = W[, x])}) # make a grm of all snps (in a form that greml likes)
dim(GB[[1]]) # 140 x 140

rownames(X) <- names(pheno)

CVfolds <- c(2, 4, 6, 8, 16, 32)
ValidationCorrMat <- matrix(nr=50, nc=length(CVfolds))
testsize <-  print(round(140/CVfolds, 0))
traintestsize <- paste0(140-testsize, ':', testsize)

set.seed(2)

for(i in 1:length(testsize)){
validate <- replicate(50, sample(1:length(pheno), testsize[i]))
cvG <- greml(y=pheno, X=X, GRM=GB, validate = validate) # this function performs 'training' (calc gBLUPs) on the training data (all data - 'validate' sets, then gives predictions in the validate sets)
ValidationCorrMat[ ,i] <- cvG$accuracy$Corr  }

colnames(ValidationCorrMat) <- paste0("[", traintestsize , "]")
boxplot(ValidationCorrMat, las=1,ylab=expression('cor(obs~predicted)'), xlab='CV [train:test] (n lines)', outline=F, col=0, ylim=c(-1,1))
out <- reshape2::melt(ValidationCorrMat)
head(out)
colnames(out) <- c('iteration', 'CV[Train:Test]', 'cor')
out$`CV[Train:Test]` <- as.factor(out$`CV[Train:Test]`)
stripchart(cor ~ `CV[Train:Test]`, out, add=T, vertical=T, pch=16, method='jitter', cex=0.8)
# these genomic predictions fail miserably on held-out data.

head(out)

ggplot(out, aes(y=cor, x=`CV[Train:Test]`))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width=0.1)+
  theme_bw(base_size = 18)+
  xlab('CV [train:test] (n lines)')+
  ylab(expression('cor(obs~predicted)'))
  
HsnpDf <- data.frame('GBLUP'=GBLUPs, 'pheno'=pheno)

ggplot(HsnpDf, aes(y=GBLUPs, x=pheno))+
  geom_point()+
  theme_bw(base_size = 16)+
  xlab('rapamycin sensitivity')


p1 <- ggplot(HsnpDf, aes(y=GBLUPs, x=pheno))+
  geom_point(size=1)+
  theme_classic(base_size = 16)+
  xlab('rapamycin sensitivity')

p2 <- ggplot(out, aes(y=cor, x=`CV[Train:Test]`, fill=`CV[Train:Test]`)) +
  geom_hline(yintercept=0, color='grey30', linetype='dashed')+
  geom_boxplot(outlier.shape = NA, alpha=0.4)+
  geom_jitter(width=0.1, size=0.6)+
  theme_classic(base_size = 16)+
  xlab('')+
  ylab(expression('cor(obs~predicted)'))+
  theme(axis.text.x = element_text(size=8))

ggarrange(p1, p2, widths=c(0.6, 1))


#########################################################

save(GRM.5.5, X, pheno, fit, file='data.for.CVAT.analysis')
rm(list=ls())

############################################################
# single marker effects:
############################################################
load('data.for.CVAT.analysis')
load('W.maf05.geno20')

fit <- greml(y=pheno, X=X, GRM=list(GRM.5.5), verbose=T) 

SNPlevel <- lma(fit=fit, W=W) # [NOTE] uses a statistic based on the 'MASTOR' statistic.
SNPlevel <- as.data.frame(SNPlevel)
head(SNPlevel)
SNPlevel$FDR <- p.adjust(SNPlevel$p, method='BH')
SNPlevel[order(SNPlevel$p), ][1:10, ] # view 'top hits'

############################################################
# try clumping in plink to reduce SNPs to more independent clumps (i.e haplotypes) that are represented by lead SNPS (most-associated in the GWAS):
head(SNPlevel)
clump.these <- data.frame('SNP' = rownames(SNPlevel), 'P'=SNPlevel[ ,4])
head(clump.these)
write.table(clump.these, file='/Users/benharrison/Google Drive/Documents/rapa/GWAS/plink/clump.these.in.plink', row.names=F, quote=F)

clumps <- read.table('/Users/benharrison/Google Drive/Documents/rapa/GWAS/plink/plink.clumped', head=T, stringsAsFactors = F) 
head(clumps) # these might be helpful in adding annotations, to manhattan plot perhaps
############################################################


# qq(SNPlevel$p, las=1, main='Single Marker Rapamycin Sensitivity') # make a Q-Q plot

BP <- sapply(strsplit(rownames(SNPlevel),  split='_'), "[[", 2)
CHR <- sapply(strsplit(rownames(SNPlevel),  split='_'), "[[", 1)
CHR <- as.factor(CHR)
CHR <- factor(CHR, levels=c('X', '2L', '2R', '3L', '3R', '4'))

SNPlevel$CHR <- as.numeric(CHR)
SNPlevel$chrom <- CHR
SNPlevel$BP <- as.numeric(BP)
SNPlevel$SNP <- rownames(SNPlevel)

rm(CHR)
rm(BP)

manhattan(SNPlevel, snp='SNP', chr='CHR', bp='BP', p='p', col=c('cyan2' ,'coral'), chrlabs = levels(CHR), suggestiveline = NULL, genomewideline = NULL, ylim=c(-log10(5e-3), -log10(2e-6)), annotatePval = 0.01)

?manhattan

## add annotations to markers
load("rapaA_FB5.57")
rapaA[1:12, ]
tmp <- rownames(rapaA)
x <- as.data.frame(cbind(rapaA, tmp))
x[1:12, ]
rownames(x) =NULL
colnames(x)[5] <- 'SNP'
head(SNPlevel)
head(x)
SNPlevel <- merge(SNPlevel, x)
head(SNPlevel)

rm(x)
rm(rapaA)

# load genome coordinates
map <- read.table('gene_map_table_fb_2022_02.tsv', header=F, fill=T, sep='\t', stringsAsFactors = F)
head(map)
colnames(map) <- c('species', 'GeneName', 'FBid', 'recombMap', 'cytoMap', 'seq location')
head(map)
map <- map[ ,c(3, 6)]
head(map)
table(map=='') # blank entries for sequence locations 

map <- map[!map$`seq location`=='', ]
x <- strsplit(map$`seq location`, split=':')
head(x)
locs <- sapply(x, "[[", 2)
head(locs)
strand <- ifelse(grepl('-1', locs), 'minus', 'plus')

locs <- sapply(strsplit(locs, split='\\('), "[[", 1)
head(locs)
left <- as.numeric(sapply(strsplit(locs, split='\\..'), "[[", 1))
right <- as.numeric(sapply(strsplit(locs, split='\\..'), "[[", 2))
coordinates <- data.frame('FBid'=map$FBid, 'strand'=strand, 'left.end'=left, 'right.end'=right)
head(coordinates)

rm(locs)

SNPlevel <- merge(SNPlevel, coordinates, by='FBid', all.x=T, all.y=F)
SNPlevel[order(SNPlevel$p), ][1:20, ] # view 'top hits'

write.table(SNPlevel, file='CVAToutput/SNPlevel', quote=F)


################################################################################
### read SNPlevel results:
################################################################################
SNPlevel <- read.table('CVAToutput/SNPlevel', stringsAsFactors = F, fill=T)
SNPlevel[order(SNPlevel$p), ][1:10, ]

## add descriptions of genes from FlyBase 
snapshots <- read.table('gene_snapshots_fb_2022_02.txt', header=T, fill=T, sep='\t')
head(snapshots)
colnames(snapshots)
colnames(snapshots)[c(1, 3)] <- c('FBid', 'Name')

geneSummaries <- read.table('automated_gene_summaries.tsv', header=F, fill=T, sep='\t') 
head(geneSummaries)
colnames(geneSummaries) <- c('FBid', 'geneSummary')

head(SNPlevel)

table(unique(SNPlevel$FBid) %in% geneSummaries$FBid)
notInSummaries <- unique(SNPlevel$FBid)[!unique(SNPlevel$FBid) %in% geneSummaries$FBid]
table(unique(SNPlevel$FBid) %in% snapshots$FBid)
table(notInSummaries %in% snapshots$FBid)
table(unique(SNPlevel$FBid) %in% geneSummaries$FBid | unique(SNPlevel$FBid) %in% snapshots$FBid)
rm(notInSummaries)
SNPlevel <- merge(SNPlevel, snapshots[ ,c(1,3,5)], by='FBid', all.x=T, all.y=F)
SNPlevel <- merge(SNPlevel, geneSummaries, by='FBid', all.x=T, all.y=F)
head(SNPlevel)



rm(geneSummaries)
rm(coordinates)
rm(snapshots)

# search for a delimiter to save the data with (should not be present in the data as-is)
table(grepl(',', SNPlevel$geneSummary))  # commas won't worjk
table(grepl(',', SNPlevel$gene_snapshot_text))
 
table(grepl('_', SNPlevel$geneSummary)) # underscores won't work 
table(grepl('_', SNPlevel$gene_snapshot_text))

table(grepl(':', SNPlevel$geneSummary)) # colon won't work 
table(grepl(':', SNPlevel$gene_snapshot_text))

table(grepl(';', SNPlevel$geneSummary)) # colon won't work 
table(grepl(':', SNPlevel$gene_snapshot_text))



# plot an example association:
################################################
GBLUPs <- fit$Py # these are the 'genetic effects', the residual phenotype after the greml fit of y ~ Xb + GRM

# confirm that Py are the 'genetic effects' after removing fixed effects of wolb+INVs and random effects of GRM
summary(aov(fit$Py ~ wolbachia + In.2L.t + In.2R.NS + In.3R.P + In.3R.K + In.3R.Mo, dat)) # the fixed effects are gone
fit2 <- greml(y=GBLUPs, X=X, GRM=list(GRM.5.5), verbose=T) 
fit2$theta # and the genetic effects are gone
rm(fit2)
################################################

SNPlevel[order(SNPlevel$p), ][1:10, ] # view 'top hits'
topHit <- SNPlevel$SNP[SNPlevel$p == min(SNPlevel$p)]

topHit <- '3L_12790278_SNP' # chose your own top hit.

top10 <- SNPlevel$SNP[order(SNPlevel$p)][1:10]

par(mfrow=c(2, 5))
for(i in 1:length(top10)) {
hit <- top10[i]
hit.genotypes <- as.factor(W[ ,hit])
levels(hit.genotypes)[c(1, 3)] <- c('major', 'minor')
tmpGeno <- droplevels(hit.genotypes[hit.genotypes != 0])
tmpPheno <- GBLUPs[hit.genotypes != 0]
plot(tmpPheno  ~ tmpGeno, las=1, main=hit, xlab=NULL, ylab='GBLUP', outline=F, ylim=range(GBLUPs))
stripchart(tmpPheno  ~ tmpGeno, vertical=T, add=T, method='jitter', pch=19, col='grey20') }

save(dat, fit, W, pheno, file = 'GBLUPobjects_for_CVAT')

rm(list=ls())


############################################################

############################################################
# Genome Feature Analysis:
############################################################
# Genome features are any lists of markers assoc. with a grouping, ie:
# Genes
# GO terms
# KEGG pathways
# others? i.e. genes-coexpressed with TOR?

# use the cvat method (covariance association test, Rhode et al., 2016).  gsea( , method='cvat').  The gsea function will run nperms to assess signifcance of the cvat score 'setT'.  Default is nperm=1000, which takes a while.  Try step-wise permutation to probe the small P from the prior round(s) of permutations. 

############################################################
# Genes
############################################################
load('GBLUPobjects_for_CVAT')

# collect garbage
gc()
dir()
load(file="rapa_FB5.57.fbSets_1kbp")
round1.Sets <- fbSets1kbp 
rm(fbSets1kbp)
############################################################
set.name <- 'GENElevel_1kb' #user defined output file name based on the nature of the genome feature set
output.path <- paste0('./CVAToutput/', set.name, '.cvat')

# 4 rounds of permutation, successive deeper permutations are used to resolve the lowest P values
round1.cvat <- gsea(fit=fit, sets=round1.Sets, W=W, method = "cvat", nperm=1000, ncores = n.cores)
rownames(round1.cvat)[round1.cvat$p<0.01] -> pass_to_10K # pull sets below a P threshold 
round2.Sets <- round1.Sets[pass_to_10K]
print(paste('round 1 complete,', length(pass_to_10K), 'features sent to round 2'))

round2.cvat <- gsea(fit=fit, sets=round2.Sets, W=W, method = "cvat", nperm=10000, ncores = n.cores)
rownames(round2.cvat)[round2.cvat$p<0.001] -> pass_to_100K # pull sets below a P threshold
round3.Sets <- round1.Sets[pass_to_100K]
rm(round2.Sets)
print(paste('round 2 complete,', length(pass_to_100K), 'features sent to round 3'))

round3.cvat <- gsea(fit=fit, sets=round3.Sets, W=W, method = "cvat", nperm=100000, ncores = n.cores)
rownames(round3.cvat)[round3.cvat$p<0.0001] -> pass_to_1M # pull sets below a P threshold
round4.Sets <- round1.Sets[pass_to_1M]
rm(round3.Sets)
print(paste('round 3 complete,', length(pass_to_1M), 'features sent to round 4'))

# combine all rounds:
if(length(pass_to_1M) > 0) { round4.cvat <- gsea(fit=fit, sets=round4.Sets, W=W, method = "cvat", nperm=1000000, ncores = n.cores)
final.cvat <- rbind(round1.cvat[!rownames(round1.cvat) %in% rownames(round2.cvat), ], round2.cvat[!rownames(round2.cvat) %in% rownames(round3.cvat), ], round3.cvat[!rownames(round3.cvat) %in% rownames(round4.cvat), ], round4.cvat) } else {
  final.cvat <- rbind(round1.cvat[!rownames(round1.cvat) %in% rownames(round2.cvat), ], round2.cvat[!rownames(round2.cvat) %in% rownames(round3.cvat), ], round3.cvat)}

final.cvat$FDR <- p.adjust(final.cvat$p, method='BH') # estimate FDR

# add gene symbol:
final.cvat$symbol <- mapIds(org.Dm.eg.db, 
                            keys=row.names(final.cvat), 
                            column="SYMBOL", 
                            keytype="FLYBASE",
                            multiVals="first")

write.table(final.cvat, file=output.path) # save output to working directory
#################

# clean up priot to million permutation round
rm(round2.cvat)
rm(round2.Sets)
rm(round3.cvat)
rm(round3.Sets)
rm(round1.cvat)
gc()


############################################################
# +/-0 gene sets defined by the primary transcript coordinates:
load(file="rapa_FB5.57.fbSets_0bp")
round1.Sets <- fbSets0bp
rm(fbSets0bp)
set.name <- 'GENElevel_0kb' #user defined output file name based on the nature of the genome feature set
############################################################

output.path <- paste0('./CVAToutput/', set.name, '.cvat')

# 4 rounds of permutation, successive deeper permutations are used to resolve the lowest P values
round1.cvat <- gsea(fit=fit, sets=round1.Sets, W=W, method = "cvat", nperm=1000, ncores = n.cores)
rownames(round1.cvat)[round1.cvat$p<0.01] -> pass_to_10K # pull sets below a P threshold 
round2.Sets <- round1.Sets[pass_to_10K]
print(paste('round 1 complete,', length(pass_to_10K), 'features sent to round 2'))

round2.cvat <- gsea(fit=fit, sets=round2.Sets, W=W, method = "cvat", nperm=10000, ncores = n.cores)
rownames(round2.cvat)[round2.cvat$p<0.001] -> pass_to_100K # pull sets below a P threshold
round3.Sets <- round1.Sets[pass_to_100K]
rm(round2.Sets)
print(paste('round 2 complete,', length(pass_to_100K), 'features sent to round 3'))

round3.cvat <- gsea(fit=fit, sets=round3.Sets, W=W, method = "cvat", nperm=100000, ncores = n.cores)
rownames(round3.cvat)[round3.cvat$p<0.0001] -> pass_to_1M # pull sets below a P threshold
round4.Sets <- round1.Sets[pass_to_1M]
rm(round3.Sets)
print(paste('round 3 complete,', length(pass_to_1M), 'features sent to round 4'))

# combine all rounds:
if(length(pass_to_1M) > 0) { round4.cvat <- gsea(fit=fit, sets=round4.Sets, W=W, method = "cvat", nperm=1000000, ncores = n.cores)
final.cvat <- rbind(round1.cvat[!rownames(round1.cvat) %in% rownames(round2.cvat), ], round2.cvat[!rownames(round2.cvat) %in% rownames(round3.cvat), ], round3.cvat[!rownames(round3.cvat) %in% rownames(round4.cvat), ], round4.cvat) } else {
  final.cvat <- rbind(round1.cvat[!rownames(round1.cvat) %in% rownames(round2.cvat), ], round2.cvat[!rownames(round2.cvat) %in% rownames(round3.cvat), ], round3.cvat)}

final.cvat$FDR <- p.adjust(final.cvat$p, method='BH') # estimate FDR

# add gene symbol:
final.cvat$symbol <- mapIds(org.Dm.eg.db, 
                            keys=row.names(final.cvat), 
                            column="SYMBOL", 
                            keytype="FLYBASE",
                            multiVals="first")

write.table(final.cvat, file=output.path) # save output to working directory
#################

# clean up priot to million permutation round
rm(round2.cvat)
rm(round2.Sets)
rm(round3.cvat)
rm(round3.Sets)
rm(round1.cvat)
rm(round4.Sets)
rm(round4.cvat)
rm(round1.Sets)
gc()


############################################################
# KEGG pathways
############################################################
# ensure that working directory has a folder named 'qgg'
# define genome feature set
load(file="rapa.keggSets_1kb")
round1.Sets <- keggSets_1kb
rm(keggSets_1kb)
set.name <- 'KEGGlevel_1kb' #user defined output file name based on the nature of the genome feature set
############################################################

output.path <- paste0('./CVAToutput/', set.name, '.cvat')

# 4 rounds of permutation, successive deeper permutations are used to resolve the lowest P values
round1.cvat <- gsea(fit=fit, sets=round1.Sets, W=W, method = "cvat", nperm=1000, ncores = n.cores)
rownames(round1.cvat)[round1.cvat$p<0.01] -> pass_to_10K # pull sets below a P threshold 
round2.Sets <- round1.Sets[pass_to_10K]
print(paste('round 1 complete,', length(pass_to_10K), 'features sent to round 2'))

round2.cvat <- gsea(fit=fit, sets=round2.Sets, W=W, method = "cvat", nperm=10000, ncores = n.cores)
rownames(round2.cvat)[round2.cvat$p<0.001] -> pass_to_100K # pull sets below a P threshold
round3.Sets <- round1.Sets[pass_to_100K]
rm(round2.Sets)
print(paste('round 2 complete,', length(pass_to_100K), 'features sent to round 3'))

round3.cvat <- gsea(fit=fit, sets=round3.Sets, W=W, method = "cvat", nperm=100000, ncores = n.cores)
rownames(round3.cvat)[round3.cvat$p<0.0001] -> pass_to_1M # pull sets below a P threshold
round4.Sets <- round1.Sets[pass_to_1M]
rm(round3.Sets)
print(paste('round 3 complete,', length(pass_to_1M), 'features sent to round 4'))

# combine all rounds:
if(length(pass_to_1M) > 0) { round4.cvat <- gsea(fit=fit, sets=round4.Sets, W=W, method = "cvat", nperm=1000000, ncores = n.cores)
final.cvat <- rbind(round1.cvat[!rownames(round1.cvat) %in% rownames(round2.cvat), ], round2.cvat[!rownames(round2.cvat) %in% rownames(round3.cvat), ], round3.cvat[!rownames(round3.cvat) %in% rownames(round4.cvat), ], round4.cvat) } else {
  final.cvat <- rbind(round1.cvat[!rownames(round1.cvat) %in% rownames(round2.cvat), ], round2.cvat[!rownames(round2.cvat) %in% rownames(round3.cvat), ], round3.cvat)}

final.cvat$FDR <- p.adjust(final.cvat$p, method='BH') # estimate FDR

columns(org.Dm.eg.db)

# add gene symbol:
final.cvat$KEGG <- mapIds(org.Dm.eg.db, 
                            keys=row.names(final.cvat), 
                            column="PATH", 
                            keytype="PATH",
                            multiVals="first")

# clean up
rm(round2.cvat)
rm(round2.Sets)
rm(round3.cvat)
rm(round3.Sets)
rm(round1.cvat)
rm(round4.Sets)
rm(round4.cvat)
rm(round1.Sets)
gc()


tmp <- paste0('path:map', rownames(final.cvat))
final.cvat$pathway <- keggList("pathway")[tmp]

write.table(final.cvat, file=output.path) # save output to working directory

############################################################
load(file="rapa.keggSets_0bp")
round1.Sets <- keggSets_0bp
rm(keggSets_0bp)
set.name <- 'KEGGlevel_0bp' #user defined output file name based on the nature of the genome feature set
############################################################

output.path <- paste0('./CVAToutput/', set.name, '.cvat')

# 4 rounds of permutation, successive deeper permutations are used to resolve the lowest P values
round1.cvat <- gsea(fit=fit, sets=round1.Sets, W=W, method = "cvat", nperm=1000, ncores = n.cores)
rownames(round1.cvat)[round1.cvat$p<0.01] -> pass_to_10K # pull sets below a P threshold 
round2.Sets <- round1.Sets[pass_to_10K]
print(paste('round 1 complete,', length(pass_to_10K), 'features sent to round 2'))

round2.cvat <- gsea(fit=fit, sets=round2.Sets, W=W, method = "cvat", nperm=10000, ncores = n.cores)
rownames(round2.cvat)[round2.cvat$p<0.001] -> pass_to_100K # pull sets below a P threshold
round3.Sets <- round1.Sets[pass_to_100K]
rm(round2.Sets)
print(paste('round 2 complete,', length(pass_to_100K), 'features sent to round 3'))

round3.cvat <- gsea(fit=fit, sets=round3.Sets, W=W, method = "cvat", nperm=100000, ncores = n.cores)
rownames(round3.cvat)[round3.cvat$p<0.0001] -> pass_to_1M # pull sets below a P threshold
round4.Sets <- round1.Sets[pass_to_1M]
rm(round3.Sets)
print(paste('round 3 complete,', length(pass_to_1M), 'features sent to round 4'))

# combine all rounds:
if(length(pass_to_1M) > 0) { round4.cvat <- gsea(fit=fit, sets=round4.Sets, W=W, method = "cvat", nperm=1000000, ncores = n.cores)
final.cvat <- rbind(round1.cvat[!rownames(round1.cvat) %in% rownames(round2.cvat), ], round2.cvat[!rownames(round2.cvat) %in% rownames(round3.cvat), ], round3.cvat[!rownames(round3.cvat) %in% rownames(round4.cvat), ], round4.cvat) } else {
  final.cvat <- rbind(round1.cvat[!rownames(round1.cvat) %in% rownames(round2.cvat), ], round2.cvat[!rownames(round2.cvat) %in% rownames(round3.cvat), ], round3.cvat)}

final.cvat$FDR <- p.adjust(final.cvat$p, method='BH') # estimate FDR

columns(org.Dm.eg.db)

# add KEGG pathway names:
final.cvat$KEGG <- mapIds(org.Dm.eg.db, 
                          keys=row.names(final.cvat), 
                          column="PATH", 
                          keytype="PATH",
                          multiVals="first")

tmp <- paste0('path:map', rownames(final.cvat))
final.cvat$pathway <- keggList("pathway")[tmp]

write.table(final.cvat, file=output.path) # save output to working directory

# clean up
rm(round2.cvat)
rm(round2.Sets)
rm(round3.cvat)
rm(round3.Sets)
rm(round1.cvat)
rm(round4.Sets)
rm(round4.cvat)
rm(round1.Sets)
gc()

############################################################



############################################################
# GO terms
############################################################
# ensure that working directory has a folder named 'qgg'
# define genome feature set
load(file="rapa.goSets_1kb")
round1.Sets <- goSets_1kb
rm(goSets_1kb)
set.name <- 'GOlevel_1kb' #user defined output file name based on the nature of the genome feature set
############################################################

output.path <- paste0('./CVAToutput/', set.name, '.cvat')

# 4 rounds of permutation, successive deeper permutations are used to resolve the lowest P values
round1.cvat <- gsea(fit=fit, sets=round1.Sets, W=W, method = "cvat", nperm=1000, ncores = n.cores)
rownames(round1.cvat)[round1.cvat$p<0.01] -> pass_to_10K # pull sets below a P threshold 
round2.Sets <- round1.Sets[pass_to_10K]
print(paste('round 1 complete,', length(pass_to_10K), 'features sent to round 2'))

round2.cvat <- gsea(fit=fit, sets=round2.Sets, W=W, method = "cvat", nperm=10000, ncores = n.cores)
rownames(round2.cvat)[round2.cvat$p<0.001] -> pass_to_100K # pull sets below a P threshold
round3.Sets <- round1.Sets[pass_to_100K]
rm(round2.Sets)
print(paste('round 2 complete,', length(pass_to_100K), 'features sent to round 3'))

round3.cvat <- gsea(fit=fit, sets=round3.Sets, W=W, method = "cvat", nperm=100000, ncores = n.cores)
rownames(round3.cvat)[round3.cvat$p<0.0001] -> pass_to_1M # pull sets below a P threshold
round4.Sets <- round1.Sets[pass_to_1M]
rm(round3.Sets)
print(paste('round 3 complete,', length(pass_to_1M), 'features sent to round 4'))

# combine all rounds:
if(length(pass_to_1M) > 0) { round4.cvat <- gsea(fit=fit, sets=round4.Sets, W=W, method = "cvat", nperm=1000000, ncores = n.cores)
final.cvat <- rbind(round1.cvat[!rownames(round1.cvat) %in% rownames(round2.cvat), ], round2.cvat[!rownames(round2.cvat) %in% rownames(round3.cvat), ], round3.cvat[!rownames(round3.cvat) %in% rownames(round4.cvat), ], round4.cvat) } else {
  final.cvat <- rbind(round1.cvat[!rownames(round1.cvat) %in% rownames(round2.cvat), ], round2.cvat[!rownames(round2.cvat) %in% rownames(round3.cvat), ], round3.cvat)}

final.cvat$FDR <- p.adjust(final.cvat$p, method='BH') # estimate FDR

write.table(final.cvat, file=output.path) # save output to working directory

head(final.cvat)

columns(org.Dm.eg.db)

# add gene symbol:
final.cvat$GO <- mapIds(org.Dm.eg.db, 
                          keys=row.names(final.cvat), 
                          column="GO", 
                          keytype="GO",
                          multiVals="first")


library(GO.db)
goterms <- Term(GOTERM)

final.cvat$GO.name <- goterms[rownames(final.cvat)]

write.table(final.cvat, file=output.path) # save output to working directory

# clean up
rm(goterms)
rm(round2.cvat)
rm(round2.Sets)
rm(round3.cvat)
rm(round3.Sets)
rm(round1.cvat)
rm(round4.Sets)
rm(round4.cvat)
rm(round1.Sets)
gc()


############################################################
# make sortable html tables of the results:
############################################################
files <- as.character(dir('CVAToutput/'))
files <- files[files != 'Icon\r']

for(i in 1:length(files)) {
  f <- paste0('CVAToutput/', files[i])
  d <- read.table(f)
  d <- d[order(d$p), ]
  d[ ,c(1, 3,4)] <- apply(d[ ,c(1, 3,4)], 2, function(x) round(x, 5)) # round for the table
  d2 <- datatable(d[1:200,])
  htmlwidgets::saveWidget(d2, paste0(files[i], '.html')) }

############################################################

############################################################
# reduce redundancy among GO terms:  
############################################################
library(rrvgo)

# 1. load P values
# 2. reduce terms
dir('CVAToutput')
go0 <- read.table('CVAToutput/GOlevel_0bp.cvat')
head(go0)

bpSimMat <- calculateSimMatrix(go0$GO, orgdb="org.Dm.eg.db", ont="BP", method="Rel") # subset of GO terms re: Biological Process
mfSimMat <- calculateSimMatrix(go0$GO, orgdb="org.Dm.eg.db", ont="MF", method="Rel") # subset of GO terms re: Molecular Function
ccSimMat <- calculateSimMatrix(go0$GO, orgdb="org.Dm.eg.db", ont="CC", method="Rel") # subset of GO terms re: Cell Compartment

dim(bpSimMat)
dim(mfSimMat)
dim(ccSimMat)
# [NOTE]: there is no intersection among GO terms, GO terms are unique to each category

# [NOTE]: all GO terms in all categories are among the GO terms tested by CVAT
table(go0$GO %in% colnames(bpSimMat)) # the Fs may be those not in the BP (biological process category)
table(go0$GO %in% colnames(mfSimMat)) # the Fs may be those not in the BP (biological process category)
table(go0$GO %in% colnames(ccSimMat)) # the Fs may be those not in the BP (biological process category)

scores <- setNames(-log10(go0$FDR), go0$GO)
reducedTerms <- reduceSimMatrix(bpSimMat, scores, threshold=0.7, orgdb="org.Dm.eg.db")

head(reducedTerms)
reducedTerms <- aggregate(.~parent, reducedTerms, max) # this will reduce the go list to the highest scoring (lowest P) term within each parent term (this includes the parent terms as well)
reducedTerms$score <- as.numeric(reducedTerms$score)
reducedTerms$FDR <- 10^-reducedTerms$score # get back FDR
head(reducedTerms)
hits <- reducedTerms[reducedTerms$FDR <=0.1, ] # subset by threshold for significance
head(hits) #'hits' is a table summarizing the signifcant Go Terms and their less redundant representative 'parentTerm'


