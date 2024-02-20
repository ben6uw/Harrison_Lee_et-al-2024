####################################################
# does the effect of rapamycin on larvae resemble that of starvation? 
####################################################

library(ComplexHeatmap)
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

rm(list=ls())
gc()
dev.off()

setwd("~/Google Drive/My Drive/Documents/rapa/metabolomics/starving_larva_Jouandin_Science_2022")
dir()

# metabolome data
dat <- read.table('raw.data.csv', sep=',', head=T) # data from Jouandin et al
head(dat)

# sample information
# from file: 'Sample collection Sheet.xlsx', which contained the note: Larvae were raised on 20g/L  yeast 8% cornmeal 5%sugar and processed near 96h AEL						
meta <- read.table('metadata.csv', sep=',', head=T, stringsAsFactors = T)
head(meta)

# try to understand experiment: 2 variables confounded with starvation time were noticed:
table(meta$Stage, meta$Starvation.time..PBS..kimwipe.) # [NOTE]: in this study, time under starvation is confounded with time since egg laying (AEL).

table(meta$Mass.Spec.ID, meta$Starvation.time..PBS..kimwipe.) # sample IDs, presumably for the LCMS.  [NOTE]: if these ran in order, then starvation time (and development time) are confounded with run order.

dat$Mass.Spec.ID <- as.numeric(sapply(strsplit(dat$Sample.Name, split='J'), '[[', 2))
table(dat$Mass.Spec.ID %in% meta$Mass.Spec.ID)

dat <- merge(meta, dat)
dat[1:6,1:19]

# try to understand experiment:
table(dat$Stage, dat$Starvation.time..PBS..kimwipe.)
table(dat$Starvation.time..PBS..kimwipe.)
table(dat$Time)
table(dat$Time, dat$Starvation.time..PBS..kimwipe.)

colnames(dat)[1:14]
par(mfrow=c(2,2))
plot(Larval.weight ~ n.larvae, dat, pch=19)
plot(Larval.weight ~ Sample.weight, dat, pch=19)
pairs(dat[ ,c(7:9, 13)], pch=19, las=1)


colnames(dat)[1:14]
mat <- as.matrix(dat[ ,c(14:ncol(dat))])
mat[mat == 'N/A'] <- NA
mat <- apply(mat, 2, as.numeric)
boxplot(mat)
boxplot(log(mat))
mat <- log(mat)

table(colSums(is.na(mat))) #248 complete metabolites, 10 missign in one sample

plot(colSums(is.na(mat)), colMeans(mat, na.rm=T), las=1, ylab='mean log(metabolite) level', xlab='n missing data', main='missingness at limit of detection')
table(colSums(is.na(mat)) <2)

mat <- mat[, colSums(is.na(mat)) <2] # remove mzs missing in more than one sample
imp <- impute::impute.knn(mat)$data # of the 10 imputed mzs, 4 of them: betaine, choline, creatine, thiamine, are among the 84 metabolties that will be matched with the rapa study (see below)
table(is.na(imp))
mat <- imp
rm(imp)

rownames(mat) <- rownames(dat)

boxplot(mat)
boxplot(t(mat))
boxplot(scale(t(mat)), las=1) # centered and scaled by sample

scat <- t(scale(t(mat)))

dat[1:6,1:15]
dat <- cbind(dat[ ,1:13], scat) # down to 251 complete metabolites

save(dat, mat, scat, meta, file='processed_dat')
################################################################################
rm(list=ls())


################################################################################
load('processed_dat')
# look for effects of covariates:
dat[1:6,1:15]

mzs <- colnames(mat)
covarP <- matrix(nr=length(mzs), nc=3)

for(i in 1:length(mzs)) {
s <- summary(lm(dat[ ,mzs[i]] ~ Sample.weight + n.larvae + Larval.weight, dat))
covarP[i,] <- s$coefficients[2:4, 4] }

rownames(covarP) <- mzs
colnames(covarP) <- c('Sample.weight', 'n.larvae', 'massPerLarva')
covarFDR <- apply(covarP, 2, function(x) p.adjust(x, 'BH'))

colSums(covarFDR < 0.05) # many metabolites affected by n.larvae as well as mass per Larva

par(mfrow=c(2,3))
moi <- 'uracil'
plot(lm(dat[ ,moi] ~ Sample.weight + n.larvae, dat)$residuals ~ dat$Larval.weight, pch=19, ylab=moi)
plot(lm(dat[ ,moi] ~ Sample.weight + Larval.weight, dat)$residuals ~ dat$n.larvae, pch=19, ylab=moi)
plot(lm(dat[ ,moi] ~ Larval.weight + n.larvae, dat)$residuals ~ dat$Sample.weight, pch=19, ylab=moi)


scat[1:4,1:4]
pca <- prcomp(scat, scale=T)$x
plot(pca[ ,1:2], col=as.factor(dat$Time), pch=19)


dat[1:5,1:19]
pca <- cbind(dat[ ,1:13], pca)

head(pca)
pairs(pca[ ,c(4, 8:10, 13:19)], las=1, pch=16)


# about to load mzs from the rapa study, rename the mzs for the current, starvation study first, so that 'mzs' isn't written over
starv.mzs <- mzs

load('~/Google Drive/My Drive/Documents/rapa/metabolomics/larva_metabolomics/data/processed.data')

rapa.dat <- dat
rapa.mzs <- mzs

rm(dat)
rm(mz.info)

rapa.mzs
starv.mzs
table(rapa.mzs %in% starv.mzs) # only 8 names initially intersect
rapa.mzs <- toupper(rapa.mzs)
starv.mzs <- toupper(starv.mzs)
table(rapa.mzs %in% starv.mzs) # 48 after uppering the case


# some starvation mzs are complementary ions (posi and nega). Find them, scale them and take their mean
tmp <- starv.mzs

table(table(tmp))

tmp

table(grepl('.POSI', tmp))
table(grepl('.NEGA', tmp))
table(grepl('_POSI', tmp))
table(grepl('_NEGA', tmp))

tmp <- gsub('.NEGA', '', tmp)
tmp <- gsub('.POSI', '', tmp)
tmp <- gsub('_NEGA', '', tmp)
tmp <- gsub('_POSI', '', tmp)

table(table(tmp)) # found 9 mzs with both ions

# scale the mzs, take the mean of the duplicate ions:
x <- scat
colnames(x) <- tmp
boxplot(x)
boxplot(scale(x))
x <- scale(x)
cn = colnames(x)
scat2 <- sapply(unique(cn), function(g) rowMeans(x[ ,cn==g, drop=F])) # this gives the mean of each of the duplicate ions
rm(scat)
rm(cn)
rm(x)

starv.mzs <- colnames(scat2)

table(rapa.mzs %in% starv.mzs) # 50 now
starv.mzs <- gsub('\\.', '-', starv.mzs)

x <- rapa.mzs[!rapa.mzs %in% starv.mzs]
y <- starv.mzs[!starv.mzs %in% rapa.mzs]

x[order(x)]
y[order(y)]

# the following code was determine, by hand, to convert metabolite names so that the same metabolite matches across studies 
starv.mzs <- gsub('-DL', '', starv.mzs)
starv.mzs <- gsub('-NEGA', '', starv.mzs)
starv.mzs <- gsub('D-', '', starv.mzs)
starv.mzs <- gsub('DL-', '', starv.mzs)
starv.mzs <- gsub('-_NEGA', '', starv.mzs)
starv.mzs <- gsub('-DL', '', starv.mzs)
starv.mzs <- gsub('-_NEGA', '', starv.mzs)
starv.mzs <- gsub('_NEGA', '', starv.mzs)
starv.mzs <- gsub('-_POSI', '', starv.mzs)
starv.mzs <- gsub('-POSI', '', starv.mzs)
starv.mzs <- gsub('_POSI', '', starv.mzs)
starv.mzs <- gsub('-ACID', ' ACID', starv.mzs)
starv.mzs <- gsub('-ACETYL-', '-AC-', starv.mzs)
starv.mzs <- gsub('KYNURENINE', 'L-KYNURENINE', starv.mzs)
starv.mzs <- gsub('A-KETOGLUTARATE', 'ALPHA-KETOGLUTARIC ACID', starv.mzs)
starv.mzs <- gsub('X1-METHYL-HISTIDINE', '1/3-METHYLHISTIDINE', starv.mzs)
starv.mzs <- gsub('X1-METHYLADENOSINE', '1-METHYLADENOSINE', starv.mzs)
starv.mzs <- gsub('LEUCINE-ISOLEUCINE', 'LEUCINE /D-NORLEUCINE', starv.mzs)
starv.mzs <- gsub('GLUCOSE-1-PHOSPHATE', 'G1P/F1P/F6P', starv.mzs)
starv.mzs <- gsub('GLUCOSE-6-PHOSPHATE', 'G6P', starv.mzs)
starv.mzs <- gsub('THYMINE', 'THIAMINE', starv.mzs)
starv.mzs <- gsub('GLYCERALDEHDYE-3-PHOSPHATE', 'DHAP/D-GA3P', starv.mzs)
starv.mzs <- gsub('CHOLESTERYL-SULFATE', 'CHOLESTERYL SULFATE', starv.mzs)
starv.mzs <- gsub('GLUTAMATE', 'GLUTAMIC ACID', starv.mzs)
starv.mzs <- gsub('SEDOHEPTULOSE-1-7-PHOSPHATE', 'SEDOHEPTULOSE 7-PHOSPHATE', starv.mzs)
starv.mzs <- gsub('METHIONINE-SULFOXIDE', 'METHIONINE SULFOXIDE', starv.mzs)
starv.mzs <- gsub('X2-HYDROXYGLUTERATE', '2-HYDROXYGLUTARATE', starv.mzs)
starv.mzs <- gsub('SN-GLYCEROL-3-PHOSPHATE', 'GLYCERO3-P', starv.mzs)
starv.mzs <- gsub('X3-PHOSPHOGLYCERATE', '2/3-PHOSPHOGLYCERIC ACID', starv.mzs)
starv.mzs <- gsub('METHIONINE-SULFOXIDE', 'METHIONINE SULFOXIDE', starv.mzs)
starv.mzs <- gsub('PIPECOLIC ACID', 'PIPECOLATE', starv.mzs)
starv.mzs <- gsub('URIC ACID', 'URATE', starv.mzs)
starv.mzs <- gsub('N-AC-GLUTAMIC ACID', 'N-AC-GLUTAMATE', starv.mzs)
starv.mzs <- gsub('ASPARTATE', 'ASPARTIC ACID', starv.mzs)
starv.mzs <- gsub('OXALOACETATE', 'OXALACETATE', starv.mzs)
starv.mzs <- gsub("RIBOSE 5-PHOSPHATE", 'RIBOSE-5-P', starv.mzs)
starv.mzs <- gsub('X3-HYDROXYBUTERATE', '3HBA', starv.mzs)
starv.mzs <- gsub('S-ADENOSYL-L-METHIONINE', 'S-ADENOSYLMETHIONINE (SAM)', starv.mzs)
starv.mzs <- gsub('ACETYLLYSINE', 'N6-ACETYL-LYSINE', starv.mzs)
starv.mzs <- gsub('N-AC-GLUTAMINE', 'N-AC-L-GLUTAMINE', starv.mzs)
starv.mzs <- gsub('GLYCERO3-P', 'GLYCEROL-3-P', starv.mzs)
starv.mzs <- gsub('S-ADENOSYL-L-HOMOCYSTEINE', 'SAH', starv.mzs)

table(rapa.mzs %in% starv.mzs) # 84 mathcing metabolites now 

x <- rapa.mzs[!rapa.mzs %in% starv.mzs]
y <- starv.mzs[!starv.mzs %in% rapa.mzs]
x[order(x)]
y[order(y)]


# reanalyze the starvation data using only the mzs that intersect with the rapa mzs:
load('processed_dat')
mzs <- colnames(scat2)
colnames(scat2) <- starv.mzs # add the new names
scat2 <- scat2[ ,starv.mzs %in% rapa.mzs] # reduce scat to only those mzs that are also measured in the rapa study

scat2[1:4,1:6]
dat[1:4,1:15]

dat <- cbind(dat[ ,1:13], scat2)
dat[1:6,1:15]

scat2[1:4,1:4]
pca <- prcomp(scat2)$x

p <- prcomp(scat2)
p$sdev/sum(p$sdev) # proportionof variation by PC
plot(p)

dat[1:5,1:19]
pca <- cbind(dat[ ,1:13], pca)

head(pca)
pairs(pca[ ,c(4, 8:10, 13:19)], las=1, pch=16, col=as.factor(dat$Time))

pairs(pca[ ,c(13:17)], las=1, pch=16, col=as.factor(dat$Time))

z <- prcomp(scat2) # get the loadings from PC1 to use on the rapa data
starv.loadings <- z$rotation

table(toupper(colnames(rapa.dat)) %in% starv.mzs)
rapamzdat <- rapa.dat[ ,toupper(colnames(rapa.dat)) %in% starv.mzs]

colnames(rapamzdat) <- toupper(colnames(rapamzdat))
table(colnames(rapamzdat) %in% rownames(starv.loadings))

starv.loadings[ ,1]

rapamzdat <- rapamzdat[ ,rownames(starv.loadings)] # align the mz data 

boxplot(rapamzdat)
boxplot(scale(rapamzdat))

rapa.dat$rapa_over_starvPC1 <- apply(scale(rapamzdat), 1, function(x) sum(x * starv.loadings[ ,'PC1']))
rapa.dat$rapa_over_starvPC2 <- apply(scale(rapamzdat), 1, function(x) sum(x * starv.loadings[ ,'PC2']))

rapa.dat$rapa_over_starvPC1
rapa.dat$rapa_over_starvPC2

treat <- ifelse(rapa.dat$treatment == 'R', 'rapa', '(-)')
rapa.dat$phenotype.treat <- as.factor(paste(rapa.dat$phenotype, treat))
levels(rapa.dat$phenotype.treat)


save(file='data processed for cross-study analysis', pca, rapa.dat)
rm(list=ls())

################################################################
load('data processed for cross-study analysis')

################################################################
# non-linear model of PC1 ~ Time

lmod <- lm(PC1 ~ Time, pca)
summary(lm(PC1 ~ Time, pca))
abline(lm(PC1 ~ Time, pca), col=2)

y <- pca$PC1
x <- pca$Time
lmod <- nls(y~a*x-c, start=list(a=1, c=-10))
nlmod <- nls(y~a*x/(b+x)-c, start=list(a=1, b=1, c=-10))

cor(pca$PC1, predict(lmod))^2
cor(pca$PC1, predict(nlmod))^2

cor.test(pca$PC1, predict(nlmod))

summary(nlmod)
AIC(lmod, nlmod)
anova(nlmod, lmod) # non-linear model is significantly better than linear model.

pars <- nlmod$m$getPars()

############################################
# plotting:
############################################

summary(lmer(rapa_over_starvPC1 ~ treatment * phenotype + (1|line), rapa.dat))


Fig4A <- ggplot(pca, aes(x=Time, y=PC1))+
  geom_point()+
  theme_bw() +
  xlab('Time (h)') +
  theme(text = element_text(size = 14))+
stat_function(fun=function(x) pars['a'] * x / (pars['b'] + x) - pars['c'], col=2)


means <- aggregate(rapa_over_starvPC1 ~ line + treatment + phenotype, rapa.dat, mean)

levels(means$phenotype) <- c('resistant', 'sensitive')
levels(means$treatment) <- c('control', 'rapamycin')

Fig4B <- ggplot(means, aes(x=treatment, y=rapa_over_starvPC1, color=treatment, group=line))+
  geom_line(color='grey')+
  geom_point(size=2)+
  labs(fill='treatment') +
  ylab(bquote(PC[starvation]))+
  labs(x='treatment')+
  facet_wrap(~phenotype)+
  theme_bw() +
  theme(text = element_text(size = 14), axis.text.x=element_blank(), panel.spacing = unit(0.2, "lines"), legend.position="none", strip.background = element_blank())    

ggpubr::ggarrange(Fig4A, Fig4B, nrow=1)


#
# add analysis from food dilution experiment:

load('/Users/ben/Library/CloudStorage/GoogleDrive-ben6@uw.edu/My Drive/Documents/rapa/starvation_rapa_devtime/dataForFigure_food_dilution')
head(tmp)
foodDilution<- aggregate(day ~ food * rapa * genotype, tmp, mean)
head(foodDilution)
P <- tmp$phenotype
names(P) <- tmp$genotype

foodDilution$phenotype <- P[as.character(foodDilution$genotype)]
table(foodDilution$phenotype, foodDilution$genotype)
foodDilution$treatment <- foodDilution$rapa
foodDilution$food <- plyr::revalue(foodDilution$food, c("100"="100% food", '75'='75% food', '50'='50% food'))

foodDilution$food_genotype <- paste(foodDilution$food, foodDilution$genotype)       
foodDilution$food_treatment <- paste(foodDilution$food, foodDilution$treatment)   
levels(foodDilution$phenotype) <- c('resistant', 'sensitive')


Fig4A <- ggplot(pca, aes(x=Time, y=PC1))+
  geom_point(color='grey30')+
  theme_bw() +
  xlab('Time (h)') +
  theme(text = element_text(size = 12))+
  stat_function(fun=function(x) pars['a'] * x / (pars['b'] + x) - pars['c'], col=2)


Fig4B <- ggplot(means, aes(x=treatment, y=rapa_over_starvPC1, color=treatment, group=line))+
  geom_line(color='grey')+
  geom_point(size=1.5)+
  labs(fill='treatment') +
  ylab(bquote(PC[starvation]))+
  labs(x='treatment')+
  facet_wrap(~phenotype)+
  theme_bw() +
  theme(text = element_text(size = 12), axis.text.x=element_blank(), panel.spacing = unit(0.2, "lines"), legend.position="none", strip.background = element_blank()) 

Fig4C <- ggplot(foodDilution, aes(y=day, x=treatment, color=treatment, group=food_genotype))+
  geom_line(color='grey')+
  geom_point(size=1.5)+
  labs(fill='food_treatment') +
  labs(x='treatment', y='development time (day)')+
  facet_grid(phenotype~food)+
  coord_cartesian(ylim = c(7, 20))+
  theme_bw() +
  theme(text = element_text(size = 12), axis.text.x=element_blank(), panel.spacing = unit(0.2, "lines"), strip.background = element_blank())  

library(ggpubr)
ggarrange(Fig4A, Fig4B, Fig4C)

ggarrange(ggarrange(Fig4A, Fig4B, ncol = 2, labels = c("A", "B")), Fig4C, nrow=2, labels = c("","","C"))





