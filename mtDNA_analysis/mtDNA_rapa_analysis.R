library(plyr)
library(corrplot)

setwd( "/Users/ben/Google_Drive/Documents/rapa")

read.table("pheno", header=T) -> pheno
head(pheno)
read.table('/Users/ben/Google_Drive/Applications/plink_mac/DGRP/rapa/covars', header = T) -> covars
head(covars)

merge(pheno, covars, by='FID') -> dat
head(dat)
dat[ ,c(1, 3, 5:ncol(dat))] -> dat
names(dat)[2] <- 'pheno'

head(dat)


# several sources of mtDNA genotype data:

# Richardson et al., 2012 https://doi.org/10.1371/journal.pgen.1003129
read.table('~/Google_Drive/Documents/generic.GWAS.analysis/mt.variants_Richardson_2012', header=T, stringsAsFactors = T) -> Richvars
Richvars[1:4, 1:4]
str(Richvars)


# mtDNA variants in 169 DGRP lines [Bevers, Litovchenko, Kapopoulou. et al. Mitochondrial haplotypes affect metabolic phenotypes in the Drosophila Genetic Reference Panel. Nat Metab 1, 1226–1242 (2019). https://doi.org/10.1038/s42255-019-0147-3]

# Haplotype Defining Variants (HDV, from Bevers et al., 2019; supplemental Table 12): these are 9 markers that are in LD with mt haplotypes 
HDV <- read.table('~/Google_Drive/Documents/generic.GWAS.analysis/mt.haplotype.markers_Bevers.2019.txt', header=T, stringsAsFactors = F)
head(HDV)
plyr::ldply(HDV, function(x) sum(x==1, na.rm=T)) # summary of haplotype frequencies

# mtDNA variants
Bevars <- read.table('~/Google_Drive/Documents/generic.GWAS.analysis/mt.variants_Bevers_2019', header=T, stringsAsFactors = T)
Bevars[1:4, 1:4]
# Genotype codes: 0 = reference, 1 = homozygote alternate, 2 = heterozygote alternate, 3 = homozygote for the second alternate 4 = heterozygote for the second alternate. 

# snp-info
snp.info <- read.table('~/Google_Drive/Documents/generic.GWAS.analysis/mt.variant_info_Bevers_2019', header=T, stringsAsFactors = F, sep= '\t') 
head(snp.info)

table(as.matrix(Bevars)) # Genotype codes: 0 = reference, 1 = homozygote alternate, 2 = heterozygote alternate, 3 = homozygote for the second alternate 4 = heterozygote for the second alternate. 
# how to handle the various genotypes?  a few ideas:
# 1. just to 'go for it', using the genotypes as they are.
# 2. convert the heterozygous (heteroplasmic) genotypes (codes= 2 and 4) to NA.  This seems a better option becuase we do not know the state of these nucletides in our DGRP lines when they were phenotyped.

# regardless, make sure that the allele codes are converted into a format that is compatible with the tool you are using.  coding them as factors is probably the most appropriate.  AVOID keeping them as numeric or integers, otherwise R will think that they are numerically-meaningful i.e. 1 is lower than 2, which is lower than 3, but not as low as 1 (a ghaslty scenario).

# remember, the goal is to have fun!

# option 2. convert 'heterozygous' genotype codes to NA:
as.matrix(Bevars) -> Bevars
table(Bevars)
Bevars[Bevars==2] = NA
Bevars[Bevars==4] = NA
table(Bevars) #[NOTE]: NAs don't show up in table()

head(dat)
table(dat$FID %in% rownames(Bevars)) # 130 lines in the rapa study have mtDNA data in Bevers data
table(dat$FID %in% rownames(Richvars)) # 123 lines in the rapa study have mtDNA data in Richardson data
table(rownames(Richvars) %in% rownames(Bevars)) # 134 Richards lines in Bevers
table(rownames(Bevars) %in% rownames(Richvars)) # 134 Bevers lines in Richards

nrow(Bevars)
nrow(Richvars) # more lines in Richards, but there are more 'rapa' lines in Bevars


HDV[1:4, 1:4]
tmp <- as.data.frame(HDV)
x <- cbind(dat, tmp[dat$FID, ])
x[1:4,1:18]
str(x)

plyr::ldply(x, function(c) sum(c==1, na.rm=T))# no invariant mtDNA haplotypes in these data

x[c(3:8,14:ncol(x))] <- lapply(x[c(3:8,14:ncol(x))], factor) 
str(x)
summary(lm(pheno ~ . -FID, x)) 
plot(pheno ~ In.2R.NS, x)

o <- summary(lm(pheno ~ chrM.791 + chrM.1716 + chrM.2071 + chrM.3892 + chrM.6872 + chrM.7424 + chrM.11128 + chrM.11416 + chrM.14666, x)) # simple model, without inversions, wolb and geno PCs

mod <- lm(pheno ~ . -FID, dat)
residual.pheno <- resid(mod)
names(residual.pheno) <- dat$FID
residual.pheno

x <- cbind(residual.pheno, tmp[names(residual.pheno), ])
x[1:4, ]
str(x)

x[2:ncol(x)] <- lapply(x[2:ncol(x)], factor) 
str(x)
summary(lm(residual.pheno ~ ., x))

par(mfrow=c(1,1))
plot(residual.pheno ~ chrM.6872, x)
plot(x$residual.pheno ~ x$chrM.6872, main='mtDNA')
stripchart(x$residual.pheno ~ x$chrM.6872, pch=16, add=T, vertical=T, method='jitter')
# not much to write home about.

snp.info[snp.info$ID == 'chrM:6872', ]
head(snp.info)

Bevars[1:4,1:4]
table(is.na(Bevars))
table(colSums(is.na(Bevars)) == 0)
Bevars[ ,colSums(is.na(Bevars)) == 0] -> B

p <- prcomp(B) 
plot(p$x, las=1, main='mtDNA PCA') # looks like three 'haplotypes' (mitotypes?)
dim(p$x)
pairs(p$x[ ,1:5])

mitotypes <- p$x[ ,1]+p$x[ ,2] # this vector could define three mitotypes
hist(mitotypes)
table(is.na(mitotypes))
bins <- c(-4, -2, 1, 4)
mtypes <- cut(mitotypes, bins)
table(mtypes)
table(is.na(mtypes))

names(mtypes) <- names(mitotypes)
head(dat)
x <- cbind(dat, mtypes[dat$FID])
length(mtypes)
nrow(dat)
table(is.na(x$`mtypes[dat$FID]`))

## no asssociation between rapa phenotype and PC1+2 'mitotype:
x <- x[!is.na(x$`mtypes[dat$FID]`), ]
summary(lm(pheno ~ as.factor(`mtypes[dat$FID]`), x))
plot(pheno ~ as.factor(`mtypes[dat$FID]`), x, ylab='rapa phenotype', xlab="PC1+2 'mitotype'", outline=F, ylim=c(-3.3, 3), las=1)
stripchart(pheno ~ as.factor(`mtypes[dat$FID]`), x, pch=16, add=T, method='jitter', vertical=T)


# correlation among mtDNA variants
cor(B) -> c # correlation matrix of non-missing data
c[1:10,1:10]

corrplot::corrplot(c)
heatmap(c)

hist(c['chrM.6872', ], 30)

