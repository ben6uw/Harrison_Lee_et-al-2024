################################################################################
# the phenotype values for GWAS/genetic analysis:
################################################################################

library(dplyr)
library(bestNormalize)

getwd()
setwd("/Volumes/GoogleDrive/My Drive/Documents/rapa/paper/code")
dir()

#read raw developmental time data, use all lines, even those later shown to be incorrect genotypes, to perform normalization and batch correction - leaving all lines in for this will allow more accurate corrections (i.e. these data will more accurately reflect variation by environment b/w batches)
read.csv("fDEVrapa_contains impostors.csv", header = TRUE, stringsAsFactors = FALSE) -> rw.data

# load genotype-verified line IDs:
genotyped <- read.csv("fDEVrapa_COMPILED_GENOTYPED.csv", header = TRUE, stringsAsFactors = FALSE) 
head(genotyped)
validated_lines <- unique(genotyped$Genotype) 
validated_lines <-validated_lines[!validated_lines %in% c("Canton-S", "W1118")]# use these later to trim GWAS down to only lines of verified identity

#convert some variables character
as.character(rw.data$Genotype) -> rw.data$Genotype
as.character(rw.data$Replicate) -> rw.data$Replicate
as.character(rw.data$Day) -> rw.data$Day
as.character(rw.data$Time) -> rw.data$Time

head(rw.data)

#aggregate replicates by Experiment, Treatment, Day, Time 
aggregate(Pupae~Experiment+Genotype+Treatment+Day+Time, data = rw.data, FUN = sum) -> rw.agg 

#compute Devtime
as.numeric(rw.agg$Day) -> rw.agg$Day
as.numeric(rw.agg$Time) -> rw.agg$Time
rw.agg$Devtime <- ((rw.agg$Time / 2400) + rw.agg$Day)

head(rw.agg)

#expand data frame by pupae
rw.final <- rw.agg[rep(seq_len(nrow(rw.agg)), times = rw.agg$Pupae), c(1:3,7)]

head(rw.final)

table(rw.final$Genotype, rw.final$Treatment) # sample sizes (n pupa) across the screen
 
#step 5 compute DEVOctrl and DEVOrapa
rw.sum <- rw.final %>% group_by(Experiment, Genotype, Treatment) %>% summarise_at(vars(Devtime), list(mean = mean, n = length, sd = sd))

#step 6 taking difference between DEVOrapa and DEVOctrl 
rw.sum$ID <- paste(rw.sum$Experiment, rw.sum$Genotype, sep = "-")
rw.sum[rw.sum$Treatment == "EtOH",] -> ctrl
rw.sum[rw.sum$Treatment == "Rapa",] -> rapa
merge(ctrl, rapa, by = "ID") -> tmp
head(tmp)
tmp2 <- tmp[, c("ID", "Experiment.x", "Genotype.x", "mean.x", "n.x", "sd.x", "mean.y", "n.y", "sd.y")]
colnames(tmp2) <-c("ID", "Experiment", "Genotype",  "mean.ctrl", "n.ctrl", "sd.ctrl", "mean.rapa", "n.rapa","sd.rapa")
tmp2$diff <- tmp2$mean.rapa - tmp2$mean.ctrl

head(tmp2)

tmp2$diff.pooledSD <- sqrt((tmp2$sd.ctrl^2 * (tmp2$n.ctrl-1) + tmp2$sd.rapa^2 * (tmp2$n.rapa-1))/(tmp2$n.ctrl+tmp2$n.rapa-2)) # pooled SD
tmp2$diff.pooledSE <- sqrt(tmp2$diff.pooledSD * sqrt((1/tmp2$n.ctrl) + (1/tmp2$n.rapa))) # pooled SE

diff.sum <- tmp2[, c("ID", "Experiment", "Genotype", "diff", "diff.pooledSD", 'diff.pooledSE')]
head(diff.sum)


# boxcox transformation
min(diff.sum$diff)
bestNormalize::boxcox(diff.sum$diff + 1) -> x
x$x.t -> diff.sum$bc.diff
hist(diff.sum$diff)
hist(diff.sum$bc.diff)

head(diff.sum)

# calculate mean development time differences of the 4 repeated 'batch control lines' within each batch
key <- c("W1118", "Canton-S", "45", "321") # the batch control lines
x <- diff.sum[diff.sum$Genotype %in% key, ]

diff.sum <- diff.sum %>% group_by(Experiment) %>% filter(Genotype %in% key) %>% summarise_at(vars(bc.diff), list(batch.ctrl.mean = mean)) %>% inner_join(diff.sum)

# and the standard deviation of all lines in each batch
DEVO <- diff.sum %>% group_by(Experiment) %>% summarise_at(vars(bc.diff), list(batch.sd = sd)) %>% inner_join(diff.sum)
head(diff.sum)

head(DEVO)
# calculate centered/scaled/center&scaled phenotype
DEVO$centered <- DEVO$bc.diff - DEVO$batch.ctrl.mean
DEVO$cs <- DEVO$centered / DEVO$batch.sd

# aggregate by genotype
GWAS <- as.data.frame(DEVO %>% group_by(Genotype) %>% summarise_at(vars(cs), list(GWAS = mean)))
head(GWAS)
GWAS <- GWAS[GWAS$Genotype %in% validated_lines, ] # trim GWAS data down to those of validated lines
GWAS$FID <- paste0('line_', GWAS$Genotype)
GWAS$IID <- GWAS$FID
GWAS <- GWAS[ ,c(3:4,2)]
names(GWAS)[3] <- 'GWAS_pheno'
write.table(GWAS, "GWASinput.txt", row.names = F, quote=F)

################################################################################
rm(list=ls())
################################################################################
