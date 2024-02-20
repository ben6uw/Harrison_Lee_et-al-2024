#####################################################################################
#####################################################################################
###    fDEV RAPA combined pupation analysis 
###    MBL analysis 10222021 - using genotyped final dataset
#######################################################################################
######################################################################################

###load packages
library(outliers)
library(Hmisc)
library(ggplot2)
library(smatr)
library(dplyr)

rm(list=ls())
##############################################################################################################
getwd()
###  set working directory
setwd("/Users/ben/Library/CloudStorage/GoogleDrive-ben6@uw.edu/My Drive/Documents/rapa/paper/code")
dir()
### read the raw data, all lines validated by PCR genotyping:
raw <- read.csv("fDEVrapa_COMPILED_GENOTYPED.csv") 

data <- raw[,c("Exp_Genotype", "Genotype", "Time", "Day", "Replicate","Treatment","Pupae")]

###################################################
### Summary table by genotype&treatment
###################################################
#combine vials 
aggregate(Pupae~Treatment+Time+Day+Genotype, data = data, FUN = sum) -> data.combined

###alternative way to expand rows by number of pupae per time point
devtime <- data.combined[rep(seq_len(nrow(data.combined)), times = data.combined$Pupae), 1:4]
devtime$Realtime = devtime$Time/2400 +devtime$Day

#extract information
se <- function(x) sd(x)/sqrt(length(x))
head(devtime)
stats <- devtime %>% group_by(Genotype, Treatment) %>% summarise_at(vars(Realtime), list(mean = mean, median = median, var=var, n = length, sd = sd, se = se))

head(stats)
range(stats$n)
hist(stats$n [!stats$Genotype %in% c('45', '321', 'W1118', 'Canton-S')]) # sample size distribution for all gentypes (except for the 4 genotypes run in each batch)
mean(stats$n [!stats$Genotype %in% c('45', '321', 'W1118', 'Canton-S')]) # sample size distribution for all genotypes (except for the 4 genotypes run in each batch)
range(stats$n [!stats$Genotype %in% c('45', '321', 'W1118', 'Canton-S')]) # sample size distribution for all genotypes (except for the 4 genotypes run in each batch)

# look at relationship between n (number of pupae counted) and pupation time.  Were the most-delayed lines less likely to form pupae?
lattice::xyplot(mean ~ n | Treatment, stats[!stats$Genotype %in% c('45', '321', 'W1118', 'Canton-S'), ], pch=19, xlab='pupae counted per genotype', ylab = 'mean development time (days)') # seems not, the low-n lines are not strongly biased toward the most delayed

#make a data frame out of stats_vial
as.data.frame(stats) -> fDEVgenotype
#make a csv file
write.table(fDEVgenotype,file="fDEV DGRP RAPA genotype summary.csv",sep=",",row.names=F)

####################################################################################################
######################################################################################################
##### difference in mean pupation time
##################################
stats[stats$Treatment == "EtOH", ] -> ctrl
stats[stats$Treatment == "Rapa", ] -> trtm
diff <- merge(ctrl, trtm, by="Genotype")

diff$diff <- diff$mean.y - diff$mean.x
hist(diff$diff)

head(diff)

######################################################################################################
# average pupation time difference for lines used in the larva size experiment (Bill Young)
mean(diff$diff[diff$Genotype %in% c('181', '348', '517', '57')]) # sensitive lines
sd(diff$diff[diff$Genotype %in% c('181', '348', '517', '57')]) # sensitive lines

mean(diff$diff[diff$Genotype %in% c('176' ,'535')]) # resistant lines
sd(diff$diff[diff$Genotype %in% c('176' ,'535')]) # resistant lines

# average pupation time difference for lines used in the larva metabolome sampling (Mitchell Lee)
R <- c('142', '153', '176', '324', '383', '441', '535') # resistant lines 
S <- c('287', '318', '348', '517', '57', '796', '822') # sensitive lines

diff$diff[diff$Genotype %in% S]
mean(diff$diff[diff$Genotype %in% S]) # sensitive lines
sd(diff$diff[diff$Genotype %in% S]) # sensitive lines
range(diff$diff[diff$Genotype %in% S])

mean(diff$diff[diff$Genotype %in% R]) # resistant lines
sd(diff$diff[diff$Genotype %in% R]) # resistant lines
range(diff$diff[diff$Genotype %in% R])

# for Fig1 stats
diff[diff$Genotype=='181', ]
diff[diff$Genotype=='176', ]

######################################################################################################

# pooled variance estimates
#pooled standard deviation
pooledSD <- sqrt((diff$sd.x^2 * (diff$n.x-1) + diff$sd.y^2 * (diff$n.y-1))/(diff$n.x + diff$n.y-2)) # pooled SD
diff$pooledSD <- pooledSD

#pooled standard error
pooledSE <- sqrt(diff$pooledSD * sqrt((1/diff$n.x) + (1/diff$n.y)))
diff$pooledSE <- pooledSE

#ranked line plot 
tmp <- diff[order(diff$diff), ]

X= which(tmp$Genotype %in% c('Canton-S', 'W1118'))
Y= tmp$diff[tmp$Genotype %in% c('Canton-S', 'W1118')]

plot(tmp$diff, xaxt="n", xlab = "", ylim=c(-1.5,7.6), cex =0.8, ylab = "", las = 1, pch=20, col=1, cex.axis=1.5)
mtext('ranked lines (n=142)', side=1, line=1, at=(142/2), cex=1.4)
mtext("mean developmental delay (days +/- pooled SD)", side=2, line=2.5, at=3.5, cex=1.2)
arrows(1:142, tmp$diff -tmp$pooledSD, 1:142, tmp$diff + tmp$pooledSD, length=0, angle=0, code=3, col='grey30')
points(tmp$diff, cex =0.5, pch=20)
arrows(x=X, y=Y-1.1, x1=X, y1=Y-0.5, length=0.1)
text(x=X, y=Y-1.3, labels=c('Canton-S', 'W1118'), cex=1)

head(tmp)

tmp$Genotype <- factor(tmp$Genotype, levels = tmp$Genotype)

ggplot(tmp, aes(y=diff, x=Genotype)) +
  geom_point(size=0.5)+
  theme_classic(base_size = 14)+
  xlab('genotype')+
  ylab('mean developmental delay\n(days +/-se)')+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y = element_text(size = 14))+
  geom_errorbar(aes(ymin=diff-pooledSE, ymax=diff+pooledSE), width=0, size=0.2)+
  geom_segment(aes(x = X[1], y = Y[1]-1.5, xend = X[1], yend = Y[1]-0.5), lwd=0.1, arrow = arrow(length = unit(0.2, "cm")))+
  annotate(geom="text", x = X[1], y = Y[1]-2, label="Canton-S", color="grey20")+
  geom_segment(aes(x = X[2], y = Y[2]-1.5, xend = X[2], yend = Y[2]-0.5), lwd=0.1, arrow = arrow(length = unit(0.2, "cm")))+
  annotate(geom="text", x = X[2], y = Y[2]-2, label="W1118", color="grey20")

head(tmp)
mean(tmp$diff)
sd(tmp$diff)

range(tmp$diff)/mean(tmp$mean.x) # range of developmental delay as a proportion of the control developmental time. 




#############################################################################################################
###   y phenotype x control development time:
############################################################################################################

test <- psych::corr.test(diff$mean.x, diff$diff, method='spearman')
test$r # spearman rho (if method='spearman' used above)
test$p # P

head(tmp)

ggplot(tmp, aes(y=diff, x=mean.x)) +
  geom_point()+
  theme_classic(base_size = 14)+
  xlab('development time - control (days)')+
  ylab('mean developmental delay (days)')


#############################################################################################################
###   rapamycin dose-response
############################################################################################################

library(dplyr)
library(tidyr)
library(tidyverse)
library(lubridate)
library(splitstackshape)
library(lme4)
library(lmerTest)
library(ggplot2)


rm(list=ls())

# pupa counts for vials started with eggs on 4/27/23 places 12:30pm-3pm (pg 106, Ben's notebook IV)
setwd('/Users/ben/Library/CloudStorage/GoogleDrive-ben6@uw.edu/My Drive/Documents/rapa/dose_response')
dir()

key <- read.csv('key.csv', head=T) 
head(key)
colnames(key)[4] <- 'rapamycin_nmol' # raelized that the amount of rapa was in nmol, not nmol; replace labels here

table(key$stripe, key$rapamycin_nmol) # vials had stripes on the side indicating the rapa treatment level
table(key$genotype, key$rapamycin_nmol) # two vials per genotype per rapa condition

dat <- readxl::read_excel('data_May_2023.xlsx')
head(dat)
table(dat[ ,1]==key[ ,1])

dat$vial <- rownames(dat)
key$vial <- rownames(key)

dat[is.na(dat)] = 0

# convert census times to elapsed time in days
censuses <- colnames(dat)[-c(1)]

day.o.the.month <- sapply(strsplit(censuses, split='_'), '[[', 2)

# eggs deposited on 4/27 from 12:30pm-3pm, so day 5/4 (noon) was 'day 7':
day <- 3 + as.numeric(day.o.the.month)
# time centered on noon (i.e. 11am would be -1 h from day x 24h)
time_of_day <- as.data.frame(print(sapply(strsplit(censuses, split='_'), '[[', 3)))
names(time_of_day) <- 'time'
time_of_day <- time_of_day %>% separate(time, into = c('clocktime', 'ampm'), sep = -2, convert = TRUE)
head(time_of_day)

time_of_day$clocktime[!grepl(':', time_of_day$clocktime)] <- paste0(time_of_day$clocktime[!grepl(':', time_of_day$clocktime)], ':00')
time_of_day$clocktime    
time_of_day$clocktime <- hm(time_of_day$clocktime)      
time_of_day$clocktime
hour(time_of_day$clocktime) <- hour(time_of_day$clocktime) + ifelse(time_of_day$ampm == 'pm', 12, 0)
hour <- hour(time_of_day$clocktime) + minute(time_of_day$clocktime)/60
hour <- hour-12 # this makes the noon deposition of eggs time=0
day <- day + hour/24 # calculate numeric time in days 
names(day) <- censuses

dat <- merge(key, dat)
head(dat)

long <- pivot_longer(dat, col=5:ncol(dat), values_to='pupae', names_to='census')
long$day <- day[long$census]
head(long)

long$genotype <- as.factor(long$genotype)
long$rapamycin_nmol <- as.factor(long$rapamycin_nmol)
long <- expandRows(long, 'pupae')
head(long) # each row is a pupae
long$nmol <- as.numeric(as.character(long$rapamycin_nmol))

summary(lm(day ~ rapamycin_nmol * genotype, long))

sensitivity <- read.table('~/Library/CloudStorage/GoogleDrive-ben6@uw.edu/My Drive/Documents/rapa/lines.txt', head=T)
sensitivity$phenotype <- as.factor(sensitivity$phenotype)
levels(sensitivity$phenotype) <- c('resistant', 'sensitive')
long$genotype <- gsub('Ral', 'line', long$genotype)

table(long$genotype %in% sensitivity$line) # all lines accounted for
head(long)

s <- merge(long, sensitivity, by.x='genotype', by.y='line') 
s$phenotype <- as.factor(s$phenotype)
head(s)

summary(lmer(day ~ rapamycin_nmol * phenotype + (1|genotype), s))
summary(lmer(day ~ nmol * phenotype + (1|genotype), s)) # rapa dose as a numeric
summary(lm(day ~ nmol * genotype, s)) # rapa dose as a numeric
summary(lm(day ~ rapamycin_nmol, s[s$genotype == 'line_383', ]))

plot(day ~ nmol, s[s$genotype == 'line_383', ], pch=16)
abline(lm(day ~ nmol, s[s$genotype == 'line_383', ]))
summary(lm(day ~ nmol, s[s$genotype == 'line_383', ]))
hist(residuals(lm(day ~ nmol, s[s$genotype == 'line_383', ])))

n <- print(s %>% group_by(genotype, rapamycin_nmol) %>%tally())
means <- as.data.frame(s %>% group_by(genotype, rapamycin_nmol)  %>% summarise_at(vars("day"), funs(mean = mean, sd = sd)))
head(means)
means <- merge(means, n)
means$se <- means$sd/sqrt(means$n)
means <- merge(means, sensitivity, by.x='genotype', by.y='line') 

head(means)
str(means)

head(s)
str(s)

#######################################
# wilcox comparison to control:

treats <- levels(s$rapamycin_nmol)[-c(1)]
lines <- unique(s$genotype)

pMat <- matrix(nr=length(lines), nc=length(treats))

for(g in 1:length(lines)){
  tmp <- s[s$genotype==lines[g], ]
  for(i in 1:length(treats)) {
    tryCatch({
      pMat[g, i] <- wilcox.test(tmp$day[tmp$rapamycin_nmol=='0'], tmp$day[tmp$rapamycin_nmol==treats[i]], alternative = c("less"))$p.value}, error=function(e){}) 
  } }

colnames(pMat) <- treats
rownames(pMat) <- lines
pMat

bonferroniMat <- t(apply(pMat, 1, function(x) p.adjust(x, 'bonferroni')))
bonferroniMat <=0.05 # use these data to mark the plot

means <- as.data.frame(s %>% group_by(genotype, rapamycin_nmol)  %>% summarise_at(vars("day"), funs(mean = mean, sd = sd)))
head(means)
means <- merge(means, n)
means$se <- means$sd/sqrt(means$n)
means <- merge(means, sensitivity, by.x='genotype', by.y='line') 
head(means)

bon <- as.data.frame(bonferroniMat)
bon$genotype <- rownames(bon)
head(bon)
bon <- bon %>% pivot_longer(!genotype, names_to = 'rapamycin_nmol', values_to = 'bonferroni')
head(bon)
means <- merge(bon, means, all.y=T)
head(means)
means$bonferroni <=0.05

dose.response.plot.S1 <- ggplot(means, aes(y=mean, x=rapamycin_nmol, group=genotype, color=phenotype))+
  geom_point()+
  geom_line() +
  geom_errorbar(aes(ymin= mean - se, ymax= mean + se), width=.2, position=position_dodge(0.05)) +
  facet_wrap(~phenotype) +
  ylab('mean development time (days +/-se)')+
  xlab('rapamycin (nmol)') +
  theme_classic() +
  theme(text = element_text(size = 18))

dose.response.plot.S1 


# divide each genotype by its mean development time in 0 rapa condition:
control.mean <- means$mean[means$rapamycin_nmol=='0']
names(control.mean) <- means$genotype[means$rapamycin_nmol=='0']
control.mean
means$control.mean <- control.mean[means$genotype]
means$relative.development <- means$mean - means$control.mean

means$rapamycin_nmol <- factor(means$rapamycin_nmol, levels=c(0, 2, 4, 8, 16))

ggplot(means, aes(y=relative.development, x=rapamycin_nmol, group=genotype, color=phenotype))+
  geom_point()+
  geom_line() +
  geom_errorbar(aes(ymin= relative.development - se, ymax= relative.development + se), width=.2, position=position_dodge(0.05)) +
  facet_wrap(~phenotype) +
  ylab('relative developmental delay (days +/-se)')+
  xlab('rapamycin (nmol)') +
  theme_classic() +
  theme(text = element_text(size = 16))

means$sigshape <- ifelse(means$bonferroni <=0.05, 19, 16)
means$sigfill <- ifelse(means$bonferroni <=0.05, 'grey30', 'white')
means$sigshape[is.na(means$sigshape)] <- 16 # for the points at rapa=0
means$sigfill[is.na(means$sigfill)] <- 'white' # for the points at rapa=0


ggplot(means, aes(y=relative.development, x=rapamycin_nmol, group=genotype, color=phenotype))+
  geom_line() +
  geom_errorbar(aes(ymin= relative.development - se, ymax= relative.development + se), width=0.1) +
  geom_point(shape=21, size=2, fill=means$sigfill)+
  ylab('relative developmental delay (days +/-se)')+
  xlab('rapamycin (nmol)') +
  theme_classic() +
  theme(text = element_text(size = 16))


dose.response.1.plot <- ggplot(means, aes(y=relative.development, x=rapamycin_nmol, group=genotype, color=phenotype))+
  geom_line() +
  geom_errorbar(aes(ymin= relative.development - se, ymax= relative.development + se), width=0.1) +
  geom_point(shape=21, size=2, fill=means$sigfill)+
  ylab('relative developmental delay (days +/-se)')+
  xlab('rapamycin (nmol)') +
  theme_classic() +
  theme(text = element_text(size = 16))

save(dose.response.1.plot, dose.response.plot.S1, file='dose.response.1.plots')


rm(list=ls())

#############################################################
# 2nd experiment, testing 0, 2, 8, 16, 128 and 1024 nmol Rapa (Ben's Notebook IV, pg 137):
#############################################################
dir('dose_response_2')
key <- read.csv('dose_response_2/key.csv')
cens <- read.csv("dose_response_2/pupa_counts_dose_response_2_late pupae removed.csv", nrows=5) # To avoid counting the 2nd generation pupae, I removed any pupae that were counted 10 days after the first pupae in a vial were recorded, this eliminated ~10 pupae (the un-deleted data are in: "dose_response_2/pupa_counts_dose_response_2.csv")
n <- cens[ ,1]
cens <- as.data.frame(t(cens[ ,-c(1)]))
colnames(cens) <- n
rm(n)


cens$elapsed_day <- as.numeric(cens$elapsed_day)

dat <- read.csv("dose_response_2/pupa_counts_dose_response_2_late pupae removed.csv", skip=5) # To avoid counting the 2nd generation pupae, I removed any pupae that were counted 10 days after the first pupae in a vial were recorded, this eliminated ~10 pupae (the un-deleted data are in: "dose_response_2/pupa_counts_dose_response_2.csv")
head(dat)
dat <- (dat[-c(1), ])
dat[is.na(dat)]=0
colnames(dat)[1] <- 'vial'

dat <- as.data.frame(pivot_longer(dat, cols=c(2:ncol(dat)), names_to='day', values_to='pupae'))
dat$day <- as.numeric(gsub('X', '', dat$day))

dat <- splitstackshape::expandRows(dat, count='pupae')
dat <- merge(key, dat)
dat$rapamycin_nmol <- as.factor(dat$rapamycin_nmol)

sensitivity <- read.table('~/Library/CloudStorage/GoogleDrive-ben6@uw.edu/My Drive/Documents/rapa/lines.txt', head=T)
head(sensitivity)
sensitivity$phenotype <- as.factor(sensitivity$phenotype)
levels(sensitivity$phenotype) <- c('resistant', 'sensitive')
sensitivity$line <- gsub('line', 'Ral', sensitivity$line)

dat <- merge(dat, sensitivity, by.x='genotype', by.y='line') 
head(dat)


ggplot(dat, aes(y=day, x=rapamycin_nmol, group=rapamycin_nmol, fill=phenotype))+
  geom_violin()+
  theme_classic() +
  facet_wrap(~genotype)+
  theme(text = element_text(size = 18))


table(dat$genotype, dat$rapamycin_nmol) # summary of sample sizes, pupae per condition
n <- dat %>% group_by(genotype, rapamycin_nmol) %>%tally()

head(dat)

#######################################
# wilcox comparison to control:
treats <- levels(dat$rapamycin_nmol)[-c(1)]
lines <- unique(dat$genotype)

pMat <- matrix(nr=length(lines), nc=length(treats))

for(g in 1:length(lines)){
  tmp <- dat[dat$genotype==lines[g], ]
  for(i in 1:length(treats)) {
    tryCatch({
      pMat[g, i] <- wilcox.test(tmp$day[tmp$rapamycin_nmol=='0'], tmp$day[tmp$rapamycin_nmol==treats[i]], alternative = c("less"))$p.value}, error=function(e){}) 
  } }

colnames(pMat) <- treats
rownames(pMat) <- lines
pMat

bonferroniMat <- t(apply(pMat, 1, function(x) p.adjust(x, 'bonferroni')))
bonferroniMat <=0.05 # use these data to mark the plot

means <- as.data.frame(dat %>% group_by(genotype, rapamycin_nmol)  %>% summarise_at(vars("day"), funs(mean = mean, sd = sd)))
head(means)
means <- merge(means, n)
means$se <- means$sd/sqrt(means$n)
means <- merge(means, sensitivity, by.x='genotype', by.y='line') 
head(means)

bon <- as.data.frame(bonferroniMat)
bon$genotype <- rownames(bon)
head(bon)
bon <- bon %>% pivot_longer(!genotype, names_to = 'rapamycin_nmol', values_to = 'bonferroni')
head(bon)
means <- merge(bon, means, all.y=T)
head(means)
means$bonferroni <=0.05

means$rapamycin_nmol <- factor(means$rapamycin_nmol, levels=c('0', '2', '16', '128' ,'1024'))

means$sigshape <- ifelse(means$bonferroni <=0.05, 19, 16)
means$sigfill <- ifelse(means$bonferroni <=0.05, 'grey30', 'white')
means$sigshape[is.na(means$sigshape)] <- 16 # for the points at rapa=0
means$sigfill[is.na(means$sigfill)] <- 'white' # for the points at rapa=0

ggplot(means, aes(y=mean, x=rapamycin_nmol, group=genotype, color=phenotype))+
  geom_line()+
  geom_errorbar(aes(ymin= mean - se, ymax= mean + se), width=.2, position=position_dodge(0.05)) +
  geom_point(shape=21, size=3, fill=means$sigfill)+
  facet_wrap(~phenotype) +
  ylab('mean development time (days +/-se)')+
  xlab('rapamycin (nmol)') +
  theme_bw() +
  theme(text = element_text(size = 18), legend.position="none")


load('dose.response.1.plots')

dose.response.plot.S2 <- ggplot(means, aes(y=mean, x=rapamycin_nmol, group=genotype, color=phenotype))+
  geom_point(shape=19)+
  geom_line() +
  geom_errorbar(aes(ymin= mean - se, ymax= mean + se), width=.2, position=position_dodge(0.05)) +
  facet_wrap(~phenotype) +
  ylab('mean development time (days +/-se)')+
  xlab('rapamycin (nmol)') +
  theme_classic() +
  theme(text = element_text(size = 16))


ggpubr::ggarrange(dose.response.plot.S1+ theme(text = element_text(size = 14))
                  , dose.response.plot.S2+ theme(text = element_text(size = 14)), 
                  common.legend = T, legend='right')


# divide each genotype by its mean development time in 0 rapa condition:
control.mean <- means$mean[means$rapamycin_nmol=='0']
names(control.mean) <- means$genotype[means$rapamycin_nmol=='0']
control.mean
means$control.mean <- control.mean[means$genotype]
means$relative.development <- means$mean - means$control.mean

ggplot(means, aes(y=relative.development, x=rapamycin_nmol, group=genotype, color=phenotype))+
  geom_point()+
  geom_line() +
  geom_errorbar(aes(ymin= relative.development - se, ymax= relative.development + se), width=.2, position=position_dodge(0.05)) +
  facet_wrap(~phenotype) +
  ylab('relative developmental delay (days +/-se)')+
  xlab('rapamycin (nmol)') +
  theme_classic() +
  theme(text = element_text(size = 16))


dose.response.2.plot <- ggplot(means, aes(y=relative.development, x=rapamycin_nmol, group=genotype, color=phenotype))+
  geom_line() +
  geom_errorbar(aes(ymin= relative.development - se, ymax= relative.development + se), width=0.1) +
  geom_point(shape=21, size=2, fill=means$sigfill)+
  ylab('relative developmental delay (days +/-se)')+
  xlab('rapamycin (nmol)') +
  theme_classic() +
  theme(text = element_text(size = 16))


ggpubr::ggarrange(dose.response.1.plot, dose.response.2.plot, common.legend = T, legend='right')


ggplot(means, aes(y=relative.development, x=rapamycin_nmol, group=genotype, color=phenotype))+
  geom_line() +
  geom_errorbar(aes(ymin= relative.development - se, ymax= relative.development + se), width=0.1) +
  geom_point(size=2)+
  ylab('relative developmental delay (days +/-se)')+
  xlab('rapamycin (nmol)') +
  theme_classic() +
  theme(text = element_text(size = 16))


rm(list=ls())


################################################################################
# larval size and growth rate analysis
################################################################################

# analyzing Bill Young's larva size data
# Bill had taken images of larvae that Alia had used in one of her development timing/staging experiments in 2021
# Bill then used ImageJ to measure the area of the larvae

setwd("/Users/ben/Library/CloudStorage/GoogleDrive-ben6@uw.edu/My Drive/Documents/rapa/larva size_Bill")
dat <- read.table('Larvae Size Data.txt', head=T, stringsAsFactors = T)

head(dat)
str(dat)
dat$Genotype <- as.factor(paste0('Ral_', dat$Genotype))
dat$Day <- as.factor(dat$Day)
dat$Area <- dat$area_sq_mm

table(dat$Sensitivity, dat$Genotype) # 2 resistant strains, 4 sensitive strains
table(dat$Genotype, dat$Treatment, dat$Day) # n per condition, day and genotype
range(table(dat$Genotype, dat$Treatment, dat$Day)) # range of n
mean(table(dat$Genotype, dat$Treatment, dat$Day)) # mean of n


levels(dat$Genotype)
dat$Genotype <- factor(dat$Genotype, levels = c("Ral_181", "Ral_348", "Ral_517", "Ral_57", "Ral_176", "Ral_535")) # arrange genotypes by resistance phenotype

library(ggplot2)

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


head(dat)

ggplot(dat, aes(x = Day, y = Area, group=Treatment, color = Treatment)) +
  geom_jitter(width = 0.1, size=0.5) +
  geom_smooth(method='lm', se=F)+
  theme_classic(base_size = 14) +
  labs(y = bquote('area'~(mm^2)), x = "time (day)") +
  scale_color_manual(values=c("grey60", cbbPalette[2])) +
  facet_grid(~ Genotype) 

ggplot(subset(dat,  Genotype %in% c('Ral_535', 'Ral_348')), aes(x = Day, y = Area, group=Treatment, color = Treatment)) +
  geom_jitter(width = 0.1, size=0.5) +
  geom_smooth(method='lm', se=F)+
  theme_classic(base_size = 14) +
  labs(y = bquote('area'~(mm^2)), x = "time (day)") +
  scale_color_manual(values=c("grey60", cbbPalette[2])) +
  facet_grid(~ Genotype) 


################################################################################
################################################################################
