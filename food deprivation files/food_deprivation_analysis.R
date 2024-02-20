########################################################################
# does food deprevation modify the rapamycin-sensitivity phenotype?
########################################################################
library(dplyr)
library(plyr)
library(tidyr)
library(splitstackshape)
library(lme4)
library(lmerTest)
library(ggplot2)

rm(list=ls())

# pupa counts for vials started with eggs on 3/21/23 (pg 92-93, Ben's notebook IV)
setwd('/Users/ben/Library/CloudStorage/GoogleDrive-ben6@uw.edu/My Drive/Documents/rapa/starvation_rapa_devtime')
dir()

key <- read.csv('vial_key.csv', head=T) 
head(key)

dat <- read.csv('dat.csv', head=T)
head(dat)

dat$vial <- rownames(dat)
key$vial <- rownames(key)
dat[is.na(dat)] = 0

censuses <- colnames(dat)[-c(1:2)]
censuses
day <- c(7:17, 19, 21, 23, 27)
day

# could do hours, centered on noon (i.e. 11am would be -1 h from day x 24h)
time_of_day <- print(sapply(strsplit(censuses, split='_'), '[[', 3))
time_of_day <- c(10, 9.5, 10, 12, 8, 11, 10, 14, 12, 14, 13.5, 10, 10, 14, 10) # manual times
time_of_day <- time_of_day-12 
time_of_day

hour <- (day*24) + time_of_day
day <- hour/24
day
names(day) <- censuses


dat <- merge(key, dat)
dat[1:4,1:12]

long <- pivot_longer(dat, col=9:ncol(dat), values_to='pupae', names_to='census')
head(long)
long <- long[ ,-c(2:4)]
head(long)
long$day <- day[long$census]

str(long)

long$genotype <- as.factor(long$genotype)
long$rapa <- as.factor(long$rapa)
long$food <- as.factor(long$food)

long$food <- factor(long$food, levels=c(100, 75, 50))
long$food

head(long)

long <- expandRows(long, 'pupae')
head(long) # each row is a pupae

save(long, dat, file='processed.data')
############################################################################################


############################################################################################
rm(list=ls())
load('processed.data')

table(long$genotype, long$food, long$rapa)

apply(table(long$genotype, long$rapa), 1, min) # limiting sample size across rapa condition, eliminate those with <10 observations
lines <- rownames(table(long$genotype, long$rapa))
lines.to.keep <- lines[apply(table(long$genotype, long$rapa), 1, min) >=10]

long <- long[long$genotype %in% lines.to.keep, ]
long$genotype <- droplevels(long$genotype)

table(long$genotype, long$food, long$rapa)
summary(lm(day ~ rapa * food * genotype, long))

sensitivity <- read.table('~/Library/CloudStorage/GoogleDrive-ben6@uw.edu/My Drive/Documents/rapa/lines.txt', head=T)
head(sensitivity)

long$genotype <- paste0('line_', long$genotype)

table(long$genotype %in% sensitivity$line) # all lines accounted for

s <- merge(long, sensitivity, by.x='genotype', by.y='line') 
head(s)
table(s$phenotype)

s$phenotype <- as.factor(s$phenotype)

summary(lmer(day ~ rapa * food * phenotype + (1|genotype), s)) # check out the fixed effects!  food dilution slows them down, but doesn't modify the rapa phenotype!


Sp <- ggplot(subset(s, phenotype=='S'), aes(y=day, x=food, fill=rapa))+
  geom_violin()+
  facet_grid(~food)+
  facet_wrap(~genotype) +
  theme_classic()


Rp <- ggplot(subset(s, phenotype=='R'), aes(y=day, x=food, fill=rapa))+
  geom_violin()+
  facet_grid(~food)+
  facet_wrap(~genotype) +
  theme_classic()

ggpubr::ggarrange(Sp, Rp, labels = c("sensitive", "resistant"), ncol = 2)



## for plotting with a cutoff, may just want to manually chop down the data 
table(long$genotype, long$food, long$rapa)
tmp <- s[!s$genotype %in% c('line_318', 'line_439', 'line_776', 'line_287'), ] # remove lines with few observations
tmp$genotype <- as.factor(tmp$genotype)
table(tmp$genotype)
tmp$genotype <- droplevels(tmp$genotype)

summary(lmer(day ~ rapa * food * phenotype + (1|genotype), tmp))

table(tmp$phenotype, tmp$genotype)
tmp$rapa
levels(tmp$rapa) <- c("control", "rapamycin")
summary(mod <- lmer(day ~ rapa * food * phenotype + (1|genotype), tmp))

library(broom)
library(jtools)
summ(mod, digits = 3) 

# plot 
save(tmp, file='dataForFigure_food_dilution')
rm(list=ls())


load('/Users/ben/Library/CloudStorage/GoogleDrive-ben6@uw.edu/My Drive/Documents/rapa/starvation_rapa_devtime/dataForFigure_food_dilution')
head(tmp)
means <- aggregate(day ~ food * rapa * genotype, tmp, mean)
head(means)
P <- tmp$phenotype
names(P) <- tmp$genotype
P

means$phenotype <- P[as.character(means$genotype)]
table(means$phenotype, means$genotype)

means$treatment <- means$rapa
means$food <- plyr::revalue(means$food, c("100"="100% food", '75'='75% food', '50'='50% food'))


ggplot(means, aes(y=day, x=treatment, color=treatment, group=genotype))+
  geom_line(color='grey')+
  geom_point(size=2)+
  labs(fill='treatment') +
  labs(x='treatment', y='development time (day)')+
  facet_wrap(~phenotype ~food)+
  theme_bw()+
  theme(text = element_text(size = 14), axis.text.x=element_blank())

#############################################################################################


