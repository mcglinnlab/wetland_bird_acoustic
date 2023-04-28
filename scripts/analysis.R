library(ggplot2)
library(dplyr)
library(tidyr)
#library(vegan)
#library(readr)
#library(easyCODA)
#library(tidyverse)
#library(here)

#Read in data

sumcomm <- read.csv('./data/clean_data/sumcomm.csv') # comm matrix when all detections are totaled
comm <- read.csv('./data/clean_data/comm.csv')       # comm matrix when detections are converted into proportions
dat_sub <- read.csv('./data/clean_data/subset_data.csv') #the cleaned acoustic data with confidence levels and times

# look at species ranks based upon number of occurrences
ac_abundance <- sort(colSums(comm[ , -1] > 0), dec = T)
#write.csv(ac_abundance, file = "./data/clean_data/ac_abundance.csv", row.names = TRUE)

colnames(comm)
comm_wetland_id <- comm$wetland_id
#date_time <- sapply(strsplit(row.names(comm), '_'), function(x) x[2])
habitat <- ifelse(substring(comm_wetland_id, 1, 1) == 'U', 'upland', 'wetland')
property <- ifelse(substring(comm_wetland_id, 1, 1) == 'H', 'Halidon', 'Stono')

# need to drop upland sites

# read in point count data
dat_pc <- read.csv('https://raw.githubusercontent.com/mcglinnlab/wetland_birds/main/data/filtered_data/clean_bird_dat.csv')
#note: in dat_pc from Jackson that wetland_id is equivalent to what I've called
# site_id; therefore let's rename wetland_id to site_id for simplicity
#jackson's data is subset from 0y-25 m rather than previous total which was 50m.
head(dat_pc)
comm_pc <- dat_pc[ , c('wetland_id', names(dat_pc)[12:66])]
comm_pc[1:5, 1:5]
#Bachman's sparrow code is incorrect in ptct data. Switch BASP to BACS
colnames(comm_pc)[21] <- "BACS"
comm_pc = subset(comm_pc, select = -c(BGNN))
comm_pc$BASP
comm_pc$BCSP
comm_pc$BACS = comm_pc$BACS + comm_pc$BCSP
comm_pc$BHNU = 

  
#drop comm_pc$BCSP
comm_pc = subset(comm_pc, select = -c(BCSP))


#look at species ranks based upon number of occurrences
pc_abundance <- sort(colSums(comm_pc[,-1] > 0), dec = T)
#write.csv(pc_abundance, file = "./data/clean_data/pc_abundance.csv", row.names = TRUE)

#create sum ptct
sumcomm_ptct <- comm_pc %>% group_by(wetland_id) %>% summarise_each(funs(sum))

#transform into data frame
comm_ptct <- as.data.frame(comm_pc)

#compress point counts and acoustic into a single value per sample, tapply or aggregate
#merge & get 2 NOCA's- one for PC & one for acoustic
comm_ptct <- aggregate(comm_ptct, by = list(comm_ptct$wetland_id),
                       function(x) sum(x > 0) / length(x))

# drop wetland id and rename group 1 to wetland id
names(comm_ptct)
comm_ptct <- subset(comm_ptct, select = -wetland_id)
names(comm_ptct) <- c("wetland_id", names(comm_ptct) [-1])
names(comm_ptct)
names(sumcomm) <- c("wetland_id", names(sumcomm) [-1])

#pivot data frame to long format

comm_lg <- pivot_longer(comm, cols = !wetland_id)
comm_lg_ptct <- pivot_longer(comm_ptct, cols = !wetland_id)

#add column to dataframes indicating methods
comm_lg_ptct$method <- "ptct"
comm_lg$method <- "acoustic"
sumcomm$method <- "acoustic"
sumcomm_ptct$method <- "ptct"

#row append community matrices
comm_both <- rbind(comm_lg, comm_lg_ptct)

# drop sites not found in both sampling methods 
acoustic_sites <- unique(comm_both$wetland_id[comm_both$method == 'acoustic'])
sites_with_both <- acoustic_sites[acoustic_sites %in% unique(comm_both$wetland_id[comm_both$method == 'ptct'])]
sites_with_both

# drop upland sites
sites_with_both <- sites_with_both[!grepl('UP', sites_with_both)]

#create combined community matrix
comm_both <- subset(comm_both, wetland_id %in% sites_with_both)

# remove value of zero
comm_both <- subset(comm_both, value > 0)

#create occupancy matrix
sp_occ <- pivot_wider(comm_both, names_from = 'method', values_from = value)
sp_occ$acoustic <- ifelse(is.na(sp_occ$acoustic), 0, sp_occ$acoustic)
sp_occ$ptct <- ifelse(is.na(sp_occ$ptct), 0, sp_occ$ptct)

#plot
plot(acoustic ~ ptct, data = sp_occ, type = 'n')
with(sp_occ, text(ptct, acoustic, labels = name, cex = 0.5))
abline(a = 0 , b= 1)

#aggregate across methods
sp_occ_agg <- sp_occ %>% 
    group_by(name) %>% 
    summarize(acoustic = mean(acoustic), ptct = mean(ptct))

plot(acoustic ~ ptct, data = sp_occ_agg, type = 'n',
     xlab = 'Point Count Average Species Occupancy', 
     ylab = 'Passive Acoustic Average Species Occupancy')
with(sp_occ_agg, text(ptct, acoustic, labels = name, cex = 1))
abline(a = 0 , b= 1)

plot(I(sp_occ_agg$ptct + 0.01), I(sp_occ_agg$acoustic + 0.01), type = 'n', log = 'xy',
     xlab = 'Point Count Average Species Occupancy', 
     ylab = 'Passive Acoustic Average Species Occupancy')
abline(a = 0 , b= 1)
with(sp_occ_agg, text(I(ptct + 0.01), I(acoustic + 0.01), labels = name, cex = 0.75))

# mockup of barplot drafts
tmp <- pivot_longer(sp_occ_agg, cols = 2:3, names_to = 'method')
#stacked
barplot(value ~ method + name, data = tmp)
#side-by-side
barplot(value ~ method + name, data = tmp, beside = TRUE)
# subset to species with largest differences
occ_diff <- with(sp_occ_agg, abs(ptct - acoustic))
top_sp <- sp_occ_agg$name[order(occ_diff, decreasing = T)][1:10]
top_sp
barplot(value ~ method + name, data = tmp, subset = name %in% top_sp,
        beside = TRUE)

# species level occupancy analysis --------------------------------------

# need CI for occ estimates
boot_occ <- function(comm, grp, return_S = FALSE) {
  grp_lvs <- unique(grp)
  row_indices <- NULL
  for (i in seq_along(grp_lvs)) { 
    row_indices <- c(row_indices, 
                     sample(which(grp == grp_lvs[i]), replace = TRUE))
  }
  # rand set of comm
  occ <- aggregate(comm[row_indices, ], list(grp), function(x) sum(x > 0) / length(x))
  if (return_S)
    rowSums(occ[ , -1])
  else {
    rownames(occ) <- occ$Group.1
    occ <- occ[ , -1]
    as.matrix(occ)
  }
}

boot_raw_S <- function(comm, grp) {
  grp_lvs <- unique(grp)
  row_indices <- NULL
  for (i in seq_along(grp_lvs)) { 
    row_indices <- c(row_indices, 
                     sample(which(grp == grp_lvs[i]), replace = TRUE))
  }
  # rand set of comm
  S <- rowSums(comm[row_indices, ] > 0)
  Savg <- tapply(S, list(grp), mean)
  Savg
}


boot_sp_occ <- function(comm, grp, return_S = FALSE) {
  grp_lvs <- unique(grp)
  row_indices <- NULL
  for (i in seq_along(grp_lvs)) { 
    row_indices <- c(row_indices, 
                     sample(which(grp == grp_lvs[i]), replace = TRUE))
  }
  # rand set of comm
  occ <- tapply(comm[row_indices], list(grp), function(x) sum(x > 0) / length(x))
  if (return_S)
    sum(occ)
  else {
    t(occ)
  }
}


# demo with acoustic data
# you have a choice randomly select 5- min pt counts or randomly select sites
# simpler to bootstrap on just randomly drawn sites so let's do that first

comm_aco <- with(dat_sub, tapply(Confidence, INDEX = list(site_date, sp_code),
                                 function(x) (sum(x) > 0)*1))
comm_aco <- ifelse(is.na(comm_aco), 0, comm_aco)
comm_aco
comm_aco_sites <- substr(rownames(comm_aco), 1, 4)

# drop upland sites
comm_aco <- subset(comm_aco, !grepl('UP', comm_aco_sites))
comm_aco <- as.data.frame(comm_aco)
comm_aco
comm_aco_sites <- substr(rownames(comm_aco), 1, 4)
comm_aco_sites

# get ptct in correct format
comm_pc <- subset(comm_pc, wetland_id %in% comm_aco_sites)
comm_pc_sites <- comm_pc$wetland_id
comm_pc <- subset(comm_pc, select = -wetland_id)
comm_pc <- (comm_pc > 0) * 1

# 
#dat_pc <- pivot_longer(data.frame(wetland_id = comm_pc_sites, comm_pc),
#                       cols = -wetland_id)
#dat_aco <- pivot_longer(data.frame(wetland_id = comm_aco_sites, comm_aco), 
#                        cols = -wetland_id)
#dat_both <- rbind(data.frame(dat_pc, method = 'ptct'),
#                  data.frame(dat_aco, method = 'acoustic'))
#dim(dat_both)
#table(dat_both$value)
#head(dat_both)
#pivot_wider(dat_both, names_from = name, values_from = value)

boot_occ(comm_aco, grp = comm_aco_sites)
boot_occ(comm_pc, grp = comm_pc_sites)

aco_boots <- replicate(1e3, boot_occ(comm_aco, grp = comm_aco_sites), simplify = FALSE)
aco_sp_q <- apply(sapply(aco_boots, function(x) apply(x, 2, mean)), 1, quantile, c(0.025, 0.5, 0.975))
aco_sp_a <- apply(sapply(aco_boots, function(x) apply(x, 2, mean)), 1, mean)

pc_boots<- replicate(1e3, boot_occ(comm_pc, grp = comm_pc_sites), simplify = FALSE)
pc_sp_q <- apply(sapply(pc_boots, function(x) apply(x, 2, mean)), 1, quantile, c(0.025, 0.5, 0.975))
pc_sp_a <- apply(sapply(pc_boots, function(x) apply(x, 2, mean)), 1, mean)

aco_sp_q
pc_sp_q

aco_sp_a
pc_sp_a

all_sp_list <- sort(unique(c(names(aco_sp_a), names(pc_sp_a))))
all_sp_list

occ_sum <- data.frame()
for (i in seq_along(all_sp_list)) { 
    occ <- lo <- me <- hi <- c(0, 0)
    if (any(names(aco_sp_a) == all_sp_list[i])) {
        occ[1] <- aco_sp_a[names(aco_sp_a) == all_sp_list[i]]
        lo[1] <- aco_sp_q[1, colnames(aco_sp_q) == all_sp_list[i]]
        me[1] <- aco_sp_q[2, colnames(aco_sp_q) == all_sp_list[i]]
        hi[1] <- aco_sp_q[3, colnames(aco_sp_q) == all_sp_list[i]]
    }
    if (any(names(pc_sp_a) == all_sp_list[i])) {
        occ[2] <- pc_sp_a[names(pc_sp_a) == all_sp_list[i]]
        lo[2] <- pc_sp_q[1, colnames(pc_sp_q) == all_sp_list[i]]
        me[2] <- pc_sp_q[2, colnames(pc_sp_q) == all_sp_list[i]]
        hi[2] <- pc_sp_q[3, colnames(pc_sp_q) == all_sp_list[i]]
    }
    # examine if intervals do not overlap
    occ_diff <- abs(occ[1] - occ[2])
    sig_diff <- lo[1] > hi[2] | lo[2] > hi[1]
    
    occ_sum <- rbind(occ_sum, 
                     data.frame(name = all_sp_list[i],
                                occ, lo, me, hi, method = c('acoustic', 'ptct'),
                                occ_diff, sig_diff))
}

# which species have sig diffs?
tmp <- subset(occ_sum, sig_diff == TRUE)
tmp <- tmp[order(tmp$occ_diff, decreasing = TRUE), ]
tmp$name <- factor(tmp$name, levels = unique(tmp$name))

ggplot(tmp, aes(x = name, y = me)) +
  geom_point(aes(col = method)) + 
  geom_errorbar(aes(ymin = lo, ymax = hi, col = method), width = 0.2) + 
  xlab('Species Code') + 
  ylab('Species Site Occupancy') +
  scale_color_discrete(labels=c('Passive', 'Active')) +
  theme_minimal()
ggsave('./figs/species_occupancy_with_quantiles.pdf', width = 11.3, height = 6.15)

occ_sum_wide <- pivot_wider(occ_sum, names_from = method, values_from = c(occ, lo, me, hi))
names(occ_sum_wide)

pdf('./figs/species_occupacny_one_to_one.pdf')
plot(occ_sum_wide$`me_point count`, occ_sum_wide$me_passive, type='n',
     xlab = 'Active Count Median Species Occupancy',
     ylab = 'Passive Count Median Species Occupancy')
abline(a = 0, b=1)
with(occ_sum_wide, 
     text(`me_point count`, me_passive, labels = name, col = 'grey', cex = 0.75))
with(subset(occ_sum_wide, sig_diff == TRUE), 
     text(`me_point count`, me_passive, labels = name, cex = 0.75))
dev.off()




occ_sum <- data.frame()
for (i in seq_along(all_sp_list)) { 
  occ <- c(0, 0)
  if (any(names(aco_sp_a) == all_sp_list[i]))
    occ[1] <- aco_sp_a[names(aco_sp_a) == all_sp_list[i]]
  if (any(names(pc_sp_a) == all_sp_list[i]))
    occ[2] <- pc_sp_a[names(pc_sp_a) == all_sp_list[i]]
  occ_sum <- rbind(occ_sum, 
                   data.frame(name = all_sp_list[i],
                              acoustic = occ[1], ptct = occ[2]))
}

plot(occ_sum$ptct, occ_sum$acoustic)
abline(a = 0, b=1)
cor(occ_sum$ptct, occ_sum$acoustic)


ggplot(occ_sum, aes(x = name, y = me)) + 
  geom_bar(stat = 'identity', aes(fill = method), position=position_dodge())


# community level diversity analysis --------------------------------------

# first just use pres-absence and do bootstrapping

boot_raw_S(comm_aco, comm_aco_sites)


plot( boot_raw_S(comm_pc, comm_pc_sites), 
     boot_raw_S(comm_aco, comm_aco_sites), ylim = c(0, 5), xlim = c(0,5))
replicate(100, points(boot_raw_S(comm_pc, comm_pc_sites), 
     boot_raw_S(comm_aco, comm_aco_sites)))

aco_div_boots <- replicate(1e3, boot_raw_S(comm_aco, comm_aco_sites))
pc_div_boots <- replicate(1e3, boot_raw_S(comm_pc, comm_pc_sites))

# quantiles of site level richness values
aco_div_q <- apply(aco_div_boots, 1, quantile, c(0.025, 0.5, 0.975))
pc_div_q <- apply(pc_div_boots, 1, quantile, c(0.025, 0.5, 0.975))

# quantiles of average richness values
aco_div_q <- quantile(colMeans(aco_div_boots), c(0.025, 0.5, 0.975))
pc_div_q <- quantile(colMeans(pc_div_boots), c(0.025, 0.5, 0.975))

aco_div_q
pc_div_q
# not sure which of the two above is best

plot(pc_div_q[2, ], aco_div_q[2, ])

# ignore site level information
aco_div_boots <- replicate(1e3, boot_raw_S(comm_aco, rep(1, nrow(comm_aco))))
pc_div_boots <- replicate(1e3, boot_raw_S(comm_pc, rep(1, nrow(comm_pc))))

aco_div_q2 <- quantile(aco_div_boots, c(0.025, 0.5, 0.975))
pc_div_q2 <- quantile(pc_div_boots, c(0.025, 0.5, 0.975))

aco_div_q2
#2.5%      50%    97.5% 
#2.113260 2.240331 2.378522 
pc_div_q2
#2.5%      50%    97.5% 
#2.611111 2.930556 3.222222 

# both bootsrap procedures are giving same take home message
# the richness is lower in the acoustic dataset at a 5 min interva

# then repeat but using occupancy estimates

dim(aco_boots)

# site level
apply(sapply(aco_boots, rowSums), 1, quantile, 
      c(0.025, 0.5, 0.975))

apply(sapply(pc_boots, rowSums), 1, quantile, 
      c(0.025, 0.5, 0.975))

# study level but still at site grain
aco_div_occ <- quantile(colMeans(sapply(aco_boots, rowSums)),
                        c(0.025, 0.5, 0.975))
pc_div_occ <- quantile(colMeans(sapply(pc_boots, rowSums)),
                        c(0.025, 0.5, 0.975))

aco_div_occ
pc_div_occ
#2.5%      50%    97.5% 
#2.089535 2.205645 2.332023 

#2.5%      50%    97.5% 
#2.511111 2.766667 3.033611 

#same basic pattern as with raw richness

# make some graphs of richness

# from raw occurrences

pdf('./figs/species_richness_comparison_from_raw_occurrences.pdf')
bp <- barplot(c(aco_div_q[2], pc_div_q[2]), 
              names =  c('passive acoustic', 'point count'),
              ylab = 'Median Species Richness',
              ylim = c(0, 3.5),
              width = 0.35, xlim = c(0, 1))
arrows(bp[1], aco_div_q[1], bp[1], aco_div_q[3], 
       angle = 90, code =3, length = 0.2)
arrows(bp[2], pc_div_q[1], bp[2], pc_div_q[3], 
       angle = 90, code =3, length = 0.2)
dev.off()

# from occupancy sums

bp <- barplot(c(aco_div_occ[2], pc_div_occ[2]), 
        names =  c('passive acoustic', 'point count'),
        ylim = c(0, 3.5),
        width = 0.35, xlim = c(0, 1))
arrows(bp[1], aco_div_occ[1], bp[1], aco_div_occ[3], 
       angle = 90, code =3, length = 0.2)
arrows(bp[2], pc_div_occ[1], bp[2], pc_div_occ[3], 
       angle = 90, code =3, length = 0.2)


# the two above barplots are identical
# so it doesn't matter if you use raw occur and compute richness
# or you first convert to occupancy and then sum the occupancies

# now to examine beta diversity

aco_boots[[1]][1:5, 1:10]

calc_beta <- function(occ_comm) {
  mean_sp_occs <- colMeans(occ_comm)
  gamma <- sum(mean_sp_occs > 0)
  1 / mean(mean_sp_occs)
}
  

aco_beta_q <- quantile(sapply(aco_boots, calc_beta), c(0.025, 0.5, 0.975))
pc_beta_q <- quantile(sapply(pc_boots, calc_beta), c(0.025, 0.5, 0.975))

pdf('./figs/whittakers_beta_comparison.pdf')
bp <- barplot(c(aco_beta_q[2], pc_beta_q[2]), 
              names =  c('passive acoustic', 'point count'),
              ylab = expression('Species Turnover ('*beta[w]*')'),
              ylim = c(1, 30),
              width = 0.35, xlim = c(0, 1))
arrows(bp[1], aco_beta_q[1], bp[1], aco_beta_q[3], 
       angle = 90, code =3, length = 0.2)
arrows(bp[2], pc_beta_q[1], bp[2], pc_beta_q[3], 
       angle = 90, code =3, length = 0.2)
dev.off()


#color mimics separately from others on 1:1 plot. Maybe raptors too?
#NOMO, BHCO, BRTH

##TABLE: rank which species are the most different. sp_occ_agg with largest diff.
difference <- sp_occ_agg$acoustic - sp_occ_agg$ptct
difference_clean <- cbind(sp_occ_agg$name, difference)
colnames(difference_clean) <- c("name", "difference")
class(difference_clean)
write.csv(difference_clean, file = './data/Spec_Occ_Difference.csv', col.names = TRUE,
          row.names = FALSE)

#aggregate species occurrence 
S_acou <- with(sp_occ, tapply(acoustic, list(wetland_id), sum))
S_ptct <- with(sp_occ, tapply(ptct, list(wetland_id), sum))

#compute variance across methods (p*q, q= 1-p) ##- sig correlation between ptct
#sr and acoustic richness???
V_acou <- with(sp_occ, tapply(acoustic, list(wetland_id), 
                              function(x) sum(x * (1 - x))))
V_ptct <- with(sp_occ, tapply(ptct, list(wetland_id),
                              function(x) sum(x * (1 - x))))

#fix variance 

#compute correlation
correlation <- cor.test(S_ptct, S_acou, method=c("pearson"))

#-0.1 cor coefficient indicates no association between the variables
par(mfrow=c(1,2))
plot(S_ptct, S_acou)
abline(a = 0, b=1)
plot(V_ptct, V_acou)
abline(a = 0, b=1)

#Which species were only detected acoustic/ptct - which fit best? Filter out 
#non-songbird

library(readr)
tmp <- read_fwf('SpeciesList.txt', skip = 14,
                fwf_widths(c(7, 6, 51, 51, 51, 51, 51, 51, 50)))
hdr <- read_fwf('SpeciesList.txt', skip = 12, n_max = 1)

names(tmp) <- hdr[1, ]
head(tmp)
tmp$AOU <- as.integer(tmp$AOU)

write.csv(tmp, 'SpeciesList.csv', row.names = FALSE)

#pair AOU codes to species name codes
species_aou <- read.csv('./data/SpeciesList_w:AOU.csv')
sp_occ_both <- merge(sp_occ_agg, species_aou)

#filter to get diurnal land birds only

filter_diurnal <-
    sp_occ_both %>% filter(AOU > 2880) %>%
    
    filter(AOU < 3650 | AOU > 3810) %>%
    
    
    filter(AOU < 3900 | AOU > 3910) %>%
    
    filter(AOU < 4160 | AOU > 4210) %>%
    
    filter(AOU != 7010) %>% 
    print

diurnal_agg <- filter_diurnal %>% 
    group_by(name) %>% 
    summarize(acoustic = mean(acoustic), ptct = mean(ptct))

#plot 1:1 graph with only diurnal birds
plot(acoustic ~ ptct, data = filter_diurnal, type = 'n', main = 'Mean Species Occupancy By Method')
with(filter_diurnal, text(ptct, acoustic, labels = name, cex = 0.5, col = cols))
abline(a = 0 , b= 1)

#find unique species detected by one method & not the other
name_acou <- unique(comm$name)
length(name_acou)
name_ptct <- unique(comm_ptct$name)
length(name_ptct)

#create unique names for for loop
uniq_sp <- unique(comm_both$name)
length(uniq_sp)
uniq_method <- unique(comm_both$method)
uniq_site <- unique(comm_both$wetland_id)
output <- NULL

#for loop: if 'species name' exists in 1 method and not the other, print 'species name'.
for(i in seq_along(uniq_sp)){
    for(j in seq_along(uniq_site)){
            temp <- subset(comm_both, wetland_id == uniq_site[j] & 
                              name == uniq_sp[i])
            if(any(temp$value > 0)) {
                if(any(temp$value == 0)){
                 output <- rbind(output, subset(temp, value == 0)) 
                }
            }
    }
}

#remove seq along site loop & sum across sites to get species difference between methods
#disregarding site.

for(i in seq_along(uniq_sp)){
         hold <- subset(comm_both, name == uniq_sp[i])
            if(any(hold$value > 0)) {
                if(any(hold$value == 0)){
                sun <- rbind(hold, subset(hold, value == 0)) 
                }
         }
}

#fix comm_both, push to git, Dan will fix for loop & try to implement conf intervals.

#Q1: SR comparison: sum diversity anova

#sr_ptct <-comm_pc %>%
    #group_by(wetland_id) %>%
   # summarise(species.richness=n()) %>%
    #arrange(-species.richness)
#sr_ptct 

#sr_acou <- comm %>%
    #group_by(wetland_id) %>%
    #summarise(species.richness=n()) %>%
    #arrange(-species.richness)
#sr_acou

#fig 1 SR analysis: summed occupancy 
boxplot(S_acou, S_ptct)

#fig 2 species by species comparison

#fig 3 descriptive nature of why are these differing? could be table of species
#traits (large body, quiet singing). Chi sq tests? We have 2 numbers are they more \
#diff than expected?

boxplot(value ~ name, data = comm_both)
boxplot(value ~ method, data= comm_both)

#null: no SR difference between methods. determine if two categorical variables
#have a significant correlation between them.
library(dplyr)
sumcomm_pc <- comm_pc %>% group_by(wetland_id) %>% summarise_each(funs(sum))
sumcomm_pc$method <- "ptct"
length(sumcomm_pc) <- length(sumcomm)
as.data.frame(sum)
occ_both <- rbind(sumcomm_pc, sumcomm)
chisq.test(comm_both$value)


#3rd analysis: ptct upland vs lowland, acoustic up vs low, mean combined of the two
#if acoustic can be brought into check more cleanly with standard set of bird filters: 
# papers by ethan white and allen hurlburt BBL (2005 rangemaps)

#pattern across ptct vs acoustic (synthetic analysis combines these two points: 
# avg occupancy for each species)


#Q2:
#compare across method - sum diversity rather than mean diversity 
boxplot(comm_both$value~ comm_both$method, 
        main = "Average Species Richness by Monitoring Method", 
        xlab = "Monitoring Method", ylab = "Average Species Richness")
aov_one.way <- aov(comm_both$value ~ comm_both$method)
summary(aov_one.way)
aov_two.way <- aov(comm_both$value ~ comm_both$method + comm_both$habitat)
summary(aov_two.way)

sum(comm$value)
sum(comm_ptct$value)

AIC(aov_one.way)
AIC(aov_two.way)

tukey.one.way <- TukeyHSD(aov_one.way)
tukey.two.way <- TukeyHSD(aov_two.way)

tukey.one.way
tukey.two.way

#compare across wetland ID
boxplot(comm_both$value ~ comm_both$wetland_id)
summary(aov(comm_both$value ~ comm_both$wetland_id))

#compare across habitat
comm_both$habitat <- ifelse(substring(comm_both$wetland_id, 1, 1) == 'U', 'upland', 'wetland')

boxplot(comm_both$value ~ comm_both$habitat)
summary(aov(comm_both$value ~ comm_both$habitat))

#add dat_pc variables to occupancy matrix (block, veg data, etc.)
block <- dat_pc$block
comm_both <- tapply(comm_both, block) ###!
block <- ifelse(substring(dat_pc$wetland_id) == comm_both$wetland_id, 
                paste(dat_pc$block))
block

#rda to assess the variation among habitat & method

bird_rda <- rda(comm_both$value ~ method + habitat, 
               comm_both)
bird_rda
RsquareAdj(bird_rda)

plot(bird_rda)


#need to fuse dat_pc data to comm_both compare across both comm matrixes? Block, tree_sr, 
# tree_dbh, canopy_cover (Use wetland_id as unique id).
comm_both <- tapply() ###!

#next: Dig into details & see what, why, where. Test upland vs wetlands to see diff
#between habitat type - open vs veg, vary depending on ptct radius or conf level.
#study scale shows that methods are slightly redundant but for most part they are
#detecting diff community comp. 
#(Should we remove all species but songbirds? One more layer of processing & see if 
#that corrects the axes).
#which species did Jackson miss that I picked up? Can go back & validate by listening.
#remove the species w 0 observations
#NEED A PLOT THAT SHOWS DIFF BETWEEN ACOUSTIC SR & POINT COUNT SR
# diversity accumulation; point counts miss rare things that come through (rare tail)
# separate analysis pathway. 

#examine the redundancy vs complementary info in acoustic methods.
