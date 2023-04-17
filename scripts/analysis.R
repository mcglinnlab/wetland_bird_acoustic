library(vegan)
library(readr)
library(dplyr)
library(easyCODA)

#Read in data

comm <- read.csv('./data/clean_data/comm.csv')
dat_sub <- read.csv('./data/clean_data/subset_data.csv')

comm$CHNA
comm$BGGN

#Chuck-will's widow - NIGHTIME
write.csv(row.names(comm[comm$CHNA > 1,]), file = './data/CHNA_detections.csv')

# distribution of species richness
hist(rowSums(comm > 0), xlab= "# of Observations within 5 min interval",
     main = "Distribution of Species Richness")

#subsample out in equal amounts of upland and wetlands or statistically model for
#the temporal and spatial autocorrelation. pseudo-replication.

# look at species ranks based upon number of occurrences
sort(colSums(comm > 0), dec = T)

row.names(comm)
habitat <- ifelse(substring(comm$wetland_id, 1, 1) == 'U', 'upland', 'wetland')
#need to assign Halidon vs stono if we care to compare
#property <- ifelse(substring(site_id, 1, 1) == 'H', 'Halidon', 'Stono')

table(habitat)

sr <- rowSums(comm > 0)
plot(density(sr))
sr_wetland_id <- tapply(sr, wetland_id, mean)
srlog <- log(sr)

boxplot(srlog ~ habitat, ylab = "log(Species Richness)", xlab="Habitat", 
        main="Species Richness by Habitat")
boxplot(sr~wetland_id, ylab = "Species Richness", xlab="Wetland ID", 
        main="Species Richness per Site")
boxplot(sr~property_wet, ylab = "Species Richness", xlab="Wetland Sites by Property",
        main="Wetland SR by Property")

hist(sr, xlab="Species Richness", main = "Histogram of Species Richness")

#separate comm matrix by habitat type
library(dplyr)
srhab <- cbind(sr,habitat)
srhab <- as.data.frame(srhab)
wetsr <- srhab %>% filter(habitat=='wetland', na.rm=TRUE)
upsr <- srhab %>% filter(habitat=='upland', na.rm=TRUE)

dat_pc <- read.csv('https://raw.githubusercontent.com/mcglinnlab/wetland_birds/main/data/filtered_data/clean_bird_dat.csv')
#note: in dat_pc from Jackson that wetland_id is equivalent to what I've called
# site_id; therefore let's rename wetland_id to site_id for simplicity
#jackson's data is subset from 0y-25 m rather than previous total which was 50m.
head(dat_pc)
comm_pc <- dat_pc[ , c('wetland_id', names(dat_pc)[12:66])]
comm_pc[1:5, 1:5]
sr_pc <- rowSums(comm_pc[ , -(1)] > 0)
plot(density(sr_pc))

sr_avg <- tapply(sr_pc, comm_pc$wetland, mean)
sr_avg
plot(sr_avg)

#transform into data frame
comm_ptct <- as.data.frame(comm_pc)

sr_avg_df <- data.frame(wetland_id = names(sr_avg), S = sr_avg)
sr_wetland_id_df <- data.frame(wetland_id = names(sr_wetland_id), S = sr_wetland_id)

sr_merge <- merge(sr_avg_df, sr_wetland_id_df, by = 'wetland_id', all = TRUE)

plot(S.y ~ S.x, data = sr_merge)
abline(a = 0 , b=1)

#compress point counts and acoustic into a single value per sample, tapply or aggregate
#merge & get 2 NOCA's- one for PC & one for acoustic
comm_ptct <- aggregate(comm_ptct, by = list(comm_ptct$wetland_id), function(x) sum(x > 0) / length(x))
# drop wetland id and rename group 1 to wetland id
names(comm_ptct)
comm_ptct = subset(comm_ptct, select = -wetland_id)
names(comm_ptct)= c("wetland_id", names(comm_ptct) [-1])
names(comm_ptct)
## wetland id pivots wrong.

#pivot data frame to long format
library(tidyr)
comm <- pivot_longer(comm, cols = !wetland_id)
comm_ptct <- pivot_longer(comm_ptct, cols = !wetland_id)

#add column to dataframes indicating methods
comm_ptct$method <- "ptct"
comm$method <- "acoustic"

#row append community matrices
comm_both <- rbind(comm, comm_ptct)

# drop sites not found in both sampling methods 
acoustic_sites <- unique(comm_both$wetland_id[comm_both$method == 'acoustic'])
sites_with_both <- acoustic_sites[acoustic_sites %in% unique(comm_both$wetland_id[comm_both$method == 'ptct'])]
sites_with_both

comm_both <- subset(comm_both, wetland_id %in% sites_with_both)

sp_occ <- pivot_wider(comm_both, names_from = 'method', values_from = value)
sp_occ$acoustic <- ifelse(is.na(sp_occ$acoustic), 0, sp_occ$acoustic)
sp_occ$ptct <- ifelse(is.na(sp_occ$ptct), 0, sp_occ$ptct)

plot(acoustic ~ ptct, data = sp_occ, type = 'n')
with(sp_occ, text(ptct, acoustic, labels = name, cex = 0.5))
abline(a = 0 , b= 1)

sp_occ_agg <- sp_occ %>% 
    group_by(name) %>% 
    summarize(acoustic = mean(acoustic), ptct = mean(ptct))

plot(acoustic ~ ptct, data = sp_occ_agg, type = 'n')
with(sp_occ_agg, text(ptct, acoustic, labels = name, cex = 0.5))
abline(a = 0 , b= 1)

#compute sr 
S_acou <- with(sp_occ, tapply(acoustic, list(wetland_id), sum))
S_ptct <- with(sp_occ, tapply(ptct, list(wetland_id), sum))

#compute variance across methods (p*q, q= 1-p) ##- sig correlation between ptct
#sr and acoustic richness???
V_acou <- with(sp_occ, tapply(acoustic, list(wetland_id), 
                              function(x) sum(x * (1 - x))))
V_ptct <- with(sp_occ, tapply(ptct, list(wetland_id),
                              function(x) sum(x * (1 - x))))

par(mfrow=c(1,2))
plot(S_ptct, S_acou)
abline(a = 0, b=1)
plot(V_ptct, V_acou)
abline(a = 0, b=1)

#Which species were only detected acoustic/ptct - which fit best? Filter out 
#non-songbird

#find unique species detected by one method & not the other
library(dplyr)
num_species <- comm_both %>% count(name)
#77 unique species names
unique_species <- unique(comm_both$name[comm_both$method])
as.data.frame(acoustic_species) <- list(comm_both$name[comm_both$method == 'acoustic'])
as.data.frame(ptct_species) <- list(comm_both$name[comm_both$method == 'ptct'])
unique_species <- setdiff(acoustic_species, ptct_species)
?setdiff


#Q1: SR comparison: sum diversity aov

#fig 1 SR comparison fig 2 species by species comparison
#fig 3 descriptive nature of why are these differing? could be table of species
#traits (large body, quiet singing). QI sq tests? We have 2 numbers are they more \
#diff than expected?

#3rd analysis: ptct upland vs lowland, acoustic up vs low, mean combined of the two
#if acoustic can be brought into check more cleanly with standard set of bird filters: 
# papers by ethan white and allen hurlburt BBL (2005 rangemaps)

#pattern across ptct vs acoustic (synthetic analysis combines these two points: 
# avg occupancy for each species)


#Q2:
#compare across method - sum diversity rather than mean diversity 
boxplot(comm_both$value~ comm_both$method)
summary(aov(comm_both$value ~ comm_both$method))

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

bird_rda <- rda(comm_both$value ~ habitat + method, 
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
