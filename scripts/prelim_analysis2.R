library(vegan)
library(readr)
library(dplyr)
library(easyCODA)

comm$CHNA

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
sr_site_id <- tapply(sr, site_id, mean)
srlog <- log(sr)

boxplot(srlog ~ habitat, ylab = "log(Species Richness)", xlab="Habitat", 
        main="Species Richness by Habitat")
boxplot(sr~site_id, ylab = "Species Richness", xlab="Site ID", 
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
sr_site_id_df <- data.frame(wetland_id = names(sr_site_id), S = sr_site_id)

sr_merge <- merge(sr_avg_df, sr_site_id_df, by = 'wetland_id', all = TRUE)

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

S_acou <- with(sp_occ, tapply(acoustic, list(wetland_id), sum))
S_ptct <- with(sp_occ, tapply(ptct, list(wetland_id), sum))

V_acou <- with(sp_occ, tapply(acoustic, list(wetland_id), 
                              function(x) sum(x * (1 - x))))
V_ptct <- with(sp_occ, tapply(ptct, list(wetland_id),
                              function(x) sum(x * (1 - x))))

par(mfrow=c(1,2))
plot(S_ptct, S_acou)
abline(a = 0, b=1)
plot(V_ptct, V_acou)
abline(a = 0, b=1)

sr_both <- rowSums(comm_both > 0)
sr_both_avg <- tapply(sr_both, comm_both$wetland, mean)
length(sr_both)
length(comm_both)
plot(sr_both_avg)

boxplot(comm_both$value~ comm_both$method)

upland <- comm_both$wetland_id
boxplot(comm_both$value ~ comm_both$wetland_id)

comm_both$habitat <- ifelse(substring(comm_both$wetland_id, 1, 1) == 'U', 'upland', 'wetland')

boxplot(comm_both$value ~ comm_both$habitat)

#next: Dig into details & see what, why, where. Test upland vs wetlands to see diff
#between habitat type - open vs veg, vary depending on ptct radius or conf level.
#study scale shows that methods are slightly redundant but for most part they are
#detecting diff community comp. 
#(Should we remove all species but songbirds? One more layer of processing & see if 
#that corrects the axes).
#which species did Jackson miss that I picked up? Go back & validate by listening.

#NEED A PLOT THAT SHOWS DIFF BETWEEN ACOUSTIC SR & POINT COUNT SR
# diversity accumulation; point counts miss rare things that come through (rare tail)
# separate analysis pathway. 

#build occupancy matrix for both mine and jackson's species
#proceed as if we dropped the wonky guys 

#examine the redundancy vs complementary info in acoustic methods.
