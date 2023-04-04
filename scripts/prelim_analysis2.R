library(vegan)
library(readr)
library(dplyr)
library(easyCODA)

make_sp_codes <- function(x) {
    sp_genus <- toupper(substring(x, 1, 2))
    sp_sp    <- toupper(sapply(strsplit(x, ' ', fixed = TRUE),
                               function(x)substring(x[2], 1, 2)))
    sp_code <- paste(sp_genus, sp_sp, sep='')
    sp_code
}

# read in data and reformat -----------
dat <- read.csv('./data/compiled_bird_audio_survey_new1.csv')
#use this once formatting has been updated.
dat <- read.csv('./data/compiled_bird_audio_survey_master.csv')

dat$names <- make_sp_codes(dat$Common.name)

dat$site_time <- paste(dat$site_id, dat$date_time, sep ='_')

hab <- ifelse(grepl('U', dat$site_id), 'upland', 'wetland')
dat <- cbind(dat, hab)
#within 5 min interval, on avg seeing 5 species in upland and 3/4 in wetland with a
#lot of outliers.

# subset any confidence less than 0.5
dat_sub <- subset(dat, Confidence > 0.5)

dat_sub$date_time <- as.POSIXlt(dat_sub$date_time, tz = 'EST',
                                format = "%m/%d/%Y %H:%M")
dat_sub$start_time <-  as.POSIXlt(dat_sub$start_time, tz = 'EST',
                                  format = "%m/%d/%Y %H:%M")

# fix times of a few sites where clock was wrong
# HH01 needs 4 hours subtracted
dat_sub$start_time[dat_sub$site_id == "HH01"] <- dat_sub$start_time[dat_sub$site_id == "HH01"] - 4*60*60
# HH02 needs ? hours subtracted
dat_sub$start_time[dat_sub$site_id == "HH02"] <- dat_sub$start_time[dat_sub$site_id == "HH02"] - 4*60*60
# HH04 needs ? hours subtracted
dat_sub$start_time[dat_sub$site_id == "HH04"] <- dat_sub$start_time[dat_sub$site_id == "HH04"] - 4*60*60
# HH05 needs ? hours subtracted
dat_sub$start_time[dat_sub$site_id == "HH05"] <- dat_sub$start_time[dat_sub$site_id == "HH05"] - 4*60*60
# SP02 needs ? hours subtracted
dat_sub$start_time[dat_sub$site_id == "SP02"] <- dat_sub$start_time[dat_sub$site_id == "SP02"] - 4*60*60
# SP09 needs 4 hours subtracted
dat_sub$start_time[dat_sub$site_id == "SP09"] <- dat_sub$start_time[dat_sub$site_id == "SP09"] - 4*60*60

# sunup
dawn_times <- as.POSIXlt(paste(dat_sub$date, "6:00"), tz = 'EST',
                         format = "%m/%d/%Y %H:%M")
# survey end time
end_times <- as.POSIXlt(paste(dat_sub$date, "9:00"), tz = 'EST',
                        format = "%m/%d/%Y %H:%M")
# subset dat_sub to dawn chorus times 
dat_sub <- dat_sub[dat_sub$start_time > dawn_times & dat_sub$start_time < end_times, ]

# now adjust start times so that they are in 5 min intervals
for (i in 1:nrow(dat_sub)) {
    if ((dat_sub$start_time[i]  - dat_sub$date_time[i]) > 4) {
        dat_sub$date_time[i] <- dat_sub$date_time[i] + 5*60
    }
}

# make unique id that has site_id and date_time
dat_sub$site_date <- with(dat_sub, paste(site_id, date_time, sep = '_'))

# make confidence matrix - take the max confidence across all observations in that
# site at that date
comm <- with(dat_sub, tapply(Confidence, list(site_date, names),
                             function(x) any(x > 0) * 1))
comm <- ifelse(is.na(comm), 0, comm)
sum(comm)

comm[1:5, 1:5]

wetland_id <- sapply(strsplit(row.names(comm), split = '_', fixed = TRUE), function(x) x[1])
comm <- aggregate(comm, by = list(wetland_id), function(x) sum(x) / length(x))
comm[1:5, 1:5]
summary(comm)
names(comm)
#rename group 1 to wetland_id
names(comm) = c("wetland_id", names(comm) [-1])

# convert to data.frame
comm <- as.data.frame(comm)

comm$CHNA

#Chuck-will's widow - NIGHTIME
write.csv(row.names(comm[comm$CHNA > 1,]), file = './data/CHNA_detections.csv')

N <- rowSums(comm)
S <- rowSums(comm > 0)
comm[N == 83, ]
cbind(unique(dat$Common.name), make_sp_codes(unique(dat$Common.name)))

plot(density(S))

# distribution of species richness
hist(rowSums(comm > 0), xlab= "# of Observations within 5 min interval",
     main = "Distribution of Species Richness")

#subsample out in equal amounts of upland and wetlands or statistically model for
#the temporal and spatial autocorrelation. pseudo-replication.

# look at species ranks based upon number of occurrences
sort(colSums(comm > 0), dec = T)

row.names(comm)
site_id <- sapply(strsplit(row.names(comm), '_'), function(x) x[1])
date_time <- sapply(strsplit(row.names(comm), '_'), function(x) x[2])
habitat <- ifelse(substring(site_id, 1, 1) == 'U', 'upland', 'wetland')
property_wet <- ifelse(substring(site_id, 1, 1) == 'H', 'Halidon', 'Stono')

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

boxplot(comm_both$GRCR ~ comm_both$method)

#NEED A PLOT THAT SHOWS DIFF BETWEEN ACOUSTIC SR & POINT COUNT SR
# diversity accumulation; point counts miss rare things that come through (rare tail)
# separate analysis pathway. 

#build occupancy matrix for both mine and jackson's species
#proceed as if we dropped the wonky guys 

#examine the redundancy vs complementary info in acoustic methods.
