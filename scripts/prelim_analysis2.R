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

# fix times on a few sites where clock was wrong
# HH01 needs 4 hours subtracted
dat_sub$start_time[dat_sub$site_id == "HH01"] <- dat_sub$start_time[dat_sub$site_id == "HH01"] - 4*60*60
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

comm <- with(dat_sub, tapply(Confidence, list(site_date, names),
                        function(x) length(x)))
comm <- ifelse(is.na(comm), 0, comm)
sum(comm)

comm[, 1:5]
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

names(comm)

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
comm_pc <- dat_pc[ , c('wetland_id', 'date.x', names(dat_pc)[12:66])]
comm_pc[1:5, 1:5]
sr_pc <- rowSums(comm_pc[ , -(1:2)] > 0)
plot(density(sr_pc))

sr_avg <- tapply(sr_pc, comm_pc$wetland, mean)
sr_avg
plot(sr_avg)


#transform into data frame
comm_ptct <- as.data.frame(comm_ptct)

sr_avg_df <- data.frame(wetland_id = names(sr_avg), S = sr_avg)
sr_site_id_df <- data.frame(wetland_id = names(sr_site_id), S = sr_site_id)

sr_merge <- merge(sr_avg_df, sr_site_id_df, by = 'wetland_id', all = TRUE)

plot(S.y ~ S.x, data = sr_merge)
abline(a = 0 , b=1)

#add column to dataframes indicating methods
comm_ptct %>% mutate(method = "ptct")
comm %>% mutate(method = "acoustic")

#merge community matrixes
comm_both <- merge(comm, comm_ptct, by = [1, ] all = TRUE)
sr_both <- rowSums(comm_both > 0)
#subset time column to get the dawn chorus comparisons between my data & jackson's.

#NEED A PLOT THAT SHOWS DIFF BETWEEN ACOUSTIC SR & POINT COUNT

boxplot(comm_both$CHNA ~ comm_both$method)

#Then redefine 10 min time units and split into 2. Every start time go 5, and that
#becomes unit.



# diversity accumulation; point counts miss rare things that come through (rare tail)
# separate analysis pathway. 



