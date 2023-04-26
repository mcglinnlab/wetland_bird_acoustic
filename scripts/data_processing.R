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
dat <- read.csv('./data/compiled_bird_audio_survey_master.csv')

dat$names <- make_sp_codes(dat$Common.name)

dat$site_time <- paste(dat$site_id, dat$date_time, sep ='_')

hab <- ifelse(grepl('U', dat$site_id), 'upland', 'wetland')
dat <- cbind(dat, hab)

acou_splist <- unique(dat$Scientific.name)

#import BBL species codes & merge
bbl_sp_codes <- read.csv("data/IBP-AOS-LIST22.csv")
all(acou_splist %in% bbl_sp_codes$SCINAME)
acou_spcode <- bbl_sp_codes$SPEC[match(acou_splist, bbl_sp_codes$SCINAME)]
acou_code <- data.frame(Scientific.name = acou_splist, 
                        sp_code = acou_spcode)
head(acou_code)

dat <- merge(dat, acou_code, by = "Scientific.name")
head(dat)

# subset any confidence less than 0.5
dat_sub <- subset(dat, Confidence > 0.5)

dat_sub$date_time <- as.POSIXlt(dat_sub$date_time, tz = 'EST',
                                format = "%m/%d/%Y %H:%M")
dat_sub$start_time <-  as.POSIXlt(dat_sub$start_time, tz = 'EST',
                                  format = "%m/%d/%Y %H:%M")

# now adjust start times so that they are in 5 min intervals
for (i in 1:nrow(dat_sub)) {
    if ((dat_sub$start_time[i]  - dat_sub$date_time[i]) > 4) {
        dat_sub$date_time[i] <- dat_sub$date_time[i] + 5*60
    }
}

# fix times of a few sites where clock was wrong
# HH01 needs 4 hours subtracted
dat_sub$start_time[dat_sub$site_id == "HH01"] <- dat_sub$start_time[dat_sub$site_id == "HH01"] - 4*60*60
# HH02 needs 4 hours subtracted
dat_sub$start_time[dat_sub$site_id == "HH02"] <- dat_sub$start_time[dat_sub$site_id == "HH02"] - 4*60*60
# HH04 needs 4 hours subtracted
dat_sub$start_time[dat_sub$site_id == "HH04"] <- dat_sub$start_time[dat_sub$site_id == "HH04"] - 4*60*60
# HH05 needs 4 hours subtracted
dat_sub$start_time[dat_sub$site_id == "HH05"] <- dat_sub$start_time[dat_sub$site_id == "HH05"] - 4*60*60
#HH14 needs 4 hours subtracted
dat_sub$start_time[dat_sub$site_id == "HH14"] <- dat_sub$start_time[dat_sub$site_id == "HH14"] - 4*60*60
#HH17 needs 4 hours subtracted
dat_sub$start_time[dat_sub$site_id == "HH17"] <- dat_sub$start_time[dat_sub$site_id == "HH17"] - 4*60*60
#HH45 needs 4 hours subtracted
dat_sub$start_time[dat_sub$site_id == "HH45"] <- dat_sub$start_time[dat_sub$site_id == "HH45"] - 4*60*60
#HH46 needs 4 hours subtracted
dat_sub$start_time[dat_sub$site_id == "HH46"] <- dat_sub$start_time[dat_sub$site_id == "HH46"] - 4*60*60
#HH48 needs 4 hours subtracted
dat_sub$start_time[dat_sub$site_id == "HH48"] <- dat_sub$start_time[dat_sub$site_id == "HH48"] - 4*60*60
#SP01 needs 4 hours subtracted
dat_sub$start_time[dat_sub$site_id == "SP01"] <- dat_sub$start_time[dat_sub$site_id == "SP01"] - 4*60*60
# SP02 needs 4 hours added! (real date is 6/14/22)
dat_sub$start_time[dat_sub$site_id == "SP02"] <- dat_sub$start_time[dat_sub$site_id == "SP02"] + 4*60*60
# SP09 needs 4 hours subtracted
dat_sub$start_time[dat_sub$site_id == "SP09"] <- dat_sub$start_time[dat_sub$site_id == "SP09"] - 4*60*60
#SP05 needs 4 hours subtracted
dat_sub$start_time[dat_sub$site_id == "SP05"] <- dat_sub$start_time[dat_sub$site_id == "SP05"] - 4*60*60
#UP04 needs 4 hours subtracted
dat_sub$start_time[dat_sub$site_id == "UP04"] <- dat_sub$start_time[dat_sub$site_id == "UP04"] - 4*60*60

#fix date time
# HH01 needs 4 hours subtracted
dat_sub$date_time[dat_sub$site_id == "HH01"] <- dat_sub$date_time[dat_sub$site_id == "HH01"] - 4*60*60
# HH02 needs 4 hours subtracted
dat_sub$date_time[dat_sub$site_id == "HH02"] <- dat_sub$date_time[dat_sub$site_id == "HH02"] - 4*60*60
# HH04 needs 4 hours subtracted
dat_sub$date_time[dat_sub$site_id == "HH04"] <- dat_sub$date_time[dat_sub$site_id == "HH04"] - 4*60*60
# HH05 needs 4 hours subtracted
dat_sub$date_time[dat_sub$site_id == "HH05"] <- dat_sub$date_time[dat_sub$site_id == "HH05"] - 4*60*60
#HH14 needs 4 hours subtracted
dat_sub$start_time[dat_sub$site_id == "HH14"] <- dat_sub$date_time[dat_sub$site_id == "HH14"] - 4*60*60
#HH17 needs 4 hours subtracted
dat_sub$start_time[dat_sub$site_id == "HH17"] <- dat_sub$date_time[dat_sub$site_id == "HH17"] - 4*60*60
#HH45 needs 4 hours subtracted
dat_sub$start_time[dat_sub$site_id == "HH45"] <- dat_sub$date_time[dat_sub$site_id == "HH45"] - 4*60*60
#HH46 needs 4 hours subtracted
dat_sub$start_time[dat_sub$site_id == "HH46"] <- dat_sub$date_time[dat_sub$site_id == "HH46"] - 4*60*60
#HH48 needs 4 hours subtracted
dat_sub$start_time[dat_sub$site_id == "HH48"] <- dat_sub$date_time[dat_sub$site_id == "HH48"] - 4*60*60
#SP01 needs 4 hours subtracted
dat_sub$start_time[dat_sub$site_id == "SP01"] <- dat_sub$date_time[dat_sub$site_id == "SP01"] - 4*60*60
# SP02 needs 4 hours added (real date is 6/14/22)
dat_sub$date_time[dat_sub$site_id == "SP02"] <- dat_sub$date_time[dat_sub$site_id == "SP02"] + 4*60*60
# SP09 needs 4 hours subtracted
dat_sub$date_time[dat_sub$site_id == "SP09"] <- dat_sub$date_time[dat_sub$site_id == "SP09"] - 4*60*60
#SP05 needs 4 hours subtracted
dat_sub$start_time[dat_sub$site_id == "SP05"] <- dat_sub$date_time[dat_sub$site_id == "SP05"] - 4*60*60
#UP04 needs 4 hours subtracted
dat_sub$start_time[dat_sub$site_id == "UP04"] <- dat_sub$date_time[dat_sub$site_id == "UP04"] - 4*60*60

# sunup
dawn_times <- as.POSIXlt(paste(dat_sub$date, "6:00"), tz = 'EST',
                         format = "%m/%d/%Y %H:%M")
# survey end time
end_times <- as.POSIXlt(paste(dat_sub$date, "9:00"), tz = 'EST',
                        format = "%m/%d/%Y %H:%M")

# subset dat_sub to dawn chorus times 
dat_sub <- dat_sub[dat_sub$start_time > dawn_times & dat_sub$start_time < end_times, ]

# make unique id that has site_id and date_time
dat_sub$site_date <- with(dat_sub, paste(site_id, date_time, sep = '_'))
comm <- with(dat_sub, tapply(Confidence, list(site_date, sp_code),
                             function(x) any(x > 0) * 1))
comm <- ifelse(is.na(comm), 0, comm)
sum(comm)

comm[, 1:5]

#avg species richness by wetland_id
wetland_id <- sapply(strsplit(row.names(comm), split = '_', fixed = TRUE), function(x) x[1])
sumcomm <- aggregate(comm, by = list(wetland_id), function(x) sum(x))
comm <- aggregate(comm, by = list(wetland_id), function(x) sum(x) / length(x))
#rename group 1 to wetland_id
names(comm) = c("wetland_id", names(comm) [-1])
names(comm)

# convert to data.frame
comm <- as.data.frame(comm)
sumcomm <- as.data.frame(sumcomm)

#export data for analysis
write.csv(sumcomm, file = "./data/clean_data/sumcomm.csv", row.names = FALSE)
write.csv(comm, file = "./data/clean_data/comm.csv", row.names = FALSE)
write.csv(dat_sub, file = "./data/clean_data/subset_data.csv", row.names = FALSE)
