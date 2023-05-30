
# read in data and reformat -----------
dat <- read.csv('./data/compiled_bird_audio_survey_master.csv')

# fix upland sites names: 
# UPred -> UP05
# UPblue -> UP06
dat$site_id[dat$site_id == "UPred"] <- "UP05"
dat$site_id[dat$site_id == "UPblue"] <- "UP06"

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

# format date fields properly
dat$date_time <- as.POSIXlt(dat$date_time, tz = 'EST',
                            format = "%m/%d/%Y %H:%M")
dat$start_time <-  as.POSIXlt(dat$start_time, tz = 'EST',
                                  format = "%m/%d/%Y %H:%M")

# import the matrix of sampling times 
pct_times <- read.csv('./data/raw_data - ptct_acoustic_sample_times.csv')
pct_times$time <- with(pct_times, paste(date, time))
pct_times$time <- as.POSIXlt(pct_times$time, tz = 'EST',
                             format = "%m/%d/%Y %H:%M")
pct_times$start_time <- as.POSIXlt(pct_times$start_time, tz = 'EST',
                                   format = "%m/%d/%Y %H:%M")
head(pct_times)

# fix times 
fix_time <- function(wrong_times, ref_time) {
  first_wrong_time <- min(wrong_times)
  time_deviation <- first_wrong_time - ref_time
  corrected_times <- wrong_times - time_deviation
  return(corrected_times)
}

site_ids <- unique(dat$site_id)
sort(site_ids)

# need to double check if SP08 and UP05 need time corrections

for(i in seq_along(site_ids)) { 
  if (site_ids[i] %in% c("SP01", "SP08", "SP09", "UP05")) {
    # 4 hours subtracted
    dat$start_time[dat$site_id == site_ids[i]] <- dat$start_time[dat$site_id == site_ids[i]] - 4*60*60
    dat$date_time[dat$site_id == site_ids[i]] <- dat$date_time[dat$site_id == site_ids[i]] - 4*60*60
  } else if (site_ids[i] == "SP02") {
      # SP02 needs 4 hours added! (real date is 6/14/22)
      dat$start_time[dat$site_id == "SP02"] <- dat$start_time[dat$site_id == "SP02"] + 4*60*60
      dat$date_time[dat$site_id == "SP02"] <- dat$date_time[dat$site_id == "SP02"] + 4*60*60
  } else { 
     dat$start_time[dat$site_id == site_ids[i]] <- fix_time(dat$start_time[dat$site_id == site_ids[i]],
                                                            pct_times$time[pct_times$wetland_id == site_ids[i]])
     dat$date_time[dat$site_id == site_ids[i]] <- fix_time(dat$date_time[dat$site_id == site_ids[i]],
                                                           pct_times$time[pct_times$wetland_id == site_ids[i]])
   }
}

# survey start time 
# must only reference time of day and not date

start_hour <- 7
end_hour <- 8

rec_hours <- as.numeric(format(dat$start_time, format = "%H"))

# subset dat to hours between the start_hour and end_hour
dat <- dat[rec_hours >= start_hour & rec_hours < end_hour, ]

# now adjust start times so that they are in 5 min intervals
for (i in 1:nrow(dat)) {
  if ((dat$start_time[i]  - dat$date_time[i]) > 4) {
    dat$date_time[i] <- dat$date_time[i] + 5*60
  }
}

# make unique id that has site_id and date_time
dat$site_date <- with(dat, paste(site_id, date_time, sep = '_'))
comm <- with(dat, tapply(Confidence, list(site_date, sp_code),
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
write.csv(dat, file = "./data/clean_data/dat.csv", row.names = FALSE)
