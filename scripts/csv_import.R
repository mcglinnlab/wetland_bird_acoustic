# script for importing audio data csv

dirs_int <- c("HH01","HH02", "HH04")
file_names <- dir(dirs_int[1], recursive = TRUE)
csv_files <- file_names[grep('.csv', file_names)]

extract_time <- function(file_name) {
    x <- strsplit(file_name, split = '_')
    y <- strsplit(x[[1]][2], '.', fixed = TRUE)
    z <- y[[1]][1]
    time_out <- gsub('-', ':', z, fixed = TRUE)
    return(time_out)
}

extract_date <- function(file_name) { 
    x <- strsplit(file_name, split = '/', fixed = TRUE)
    date_out <- x[[1]][2]
    return(date_out)
}

extract_site_id <- function(file_name) { 
    x <- strsplit(file_name, split = '/', fixed = TRUE)
    site_out <- x[[1]][1]
    return(site_out)
}

# 
dir()

dirs_int <- c("HH01","HH02", "HH04", "HH05", "HH14",
              "HH17", "HH30", "HH31", "HH33", "HH34",
              "HH34R", "HH45", "HH46", "HH48", "HH49",
              "HH50", "SP01", "SP02", "SP03", "SP05",
              "SP07", "SP08", "SP09", "SP11", "SP13",
              "SP16", "UP01", "UP03", "UP04", "UPblue", "UPred",
              "Urban 1", "Urban 2")
              

# create empty data.frame
dat <- data.frame()
# loop through site ids
for (j in seq_along(dirs_int)) {
  print(paste(j, 'of', length(dirs_int)))
  file_names <- dir(dirs_int[j], recursive = TRUE)
  # isolate only csv file names from .wav files or others
  csv_files <- file_names[grep('.csv', file_names)]
  # loop through all files for a particular site
  for (i in seq_along(csv_files)) {
    tmp_file_name <- paste(dirs_int[j], csv_files[i], sep = '/')
    tmp <-  read.csv(tmp_file_name)
    if (nrow(tmp) > 0) {
      date_of_sample <- extract_date(tmp_file_name)
      time_of_day <- extract_time(tmp_file_name)
      day_time <- paste(date_of_sample, time_of_day)
      site_id <- extract_site_id(tmp_file_name)
      tmp <- data.frame(site_id,
                        date_time = day_time,
                        date = date_of_sample,
                        start_time = as.POSIXct(day_time, tz = 'EST', format = "%Y-%m-%d %H:%M") + tmp$Start..s.,
                        end_time = as.POSIXct(day_time, tz = 'EST', format = "%Y-%m-%d %H:%M") + tmp$End..s.,
                        tmp[ , -(1:2)])
      dat <- rbind(dat, tmp)
    }
  }    
}


dim(dat)
head(dat)

write.csv(dat, file = './data/compiled_bird_audio_survey_updated.csv', row.names = FALSE)
