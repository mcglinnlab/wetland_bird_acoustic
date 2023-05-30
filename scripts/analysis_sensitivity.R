library(ggplot2)
library(dplyr)
library(tidyr)
library(mobr)

source('./scripts/bootstrap_functions.R')

#Read in data

sumcomm <- read.csv('./data/clean_data/sumcomm.csv') # comm matrix when all detections are totaled
comm <- read.csv('./data/clean_data/comm.csv')       # comm matrix when detections are converted into proportions
dat <- read.csv('./data/clean_data/dat.csv') #the cleaned acoustic data with confidence levels and times

# consider dropping SP16 because only has 168 observations
#dat <- subset(dat, site_id != "SP16")
# consider dropping UP05 only has 45 observations
#dat <- subset(dat, site_id != "UP05")

pdf('./figs/condience_histogram.pdf')
hist(dat$Confidence, xlab = 'Species Detection Confidence',
     breaks = seq(0, 1, .1))
dev.off()


# read in point count data
dat_pc <- read.csv('https://raw.githubusercontent.com/mcglinnlab/wetland_birds/main/data/filtered_data/clean_bird_dat.csv')
comm_pc <- dat_pc[ , c('wetland_id', names(dat_pc)[12:66])]
#Bachman's sparrow code is incorrect in ptct data. Switch BASP to BACS
colnames(comm_pc)[21] <- "BACS"
comm_pc = subset(comm_pc, select = -c(BGNN))
comm_pc$BASP
comm_pc$BCSP
comm_pc$BACS = comm_pc$BACS + comm_pc$BCSP
#drop comm_pc$BCSP
comm_pc = subset(comm_pc, select = -c(BCSP))

#create sum ptct
sumcomm_ptct <- comm_pc %>% group_by(wetland_id) %>% summarise_each(funs(sum))

comm_ptct <- as.data.frame(comm_pc)

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

acoustic_sites <- unique(comm_both$wetland_id[comm_both$method == 'acoustic'])
sites_with_both <- acoustic_sites[acoustic_sites %in% unique(comm_both$wetland_id[comm_both$method == 'ptct'])]
sites_with_both

# Consider dropping upland sites
#sites_with_both <- sites_with_both[!grepl('UP', sites_with_both)]

#create combined community matrix
comm_both <- subset(comm_both, wetland_id %in% sites_with_both)

# remove value of zero
comm_both <- subset(comm_both, value > 0)

#create occupancy matrix
sp_occ <- pivot_wider(comm_both, names_from = 'method', values_from = value)
sp_occ$acoustic <- ifelse(is.na(sp_occ$acoustic), 0, sp_occ$acoustic)
sp_occ$ptct <- ifelse(is.na(sp_occ$ptct), 0, sp_occ$ptct)


# species level occupancy analysis --------------------------------------

occ_results <- data.frame()
confi_vals <- seq(0.1, 0.9, 0.1)
nboots <- 1e2
for(i in seq_along(confi_vals)) { 
  # subset based upon confidence lv
  dat_sub <- subset(dat, Confidence >= confi_vals[i])
  
  comm_aco <- with(dat_sub, tapply(Confidence, INDEX = list(site_date, sp_code),
                                   function(x) (sum(x) > 0)*1))
  comm_aco <- ifelse(is.na(comm_aco), 0, comm_aco)
  comm_aco_sites <- substr(rownames(comm_aco), 1, 4)
  
  # drop upland sites
  #comm_aco <- subset(comm_aco, !grepl('UP', comm_aco_sites))
  #comm_aco_sites <- substr(rownames(comm_aco), 1, 4)

  comm_aco <- as.data.frame(comm_aco)
  
  aco_boots <- replicate(nboots, boot_occ(comm_aco, grp = comm_aco_sites), simplify = FALSE)
  aco_sp_q <- apply(sapply(aco_boots, function(x) apply(x, 2, mean)), 1,
                    quantile, c(0.025, 0.5, 0.975))
  aco_sp_a <- apply(sapply(aco_boots, function(x) apply(x, 2, mean)), 1, mean)
  
  occ_results <- rbind(occ_results, 
                       data.frame(sp = names(aco_sp_a), occ_avg = aco_sp_a, 
                                  occ_lo = aco_sp_q[1, ],
                                  occ_me = aco_sp_q[2, ],
                                  occ_hi = aco_sp_q[3, ], 
                                  occ_confi = confi_vals[i]))
  

}

# run empirical results - these are given occ_confi of 100% or 1.
# get ptct in correct format
comm_pc <- subset(comm_pc, wetland_id %in% comm_aco_sites)
comm_pc_sites <- comm_pc$wetland_id
comm_pc <- subset(comm_pc, select = -wetland_id)
comm_pc <- (comm_pc > 0) * 1
# run pc analysis
pc_boots<- replicate(nboots, boot_occ(comm_pc, grp = comm_pc_sites), simplify = FALSE)
pc_sp_q <- apply(sapply(pc_boots, function(x) apply(x, 2, mean)), 1,
                 quantile, c(0.025, 0.5, 0.975))
pc_sp_a <- apply(sapply(pc_boots, function(x) apply(x, 2, mean)), 1, mean)
occ_results <- rbind(occ_results, 
                    data.frame(sp = names(pc_sp_a), occ_avg = pc_sp_a, 
                               occ_lo = pc_sp_q[1, ],
                               occ_me = pc_sp_q[2, ],
                               occ_hi = pc_sp_q[3, ], 
                               occ_confi = 1))

write.csv(occ_results, file = './results/occ_results.csv', row.names = FALSE)

# need to restructure data to a wide format for plotting against 1:1
occ_wide_me <- pivot_wider(occ_results[ , c('sp','occ_me','occ_confi')], values_from = occ_avg,
                        names_from = occ_confi, values_fill = 0)

occ_wide_lo <- pivot_wider(occ_results[ , c('sp','occ_lo','occ_confi')], values_from = occ_avg,
                        names_from = occ_confi, values_fill = 0)
occ_wide_hi <- pivot_wider(occ_results[ , c('sp','occ_hi','occ_confi')], values_from = occ_avg,
                           names_from = occ_confi, values_fill = 0)
head(occ_wide)


cols <- rev(terrain.colors(10))

plot(`1` ~ `0.1`, data = occ_wide, type = 'n')
abline(a=0, b=1)
for(i in seq_along(confi_vals)) { 
  points(as.matrix(occ_wide_me[ , 11]), 
         as.matrix(occ_wide_me[ , i + 1]), col = cols[i+1], pch = 19)
}

# average alpha richness
plot(c(confi_vals, 1), colSums(occ_wide[ , -1]), type = 'p', axes = FALSE, 
     xlab = 'Confidence Value', ylab = 'Average Richness')
axis(side = 1, labels = c(confi_vals, 'Human'), at = c(confi_vals, 1))
axis(side = 2)
# add error bars
arrows()

# gamma richness
plot(c(confi_vals, 1),  colSums(occ_wide[ , -1] > 0), type = 'h')


