library(dplyr)

# read in merlin data
merlin <- read.csv('./data/raw_data - merlin_audio_year2.csv')
dim(merlin)
# replace NA values with zeros
merlin[ , -(1:2)] <- ifelse(is.na(as.matrix(merlin[ , -(1:2)])),
                            0, as.matrix(merlin[ , -(1:2)]))

merlin$NOBO

# read in ptct data
dat_pc <- read.csv('https://raw.githubusercontent.com/mcglinnlab/wetland_birds/main/data/filtered_data/clean_bird_dat.csv')
#note: in dat_pc from Jackson that wetland_id is equivalent to what I've called
# site_id; therefore let's rename wetland_id to site_id for simplicity
#jackson's data is subset from 0y-25 m rather than previous total which was 50m.
head(dat_pc)
comm_pc <- dat_pc[ , c('wetland_id', 'date.x', names(dat_pc)[12:66])]
comm_pc[1:5, 1:5]
# fix name typo
comm_pc$BGGN <- comm_pc$BGGN + comm_pc$BGNN
comm_pc = subset(comm_pc, select = -c(BGNN))
#Bachman's sparrow code is incorrect in ptct data. Switch BASP to BACS
names(comm_pc)[22] <- "BACS"
comm_pc$BACS
comm_pc$BCSP
comm_pc$BACS <- comm_pc$BACS + comm_pc$BCSP
#drop comm_pc$BCSP
comm_pc <- subset(comm_pc, select = -c(BCSP))

# get everything matched up properly for comparison
names(comm_pc)[1] <- 'site'
names(comm_pc)[2] <- 'date'

#date field fix formatting
comm_pc$date <- as.Date(comm_pc$date, format = "%m/%d/%Y")
bad_dates <- which(sapply(comm_pc$date, nchar) == 7)
bad_yrs <- comm_pc$date[bad_dates]
comm_pc$date[bad_dates] <- as.Date(paste("2022-", 
                                         substr(comm_pc$date[bad_dates], 4, 8),
                                         sep = ''))
merlin$date <- as.Date(merlin$date, format = "%m/%d/%Y")

subset(comm_pc, site == "SP03")

comm_pc$site
merlin$site

names(merlin)
# pivot datasets to long for merging

merlin_lg <- pivot_longer(merlin, cols = -(1:2))
pc_lg <- pivot_longer(comm_pc, cols = -(1:2))

merlin_lg[ , 1:2] <- merlin_lg[ , 2:1]
names(merlin_lg)[1:2] <- c('site', 'date')
merlin_lg[1:5, ]

length(unique(pc_lg$site))
length(unique(merlin_lg$site))

names(merlin_lg)[4] <- 'pa'
names(pc_lg)[4] <- 'pc'

head(merlin_lg)
head(pc_lg)


# remove all unneeded zeros
merlin_lg <- subset(merlin_lg, pa > 0)
pc_lg <- subset(pc_lg, pc > 0)

# create merged site date field for filtering
merlin_lg$site_date <- with(merlin_lg, paste(site, date, sep='-'))
pc_lg$site_date <- with(pc_lg, paste(site, date, sep='-'))

merlin_site_dates <- unique(merlin_lg$site_date)

# ready for merging
dat <- merge(merlin_lg, subset(pc_lg, site_date %in% merlin_site_dates),
             all = TRUE)

head(dat)
tail(dat)

dat[1:100, ]

sum(is.na(dat$pa))
sum(is.na(dat$pc))

dat$pa <- ifelse(is.na(dat$pa), 0, dat$pa)
dat$pc <- ifelse(is.na(dat$pc), 0, dat$pc)
dat$pc <- (dat$pc > 0)*1

# now analysis can be carried out
plot(jitter(pa) ~ jitter(pc), data = dat)

# so most of the species occurrences are only detected via PAM
# next biggest group is species ocurrences detected by both
# only a few species occurrences were only detected in pc but not PAM

comm_pa <- with(dat, tapply(pa, list(site_date, name), sum))
comm_pc <- with(dat, tapply(pc, list(site_date, name), sum))

comm_pa <- ifelse(is.na(comm_pa), 0, comm_pa)
comm_pc <- ifelse(is.na(comm_pc), 0, comm_pc)

# species occupancy
occ_pa <- colSums(comm_pa) / nrow(comm_pa)
occ_pc <- colSums(comm_pc) / nrow(comm_pc)

pdf('./figs/occupacny_pc_vs_merlin.pdf')
plot(occ_pc, occ_pa, type = 'n')
abline(a = 0, b=1)
text(occ_pc, occ_pa, labels = names(occ_pc), cex = 0.75)
dev.off()

# species richness
S_pa <- rowSums(comm_pa)
S_pc <- rowSums(comm_pc)
S <- c(S_pa, S_pc)

# looks like no relationship
plot(jitter(S_pc), S_pa)

method <- factor(c(rep('pa', length(S_pa)), rep('pc', length(S_pc))))
mod <- glm(S ~ method, family = 'poisson')
summary(mod)
S_pred <- predict(mod, newdata = data.frame(method = c('pc', 'pa')), 
                  se = TRUE, type = 'response')
exp(confint(mod, type = 'response'))

pdf('./figs/species_richness_pc_vs_merlin.pdf')
bp <- barplot(S_pred$fit, names.arg = c('Point Count', 'Passive Acoustics'),
        ylim = c(0, 15), ylab = 'Species Richness')
arrows(bp[1], S_pred$fit[1] + 1.96 * S_pred$se.fit[1], 
       bp[1], S_pred$fit[1] - 1.96 * S_pred$se.fit[1], 
       code = 3, angle = 90)
arrows(bp[2], S_pred$fit[2] + 1.96 * S_pred$se.fit[2], 
       bp[2], S_pred$fit[2] - 1.96 * S_pred$se.fit[2], 
       code = 3, angle = 90)
dev.off()





