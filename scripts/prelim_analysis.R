#Read in data

sumcomm <- read.csv('./data/clean_data/sumcomm.csv') # comm matrix when all detections are totaled
comm <- read.csv('./data/clean_data/comm.csv')       # comm matrix when detections are converted into proportions
dat_sub <- read.csv('./data/clean_data/subset_data.csv') #the cleaned acoustic data with confidence levels and times

comm$CWWI
comm$BGGN

#Chuck-will's widow - NIGHTIME
#write.csv(row.names(comm[comm$CWWI > 1,]), file = './data/CHNA_detections.csv')

# distribution of species richness
hist(rowSums(comm > 0), xlab= "# of Species detected across all recordings at site",
     main = "Distribution of Species Richness")
sr <- rowSums(comm > 0)
plot(density(sr))
sr_site_id <- tapply(sr, wetland_id, mean)
srlog <- log(sr)

boxplot(srlog ~ habitat, ylab = "log(Species Richness)", xlab="Habitat", 
        main="Species Richness by Habitat")
boxplot(sr~wetland_id, ylab = "Species Richness", xlab="Site ID", 
        main="Species Richness per Site")
boxplot(sr~property, ylab = "Species Richness", xlab="Wetland Sites by Property",
        main="Wetland SR by Property")

hist(sr, xlab="Species Richness", main = "Histogram of Species Richness")


sr <- rowSums(comm > 0)
plot(density(sr))
sr_site_id <- tapply(sr, wetland_id, mean)
srlog <- log(sr)

boxplot(srlog ~ habitat, ylab = "log(Species Richness)", xlab="Habitat", 
        main="Species Richness by Habitat")
boxplot(sr~wetland_id, ylab = "Species Richness", xlab="Site ID", 
        main="Species Richness per Site")
boxplot(sr~property, ylab = "Species Richness", xlab="Wetland Sites by Property",
        main="Wetland SR by Property")

hist(sr, xlab="Species Richness", main = "Histogram of Species Richness")










qqnorm(sr)
qqline(sr,col=2)

#Try to transform the data with sqrt
sr_sqrt <- sqrt(sr)
qqnorm(sr_sqrt)
qqline(sr_sqrt,col=2)
hist(sr_sqrt)
plotNormalHistogram(sr_sqrt, xlab="sr_sqrt", main = "Histogram of sr_sqrt")

#Transform data w cube root
sr_cube<- sign(sr) * abs(sr)^(1/3)
qqnorm(sr_cube)
qqline(sr_cube,col=2)
hist(sr_cube)
plotNormalHistogram(sr_cube,xlab="sr_cuve", main = "Histogram of sr_cube")

#Transform data w log
srlog<- log(sr)
qqnorm(srlog)
qqline(srlog,col=2)
plotNormalHistogram(srlog, xlab="srlog", main = "Histogram of srlog")

ks.test(sr,y='pnorm',alternative = "two.sided")
#Asymptotic one-sample Kolmogorov-Smirnov test
#data:  sr
#D = 0.84134, p-value < 2.2e-16

kruskal.test(hab,Common.name,data=dat)
#Kruskal-Wallis rank sum test
#data:  hab and Common.name
#Kruskal-Wallis chi-squared = 743.49, df = 116, p-value < 2.2e-16

diversity(sr, dat, equalize.groups = TRUE)

# random sampling with function sample
random1<-comm[sample(1:nrow(comm), size = 300 ), ]
rsr<-rowSums(random1>0)
rhab<- sample(1:nrow(comm), size = 300)
boxplot(rsr~rhab)

bird_rda <- rda(comm ~ habitat)
anova(bird_rda)
PLOT.RDA(bird_rda, map = "asymmetric")

t.test(srlog~habitat)

summary(aov(srlog~habitat))
#F-value is 65.686 and p-value is 0.001, so we can reject the null hypothesis.
#species richness and shannon index. Throw away some data so that sample sizes match. 
#could try resampling approach. Resample randomly to create same size data.

#filter soundwave files so only sounds above certain decimal are kept? 
#install.packages("warbleR")
#install.packages("Rraven")
#does it match species composition? Bring in jackson's data and compare the community
#matrix. Level of diversity: species ID & compositional diff btwn the two. 
#Merge into single comm matrix. With sample column acoustic vs human.
#convert data into 5 min intervals. Add veg data and area/perimeters.
#Canopy cover variables showed signal in jackson's results. Distance from nearest
#tree and wet area.

#sr avg to occurence stat