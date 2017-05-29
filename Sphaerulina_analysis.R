###### Pathogen without borders: Sphaerulina populicola is genetically diverse 
# and spreads freely through the host tree’s range

### Fiona Paul, Imke Schmitt, Sunil Mundra, Miklós Bálint

#####################################################
# Contents
# 1. Import dataset and data cleaning
# 2. Calculate and investigate diversity
# 3. Investigate abundance of oligotypes
# 4. Visualization of sample differences with NMDS
# 5. Making a map of the balsam poplar distribution
# 6. Making the haplotype network
################################################

library(plyr) # (Wickham 2011)
library(corrplot) # (Wei and Simko 2016)
library(vegan) # (Oksanen et al. 2016)
library(MASS) # (Venables and Ripley 2002)
library(effects) # (Fox 2003)
library(lme4) # (Bates et al. 2015)
library(mvabund) # (Wang et al. 2016)
library(maps) # (Becker et al. 2016)
library(maptools) # (Bivand and Lewin-Koh 2016)
library(rgdal) # (Bivand et al. 2017)
library(scales) # (Wickham 2016)
library(raster) # (Hijmans 2015)
library(mapplots) # (Gerritsen 2014)
library(ape) # (Paradis et al. 2004)
library(pegas) # (Paradis 2010)
library(dplyr) # (Wickham and Francois 2016)
library(ade4) # (Dray and Dufour 2007)

save.image("Sphaerulina_oligotypes.image.RData")
load("Sphaerulina_oligotypes.image.RData")

##### 1. Import dataset and data cleaning #####
## Read in abundance table
oligotypeabund = read.csv(file="MATRIX-COUNT.txt", header=T, sep='', row.names = 1)
oligotypeabund = oligotypeabund[order(row.names(oligotypeabund)),]

## Sequencing depth and oligotype abundances
# Frequency distribution of number of reads per sample
readsample = apply(oligotypeabund,1,sum)
summary(readsample)
hist(readsample, main = "Histogram of summed reads per sample", 
     xlab = "Summed reads per sample")
# Frequency distribution of number of reads per oligotype
readOT = apply(oligotypeabund,2,sum)
hist(readOT, main = "Histogram of summed reads per oligotype", 
     xlab = "Summed reads per oligotype")
# Freuqeuncy distribution of maximum read number per sample
hist(apply(oligotypeabund, 1, max), main = "Histogram of maximum read number per sample", 
     xlab = "Max read number per sample", ylab = "Frequency")
# Frequency distribution of highest read count per oligotype
hist(apply(oligotypeabund, 2, max), main = "Histogram of maximum read number per oligotype",
     xlab = "Max read number per oligotype", ylab = "Frequency")

## Oligotype abundances
colors = rainbow(43)
pie(readOT, labels = names(readOT), col = colors)
pie(sqrt(readOT), labels = names(readOT), col = colors)

### Clean up negative controls
# Remove the maximum read number of a sequence variant found in a negative control from 
# every sample that contains that variant

# Blank controls
Blank = grep("^B.*", row.names(oligotypeabund))

# Negative controls
Negative = grep("^N.*", row.names(oligotypeabund))

# Maximum number of reads in any control sample
MaxControl = apply(oligotypeabund[c(Blank, Negative),], 2, max)

# Extract the highest read number of a sequence variant in a control from every sample
negcleanoligotypeabund = oligotypeabund

negcleanoligotypeabund[grep("^[a-z]", row.names(oligotypeabund)),] <- 
  sweep(oligotypeabund[grep("^[a-z]", row.names(oligotypeabund)),], 2, MaxControl, "-")

# Set negative values to 0. Warnings are because the non-numeric cells
negcleanoligotypeabund[negcleanoligotypeabund < 0] <- 0

# Remove the negative samples
negcleanoligotypeabund = negcleanoligotypeabund[grep("^[a-z]", row.names(oligotypeabund)),]
# removal of 21 negative control samples

### Filter rare observations
## Per oligotype threshold
# Set counts within an Oligotype falling below 0.1% of the highest count of an Oligotype
# to 0 to remove errors

# Define threshold
rarethreshold = apply(negcleanoligotypeabund, 2, max)*0.001

rm.small.obs = colwise(function(x){return(ifelse(x <= rarethreshold, 0, x))})

rarecleanedabund = rm.small.obs(as.data.frame(t(negcleanoligotypeabund)))
rarecleanedabund = as.data.frame(t(rarecleanedabund))
colnames(rarecleanedabund) = colnames(negcleanoligotypeabund)

## Global threshold
# Filter out observations that fall below 0.05% of globally highest observation

# Define global cutoff value
maxcount = max(rarecleanedabund)
cutoff.global = round(maxcount*0.0005) 

# Set the counts falling below cutoff to 0
globallycleanedabund = rarecleanedabund
globallycleanedabund[globallycleanedabund < cutoff.global] = 0

## Filter out oligotypes if they occur in less than 5 samples/trees
low.presence = apply(globallycleanedabund,2,function(vec) sum(vec>0))
IsFreq = low.presence > 5
prescleanedabund = globallycleanedabund[,IsFreq]

## Remove samples with no reads left
finalcleanedabund = prescleanedabund
samples.to.keep = (apply(finalcleanedabund,1,sum)) > 0
finalcleanedabund = finalcleanedabund[samples.to.keep,]
# removal of 393 samples

## Remove oligotypes with no reads left
ots.to.keep = (apply(finalcleanedabund,2,max)) > 0
finalcleanedabund = finalcleanedabund[,ots.to.keep]
# Removal of 14 oligotypes

### Remove problematic samples from the dataset
finalcleanedabund = finalcleanedabund[4:69,]
finalcleanedabund = finalcleanedabund[-c(1,25,48,52,54),]

### Sequencing depth and oligotype abundance in the cleaned dataset
# Frequency distribution of number of reads per sample
# Abundance distribution plots
cleanreadsample = apply(finalcleanedabund,1,sum)
hist(cleanreadsample,breaks = 100 ,main = "Histogram of summed reads per sample", 
     xlab = "Summed reads per sample")
# Frequency distribution of number of reads per oligotype
cleanreadOT = apply(finalcleanedabund,2,sum)
hist(cleanreadOT,breaks = 100, main = "Histogram of summed reads per oligotype", 
     xlab = "Summed reads per oligotype")
# Freuqeuncy distribution of maximum read number per sample
hist(apply(finalcleanedabund, 1, max), main = "Histogram of maximum read number per sample", 
     xlab = "Max read number per sample", ylab = "Frequency")
# Frequency distribution of highest read count per oligotype
hist(apply(finalcleanedabund, 2, max), main = "Histogram of maximum read number per oligotype",
     xlab = "Max read number per oligotype", ylab = "Frequency")
png(file="histogram.png", units="mm", height=90, width=90, 
    pointsize=10, bg="white", res=1200)

## Abundance of oligotypes
###### Oligotype pie charts ######
# colors for final oligotype plotting
my_color = c("#d896ff", "#800080", "#ee4035", "#fdf498", "#7bc043", "#0392cf", "#028900",
             "#49796b", "#602320", "#011f4b", "#000000", "#a0d6b4", "#ffbf00", "#a67c00", 
             "#ff0097", "#ff0000")
par(mar=c(0,0,0,0))
pie(cleanreadOT, labels = c("155", "159", "125", "92", "6", "15", "102", "79", "80", "156", "78", 
                            "117", "58", "105", "8", "161"), col = my_color)
pdf("piechart_oligo_freq.pdf",height = 8.27, width= 11.69, pointsize = 18 )
par(mar=c(0,0,0,0))
pie(sqrt(cleanreadOT), labels = c("155", "159", "125", "92", "6", "15", "102", "79", "80", "156", 
                                  "78", "117", "58", "105", "8", "161"), col = my_color)
dev.off()

## Print out abundance table into a csv file
write.table(finalcleanedabund, file = "cleaned_abundance_table.csv",sep = ";", 
            col.names = NA, row.names = TRUE)

##### 2. Calculate and investigate diversity (creating the diversity boxplot) #####
# load in the abundance table with reads of oligotypes per tree
abundance = read.csv(file="cleaned_abundance_table.csv", header=T, sep=';', row.names = 1)

# load in the metadata by sample with raw read numbers
metadata = read.csv(file = "metadata.csv", header = T, sep = ";", row.names = 1)

## Correlation of explanatory variables
###### Correlation plot of environmental variables ######
pdf("correlation_matrix_environment.pdf",height = 8.27, width= 11.69, pointsize = 18 )
par(mar=c(2,0,0,0), oma=c(2,0,0,0))
cor.data = cbind(metadata[,c((2:10),15,16)])
corrplot.mixed(cor(cor.data), upper = "number", lower = "circle", tl.pos = "lt", order = "hclust", 
               tl.cex= 0.75, tl.col = "black", tl.srt = 45)
dev.off()


### Calculate Hill diversities
OTHill = renyi(abundance, scale=c(0,1,2), hill=T)
# Hill 1
hill.1 = OTHill$"0"
names(hill.1) = rownames(abundance)
hist(hill.1)
plot(hill.1 ~ metadata$reads, pch = 19, ylab= "richness", xlab="number of raw reads")
shapiro.test(hill.1)
# Hill 2
hill.2 = OTHill$"1"
names(hill.2) = rownames(abundance)
hist(hill.2)
plot(hill.2 ~ metadata$reads, pch = 19)
# Hill 3
hill.3 = OTHill$"2"
names(hill.3) = rownames(abundance)
hist(hill.3)
plot(hill.3 ~ metadata$reads, pch = 19)

### GLMs for factors influencing diversity
## Hill 1
# Read number
hill1.glm.reads = glm(hill.1 ~ reads, data = metadata)
summary(hill1.glm.reads)
hill1.reads = glm.nb(hill.1 ~ reads, data = metadata)
summary(hill1.reads)
anova(hill1.glm.reads,hill1.reads)
AIC(hill1.glm.reads,hill1.reads)
# negative binomial is better than normal
plot(allEffects(hill1.reads))

# Site/Location effect
hill1.site = glm.nb(hill.1 ~ reads + site, data = metadata)
summary(hill1.site)
anova(hill1.site)
anova(hill1.reads,hill1.site)
AIC(hill1.reads,hill1.site)

# Read number and site effect plot
plot(allEffects(hill1.site))

# Replication Effect
hill1.rep = glm.nb(hill.1 ~ reads + replicates, data = metadata)
summary(hill1.rep)
anova(hill1.rep)
anova(hill1.reads,hill1.site,hill1.rep)
AIC(hill1.reads,hill1.site,hill1.rep)
plot(allEffects(hill1.rep))

# Region effect, Canada vs Alaska
hill1.reg = glm.nb(hill.1 ~ reads + Region, data = metadata)
summary(hill1.reg)
plot(allEffects(hill1.reg))

# Latitude effect
hill1.lat = glm.nb(hill.1 ~ reads + latitude, data = metadata)
summary(hill1.lat)
anova(hill1.reads,hill1.site,hill1.lat)
AIC(hill1.reads,hill1.site,hill1.lat)
plot(allEffects(hill1.lat))
hill1.lat2 = glm.nb(hill.1 ~ reads + latitude + site, data = metadata)
summary(hill1.lat2)
anova(hill1.lat2)

# Latitude, read number, site effect plot
plot(allEffects(hill1.lat2))
plot(effect("reads", hill1.lat2), main= NULL, xlab= "Sequencing depth (DNA sequences)", 
     ylab= "Hill's N0")
plot(effect("latitude", hill1.lat2), main= NULL, xlab= "Latitude (° N)", ylab= "Hill's N0")

# Landuse type effect
hill1.land = glm.nb(hill.1 ~ reads + Landuse_type, data = metadata)
summary(hill1.land)
anova(hill1.reads,hill1.site,hill1.land)
AIC(hill1.reads,hill1.site,hill1.land)
plot(allEffects(hill1.land))


## Hill 2
# Read number
hill2.glm.nb.reads = glm.nb(hill.2 ~ reads, data = metadata)
summary(hill2.glm.nb.reads)
hill2.glm = glm(hill.2 ~ reads, data = metadata)
summary(hill2.glm)
anova(hill2.glm.nb.reads,hill2.glm)
AIC(hill2.glm, hill2.glm.nb.reads)
plot(allEffects(hill2.glm))

# Site/Location effect
hill2.site = glm(hill.2 ~ reads + site, data = metadata)
summary(hill2.site)
anova(hill2.glm,hill2.site)
AIC(hill2.glm,hill2.site)
plot(allEffects(hill2.site))

# Replication Effect
hill2.rep = glm(hill.2 ~ reads + replicates, data = metadata)
summary(hill2.rep)
anova(hill2.rep)
anova(hill2.glm,hill2.site, hill2.rep)
AIC(hill2.glm,hill2.site, hill2.rep)
plot(allEffects(hill2.rep))

# Region effect, Canada vs Alaska
hill2.reg = glm(hill.2 ~ reads + Region, data = metadata)
summary(hill2.reg)
plot(allEffects(hill2.reg))

# Latitude effect
hill2.lat = glm(hill.2 ~ reads + latitude, data = metadata)
summary(hill2.lat)
plot(allEffects(hill2.lat))

# Landuse type effect
hill2.land = glm(hill.2 ~ reads + Landuse_type, data = metadata)
summary(hill2.land)
plot(allEffects(hill2.land))


## Hill 3
# Read number
hill3.glm.nb.reads = glm.nb(hill.3 ~ reads, data = metadata)
summary(hill3.glm.nb.reads)
hill3.glm = glm(hill.3 ~ reads, data = metadata)
summary(hill3.glm)
anova(hill3.glm.nb.reads,hill3.glm)
AIC(hill3.glm, hill3.glm.nb.reads)
plot(allEffects(hill3.glm))

# Site/Location effect
hill3.site = glm(hill.3 ~ reads + site, data = metadata)
summary(hill3.site)
anova(hill3.glm, hill3.site)
AIC(hill3.glm, hill3.site)
plot(allEffects(hill3.site))

# Replication Effect
hill3.rep = glm(hill.3 ~ reads + replicates, data = metadata)
summary(hill3.rep)
anova(hill3.rep)
anova(hill3.glm,hill3.site, hill3.rep)
AIC(hill3.glm,hill3.site, hill3.rep)
plot(allEffects(hill3.rep))

# Region effect, Canada vs Alaska
hill3.reg = glm(hill.3 ~ reads + Region, data = metadata)
summary(hill3.reg)
plot(allEffects(hill3.reg))

# Latitude effect
hill3.lat = glm(hill.3 ~ reads + latitude, data = metadata)
summary(hill3.lat)
plot(allEffects(hill3.lat))

# Landuse type effect
hill3.land = glm(hill.3 ~ reads + Landuse_type, data = metadata)
summary(hill3.land)
plot(allEffects(hill3.land))


###### Making the boxplots for richness and diversity in the two geographic demes ######
shannon = diversity(abundance, index = "shannon", MARGIN = 1)
pdf(file = "richness_diversity_region_boxplot.pdf", height = 6.5, width= 11.69, pointsize = 18)
par(mar=c(2,4,2,2), las=1, oma=c(2,1,1,1), mfrow=c(1,2))
boxplot(hill.1 ~ metadata$Region, ylab="Oligotype richness", fill=TRUE, col="gray")
boxplot(shannon ~ metadata$Region, ylab="Shannon diversity", fill=TRUE, col="gray")
dev.off()
###

### Mixed effects GLMs for diversity patterns
## Hill 1
# Read number, null model
hill.glmer.reads = glmer.nb(hill.1 ~ (scale(reads)) + (1|site), data = metadata)
summary(hill.glmer.reads)
anova(hill.glmer.reads)
AIC(hill.glmer.reads)
plot(allEffects(hill.glmer.reads))

# Replication number effect
hill.glmer.rep = glmer.nb(hill.1 ~ (scale(reads)) + replicates + (1|site), 
                          data= metadata)
summary(hill.glmer.rep)
anova(hill.glmer.reads, hill.glmer.rep)
plot(allEffects(hill.glmer.rep))

# Region effect
hill.glmer.reg = glmer.nb(hill.1 ~ (scale(reads)) + Region + (1|site), 
                          data = metadata)
summary(hill.glmer.reg)
anova(hill.glmer.reads,hill.glmer.reg)
plot(allEffects(hill.glmer.reg))

# Latitude effect
hill.glmer.lat = glmer.nb(hill.1 ~ (scale(reads)) + latitude + (1|site),
                          data = metadata)
summary(hill.glmer.lat)
anova(hill.glmer.reads,hill.glmer.lat)
plot(allEffects(hill.glmer.lat))

# Landuse type effect
hill.glmer.land = glmer.nb(hill.1 ~ (scale(reads)) + Landuse_type + (1|site), 
                           data = metadata)
summary(hill.glmer.land)
anova(hill.glmer.reads, hill.glmer.land)
plot(allEffects(hill.glmer.land))

# Temperature effect
hill.glmer.temp = glmer.nb(hill.1 ~ (scale(reads)) + temp_annual + (1|site), 
                           data = metadata)
summary(hill.glmer.temp)
anova(hill.glmer.reads, hill.glmer.temp)
plot(allEffects(hill.glmer.temp))

# Precipitation effect
hill.glmer.ppt = glmer.nb(hill.1 ~ (scale(reads)) + scale(prec_annual) + 
                            (1|site), data = metadata)
summary(hill.glmer.ppt)
anova(hill.glmer.reads, hill.glmer.ppt)
plot(allEffects(hill.glmer.ppt))

anova(hill.glmer.reads,hill.glmer.rep,hill.glmer.reg,hill.glmer.lat,hill.glmer.temp,hill.glmer.ppt)

## Hill 2
# Read number effect 
hill2.glmer.reads = glmer.nb(hill.2 ~ (scale(reads)) + (1|site),data = metadata)
summary(hill2.glmer.reads)
anova(hill2.glmer.reads)
AIC(hill2.glmer.reads)


## Hill 3
# Read number effect
hill3.glmer.reads = glmer.nb(hill.3 ~ (scale(reads)) + (1|site),data = metadata)
summary(hill3.glmer.reads)
anova(hill3.glmer.reads)
AIC(hill3.glmer.reads)


##### 3. Investigate abundance of oligotypes #####
# Checking the response to factors of the abundance of individual oligotypes
# Investigating abundance patterns with manyglms
abund.mva
# convert the abundance table into mvabund table
abund.mva = mvabund(abundance)
# plot of the mean-variance relationship

# Mean Variance plot for dataset
meanvar.plot(abund.mva, xlab= "Mean", ylab= "Variance", table=T)

# Read number effect
glm.readnum = manyglm(abund.mva ~ reads, family = "negative.binomial", 
                      data = metadata)
glm.readnum.anova = anova.manyglm(glm.readnum, test = "LR",nBoot=1000)
glm.readnum.anova
plot(glm.readnum)

# Site effect
glm.site = manyglm(abund.mva ~ site, family = "negative.binomial", 
                   data = metadata)
glm.site.anova = anova.manyglm(glm.site, test = "LR",nBoot=1000)
glm.site.anova

# Replication effect
glm.repl = manyglm(abund.mva ~ replicates, family = "negative.binomial",
                   data = metadata)
glm.repl.anova = anova.manyglm(glm.repl, test = "LR", nBoot = 1000)
glm.repl.anova

# Region effect
glm.region = manyglm(abund.mva ~ Region, family = "negative.binomial", 
                     data = metadata)
glm.region.anova = anova.manyglm(glm.region, test = "LR",nBoot=1000)
glm.region.anova

# Latitude effect
glm.lat = manyglm(abund.mva ~ latitude, family = "negative.binomial", 
                  data = metadata)
glm.lat.anova = anova.manyglm(glm.lat, test = "LR",nBoot=1000)
glm.lat.anova

# Landuse type effect
glm.land = manyglm(abund.mva ~ Landuse_type, data = metadata, 
                   family = "negative.binomial")
glm.land.anova = anova.manyglm(glm.land, test = "LR",nBoot=1000)
glm.land.anova

# Temperature effect
glm.temp = manyglm(abund.mva ~ temp_annual, family = "negative.binomial",
                   data = metadata)
glm.temp.anova = anova.manyglm(glm.temp, test = "LR", nBoot = 1000)
glm.temp.anova

# Precipitation effect
glm.ppt = manyglm(abund.mva ~ prec_annual, family = "negative.binomial",
                  data = metadata)
glm.ppt.anova = anova.manyglm(glm.ppt, test = "LR", nBoot = 1000)
glm.ppt.anova


##### 4. Visualization of sample differences with NMDS #####
MDS = metaMDS(abundance, distance = "bray")
MDS = metaMDS(abundance, previous = MDS)
stressplot(MDS)

##### NMDS plot with sample dots coloured by sampling location #####
stand.data=  read.csv(file="stand.csv", sep=";", header=T, row.names = 1)
attach(stand.data)
head(stand.data)
mode(stand)
stand<-as.factor(stand)
mode(stand) 
is.factor(stand)

par(mfrow=c(1,1), mar=c(4,4,2,2))
gnmds1 <- jitter(MDS$points[,1],600)
gnmds2 <- jitter(MDS$points[,2],600)

mds.df<-data.frame(gnmds1,gnmds2)

pdf("NMDS_locations.pdf",height = 8.27, width= 11.69, pointsize = 18 )
par(mar=c(5,5,2,2))
plot(MDS$points, type="n", xlab="NMDS1", ylab="NMDS2", xlim=c(-4,4), ylim= c(-6,3))
points(gnmds1[stand==1],gnmds2[stand==1],cex=1,pch=16,col="firebrick1")
points(gnmds1[stand==5],gnmds2[stand==5],cex=1,pch=16,col="brown")
points(gnmds1[stand==4],gnmds2[stand==4],cex=1,pch=16,col="chocolate")
points(gnmds1[stand==8],gnmds2[stand==8],cex=1,pch=16,col="yellow")
points(gnmds1[stand==6],gnmds2[stand==6],cex=1,pch=16,col="chartreuse")
points(gnmds1[stand==7],gnmds2[stand==7],cex=1,pch=16,col="darkgreen")
points(gnmds1[stand==2],gnmds2[stand==2],cex=1,pch=16,col="violet")
points(gnmds1[stand==12],gnmds2[stand==12],cex=1,pch=16,col="purple4")
points(gnmds1[stand==9],gnmds2[stand==9],cex=1,pch=16,col="blue")
points(gnmds1[stand==10],gnmds2[stand==10],cex=1,pch=16,col="cadetblue")
points(gnmds1[stand==11],gnmds2[stand==11],cex=1,pch=16,col="skyblue4")
points(gnmds1[stand==3],gnmds2[stand==3],cex=1,pch=16,col="darkgrey")
legend(-4, -1, c("Arctic Village","Fairbanks", "Denali N. Park","Hay River", "Fort McMurray",
                 "Grande Prairie","Boyle", "Cadxyz","Love","Melville","Portage", "Carnduff"), 
       fill=c("firebrick1","brown","chocolate","yellow","chartreuse","darkgreen",
              "violet","purple4","blue","cadetblue", "skyblue4","darkgrey"), cex = 0.7)
dev.off()

##### NMDS plot by region, i.e. Canada and Alaska #####
region.data=  read.csv(file="region.csv", sep=";", header=T, row.names = 1)
attach(region.data)
head(stand.data)
mode(reg)
reg<-as.factor(reg)
mode(reg) 
is.factor(reg)

newgnmds1 <- jitter(MDS$points[,1],300)
newgnmds2 <- jitter(MDS$points[,2],300)

mds.df<-data.frame(newgnmds1,newgnmds2)

pdf("NMDS_regions.pdf",height = 8.27, width= 11.69, pointsize = 18 )
par(mar=c(5,5,2,2))
plot(MDS$points, type="n", xlab="NMDS1", ylab="NMDS2")
points(newgnmds1[reg==2],newgnmds2[reg==2],cex=1,pch=16,col="blue")
points(newgnmds1[reg==1],newgnmds2[reg==1],cex=1,pch=16,col="red")
legend(3.5, 2, c("Alaska", "Canada"), fill = c("red","blue"))
dev.off()

##### 5. Making a map of the balsam poplar distribution #####
# Load in the balsam poplar distribution shape file
popdistr = readOGR("popubals.shp")
head(popdistr)

# Load in the sampling location coordinates
locations = read.csv("metadata.csv", header = T, sep = ";", row.names = 1)
head(locations)

# Load in the mountain range raster layer shape file
rockys = raster("alt.grd")
rockys

### Piechart preparations
# Abundances summed for locality
local_abundance_input = cbind(abundance, site = metadata$site_code)
local_abundance = aggregate(. ~ site, data = local_abundance_input, FUN = sum)

# Site coordinates
pop_coord_input = data.frame(lon = locations$longitude,
                             lat = locations$latitude,
                             site = metadata$site_code)
pop_coord = aggregate(. ~ site, data = pop_coord_input, FUN = mean)

# only the relatively large read numbers
frequent_abund = local_abundance[,2:17]

# create a pdf file map of oligotype frequency distribution across study area
pdf("poplar_distribution_pies.pdf", height = 8.27, width= 10.8, pointsize = 18)
par(mfrow=c(1,1),mar=c(4,4,2,0), oma=c(0,0,0,0))
plot(rockys, maxpixels= 750000, xlim=c(-155,-95), ylim=c(45,70), useRaster=TRUE,
     interpolate=TRUE, box=FALSE, axes=FALSE, col=gray.colors(9,start=1,end= 0.1, gamma=0.2),
     legend= F, xlab= "Longitude (° E)", ylab= "Latitude (° N)", cex.lab=1)
plot(popdistr, add=T, col=alpha("darkgreen", 0.5) , border=FALSE, xlim=c(-155,-95), ylim=c(45,70))
map(database = "world", add= TRUE, xlim = c(-155,-95), ylim = c(45,70), myborder = c(0,0))
axis(side = 1, lwd = 1, lwd.ticks = 1, tck=-0.02, cex.axis=1)
axis(side = 2, lwd = 1, lwd.ticks = 1, tck=-0.02, cex.axis=1)
for (i in 1:nrow(pop_coord)){
  add.pie(as.numeric(sqrt(frequent_abund[i,])), pop_coord[i,2], pop_coord[i,3], 
          radius=log(sum(frequent_abund[i,]))/5,
          col = my_color, labels = "")
}
dev.off()


##### 6. Making the haplotype network #####

input <- "repr_seqs_sphaerulina_aligned.fasta"
d <- ape::read.dna(input, format='fasta')
e <- dist.dna(d)
h <- pegas::haplotype(d)
h <- sort(h, what = "label")
(net <- pegas::haploNet(h))

haplo_size = log(apply(frequent_abund,2,sum))/5

pdf("haplotype_network.pdf",height = 8.27, width= 11.69, pointsize = 18 )
plot(net, scale.ratio=1, bg = my_color,
     labels = F,
     threshold = 0,
     size = haplo_size)
dev.off()
