
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
################# start script  ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

rm(list=ls())
objects()

################################# set working directory
#setwd("C:\\Eeva_jobb\\artikkeleita tekeillä\\James herbivore diversity\\analysis")
setwd('S:\\DISENTANGLE\\WP3\\ArcticHerbivoreFDPD\\Final trait classification')

#################################  get packages 

require(FactoMineR) # factorial analysis of mixed data
require(ade4) # mantel test
require(ape) # export data to biodiverse
require(RColorBrewer)#Colours to match the results figures

#################################  get data

#traits<-read.delim("herbivore_traits_sept17.txt", dec=",")
#traits<-read.delim("S:\\DISENTANGLE\\WP3\\ArcticHerbivoreFDPD\\Final trait classification\\herbivore_traits_apr18.txt", dec=".")
traits<-read.delim("herbivore_traits_apr18.txt", dec=".")
names(traits)


## select variables to dataframe that will be used
var<-c("Order", "Family", "Genus", "Species", "Binomial", 
       "body_mass", "gut_type", "group_size_summer", "group_size_winter",  "wintering_strategy",
       "Litter_clutch_size",  "Population_dynamics", "Habitat_type", "Belowground_feeding",
       "Mobility", "Human_managed",  "Diet_type",   "diet_item_forb", "diet_item_graminoid", 
       "diet_item_shrub",   "diet_item_moss"  , "diet_item_lichen")

Traits<-traits[,which(colnames(traits) %in% var)]
summary(Traits)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
################# remove species that were not part of the phylogeny analyses   ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

names(Traits)

Traits<-Traits[!Traits$Binomial=="Dicrostonyx nelsoni",]
Traits<-Traits[!Traits$Binomial=="Dicrostonyx unalascensis",]
Traits<-Traits[!Traits$Binomial=="Dicrostonyx vinogradovi",]
Traits<-Traits[!Traits$Binomial=="Dicrostonyx nunatakensis",]
Traits<-Traits[!Traits$Binomial=="Lemmus portenkoi",]
Traits<-Traits[!Traits$Binomial=="Ovis aries",]
Traits$Binomial<-droplevels(Traits$Binomial)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
################# explore data  ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


## rename and re-group some variable levels
levels(Traits$Diet_type)[which(levels(Traits$Diet_type)=="facultative_generalist")] = "fac_gen" 
levels(Traits$Diet_type)[which(levels(Traits$Diet_type)=="facultative_specialist")] = "fac_spe" 
levels(Traits$Diet_type)[which(levels(Traits$Diet_type)=="obligatory_generalist")] = "ob_gen" 
levels(Traits$Diet_type)[which(levels(Traits$Diet_type)=="obligatory_specialist")] = "ob_spe" 
levels(Traits$Diet_type)

levels(Traits$group_size_summer)[which(levels(Traits$group_size_summer)=="small_groups")] = "small_group" 
levels(Traits$group_size_summer)[which(levels(Traits$group_size_summer)=="family_group")] = "small_group" 
levels(Traits$group_size_winter)[which(levels(Traits$group_size_winter)=="small_groups")] = "small_group" 
levels(Traits$group_size_winter)[which(levels(Traits$group_size_winter)=="family_group")] = "small_group" 

## histograms of variables
options(device = "windows")
windows(width=30, height=30)

par(mfrow=c(4,4), mar=c(8,4,2,2))
hist(log(Traits$body_mass), main="Body mass (grams log)", xlab = "", ylab="")
plot(Traits$gut_type, main="Gut type", las=2)
plot(Traits$group_size_summer, main="Group size S", las=2)
plot(Traits$group_size_winter, main="Group size W", las=2)
plot(Traits$Population_dynamics, main="Popul dynamics", las=2)
plot(Traits$Habitat_type, main="Habitat type", las=2)
plot(Traits$wintering_strategy, main="Wintering strategy", las=2)
#plot(Traits$Hibernation, main="Hibernation", las=2)
hist(Traits$Litter_clutch_size, main="Litter size", xlab = "", ylab="")
plot(Traits$Belowground_feeding, main="Belowground feeding", las=2)
plot(Traits$Mobility, main="Mobility")
#plot(Traits$Human_managed, main="Human managed")
plot(Traits$Diet_type, main="Diet type", las=2)
plot(as.factor(Traits$diet_item_forb), main="Forbs in diet")
plot(as.factor(Traits$diet_item_graminoid), main="Graminoids in diet")
plot(as.factor(Traits$diet_item_shrub), main="Shrubs in diet")
plot(as.factor(Traits$diet_item_moss), main="Mosses in diet")
plot(as.factor(Traits$diet_item_lichen), main="Lichens in diet")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
################# transform variables   ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

## body mass to log(body mass); otherwise the large values will dominate enormously
Traits$body_mass_log<-log(Traits$body_mass); hist(Traits$body_mass_log)

## group size to numeric with order; gradient from small to large group size 
Traits$group_size_order<-Traits$group_size_summer
str(Traits)
levels(Traits$group_size_order)<- c(levels(Traits$group_size_order), "1", "2", "3")
Traits$group_size_order[Traits$group_size_order == "solitary"] <- "1"
Traits$group_size_order[Traits$group_size_order == "family_group"] <- "2"
Traits$group_size_order[Traits$group_size_order == "small_group"] <- "2"
Traits$group_size_order[Traits$group_size_order == "small_groups"] <- "2"
Traits$group_size_order[Traits$group_size_order == "large_group"] <- "3"
Traits$group_size_order<-droplevels(Traits$group_size_order)
Traits$group_size_summer_order<-as.numeric(Traits$group_size_order)

Traits$group_size_order<-Traits$group_size_winter
levels(Traits$group_size_order) <- c(levels(Traits$group_size_order), "1", "2", "3")
Traits$group_size_order[Traits$group_size_order == "solitary"] <- "1"
Traits$group_size_order[Traits$group_size_order == "family_group"] <- "2"
Traits$group_size_order[Traits$group_size_order == "small_group"] <- "2"
Traits$group_size_order[Traits$group_size_order == "large_group"] <- "3"
Traits$group_size_order<-droplevels(Traits$group_size_order)
Traits$group_size_winter_order<-as.numeric(Traits$group_size_order)


## population dynamics; gradient from little temporal variation ot loads of temporal variation
Traits$population_dynamics_order<-Traits$Population_dynamics
levels(Traits$population_dynamics_order) <- c(levels(Traits$population_dynamics_order), "1", "2", "3")
Traits$population_dynamics_order[Traits$population_dynamics_order == "noncyclic"] <- "1"
Traits$population_dynamics_order[Traits$population_dynamics_order == "cyclic_noncyclic"] <- "2"
Traits$population_dynamics_order[Traits$population_dynamics_order == "cyclic"] <- "3"
Traits$population_dynamics_order<-droplevels(Traits$population_dynamics_order)
Traits$population_dynamics_order<-as.numeric(Traits$population_dynamics_order)

## gut type; gradient from "unefficient" to "efficient"
Traits$gut_type_order<-Traits$gut_type
levels(Traits$gut_type_order) <- c(levels(Traits$gut_type_order), "1", "2", "3")
Traits$gut_type_order[Traits$gut_type_order == "undifferentiated"] <- "1"
Traits$gut_type_order[Traits$gut_type_order == "hindgut_fermenter"] <- "2"
Traits$gut_type_order[Traits$gut_type_order == "ruminant"] <- "3"
Traits$gut_type_order<-droplevels(Traits$gut_type_order)
Traits$gut_type_order<-as.numeric(Traits$gut_type_order)


# ## diet type; gradient from "generalist to specialist " towards increasing specialisation
Traits$diet_type_order<-Traits$Diet_type
levels(Traits$diet_type_order) <- c(levels(Traits$diet_type_order), "1", "2", "3", "4")
Traits$diet_type_order[Traits$diet_type_order == "ob_gen"] <- "1"
Traits$diet_type_order[Traits$diet_type_order == "fac_gen"] <- "2"
Traits$diet_type_order[Traits$diet_type_order == "fac_spe"] <- "3"
Traits$diet_type_order[Traits$diet_type_order == "ob_spe"] <- "4"
Traits$diet_type_order<-droplevels(Traits$diet_type_order)
Traits$diet_type_order<-as.numeric(Traits$diet_type_order)


# ## winter mode; gradient from "certainly not present " towards increasing presence
 Traits$wintering_strategy_order<-Traits$wintering_strategy
 levels(Traits$wintering_strategy_order) 

 #Below ground feeding
 levels(Traits$Belowground_feeding)[levels(Traits$Belowground_feeding)=='none']<-'no belowground feeding'

 #Mobility
 levels(Traits$Mobility)<-c('immobile','mobile')
  
str(Traits)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
################# clustering using factorial analysis of mixed data  ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

names(Traits)

var<-c( "Binomial", 
       "body_mass_log", "gut_type_order", "group_size_summer_order", "group_size_winter_order",  "wintering_strategy_order",  "Litter_clutch_size",  "population_dynamics_order", "Habitat_type", "Belowground_feeding",
       "Mobility",  "diet_type_order",   "diet_item_forb", "diet_item_graminoid", 
       "diet_item_shrub",   "diet_item_moss"  , "diet_item_lichen")

TOTO<-Traits[,which(colnames(Traits) %in% var)]
#TOTO<-na.omit(TOTO)
head(TOTO)
rownames(TOTO)<- TOTO[,1]
names(TOTO)
TOTO<-TOTO[,2:17]
Traits_famd<-TOTO ; summary(Traits_famd)

res.famd<-FAMD(Traits_famd,ncp=6) # factorial analysis of mixed data
summary(res.famd)

options(device = "windows")
windows(width=15, height=10)
par(mfrow=c(1,2))
plot(res.famd,choix="quant",las=1,title='Quantitative variables')
mtext(side=3,'(a)',adj=0,line=2,cex=2)
plot(res.famd,choix="quali",las=1,title='Quaulitative variables')
mtext(side=3,'(b)',adj=0,line=2,cex=2)


res.famd$eig
barplot(res.famd$eig$`percentage of variance`)
## how much does each variable contribute to the 
res.famd
TITO<-as.data.frame(res.famd$var$contrib)

options(device = "windows")
windows(width=15, height=15)

xaxnames<-c("Litter/clutch size", "Forbs in diet", "Graminoids in diet", "Shrubs in diet", "Mosses in diet", "Lichens in diet",
            "Body mass", "Group size summer", "Group size winter", "Population dynamics", "Gut type", "Diet type", "Habitat type", 
            "Belowground feeding", "Mobility", "Wintering strategy")

cols<-brewer.pal(7,'Blues')[2:7]
par(las=2, mar=c(14,6,2,2))
plot(TITO$Dim.1, xaxt="n", xlab="", pch=16, col=cols[6], ylim=c(0,60), ylab="variable contribution to dimensions", cex=1.5, cex.lab=1.5, cex.axis=1.5)
points(TITO$Dim.2, pch=16, col=cols[5], cex=1.5)
points(TITO$Dim.3, pch=16, col=cols[4], cex=1.5)
points(TITO$Dim.4, pch=16, col=cols[3], cex=1.5)
points(TITO$Dim.5, pch=16, col=cols[2], cex=1.5)
points(TITO$Dim.6, pch=16, col=cols[1], cex=1.5)
axis(side=1,at=1:16, labels=xaxnames, cex.axis=1.5)
legend("topleft", legend=c("dim1","dim2","dim3","dim4","dim5",'dim6' ), pch=16, col=rev(cols),
                bty="n", cex=1.5)


hc_2_unconsol<-HCPC(res.famd,method='ward',metric='euclidean',consol=F,nb.clust=-1) # Hierarchical Clustering on Principle Components
hc_2<-HCPC(res.famd,method='ward',metric='euclidean',consol=T,nb.clust=-1) # Hierarchical Clustering on Principle Components
#hc_2<-HCPC(res.famd,order=F)
#hc_2b<-HCPC(res.famd,method='average')
#hc_5<-HCPC(res.famd,nb.clust=5) # Hierarchical Clustering on Principle Components
#Identify 5 main clusters
#hc_5$desc.ind$para
#hc_5$desc.var

## plot results
plot(hc_2)
plot(hc_2, choice="tree")
plot(hc_2, choice="map")
par(mar=c(5,5,1,1))
plot(hc_2_unconsol,choice='map',title="")

##  results
hc_2
hc_2$desc.var
hc_2$desc.axes
hc_2$data.clust
hc_2$call$t

#Paragons (species closest to the cluster centres)
hc_2_unconsol$desc.ind$para

#Plot traits by FG
dataclust<-hc_2_unconsol$data.clust
dataclust$species<-rownames(dataclust)
dataclust$Clustname<-as.character(dataclust$clust)
dataclust$Clustname[dataclust$clust==1]<-'Burrowing\nhindgut-fermenters'
dataclust$Clustname[dataclust$clust==2]<-'Large-bodied\nfacultative generalists'
dataclust$Clustname[dataclust$clust==3]<-'Limnic\nmigrators'

mergeclust<-merge(Traits,dataclust,by.y='species',by.x='Binomial')
par(mfrow=c(4,4))
par(mar=c(5,10,3,0))
for(i in c(23,10)){
 boxplot(as.numeric(mergeclust[,i])~mergeclust$Clustname,main=names(mergeclust)[i],las=1)}
par(mar=c(5,12,3,3))
for(i in c(17:21,16,7,13,14,8,9,11,22,12)){
 par(las=1)
  spineplot(as.factor(mergeclust$Clustname),as.factor(mergeclust[,i]),main=names(mergeclust)[i],xlab="",ylab="")}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
################# mantel test   ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# testing the sensitivity of the clustering/tree structure for a given trait
# not done for now

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#################  export hierarchical classification as newick tree for biodiverse  ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


herbfunctree <- as.phylo(hc_2$call$t$tree) 
plot(herbfunctree)
write.tree(phy=herbfunctree, file="arcticherbivorefunctiontree_apr18.nwk")
