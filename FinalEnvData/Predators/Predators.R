
#Study region
laea<-'+proj=laea +lat_0=90 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0'

#Same as Barrio 2016
#studyregion<-raster(
#studyregion<-projectRaster(studyregion1,crs=laea)
#plot(studyregion)

#Study Region as biodiverse output
studyregion<-raster('S:\\DISENTANGLE\\WP3\\ArcticHerbivoreFDPD\\BiodiverseOutputProvisional\\Preresults_spatial_PD_SR\\PD_spatial_PD_P.tif')
crs(studyregion)<-laea

#Predator list
preds<-read.csv('S:\\DISENTANGLE\\WP3\\ArcticHerbivoreFDPD\\FinalEnvData\\list_spp_PREDATORS.csv',header=T)


# Mammals -----------------------------------------------------------------
#Terrestrial mammals
#Data from IUCN
mamm<-readOGR('S:\\DISENTANGLE\\WP3\\HerbivoreDistributions\\TERRESTRIAL_MAMMALS','TERRESTRIAL_MAMMALS')

predmamm<-(mamm[mamm$binomial%in%preds$Species,])
levels(droplevels(predmamm$binomial))

writeOGR(predmamm,'S:\\DISENTANGLE\\WP3\\ArcticHerbivoreFDPD\\FinalEnvData\\Predators','ArcticMammalPredators',driver='ESRI Shapefile')

mammproj<-spTransform(predmamm,laea)
mammproj@data<-droplevels(mammproj@data)

#Rasterize each spp
mamms<-list()
for (i in 1:length((levels(mammproj$binomial)))){
  mamms[[i]]<-rasterize(mammproj[mammproj$binomial==(levels(mammproj$binomial))[i],]
                        ,studyregion,field=1)
}
mammstack<-stack(mamms)  
names(mammstack)<-levels(mammproj$binomial)
mammstack

#Herbivore eating predators (Not polar bear)
mammstackforherbivores<-mammstack[[c(1:11,13:14)]]

mammpredrich<-calc(mammstack,sum,na.rm=T)
plot(mammpredrich)


# Birds ------------------------------------------------------------------


#Birds
#Data from BirdLife international
birds<-readOGR('S:\\DISENTANGLE\\WP3\\HerbivoreDistributions\\BOTW','BOTW')


predbird<-birds[birds$SCINAME%in%preds$Species,]
levels(droplevels(predbird$SCINAME))

##
sel<-readOGR('T:\\vm\\alle\\Bruker\\James\\Disentangle\\BOTW\\Extract','Selection_extract')
sel
levels(droplevels(sel$Scientific))
writeOGR(sel,'S:\\DISENTANGLE\\WP3\\ArcticHerbivoreFDPD\\FinalEnvData\\Predators','ArcticBirdPredators',driver='ESRI Shapefile')

s2<-spTransform(sel,laea)

birds<-list()
for (i in 1:length(levels(s2$Scientific))){
  birds[[i]]<-rasterize(s2[s2$Scientific==levels(s2$Scientific)[i],],studyregion,field=1)
}
birdstack<-stack(birds)  
names(birdstack)<-levels(s2$Scientific)

#Herbivore eating birds (not Larus or Xema)
birdstackforherbivores<-birdstack[[c(1:16,23:29)]]

birdpredrich<-calc(birdstack,sum,na.rm=T)



# All predators -----------------------------------------------------------
#Total predator richness
ArcticPredators<-stack(birdstack,mammstack)
writeRaster(ArcticPredators,'S:\\DISENTANGLE\\WP3\\ArcticHerbivoreFDPD\\FinalEnvData\\Predators\\ArcticPredators',format='GTiff',overwrite=T)
totalpredrich<-birdpredrich+mammpredrich
plot(totalpredrich)


#Total predator richness feeding on herbivores
totalpredrichforherbivores<-stack(birdstackforherbivores,mammstackforherbivores)
#HerbivoresEatingPredatorRichness
arcticpredatorsforherbivores<-calc(totalpredrichforherbivores,sum,na.rm=T)
plot(arcticpredatorsforherbivores)


arcticpredatorrichness<-stack(birdpredrich,mammpredrich,totalpredrich,arcticpredatorsforherbivores)
names(arcticpredatorrichness)<-c('AvianPredatorRichness','MammalianPredatorRichness','TotalPredatorRichness','TotalPredatorRichnessforHerbivores')
plot(arcticpredatorrichness)

writeRaster(arcticpredatorrichness,'S:\\DISENTANGLE\\WP3\\ArcticHerbivoreFDPD\\FinalEnvData\\Predators\\ArcticPredatorRichness',format='GTiff',overwrite=T)

p.strip <- list(cex=0.7, lines=2)
tiff(8,8,units='in',res=72,pointsize=1,file='S:\\DISENTANGLE\\WP3\\ArcticHerbivoreFDPD\\FinalEnvData\\Predators\\PredDists.tif')
levelplot(ArcticPredators,scales=list(draw=F),colorkey=F,par.strip.text=p.strip)+
layer(sp.points(SpatialPoints(np),col=1)) #+
#layer(sp.polygons(allarc))
dev.off()
