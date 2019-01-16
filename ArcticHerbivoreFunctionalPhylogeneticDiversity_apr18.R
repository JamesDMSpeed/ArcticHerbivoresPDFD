#Spatial phylogenetic analysis using picnate

require(raster)
require(ape)
require(picante)
require(rasterVis)
require(rgdal)
require(hexbin)
require(MuMIn)
require(nlme)

# Basic plottings ---------------------------------------------------------

#Country outlines and north pole
#Lambert azimunthual equal area
laea<-'+proj=laea +lat_0=90 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0 '
charcount <- c('NO', 'SE', 'FI','RU','CA','US','IS','GL','SJ','JP') 
allac2 <- do.call("bind", lapply(charcount, function(x)  raster::getData('GADM', country=x, level=0)))
allarc<-spTransform(allac2,CRS=crs(laea))
np<-SpatialPoints(cbind(0,0))


# Geographic data ---------------------------------------------------------

spplist<-list.files('RangeMaps',full.names=T)

herbstack<-stack(spplist)
herbstack

#Only all reindeer, and drop sheep
use1<-herbstack[[c(1:69,71:72,75:77)]]
names(use1)[3]<-'Rangifer_tarandus'
r1<-raster(use1)

sr_geo<-sum(use1,na.rm=T)
plot(sr_geo)

commdat<-getValues(use1)
#Replace NA with 0
commdat[is.na(commdat)]<-0

# Environmental data ------------------------------------------------------

#Environmental drivers...
listenvvars<-list.files('FinalEnvData/EnvironmentalDrivers',full.names=T)
envvars<-stack(listenvvars)
names(envvars)<-c('AvianPredators','MammalianPredators','TotalPredators','HerbivorePredators',
                  'NDVI','WinterMinTemp','TempRange',
                  'DistanceToCoast','HabitatType','HabitatHet',
                  'IceFreeHistory','CurrentIce','TopographicHet')


realmsP<-readOGR('CMEC regions & realms','Regions')
#Correct misspecified prime meridian of realms data
crs(realmsP)<-'+proj=eqc +lat_ts=0 +lat_0=0 +lon_0=12 +x_0=0 +y_0=0 +datum=WGS84 +pm=0 +units=m +no_defs+ellps=WGS84 +towgs84=0,0,0'
rl<-crop(spTransform(realmsP,crs(use1)),use1)
#plot(rl)
#plot(allarc,add=T,col=NA,border=2)
#realml<-rasterize(rl,use1,field='Regions')
#plot(realml)
#writeRaster(realml,'CMEC regions & realms/Regions.tif',overwrite=T)

realms<-raster('CMEC regions & realms/Regions.tif')
#levelplot(realms)+
#layer(sp.polygons(allarc,add=T))

plot(realms)
sr_geo[sr_geo<1]<-NA
ps1<-rasterToPoints(sr_geo,spatial=T)
points(ps1,cex=0.3)
ps1@data$ext<-extract(realms,ps1)
points(ps1[is.na(ps1$ext),],col=2,pch=16,cex=0.5)

#Make simpler Nearctic/Palaearctic distinction
allarc$bioregion<-allarc$NAME_ENGLISH
allarc$bioregion[allarc$NAME_ISO=='NORWAY'|allarc$NAME_ISO=='SWEDEN'|allarc$NAME_ISO=='FINLAND'|allarc$NAME_ISO=='ICELAND'|allarc$NAME_ISO=='SVALBARD AND JAN MAYEN'|allarc$NAME_ISO=='RUSSIAN FEDERATION']<-'Palaearctic'
allarc$bioregion[allarc$NAME_ISO=='CANADA'|allarc$NAME_ISO=='UNITED STATES'|allarc$NAME_ISO=='GREENLAND']<-'Nearctic'
simprealm<-rasterize(crop(allarc,use1),use1,field='bioregion')

#High Low sub Arctic
#Conservation of Arctic Flora and Fauna Working Group (2010) CAFF Map No.53 - Boundaries of the geographic area covered by the Arctic Biodiversity Assessment.
arczones<-readOGR('ABA-Boundaries','Arctic_Zones')
arczones_laea<-spTransform(arczones,laea)
plot(arczones_laea,col=arczones_laea$Zone)
arczonesR<-rasterize(arczones_laea,use1,field='Zone')

realmsrat<-ratify(mask(realms,arczonesR))
ratzoo<-levels(realmsrat)[[1]]
ratzoo$Realm<-c('Arcto-Siberian','Eurasian','North American')
levels(realmsrat)<-ratzoo

arczonesrat<-ratify(arczonesR)
ratZ<-levels(arczonesrat)[[1]]
ratZ$Zone<-c('High Arctic','Low Arctic','Sub-Arctic')
levels(arczonesrat)<-ratZ

pp1<-levelplot(realmsrat,scales=list(draw=FALSE),margin=F,main='Zoogeographic regions',col.regions=c( 'aquamarine3', 'seagreen1','cyan2'))#+
  layer(sp.polygons(allarc,add=T,lwd=0.5,col=grey(0.5)))+
  layer(sp.points(np,col=1))
pp2<-levelplot(arczonesrat,scales=list(draw=FALSE),margin=F,main='Arctic boundaries',col.regions=c( 'orchid4', 'orchid3','plum1'))#+
  layer(sp.polygons(allarc,add=T,lwd=0.5,col=grey(0.5)))+
  layer(sp.points(np,col=1))
print(pp1, split=c(1, 1, 2, 1), more=TRUE)
print(pp2, split=c(2, 1, 2, 1), more=F)
# Phylogenetic diverstiy --------------------------------------------------

#Trees
phylogeny<-read.tree('Final phylogeny/rooted.nwk')
plot(phylogeny)
phylogeny$tip.label

phylogeny$tip.label[which(!phylogeny$tip.label%in%names(use1))]

#Standardise edge lengths
phylo2<-phylogeny
#phylo2$edge.length<-phylogeny$edge.length/sum(phylogeny$edge.length)

#Use picante to trim community and phylogenetic data
phydata<-match.phylo.comm(phylo2,commdat)

write.tree(phydata$phy,file='Final phylogeny/rooted_trimeed.nwk')
#PD & SR
phydiv<-pd(phydata$comm,phydata$phy,include.root=T)

srras<-raster(use1)
srras<-setValues(srras,phydiv$SR)
srras<-mask(srras,sr_geo,maskvalue=0)

sr_p<-srras/69

pdras<-raster(use1)
pdras<-setValues(pdras,phydiv$PD)
pdras<-mask(pdras,sr_geo,maskvalue=0)
plot(pdras)
pd_p<-pdras/sum(phydata$phy$edge.length)


#Expected PD for given spp richness
expected.pd(phydata$phy)

#Standardised effect sizes

#Different null modesl
pd_es_taxalab<-ses.pd(phydata$comm,phydata$phy,null.model='taxa.labels',runs=1000)
pd_es_indswap<-ses.pd(phydata$comm,phydata$phy,null.model='independentswap',runs=1000)
pd_es_richness<-ses.pd(phydata$comm,phydata$phy,null.model='richness',runs=1000)
#pd_es_frequency<-ses.pd(phydata$comm,phydata$phy,null.model='frequency',runs=1000)
pd_es_samplepool<-ses.pd(phydata$comm,phydata$phy,null.model='sample.pool',runs=1000)
pd_es_trialswap<-ses.pd(phydata$comm,phydata$phy,null.model='trialswap',runs=1000)

pdesras<-raster(use1)

pdestaxalabras<-setValues(pdesras,pd_es_taxalab$pd.obs.p)
pdesindswapras<-setValues(pdesras,pd_es_indswap$pd.obs.p)
pdesrichras<-setValues(pdesras,pd_es_richness$pd.obs.p)
pdessampras<-setValues(pdesras,pd_es_samplepool$pd.obs.p)
pdestrialswapras<-setValues(pdesras,pd_es_trialswap$pd.obs.p)


breaks<-c(0,0.01,0.025,0.975,0.99,1)
cols<-colorRampPalette(c('yellow','gold','grey','blue','darkblue'))
levelplot(mask(stack(pdestaxalabras,pdesindswapras,pdesrichras,pdessampras,pdestrialswapras),sr_geo,maskvalue=0),
          names.attr=c('TaxaLab','Independentent swap','Richness','SamplePool','TrialSwap'),
          at=breaks,col.regions=cols,margin=F)

plot(pd_es_taxalab$pd.obs.p,pd_es_samplepool$pd.obs.p)
plot(pd_es_taxalab$pd.obs.p,pd_es_trialswap$pd.obs.p)
plot(pd_es_taxalab$pd.obs.p,pd_es_indswap$pd.obs.p)

#pd_es<-ses.pd(phydata$comm,phydata$phy,null.model='taxa.labels',runs=1000)
#write.table(pd_es,'PhylogeneticFunctionalAnalysisPicanteR/pdes_apr18.txt')
pd_es<-read.table('PhylogeneticFunctionalAnalysisPicanteR/pdes_apr18.txt')
pd_es[pd_es$ntaxa>0,]

pdesras<-raster(use1)
pdesras<-setValues(pdesras,pd_es$pd.obs.p)
pdesras<-mask(pdesras,sr_geo,maskvalue=0)

breaks<-c(0,0.01,0.025,0.975,0.99,1)
cols<-colorRampPalette(c('yellow','gold','grey','blue','darkblue'))
levelplot(pdesras,at=breaks,col.regions=cols,margin=F)

# Functional diversity ----------------------------------------------------

#FD
functree<-read.tree('Final trait classification/arcticherbivorefunctiontree_apr18.nwk')
plot(functree)
#Adjusting name to match phylogeny and range map
#Check tip label first
functree$tip.label[which(functree$tip.label=='Urocitellus_parryii')]<-'Spermophilus_parryii'
functree2<-functree
#functree2$edge.length<-functree$edge.length/sum(functree$edge.length)

funcdata<-match.phylo.comm(functree2,commdat)

funcdiv<-pd(funcdata$comm,funcdata$phy,include.root=T)

fdras<-raster(use1)
fdras<-setValues(fdras,funcdiv$PD)
fdras<-mask(fdras,sr_geo,maskvalue=0)

plot(fdras)
fd_p<-fdras/sum(funcdata$phy$edge.length)

#fd_es<-ses.pd(funcdata$comm,funcdata$phy,null.model='taxa.labels',runs=1000)
#write.table(fd_es,'PhylogeneticFunctionalAnalysisPicanteR/fdes_apr18.txt')
fd_es<-read.table('PhylogeneticFunctionalAnalysisPicanteR/fdes_apr18.txt')
fd_es[fd_es$ntaxa>0,]

fdesras<-raster(use1)
fdesras<-setValues(fdesras,fd_es$pd.obs.p)
fdesras<-mask(fdesras,sr_geo,maskvalue=0)

levelplot(fdesras,at=breaks,col.regions=cols,margin=F)

# Ratio FD:PD -------------------------------------------------------------
dfdat<-data.frame(realm=getValues(realms),sr=getValues(srras),pd=getValues(pdras),fd=getValues(fdras),sr_p=getValues(sr_p),pd_p=getValues(pd_p),fd_p=getValues(fd_p))
dfdatstack<-stack(realms,divstack2)
dfdat1<-rasterToPoints(dfdatstack)
#write.csv(dfdat1,'PhylogeneticFunctionalAnalysisPicanteR/DiversitywithCoords.csv')
tiff('PhylogeneticFunctionalAnalysisPicanteR/DiversityRatios_6pan.tif',res=150,width=6,height=9,units='in')
{
  colx<-brewer.pal(3, 'YlOrRd')
  par(mfrow=c(3,2))
par(mar=c(5,5,2,1))
with(dfdat,plot(sr_p,pd_p,xlab='Species richness \n (proportion of total)',ylab='Phylogenetic diversity \n (proportion of total)',las=1,ylim=c(0,1),type='n'))
with(dfdat[dfdat$realm==3 & dfdat$sr_p>0,],points(sr_p,pd_p,pch=15,cex=0.7,col=colx[1]))
with(dfdat[dfdat$realm==12 & dfdat$sr_p>0,],points(sr_p,pd_p,pch=16,cex=0.7,col=colx[2]))
with(dfdat[dfdat$realm==6 & dfdat$sr_p>0,],points(sr_p,pd_p,pch=17,cex=0.7,col=colx[3]))
legend('topl',pch=c(15,16,17),pt.cex=1,col=c(colx[1],colx[2],colx[3]),c('Arctico-Siberian','North American','Eurasian'))
abline(0,1)
mtext(side=3,adj=0,'(a)',line=0.5,cex=1)
with(dfdat,plot(sr_p,fd_p,xlab='Species richness \n (proportion of total)',ylab='Functional diversity \n (proportion of total)',las=1,ylim=c(0,1),type='n'))
with(dfdat[dfdat$realm==3 & dfdat$sr_p>0,],points(sr_p,fd_p,pch=15,cex=0.7,col=colx[1]))
with(dfdat[dfdat$realm==12 & dfdat$sr_p>0,],points(sr_p,fd_p,pch=16,cex=0.7,col=colx[2]))
with(dfdat[dfdat$realm==6 & dfdat$sr_p>0,],points(sr_p,fd_p,pch=17,cex=0.7,col=colx[3]))
abline(0,1)
mtext(side=3,adj=0,'(b)',line=0.5,cex=1)
with(dfdat,plot(pd_p,fd_p,xlab='Phylogenetic diversity \n (proportion of total)',ylab='Functional diversity \n (proportion of total)',las=1,xlim=c(0,0.75),ylim=c(0,1),type='n'))
with(dfdat[dfdat$realm==3 & dfdat$sr_p>0,],points(pd_p,fd_p,pch=15,cex=0.7,col=colx[1]))
with(dfdat[dfdat$realm==12 & dfdat$sr_p>0,],points(pd_p,fd_p,pch=16,cex=0.7,col=colx[2]))
with(dfdat[dfdat$realm==6 & dfdat$sr_p>0,],points(pd_p,fd_p,pch=17,cex=0.7,col=colx[3]))
abline(0,1)
mtext(side=3,adj=0,'(c)',line=0.5,cex=1)
plot(sr_p*70,fd_p/pd_p,xlab='Species richness',ylab='Functional divergence',las=1,ylim=c(0,5),type='n')
with(dfdat[dfdat$realm==3 & dfdat$sr_p>0,],points(sr_p*70,fd_p/pd_p,pch=15,cex=0.7,col=colx[1]))
with(dfdat[dfdat$realm==12 & dfdat$sr_p>0,],points(sr_p*70,fd_p/pd_p,pch=16,cex=0.7,col=colx[2]))
with(dfdat[dfdat$realm==6 & dfdat$sr_p>0,],points(sr_p*70,fd_p/pd_p,pch=17,cex=0.7,col=colx[3]))
mtext(side=3,adj=0,'(d)',line=0.5,cex=1)
with(dfdat,plot(pd_p,fd_p/pd_p,xlab='Phylogenetic diversity \n (proportion of total)',ylab='Functional divergence',las=1,ylim=c(0,5),type='n'))
with(dfdat[dfdat$realm==3 & dfdat$sr_p>0,],points(pd_p,fd_p/pd_p,pch=15,cex=0.7,col=colx[1]))
with(dfdat[dfdat$realm==12 & dfdat$sr_p>0,],points(pd_p,fd_p/pd_p,pch=16,cex=0.7,col=colx[2]))
with(dfdat[dfdat$realm==6 & dfdat$sr_p>0,],points(pd_p,fd_p/pd_p,pch=17,cex=0.7,col=colx[3]))
mtext(side=3,adj=0,'(e)',line=0.5,cex=1)
with(dfdat,plot(pd_p,fd_p/pd_p,xlab='Functional diversity \n (proportion of total)',ylab='Functional divergence',las=1,ylim=c(0,5),type='n',xlim=c(0,1)))
with(dfdat[dfdat$realm==3 & dfdat$sr_p>0,],points(fd_p,fd_p/pd_p,pch=15,cex=0.7,col=colx[1]))
with(dfdat[dfdat$realm==12 & dfdat$sr_p>0,],points(fd_p,fd_p/pd_p,pch=16,cex=0.7,col=colx[2]))
with(dfdat[dfdat$realm==6 & dfdat$sr_p>0,],points(fd_p,fd_p/pd_p,pch=17,cex=0.7,col=colx[3]))
mtext(side=3,adj=0,'(f)',line=0.5,cex=1)
dev.off()
}

tiff('PhylogeneticFunctionalAnalysisPicanteR/DiversityRatios_4pan_R1.tif',res=150,width=6,height=6,units='in')
{
  colx<-brewer.pal(3, 'YlOrRd')
  par(mfrow=c(2,2))
  par(mar=c(5,5,2,1))
  with(dfdat,plot(sr_p,pd_p,xlab='Species richness \n (proportion of total)',ylab='Phylogenetic diversity \n (proportion of total)',las=1,ylim=c(0,1),type='n'))
  with(dfdat[dfdat$realm==3 & dfdat$sr_p>0,],points(sr_p,pd_p,pch=15,cex=0.7,col=colx[1]))
  with(dfdat[dfdat$realm==12 & dfdat$sr_p>0,],points(sr_p,pd_p,pch=16,cex=0.7,col=colx[2]))
  with(dfdat[dfdat$realm==6 & dfdat$sr_p>0,],points(sr_p,pd_p,pch=17,cex=0.7,col=colx[3]))
  abline(0,1)
  mtext(side=3,adj=0,'(a)',line=0.5,cex=1)
  with(dfdat,plot(sr_p,fd_p,xlab='Species richness \n (proportion of total)',ylab='Functional diversity \n (proportion of total)',las=1,ylim=c(0,1),type='n'))
  with(dfdat[dfdat$realm==3 & dfdat$sr_p>0,],points(sr_p,fd_p,pch=15,cex=0.7,col=colx[1]))
  with(dfdat[dfdat$realm==12 & dfdat$sr_p>0,],points(sr_p,fd_p,pch=16,cex=0.7,col=colx[2]))
  with(dfdat[dfdat$realm==6 & dfdat$sr_p>0,],points(sr_p,fd_p,pch=17,cex=0.7,col=colx[3]))
  abline(0,1)
  mtext(side=3,adj=0,'(b)',line=0.5,cex=1)
  plot(pd_p,fd_p,xlab='Phylogenetic diversity \n (proportion of total)',ylab='Functional diversity \n (proportion of total)',las=1,xlim=c(0,0.75),ylim=c(0,1),type='n')
  with(dfdat[dfdat$realm==3 & dfdat$sr_p>0,],points(pd_p,fd_p,pch=15,cex=0.7,col=colx[1]))
  with(dfdat[dfdat$realm==12 & dfdat$sr_p>0,],points(pd_p,fd_p,pch=16,cex=0.7,col=colx[2]))
  with(dfdat[dfdat$realm==6 & dfdat$sr_p>0,],points(pd_p,fd_p,pch=17,cex=0.7,col=colx[3]))
  abline(0,1)
  mtext(side=3,adj=0,'(c)',line=0.5,cex=1)
  plot(sr_p*70,(1-fd_p/pd_p),xlab='Species richness',ylab='Functional convergence',las=1,ylim=c(-5,0),type='n')
  with(dfdat[dfdat$realm==3 & dfdat$sr_p>0,],points(sr_p*70,(1-fd_p/pd_p),pch=15,cex=0.7,col=colx[1]))
  with(dfdat[dfdat$realm==12 & dfdat$sr_p>0,],points(sr_p*70,(1-fd_p/pd_p),pch=16,cex=0.7,col=colx[2]))
  with(dfdat[dfdat$realm==6 & dfdat$sr_p>0,],points(sr_p*70,(1-fd_p/pd_p),pch=17,cex=0.7,col=colx[3]))
  legend('bottomr',pch=c(15,16,17),pt.cex=1,col=c(colx[1],colx[2],colx[3]),c('Arctico-Siberian','North American','Eurasian'))
      mtext(side=3,adj=0,'(d)',line=0.5,cex=1)
    dev.off()
}

#Different scaling
par(mfrow=c(2,3))
plot(dfdat$pd,dfdat$fd,main='Raw')
plot(dfdat$pd_p,dfdat$fd_p,main='Standardised by total \nbranch length')
plot(scale(dfdat$pd),scale(dfdat$fd),main='Centered and scaled')
plot(dfdat$sr,(1-dfdat$fd/dfdat$pd))
plot(dfdat$sr,(1-dfdat$fd_p/dfdat$pd_p))
plot(dfdat$sr,1-(scale(dfdat$fd)/(scale(dfdat$pd))))


#Models
with(dfdat[dfdat$realm!=9,],summary(step(lm(pd_p~sr_p*as.factor(realm)))))

lm1<-lm(dfdat$pd_p~dfdat$sr_p)
summary(lm1)
fdiv_var<-(dfdat$fd_p/dfdat$pd_p)
sr_var<-dfdat$sr_p*70
lm2<-lm(fdiv_var~sr_var)
abline(lm2)

levelplot(fdras/pdras,margin=F)
plot(divstack$Phylogenetic.diversity,divstack$Functional.diversity,las=1,xlim=c(0,1),ylim=c(0,1))
abline(0,1,col=grey(0.5))

exppd<-variance.pd(phydata$phy)
expfd<-variance.pd(funcdata$phy)
expfd$nspp<-as.numeric(rownames(expfd))
expfd

merge1<-merge(pd_es,expfd,by.x='ntaxa',by.y='nspp')


#PDFD pairwise
par(mfrow=c(1,2))
par(mar=c(5,5,2,1))
with(dfdat,plot(pd_p,fd_p,xlab='Phylogenetic diversity \n (proportion of total)',ylab='Functional diversity \n (proportion of total)',las=1,ylim=c(0,1),type='n'))
with(dfdat[dfdat$realm==3 & dfdat$sr_p>0,],points(pd_p,fd_p,pch=15,cex=0.7,col=colx[1]))
with(dfdat[dfdat$realm==12 & dfdat$sr_p>0,],points(pd_p,fd_p,pch=16,cex=0.7,col=colx[2]))
with(dfdat[dfdat$realm==6 & dfdat$sr_p>0,],points(pd_p,fd_p,pch=17,cex=0.7,col=colx[3]))
abline(0,1)
with(dfdat,plot(pd_p,fd_p/pd_p,xlab='Phylogenetic diversity \n (proportion of total)',ylab='Functional divergence',las=1,ylim=c(0,5),type='n'))
with(dfdat[dfdat$realm==3 & dfdat$sr_p>0,],points(pd_p,fd_p/pd_p,pch=15,cex=0.7,col=colx[1]))
with(dfdat[dfdat$realm==12 & dfdat$sr_p>0,],points(pd_p,fd_p/pd_p,pch=16,cex=0.7,col=colx[2]))
with(dfdat[dfdat$realm==6 & dfdat$sr_p>0,],points(pd_p,fd_p/pd_p,pch=17,cex=0.7,col=colx[3]))
with(dfdat,plot(pd_p,fd_p/pd_p,xlab='Functional diversity \n (proportion of total)',ylab='Functional divergence',las=1,ylim=c(0,5),type='n',xlim=c(0,1)))
with(dfdat[dfdat$realm==3 & dfdat$sr_p>0,],points(fd_p,fd_p/pd_p,pch=15,cex=0.7,col=colx[1]))
with(dfdat[dfdat$realm==12 & dfdat$sr_p>0,],points(fd_p,fd_p/pd_p,pch=16,cex=0.7,col=colx[2]))
with(dfdat[dfdat$realm==6 & dfdat$sr_p>0,],points(fd_p,fd_p/pd_p,pch=17,cex=0.7,col=colx[3]))




#Function to do FDPD hacked from picante ####
pdtreeratio<-function(samp,tree,tree2,runs=999,iterations=1000) #For taxa.labels null model
{pd.obs <- as.vector(pd(samp, tree)$PD)/sum(tree$edge.length)
pd.rand<-t(replicate(runs,as.vector(pd(samp, tipShuffle(tree2))$PD)/sum(tree2$edge.length)))
pd.rand.mean <- apply(X = pd.rand, MARGIN = 2, FUN = mean, 
                      na.rm = TRUE)
pd.rand.sd <- apply(X = pd.rand, MARGIN = 2, FUN = sd, na.rm = TRUE)
pd.obs.z <- (pd.obs - pd.rand.mean)/pd.rand.sd
pd.obs.rank <- apply(X = rbind(pd.obs, pd.rand), MARGIN = 2, 
                     FUN = rank)[1, ]
pd.obs.rank <- ifelse(is.na(pd.rand.mean), NA, pd.obs.rank)
data.frame(ntaxa = specnumber(samp), pd.obs, pd.rand.mean, 
           pd.rand.sd, pd.obs.rank, pd.obs.z, pd.obs.p = pd.obs.rank/(runs + 
                                                                        1), runs = runs, row.names = row.names(samp))
}

#Make a randomisation that randomises FD:PD
fdpdratiofunc<-function(samp,tree1,tree2,runs=999)
{
  fdpd.obs<-as.vector(pd(samp, tree1)$PD)/sum(tree1$edge.length)/as.vector(pd(samp, tree2)$PD)/sum(tree2$edge.length)
  fdpd.rand<-t(replicate(runs,as.vector(pd(samp, tipShuffle(tree1))$PD)/sum(tree1$edge.length))
               /as.vector(pd(samp,tipShuffle(tree2))$PD)/sum(tree2$edge.length))
  fdpd.rand.mean <- apply(X = fdpd.rand, MARGIN = 2, FUN = mean, 
                          na.rm = TRUE)
  fdpd.rand.sd <- apply(X = fdpd.rand, MARGIN = 2, FUN = sd, na.rm = TRUE)
  fdpd.obs.z <- (fdpd.obs - fdpd.rand.mean)/fdpd.rand.sd
  fdpd.obs.rank <- apply(X = rbind(fdpd.obs, fdpd.rand), MARGIN = 2, 
                         FUN = rank)[1, ]
  fdpd.obs.rank <- ifelse(is.na(fdpd.rand.mean), NA, fdpd.obs.rank)
  data.frame(ntaxa=specnumber(samp),fdpd.obs,fdpd.rand.mean,
             fdpd.rand.sd,fdpd.obs.rank,fdpd.obs.z,
             fdpd.obs.p = fdpd.obs.rank/(runs + 
                                           1), runs = runs, row.names = row.names(samp))
}
#fdpdrat<-fdpdratiofunc(phydata$comm,funcdata$phy,phydata$phy,runs=1000)
#write.table(fdpdrat,'PhylogeneticFunctionalAnalysisPicanteR/fdpd_es_apr18.txt')
fdpdrat<-read.table('PhylogeneticFunctionalAnalysisPicanteR/fdpd_es_apr18.txt')
fdpdrat[fdpdrat$ntaxa>1,]
fdpdnotsingle<-fdpdrat
fdpdnotsingle$fdpd.obs.z[fdpdnotsingle$ntaxa<=1]<-NA#Remove communities with 1spp
fdpdnotsingle$fdpd.obs.p[fdpdnotsingle$ntaxa<=1]<-NA#Remove communities with 1spp
fdpdlp<-levelplot(setValues(r1,fdpdnotsingle$fdpd.obs.z),margin=F,main='FD:PD effect size',scales=list(draw=FALSE))
diverge0(fdpdlp,'RdBu')
levelplot(setValues(r1,fdpdrat$fdpd.obs.p),at=breaks,col.regions=cols,margin=F,main='Randomisation FD:PD',scales=list(draw=FALSE))


#fdpd<-pdtreeratio(funcdata$comm,funcdata$phy,phydata$phy,runs=1000)
#write.table(fdpd,'PhylogeneticFunctionalAnalysisPicanteR/fdpd_es_apr18.txt')
fdpd<-read.table('PhylogeneticFunctionalAnalysisPicanteR/fdpd_es_apr18.txt')
fdpd[fdpd$ntaxa>1,]
fdpdras<-raster(use1)
fdpdras<-mask(setValues(fdpdras,fdpd$fdpd.obs.z),sr_geo,maskvalue=0)
levelplot(fdpdras,par.settings=YlOrRdTheme,scales=list(draw=FALSE),margin=F,main='FD:PD')



# Logistic relationship FDPD ----------------------------------------------
df2<-dfdat
df2$pd_p[df2$sr<1]<-NA
lmFDPD<-with(df2,lm(fd_p~log(pd_p)))
summary(lmFDPD)
newdat<-data.frame(pd_p=seq(min(df2$pd_p,na.rm=T),max(df2$pd_p,na.rm=T),by=0.01))
newdatpreds<-predict(lmFDPD,newdat,type='response',interval='c',level=0.99)
logfdpd<-cbind(newdat,newdatpreds)

with(dfdat,plot(pd_p,fd_p,xlab='Phylogenetic diversity \n (proportion of total)',ylab='Functional diversity \n (proportion of total)',las=1,xlim=c(0,0.75),ylim=c(0,1),type='n'))
with(dfdat[dfdat$realm==3 & dfdat$sr_p>0,],points(pd_p,fd_p,pch=15,cex=0.7,col=colx[1]))
with(dfdat[dfdat$realm==12 & dfdat$sr_p>0,],points(pd_p,fd_p,pch=16,cex=0.7,col=colx[2]))
with(dfdat[dfdat$realm==6 & dfdat$sr_p>0,],points(pd_p,fd_p,pch=17,cex=0.7,col=colx[3]))
lines(logfdpd$pd_p,logfdpd$fit)
lines(logfdpd$pd_p,logfdpd$lwr,lty=2)
lines(logfdpd$pd_p,logfdpd$upr,lty=2)


#Model with all raster cells
summary(lmFDPD)
func_con<-residuals(lmFDPD)

#divstack2$FuncPhy_logmod<-divstack2$Phylogenetic.diversity
#divstack2$FuncPhy_logmod[!is.na(divstack2$Phylogenetic.diversity)]<-func_con

func_con_ras<-fdras
func_con_ras[pdras>0 & srras>1]<-func_con
func_con_ras[srras<2]<-NA
func_conras<-mask(func_con_ras,extend(envvars$CurrentIce,pd_p),maskvalue=1,updatevalue=NA)
plot(func_conras)
lp1<-levelplot(func_conras,main='FC',margin=F)
diverge0(lp1,'RdBu')


#Randomize matrix
matlayers<-1000
rm1<-array(dim=c(nrow(phydata$comm),ncol(phydata$comm),matlayers))
for(i in 1:matlayers){
rm1[,,i]<-randomizeMatrix(phydata$comm,null.model='richness')}
dimnames(rm1)[[2]]<-colnames(phydata$comm)
pdmat<-matrix(nrow=dim(rm1)[1],ncol=matlayers)
fdmat<-matrix(nrow=dim(rm1)[1],ncol=matlayers)
for(i in 1:matlayers){
pdmat[,i]<-  pd(rm1[,,i],phydata$phy)$PD/sum(phydata$phy$edge.length)
fdmat[,i]<-  pd(rm1[,,i],funcdata$phy)$PD/sum(funcdata$phy$edge.length)
}

dfRand<-data.frame(pdrand=apply(pdmat,1,mean),fdrand=apply(fdmat,1,mean))
lmRand<-with(dfRand[dfRand$pdrand>0,],lm(fdrand~log(pdrand)))
with(dfRand[dfRand$pdrand>0,],plot(pdrand,fdrand,pch=16))
with(dfdat[dfdat$realm==3 & dfdat$sr_p>0,],points(pd_p,fd_p,pch=15,cex=0.7,col=colx[1]))
with(dfdat[dfdat$realm==12 & dfdat$sr_p>0,],points(pd_p,fd_p,pch=16,cex=0.7,col=colx[2]))
with(dfdat[dfdat$realm==6 & dfdat$sr_p>0,],points(pd_p,fd_p,pch=17,cex=0.7,col=colx[3]))
nd1<-data.frame(pdrand=seq(min(dfRand$pdrand[dfRand$pdrand>0]),max(dfRand$pdrand),by=0.001))
pd1<-as.data.frame(predict(lmRand,nd1,type='response',interval='p'))
lines(nd1$pdrand,pd1$fit)
lines(nd1$pdrand,pd1$upr,lty=2)
lines(nd1$pdrand,pd1$lwr,lty=2)

#Effect size 
#Obs-mean.rand/se.rand
resi<-matrix(nrow=1581,ncol=matlayers)
for (i in 1:matlayers){
  lmi<-lm(fdmat[,i][pdmat[,i]>0]~log(pdmat[,i][pdmat[,1]>0]))
  resi[,i]<-residuals(lmi)
}
meanresid<-apply(resi,1,mean)
sdresid<-apply(resi,1,sd)
zresid<-(residuals(lmFDPD)-meanresid)/sdresid
rankresid <- apply(X = cbind(residuals(lmFDPD), resi), MARGIN = 1, 
                     FUN = rank)[1, ]
p_resid<-rankresid/(matlayers+1)
p_resid[p_resid>0.975]
p_resid[p_resid<0.025]

funccon_es<-fdras
funccon_es[srras>0]<-zresid
funccon_es[srras==0]<-NA
funccon_es<-mask(funccon_es,extend(envvars$CurrentIce,pd_p),maskvalue=1,updatevalue=NA)

lpfces<-levelplot(funccon_es,margin=F)
diverge0(lpfces,'RdBu')

rsigfunccon<-fdras
rsigfunccon[srras>0]<-p_resid
rsigfunccon[srras==0]<-NA
rsigfunccon<-mask(rsigfunccon,extend(envvars$CurrentIce,pd_p),maskvalue=1,updatevalue=NA)
a<-rsigfunccon
a[a<0.975 & a>0.025]<-NA
a[a>0.975]<-1
a[a<0.025]<-0
rsigfuncconcells<-rasterToPolygons(a,dissolve=T)

lpfces<-levelplot(funccon_es,margin=F)+
 layer( sp.polygons(rsigfuncconcells))
diverge0(lpfces,'RdBu')

writeRaster(funccon_es,'logFuncDivEffSize.tif')
writeRaster(func_conras,'logFuncDiv.tif')


# Stacks ------------------------------------------------------------------


#Diverge function for colouring plots centered on 0
diverge0 <- function(p, ramp) {
  require(RColorBrewer)
  require(rasterVis)
  if(length(ramp)==1 && is.character(ramp) && ramp %in% 
     row.names(brewer.pal.info)) {
    ramp <- suppressWarnings(colorRampPalette(rev(brewer.pal(11, ramp))))
  } else if(length(ramp) > 1 && is.character(ramp) && all(ramp %in% colors())) {
    ramp <- colorRampPalette(ramp)
  } else if(!is.function(ramp)) 
    stop('ramp should be either the name of a RColorBrewer palette, ', 
         'a vector of colours to be interpolated, or a colorRampPalette.')
  rng <- range(p$legend[[1]]$args$key$at)
  s <- seq(-max(abs(rng)), max(abs(rng)), len=1001)
  i <- findInterval(rng[which.min(abs(rng))], s)
  zlim <- switch(which.min(abs(rng)), `1`=i:(1000+1), `2`=1:(i+1))
  p$legend[[1]]$args$key$at <- s[zlim]
  p$par.settings$regions$col <- ramp(1000)[zlim[-length(zlim)]]
  p
}

divstack<-stack(sr_p,pd_p,fd_p)
names(divstack)<-c('Species richness','Phylogenetic diversity','Functional diversity')
levelplot(divstack,par.settings=YlOrRdTheme,scales=list(draw=FALSE))#+
  layer(sp.points(np,col=1))+
  layer(sp.polygons(allarc,col='grey',lwd=0.5))

#divstack2<-stack(sr_p,pd_p,fd_p,0-func_conras)
#names(divstack2)<-c('Species richness','Phylogenetic diversity','Functional diversity','Functional convergence')
#divstack2<-mask(divstack2,sr_geo)
#divstack2m<-mask(divstack2,extend(envvars$CurrentIce,divstack2),maskvalue=1,updatevalue=NA)
#plot(divstack2m)
#writeRaster(divstack2m,'PhylogeneticFunctionalAnalysisPicanteR/DiversityStack')
divstack2<-stack('PhylogeneticFunctionalAnalysisPicanteR/DiversityStack')

my.at <- seq(0, 1, by = 0.1)
levelplot(divstack2,at=my.at,par.settings=YlOrRdTheme())
quantile(divstack2)

writeRaster(divstack2,'PhylogeneticFunctionalAnalysisPicanteR/DiversityPatterns/ArcticHerbivore',format='GTiff',bylayer=T,suffix=names(divstack2))  

divstackll<-projectRaster(divstack2,crs='+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0' )
KML(divstackll,'PhylogeneticFunctionalAnalysisPicanteR/Divstacks')

tiff('PhylogeneticFunctionalAnalysisPicanteR/DiversityMaps_R1.tif',width = 6,height=6,units='in',res=300)
p1 <- levelplot(divstack2[[1]], at=my.at,par.settings=YlOrRdTheme(layout.heights=list(top.padding=0,bottom.padding=0)),scales=list(draw=FALSE), margin=FALSE,main=list('Species richness',cex=0.8))+
  layer(sp.points(np,col=1))+
  layer(sp.polygons(allarc,lwd=0.5,col=grey(0.5)))+
  layer(sp.polygons(rl,lwd=0.4,col='blue'))+
  layer(sp.polygons(arczones_laea,lwd=0.4,lty=2,col='black'))+
  layer(panel.text(-2500000,-1800000,'AS',col='blue'))+ 
  layer(panel.text(-1500000, 2500000,'EUR',col='blue')) +
  layer(panel.text( 2500000,-1800000,'NA',col='blue')) 
p2 <- levelplot(divstack2[[2]], at=my.at, par.settings=YlOrRdTheme(layout.heights=list(top.padding=0,bottom.padding=0)),scales=list(draw=FALSE), margin=FALSE,main=list('Phylogenetic diversity',cex=0.8))+
  layer(sp.points(np,col=1))+
  layer(sp.polygons(allarc,lwd=0.5,col=grey(0.5)))
p3 <- levelplot(divstack2[[3]], at=my.at,par.settings=YlOrRdTheme(layout.heights=list(top.padding=0,bottom.padding=0)),scales=list(draw=FALSE), margin=FALSE,main=list('Functional diversity',cex=0.8))+
  layer(sp.points(np,col=1))+
  layer(sp.polygons(allarc,lwd=0.5,col=grey(0.5)))
p4 <- levelplot(divstack2[[4]],par.settings=YlOrRdTheme(layout.heights=list(top.padding=0,bottom.padding=0)),scales=list(draw=FALSE), margin=FALSE,main=list('Functional convergence',cex=0.8))+
  layer(sp.points(np,col=1))+
  layer(sp.polygons(allarc,lwd=0.5,col=grey(0.5)))

  print(p1, split=c(1, 1, 2, 2), more=TRUE)
  print(p2, split=c(2, 1, 2, 2), more=TRUE)
  print(p3, split=c(1, 2, 2, 2), more=TRUE)
  print(diverge0(p4,'RdBu'), split=c(2, 2, 2, 2))
dev.off()    


#Effect sizes and randomisation tests
#effsizestack<-stack(setValues(r1,pd_es$pd.obs.z),setValues(r1,fd_es$pd.obs.z),0-funccon_es)
#writeRaster(effsizestack,'PhylogeneticFunctionalAnalysisPicanteR/EffectSizesStack')
effsizestack<-stack('PhylogeneticFunctionalAnalysisPicanteR/EffectSizesStack')
names(effsizestack)<-c('Phylogenetic diversity_es','Functional diversity_es','Functional convergence_es')
effsizestack<-mask(effsizestack,envvars1$CurrentIce,maskvalue=1)

writeRaster(effsizestack,'PhylogeneticFunctionalAnalysisPicanteR/DiversityPatterns/ArcticHerbivore',format='GTiff',bylayer=T,suffix=names(effsizestack))

#randomstack<-stack(pdesras,fdesras,rsigfunccon)
#randomstackm<-mask(randomstack,sr_geo)
#writeRaster(randomstackm,'PhylogeneticFunctionalAnalysisPicanteR/RandomStack')
randomstack<-stack('PhylogeneticFunctionalAnalysisPicanteR/RandomStack')
randomstack<-mask(randomstack,envvars1$CurrentIce,maskvalue=1)
names(randomstack)<-c('Phylogenetic diversity_rank','Functional diversity_rank','Functional divergence_rank')
writeRaster(randomstack,'PhylogeneticFunctionalAnalysisPicanteR/DiversityPatterns/ArcticHerbivore',format='GTiff',bylayer=T,suffix=names(randomstack))  

levelplot(randomstack,at=breaks,col.regions=cols,margin=F,scales=list(draw=FALSE),main='Randomisation test')#+
# layer(sp.points(np,col=1))+
# layer(sp.polygons(allarc,lwd=0.5))

#Plot combining effect size and randomisation tests
#Make a rasterToPolygon for the randomisation tests - outlining cells statistically different
rsig<-list()
for (i in 1:3){
rsig[[i]]<-randomstack[[i]]
a<-rsig[[i]]
a[a<0.975 & a>0.025]<-NA
a[a>0.975]<-1
a[a<0.025]<-0
rsig[[i]]<-rasterToPolygons(a,dissolve=T)
}


tiff('PhylogeneticFunctionalAnalysisPicanteR/DiversityEffectSizeMaps_R1.tif',width = 8,height=3,units='in',res=300)
p.strip <- list(cex=0.8)
effsigplot<-levelplot(effsizestack,scales=list(draw=FALSE),names.attr=c('Phylogenetic diversity','Functional diversity','Functional convergence'),main='', par.strip.text=p.strip)+
 # layer(sp.polygons(allarc,lwd=0.5,col=grey(0.5)))+
    layer(sp.polygons(rsig[[panel.number()]],col='black'))+
  layer(sp.points(np,col=1))
diverge0(effsigplot,'RdBu')
dev.off()  


#Number of significant FD cells 
length(which(getValues(randomstack[[2]])>0.975))

#Effect sizes
#NA subarctic
summary(effsizestack$Functional.divergence_es[randomstack$layer.3<0.025 & realmsrat==12 & arczonesrat==3])
summary(effsizestack$Functional.divergence_es[randomstack$layer.3>0.975 & (realmsrat==3 | realmsrat==6 | (realmsrat==12 & arczonesrat<3))])

# Environmental drivers  --------------------------------------------------

#Pariwise correlations between environmental variables
pairs(envvars[[c(4,5,6,10,11,13)]])
envvars1<-extend(envvars,divstack2)
alllayers<-stack(divstack2,effsizestack,randomstack,realms,envvars1,arczonesR)
names(alllayers)[25]<-'ArcticZone'
alldata<-extract(alllayers,rasterToPoints(divstack2$Species.richness,spatial=T),sp=T)
alldatadf<-cbind(alldata@coords,alldata@data)
alldatadf$HerbivorePredators_R<-resid(lm(alldatadf$HerbivorePredators~alldatadf$NDVI))  

with(alldatadf,boxplot(WinterMinTemp~ArcticZone))#Can't use Arctic subzone as a predictor due to colinearlity with temperature

standdf<-cbind(alldatadf[,1:2],scale(alldatadf[,3:29],scale=T,center=T))
write.table(alldatadf,'PhylogeneticFunctionalAnalysisPicanteR/analysisdataframe.txt')
write.table(standdf,'PhylogeneticFunctionalAnalysisPicanteR/standardised_analysisdataframe.txt')

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

tiff('PhylogeneticFunctionalAnalysisPicanteR/EnvVarsPairs.tif',res=150,width=7,height=7,units='in')
pairs(alldatadf[,c(19,18,20,24,27,25)],upper.panel=panel.cor
      ,labels=c('Vegetation \nproductivity (NDVI)','Predator \n richness','Winter minimum\n temperature (x10)','Habitat \n heterogeneity','Topographic \n heterogeneity','Ice free period'))
dev.off()

tiff('PhylogeneticFunctionalAnalysisPicanteR/ResponseEnvPairs.tif',res=150,width=12,height=7,units='in')
par(mfrow=c(4,7))
plot(alldatadf$Species.richness~alldatadf$NDVI,ylab='Species richness',xlab='Vegetation productivity \nNDVI')
plot(alldatadf$Species.richness~alldatadf$HerbivorePredators_R,ylab='Species richness',xlab='Predator richness')
plot(alldatadf$Species.richness~alldatadf$WinterMinTemp,ylab='Species richness',xlab='Winter minimum temperature \nx10')
plot(alldatadf$Species.richness~alldatadf$HabitatHet,ylab='Species richness',xlab='Habitat heterogeneity')
plot(alldatadf$Species.richness~alldatadf$TopographicHet,ylab='Species richness',xlab='Topographic heterogeneity')
plot(alldatadf$Species.richness~alldatadf$IceFreeHistory,ylab='Species richness',xlab='Ice free history')
boxplot(alldatadf$Species.richness~alldatadf$Regions,ylab='Species richness',names=c('Eur','Arc-Sib','N Amer'),xlab='Region')
plot(alldatadf$Phylogenetic.diversity~alldatadf$NDVI,ylab='Phylogenetic diversity',xlab='Vegetation productivity\nNDVI')
plot(alldatadf$Phylogenetic.diversity~alldatadf$HerbivorePredators_R,ylab='Phylogenetic diversity',xlab='Predator richness')
plot(alldatadf$Phylogenetic.diversity~alldatadf$WinterMinTemp,ylab='Phylogenetic diversity',xlab='Winter minimum temperature \nx10')
plot(alldatadf$Phylogenetic.diversity~alldatadf$HabitatHet,ylab='Phylogenetic diversity',xlab='Habitat heterogeneity')
plot(alldatadf$Phylogenetic.diversity~alldatadf$TopographicHet,ylab='Phylogenetic diversity',xlab='Topographic heterogeneity')
plot(alldatadf$Phylogenetic.diversity~alldatadf$IceFreeHistory,ylab='Phylogenetic diversity',xlab='Ice free history')
boxplot(alldatadf$Phylogenetic.diversity~alldatadf$Regions,ylab='Phylogenetic diversity',names=c('Eur','Arc-Sib','N Amer'),xlab='Region')
plot(alldatadf$Functional.diversity~alldatadf$NDVI,ylab='Functional diversity',xlab='Vegetation productivity\nNDVI')
plot(alldatadf$Functional.diversity~alldatadf$HerbivorePredators_R,ylab='Functional diversity',xlab='Predator richness')
plot(alldatadf$Functional.diversity~alldatadf$WinterMinTemp,ylab='Functional diversity',xlab='Winter minimum temperature \nx10')
plot(alldatadf$Functional.diversity~alldatadf$HabitatHet,ylab='Functional diversity',xlab='Habitat heterogeneity')
plot(alldatadf$Functional.diversity~alldatadf$TopographicHet,ylab='Functional diversity',xlab='Topographic heterogeneity')
plot(alldatadf$Functional.diversity~alldatadf$IceFreeHistory,ylab='Functional diversity',xlab='Ice free history')
boxplot(alldatadf$Functional.diversity~alldatadf$Regions,ylab='Functional diversity',names=c('Eur','Arc-Sib','N Amer'),xlab='Region')
plot(alldatadf$Functional.Phylognetic.diversity~alldatadf$NDVI,ylab='Functional dispersion',xlab='Vegetation productivity\nNDVI')
plot(alldatadf$Functional.Phylognetic.diversity~alldatadf$HerbivorePredators_R,ylab='Functional dispersion',xlab='Predator richness')
plot(alldatadf$Functional.Phylognetic.diversity~alldatadf$WinterMinTemp,ylab='Functional dispersion',xlab='Winter minimum temperature\nx10')
plot(alldatadf$Functional.Phylognetic.diversity~alldatadf$HabitatHet,ylab='Functional dispersion',xlab='Habitat heterogeneity')
plot(alldatadf$Functional.Phylognetic.diversity~alldatadf$TopographicHet,ylab='Functional dispersion',xlab='Topographic heterogeneity')
plot(alldatadf$Functional.Phylognetic.diversity~alldatadf$IceFreeHistory,ylab='Functional dispersion',xlab='Ice free history')
boxplot(alldatadf$Functional.Phylognetic.diversity~alldatadf$Regions,ylab='Functional dispersion',names=c('Eur','Arc-Sib','N Amer'),xlab='Region')

dev.off()


# GLS ---------------------------------------------------------------------
#Run in Franklin
#ArcticHerbivore analysis

# GLS of diversity and drivers --------------------------------------------


#rm(list=ls())

require(nlme)
require(MuMIn)
fulldf<-read.table('PhylogeneticFunctionalAnalysisPicanteR/analysisdataframe.txt')

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}
pairs(fulldf[,c(19,18,20,24,27,25)],upper.panel=panel.cor
      ,labels=c('NDVI','Predator \n diversity','Winter minimum\n temperature','Habitat \n heterogeneity','Topographic \n heterogeneity','Ice free period'))

pairs(fulldf[,c(4,5,6,7,19,18,20,24,27,25)],upper.panel=panel.cor
      ,labels=c('Species richness','Phylogenetic \n diversity', 'Functional \n diversity', 'Functional \n dispersion','NDVI','Predator \n diversity','Winter minimum\n temperature','Habitat \n heterogeneity','Topographic \n heterogeneity','Ice free period'))


standdf<-read.table('~/DISENTANGLE/WP3/Arctic/AnalysisJan2018/standardised_analysisdataframe.txt')
standdf<-read.table('PhylogeneticFunctionalAnalysisPicanteR/standardised_analysisdataframe.txt')


modeldf<-standdf[!is.na(standdf$TopographicHet) & !is.na(standdf$Regions),]


#dredge with fixed correlation structure
#SR
fullgls_sr_nocor<- gls(Species.richness~NDVI+WinterMinTemp+
                         HabitatHet+TopographicHet+
                         IceFreeHistory+HerbivorePredators_R + as.factor(Regions),
                       method="ML",
                       data=modeldf)
vario0<- Variogram(fullgls_sr_nocor, form = ~modeldf$x +modeldf$y, resType = "pearson")
plot(vario0,smooth=TRUE)

fullgls_sr<- gls(Species.richness~NDVI+WinterMinTemp+
                   HabitatHet+TopographicHet+
                   IceFreeHistory+HerbivorePredators_R + as.factor(Regions),
                 correlation=corExp(form=~x+y, nugget=T), method="ML",
                 data=modeldf)
vario1 <- Variogram(fullgls_sr, form = ~x + y, resType = "pearson")
plot(vario1, smooth = TRUE, ylim = c(0, 1.2))

fullgls_sr_rat<- gls(Species.richness~NDVI+WinterMinTemp+
                       HabitatHet+TopographicHet+
                       IceFreeHistory+HerbivorePredators_R + as.factor(Regions),
                     correlation=corRatio(form=~x+y, nugget=T), method="ML",
                     data=modeldf)
vario2 <- Variogram(fullgls_sr_rat, form = ~x + y, resType = "pearson")
plot(vario2, smooth = TRUE, ylim = c(0, 1.2))

fullgls_sr_lin<- gls(Species.richness~NDVI+WinterMinTemp+
                       HabitatHet+TopographicHet+
                       IceFreeHistory+HerbivorePredators_R + as.factor(Regions),
                     correlation=corLin(form=~x+y, nugget=T), method="ML",
                     data=modeldf)
vario3 <- Variogram(fullgls_sr_lin, form = ~x + y, resType = "pearson")
plot(vario3, smooth = TRUE, ylim = c(0, 1.2))



globmodgls_sr<-gls(Species.richness~NDVI+WinterMinTemp+
                     HabitatHet+TopographicHet+
                     IceFreeHistory+HerbivorePredators_R + as.factor(Regions),
                   correlation=corExp(form=~x+y,nugget=T,value=coef(fullgls_sr$modelStruct$corStruct,unconstrained=F),fixed=T), method="ML",data=modeldf)
varioSR<-Variogram(globmodgls_sr,resType='normalized')
plot(varioSR,main='Species richness')
write.table(varioSR,'Variogram_SR.txt')

modsetcor_sr<-dredge(globmodgls_sr,trace=2)
modselcor_sr<-model.sel(modsetcor_sr)
modavgcor_sr<-model.avg(modselcor_sr)
importance(modavgcor_sr)
summary(modavgcor_sr)
write.table(importance(modavgcor_sr),'~/DISENTANGLE/WP3/Arctic/AnalysisJan2018/importance_gls_sr.txt')
write.table(summary(modavgcor_sr)$coefmat.full,'~/DISENTANGLE/WP3/Arctic/AnalysisJan2018/modcoef_gls_sr.txt')
write.table(summary(modavgcor_sr)$coefmat.subset,'~/DISENTANGLE/WP3/Arctic/AnalysisJan2018/modcoef_cond_gls_sr.txt')

#PD
fullgls_pd<-gls(Phylogenetic.diversity~NDVI+WinterMinTemp+
                  HabitatHet+TopographicHet+
                  IceFreeHistory+
                  HerbivorePredators_R + as.factor(Regions),correlation=corExp(form=~x+y, nugget=T), method="ML",data=modeldf)

globmodgls_pd<-gls(Phylogenetic.diversity~NDVI+WinterMinTemp+
                     HabitatHet+TopographicHet+
                     IceFreeHistory+HerbivorePredators_R + as.factor(Regions),
                   correlation=corExp(form=~x+y,nugget=T,value=coef(fullgls_pd$modelStruct$corStruct,unconstrained=F),fixed=T), method="ML",data=modeldf)


varioPD<- Variogram(globmodgls_pd,resType = 'normalized')
plot(varioPD,main='Phylogenetic diversity')
write.table(varioPD,'Variogram_PD.txt')

modsetcor_pd<-dredge(globmodgls_pd,trace=2)
modselcor_pd<-model.sel(modsetcor_pd)
modavgcor_pd<-model.avg(modselcor_pd)
importance(modavgcor_pd)
summary(modavgcor_pd)
write.table(importance(modavgcor_pd),'~/DISENTANGLE/WP3/Arctic/AnalysisJan2018/importance_gls_pd.txt')
write.table(summary(modavgcor_pd)$coefmat.full,'~/DISENTANGLE/WP3/Arctic/AnalysisJan2018/modcoef_gls_pd.txt')
write.table(summary(modavgcor_pd)$coefmat.subset,'~/DISENTANGLE/WP3/Arctic/AnalysisJan2018/modcoef_cond_gls_pd.txt')


#FD
fullgls_fd<-gls(Functional.diversity~NDVI+WinterMinTemp+
                  HabitatHet+TopographicHet+
                  IceFreeHistory+
                  HerbivorePredators_R + as.factor(Regions),correlation=corExp(form=~x+y, nugget=T), method="ML",data=modeldf)

globmodgls_fd<-gls(Functional.diversity~NDVI+WinterMinTemp+
                     HabitatHet+TopographicHet+
                     IceFreeHistory+HerbivorePredators_R + as.factor(Regions),
                   correlation=corExp(form=~x+y,nugget=T,value=coef(fullgls_fd$modelStruct$corStruct,unconstrained=F),fixed=T), method="ML",data=modeldf)

varioFD<- Variogram(globmodgls_fd,resType = 'normalized')
plot(varioFD,main='Functional diversity')
write.table(varioFD,'Variogram_FD.txt')

modsetcor_fd<-dredge(globmodgls_fd,trace=2)
modselcor_fd<-model.sel(modsetcor_fd)
modavgcor_fd<-model.avg(modselcor_fd)
importance(modavgcor_fd)
summary(modavgcor_fd)
write.table(importance(modavgcor_fd),'~/DISENTANGLE/WP3/Arctic/AnalysisJan2018/importance_gls_fd.txt')
write.table(summary(modavgcor_fd)$coefmat.full,'~/DISENTANGLE/WP3/Arctic/AnalysisJan2018/modcoef_gls_fd.txt')
write.table(summary(modavgcor_fd)$coefmat.subset,'~/DISENTANGLE/WP3/Arctic/AnalysisJan2018/modcoef_cond_gls_fd.txt')

#FDPD
modeldfFC<-modeldf[!is.na(modeldf$Functional.convergence),]
fullgls_fdpd<-gls(Functional.convergence~NDVI+WinterMinTemp+
                    HabitatHet+TopographicHet+
                    IceFreeHistory+
                    HerbivorePredators_R + as.factor(Regions),correlation=corExp(form=~x+y, nugget=T), method="ML",data=modeldfFC)

globmodgls_fdpd<-gls(Functional.convergence~NDVI+WinterMinTemp+
                       HabitatHet+TopographicHet+
                       IceFreeHistory+HerbivorePredators_R + as.factor(Regions),
                     correlation=corExp(form=~x+y,nugget=T,value=coef(fullgls_fdpd$modelStruct$corStruct,unconstrained=F),fixed=T), method="ML",data=modeldfFC)

varioFDPD<- Variogram(globmodgls_fdpd,resType = 'normalized')
plot(varioFDPD,main='Functional convergence')
write.table(varioFDPD,'Variogram_FDPD.txt')

modsetcor_fdpd<-dredge(globmodgls_fdpd,trace=2)
modselcor_fdpd<-model.sel(modsetcor_fdpd)
modavgcor_fdpd<-model.avg(modselcor_fdpd)
importance(modavgcor_fdpd)
summary(modavgcor_fdpd)
write.table(importance(modavgcor_fdpd),'~/DISENTANGLE/WP3/Arctic/AnalysisJan2018/importance_gls_fdpd_R1.txt')
write.table(summary(modavgcor_fdpd)$coefmat.full,'~/DISENTANGLE/WP3/Arctic/AnalysisJan2018/modcoef_gls_fdpd_R1.txt')
write.table(summary(modavgcor_fdpd)$coefmat.subset,'~/DISENTANGLE/WP3/Arctic/AnalysisJan2018/modcoef_cond_gls_fdpd_R1.txt')


#Variograms plot
tiff('PhylogeneticFunctionalAnalysisPicanteR/Variograms.tif',width=5,height=5,units='in',res=150)
V1<-plot(varioSR,main='Species richness',xlim=c(0,5000000),ylim=c(0,1.3))
V2<-plot(varioPD,main='Phylogenetic diversity',xlim=c(0,5000000),ylim=c(0,1.3))
V3<-plot(varioFD,main='Functional diversity',xlim=c(0,5000000),ylim=c(0,1.3))
V4<-plot(varioFDPD,main='Functional convergence',xlim=c(0,5000000),ylim=c(0,1.3))
print(V1, split=c(1, 1, 2, 2), more=TRUE)
print(V2, split=c(2, 1, 2, 2), more=TRUE)
print(V3, split=c(1, 2, 2, 2), more=TRUE)
print(V4, split=c(2, 2, 2, 2))
dev.off()


#GLS of residuals against species richness ####
#MDf
mdf<-data.frame(cbind(modeldf,fulldf[!is.na(standdf$TopographicHet) & !is.na(standdf$Regions),c(3,5:7)]))
mdf2<-mdf[!is.na(mdf$Functional.convergence),]
mdf2$fcsr<-residuals(lm(mdf2$Functional.convergence.1~mdf2$Species.richness.2))
fullgls_fcsr<-gls(fcsr~NDVI+WinterMinTemp+
                    HabitatHet+TopographicHet+
                    IceFreeHistory+
                    HerbivorePredators_R + as.factor(Regions),correlation=corExp(form=~x+y, nugget=T), method="ML",data=mdf2)

globmodgls_fcsr<-gls(fcsr~NDVI+WinterMinTemp+
                       HabitatHet+TopographicHet+
                       IceFreeHistory+HerbivorePredators_R + as.factor(Regions),
                     correlation=corExp(form=~x+y,nugget=T,value=coef(fullgls_fcsr$modelStruct$corStruct,unconstrained=F),fixed=T), method="ML",data=mdf2)
modsetcor_fcsr<-dredge(globmodgls_fcsr,trace=2,fixed='Species.richness')
modselcor_fcsr<-model.sel(modsetcor_fcsr)
modavgcor_fcsr<-model.avg(modselcor_fcsr)
importance(modavgcor_fcsr)
summary(modavgcor_fcsr)

 #FDSR
mdf2$fdsr1<-residuals(lm(mdf2$Functional.diversity.1~log(mdf2$Species.richness.2)))#Log for FD-SR
fullgls_fdsr<-gls(fdsr1~NDVI+WinterMinTemp+
                    HabitatHet+TopographicHet+
                    IceFreeHistory+
                    HerbivorePredators_R + as.factor(Regions),correlation=corExp(form=~x+y, nugget=T), method="ML",data=mdf2)

globmodgls_fdsr<-gls(fdsr1~NDVI+WinterMinTemp+
                       HabitatHet+TopographicHet+
                       IceFreeHistory+HerbivorePredators_R + as.factor(Regions),
                     correlation=corExp(form=~x+y,nugget=T,value=coef(fullgls_fdsr$modelStruct$corStruct,unconstrained=F),fixed=T), method="ML",data=mdf2)
modsetcor_fdsr<-dredge(globmodgls_fdsr,trace=2)
modselcor_fdsr<-model.sel(modsetcor_fdsr)
modavgcor_fdsr<-model.avg(modselcor_fdsr)
importance(modavgcor_fdsr)
summary(modavgcor_fdsr)

#PDSR
mdf2$pdsr1<-residuals(lm(mdf2$Phylogenetic.diversity.1~mdf2$Species.richness.2))#Linear for PD-SR
fullgls_pdsr<-gls(pdsr1~NDVI+WinterMinTemp+
                    HabitatHet+TopographicHet+
                    IceFreeHistory+
                    HerbivorePredators_R + as.factor(Regions),correlation=corExp(form=~x+y, nugget=T), method="ML",data=mdf2)

globmodgls_pdsr<-gls(pdsr1~NDVI+WinterMinTemp+
                       HabitatHet+TopographicHet+
                       IceFreeHistory+HerbivorePredators_R + as.factor(Regions),
                     correlation=corExp(form=~x+y,nugget=T,value=coef(fullgls_pdsr$modelStruct$corStruct,unconstrained=F),fixed=T), method="ML",data=mdf2)
modsetcor_pdsr<-dredge(globmodgls_pdsr,trace=2)
modselcor_pdsr<-model.sel(modsetcor_pdsr)
modavgcor_pdsr<-model.avg(modselcor_pdsr)
importance(modavgcor_pdsr)
summary(modavgcor_pdsr)


write.table(importance(modavgcor_pdsr),'importanceSR/importance_gls_pdsr.txt')
write.table(summary(modavgcor_pdsr)$coefmat.full,'importanceSR/modcoef_gls_pdsr.txt')
write.table(summary(modavgcor_pdsr)$coefmat.subset,'importanceSR/modcoefcond_gls_pdsr.txt')
write.table(importance(modavgcor_fdsr),'importanceSR/importance_gls_fdsr.txt')
write.table(summary(modavgcor_fdsr)$coefmat.full,'importanceSR/modcoef_gls_fdsr.txt')
write.table(summary(modavgcor_fdsr)$coefmat.subset,'importanceSR/modcoefcond_gls_fdsr.txt')
write.table(importance(modavgcor_fcsr),'importanceSR/importance_gls_fcsr.txt')
write.table(summary(modavgcor_fcsr)$coefmat.full,'importanceSR/modcoef_gls_fcsr.txt')
write.table(summary(modavgcor_fcsr)$coefmat.subset,'importanceSR/modcoefcond_gls_fcsr.txt')

#RVI&MAC PDSR ####
#Diversity
isr<-read.table('~/DISENTANGLE/WP3/Arctic/AnalysisJan2018/importance_gls_sr.txt')
ipdsr<-read.table('importanceSR/importance_gls_pdsr.txt')
ifdsr<-read.table('importanceSR/importance_gls_fdsr.txt')
ifcsr<-read.table('importanceSR/importance_gls_fcsr.txt')
#Give names
colnames(isr)<-'Species'
colnames(ipdsr)<-'PhylogeneticDiv'
colnames(ifdsr)<-'FunctionalDiv'
colnames(ifcsr)<-'FunctionalCon'

impdivsr<-cbind( isr[order(rownames(isr)),], ipdsr[order(rownames(ipdsr)),],ifdsr[order(rownames(ifdsr)),],ifcsr[order(rownames(ifcsr)),])
colnames(impdivsr)<-c('Species','Phylogenetic diversity','Functional diversity','Functional convergence')
rownames(impdivsr)<-rownames(isr)[order(rownames(isr))]
impdivsr
impdivsrsort<-impdivsr[order(-impdivsr[,1]),]
rownames(impdivsrsort)<-(c('NDVI','Predator diversity (R)','Topographic heterogeneity','Habitat heterogeneity','Ice-free history','Winter minimum temperature','Zoogeographic region (F)'))
impdivsortsrx<-impdivsrsort[c(2,1,4,3,6,5,7),]
colsImp<-brewer.pal(5,'Blues')[2:5]
par(mar=c(5,12,1,1))
barplot(t(as.matrix(impdivsortsrx)), horiz=T,las=1,xlab='Relative variable importance',col=colsImp#,col=c('darkred','red','pink4',grey(0.8))
        ,beside=T,legend.text=T,args.legend=list(y=nrow(impdivsortx)-5,x=-0.05,title='Diversity',legend=rev(c('Species','Phylogenetic','Functional','Functional:Phylogenetic'))))


#Model averaged coefficients
#Diversity 
#Full coefficents
coefsr<-read.table('~/DISENTANGLE/WP3/Arctic/AnalysisJan2018/modcoef_cond_gls_sr.txt')
coefpdsr<-read.table('importanceSR/modcoefcond_gls_pdsr.txt')
coeffdsr<-read.table('importanceSR/modcoefcond_gls_fdsr.txt')
coeffcsr<-read.table('importanceSR/modcoefcond_gls_fcsr.txt')
#Reorder rows
coefsr1<-coefsr[order(rownames(coefsr)),]
coefpdsr<-coefpdsr[order(rownames(coefpdsr)),]
coeffdsr<-coeffdsr[order(rownames(coeffdsr)),]
coeffcsr<-coeffcsr[order(rownames(coeffcsr)),]

macplot_est<-rbind(SR=coefsr[,1],PD=coefpdsr[,1],FD=coeffdsr[,1],FDPD=coeffcsr[,1])
colnames(macplot_est)<-rownames(coefsr)
macplot_est<-macplot_est[,c(6,4,7,3,8,9,5,2,3)]
colnames(macplot_est)<-c('Intercept','Predator diversity (R)','NDVI','Habitat heterogeneity','Topographic heterogeneity','Winter minimum temperature','Ice-free history','Region:Eur vs. Arc','Region: NA vs. Arc')
macplot_se<-rbind(SR=coefsr[,2],PD=coefpdsr[,2],FD=coeffdsr[,2],FDPD=coeffcsr[,2])
macplot_se<-macplot_se[,c(6,4,7,3,8,9,5,2,3)]
macplot_p<-rbind(SR=coefsr[,5],PD=coefpdsr[,5],FD=coeffdsr[,5],FDPD=coeffcsr[,5])
macplot_p<-macplot_p[,c(6,4,7,3,8,9,5,2,3)]
macplot_p[macplot_p>=0.05]<-1
macplot_p[macplot_p<0.05]<-16


tiff('PhylogeneticFunctionalAnalysisPicanteR/VarImpModAvgCoef_SR.tif',width = 8,height=5,units='in',res=150,pointsize=8)
par(mfrow=c(1,2))
par(mar=c(5,13,1,1))
#Imp
#bI<-barplot(cbind(t(as.matrix(impdivsortx)),rep(NA,times=4)), horiz=T,las=1,xlab='Relative variable importance',col=colsImp#,col=c('darkred','red','pink4',grey(0.8))
#        ,beside=T,legend.text=T,args.legend=list(cex=0.9,'topr',title='Diversity',legend=rev(c('Species','Phylogenetic','Functional','Functional divergence'))))
bI<-barplot(cbind(t(as.matrix(impdivsortsrx[,2:4])),rep(NA,times=4)), horiz=T,las=1,xlab='Relative variable importance',col=colsImp[2:4]#,col=c('darkred','red','pink4',grey(0.8))
            ,beside=T,legend.text=T,args.legend=list(cex=0.9,'topr',title='Diversity',legend=rev(c('Phylogenetic diversity','Functional diversity','Functional convergence'))),
            names.arg=c('Predator diversity (R)','Vegetation productivity \n (NDVI)','Habitat heterogeneity','Topographic heterogeneity','Climatic severity \n (winter minimum temperature)',
                        'Landscape history \n (time since glaciation)','Zoogeographic region (F)',''))
mtext('(a)',side=3,adj=0)
par(mar=c(5,2,1,1))
par(xpd=T)
bI1<-barplot(cbind(t(as.matrix(impdivsortsrx[,2:4])),rep(NA,times=3)), horiz=T,beside=T,xlim=c(-0.03,0.012),col=F,border=F,
             xlab='Model averaged coefficients',las=1,names.arg=rep(NA,times=8))#,names.arg=colnames(macplot_est[,2:ncol(macplot_est)]))
#names.arg=c(rep(NA,times=6),colnames(macplot_est[,2:ncol(macplot_est)])[7:8]))
points(macplot_est[2:4,2:ncol(macplot_est)],bI1,pch=macplot_p[2:4,2:ncol(macplot_p)],col=colsImp[2:4],lwd=2,cex=1.5) #col=c('black','orange','blue','pink4')) 
arrows(macplot_est[2:4,2:ncol(macplot_est)]+1.96*macplot_se[2:4,2:ncol(macplot_est)],bI1,
       macplot_est[2:4,2:ncol(macplot_est)]-1.96*macplot_se[2:4,2:ncol(macplot_est)],bI1,
       angle=90,length=0.05,code=3,col=colsImp[2:4])#,col=c('black','orange','blue','pink4'))
#legend('topr',pch=c(1,16),col=colsImp[4],c('P>=0.05','P<0.05'),title='Significance',pt.cex=1.5,cex=0.9)
text(-0.022,colMeans(bI1),c(rep(NA,times=6),'Zoogeographic region \n (Eurasian vs. Arctico-Siberian)','Zoogeographic region \n (N. American vs Arctico-Siberian)'),cex=0.8)
axis(side=1)
title(xlab='Model averaged coefficients')
axis(side=2,pos=0,outer=F,lwd.ticks=NA,labels=F,lty=2)
mtext('(b)',side=3,adj=0)
dev.off()

#Variograms PDSR

varioPDSR<- Variogram(globmodgls_pdsr,resType = 'normalized')
varioFDSR<- Variogram(globmodgls_fdsr,resType = 'normalized')
varioFCSR<- Variogram(globmodgls_fcsr,resType = 'normalized')
tiff('PhylogeneticFunctionalAnalysisPicanteR/Variograms_sr.tif',width=7,height=3,units='in',res=150)
Vsr1<-plot(varioPDSR,main='Phylogenetic diversity',xlim=c(0,5000000),ylim=c(0,1.3))
Vsr2<-plot(varioFDSR,main='Functional diversity',xlim=c(0,5000000),ylim=c(0,1.3))
Vsr3<-plot(varioFCSR,main='Functional convergence',xlim=c(0,5000000),ylim=c(0,1.3))
print(Vsr1, split=c(1, 1, 3, 1), more=TRUE)
print(Vsr2, split=c(2, 1, 3, 1), more=TRUE)
print(Vsr3, split=c(3, 1, 3, 1), more=TRUE)
dev.off()


#Effect size - use only cells with SR >1 ####
modeldf_es<-modeldf[!is.na(modeldf$Phylogenetic.diversity_es),]
#PD_ES
fullgls_pdes<-gls(Phylogenetic.diversity_es~NDVI+WinterMinTemp+
                    HabitatHet+TopographicHet+
                    IceFreeHistory
                  +HerbivorePredators_R + as.factor(Regions),correlation=corExp(form=~x+y, nugget=T), method="ML",data=modeldf_es)

globmodgls_pdes<-gls(Phylogenetic.diversity_es~NDVI+WinterMinTemp+
                       HabitatHet+TopographicHet+
                       IceFreeHistory+HerbivorePredators_R + as.factor(Regions),
                     correlation=corExp(form=~x+y,nugget=T,value=coef(fullgls_pdes$modelStruct$corStruct,unconstrained=F),fixed=T), method="ML",data=modeldf_es)

modsetcor_pdes<-dredge(globmodgls_pdes,trace=2)
modselcor_pdes<-model.sel(modsetcor_pdes)
modavgcor_pdes<-model.avg(modselcor_pdes)
importance(modavgcor_pdes)
summary(modavgcor_pdes)
write.table(importance(modavgcor_pdes),'~/DISENTANGLE/WP3/Arctic/AnalysisJan2018/importance_gls_pdes.txt')
write.table(summary(modavgcor_pdes)$coefmat.full,'~/DISENTANGLE/WP3/Arctic/AnalysisJan2018/modcoef_gls_pdes.txt')
write.table(summary(modavgcor_pdes)$coefmat.subset,'~/DISENTANGLE/WP3/Arctic/AnalysisJan2018/modcoef_cond_gls_pdes.txt')
#FD_ES
fullgls_fdes<-gls(Functional.diversity_es~NDVI+WinterMinTemp+
                    HabitatHet+TopographicHet+
                    IceFreeHistory+HerbivorePredators_R + as.factor(Regions),correlation=corExp(form=~x+y, nugget=T), method="ML",data=modeldf_es)

globmodgls_fdes<-gls(Functional.diversity_es~NDVI+WinterMinTemp+
                       HabitatHet+TopographicHet+
                       IceFreeHistory+HerbivorePredators_R + as.factor(Regions),
                     correlation=corExp(form=~x+y,nugget=T,value=coef(fullgls_fdes$modelStruct$corStruct,unconstrained=F),fixed=T), method="ML",data=modeldf_es)

modsetcor_fdes<-dredge(globmodgls_fdes,trace=2)
modselcor_fdes<-model.sel(modsetcor_fdes)
modavgcor_fdes<-model.avg(modselcor_fdes)
importance(modavgcor_fdes)
summary(modavgcor_fdes)
write.table(importance(modavgcor_fdes),'~/DISENTANGLE/WP3/Arctic/AnalysisJan2018/importance_gls_fdes.txt')
write.table(summary(modavgcor_fdes)$coefmat.full,'~/DISENTANGLE/WP3/Arctic/AnalysisJan2018/modcoef_gls_fdes.txt')
write.table(summary(modavgcor_fdes)$coefmat.subset,'~/DISENTANGLE/WP3/Arctic/AnalysisJan2018/modcoef_cond_gls_fdes.txt')

#FDPD_ES
fullgls_fdpdes<-gls(Functional.Phylogenetic.diversity_es~NDVI+WinterMinTemp+
                      HabitatHet+TopographicHet+
                      IceFreeHistory+HerbivorePredators_R + as.factor(Regions),correlation=corExp(form=~x+y, nugget=T), method="ML",data=modeldf_es)

globmodgls_fdpdes<-gls(Functional.Phylogenetic.diversity_es~NDVI+WinterMinTemp+
                         HabitatHet+TopographicHet+
                         IceFreeHistory+HerbivorePredators_R + as.factor(Regions),
                       correlation=corExp(form=~x+y,nugget=T,value=coef(fullgls_fdpdes$modelStruct$corStruct,unconstrained=F),fixed=T), method="ML",data=modeldf_es)

modsetcor_fdpdes<-dredge(globmodgls_fdpdes,trace=2)
modselcor_fdpdes<-model.sel(modsetcor_fdpdes)
modavgcor_fdpdes<-model.avg(modselcor_fdpdes)
importance(modavgcor_fdpdes)
summary(modavgcor_fdpdes)
write.table(importance(modavgcor_fdpdes),'~/DISENTANGLE/WP3/Arctic/AnalysisJan2018/importance_gls_fdpdes.txt')
write.table(summary(modavgcor_fdpdes)$coefmat.full,'~/DISENTANGLE/WP3/Arctic/AnalysisJan2018/modcoef_gls_fdpdes.txt')
write.table(summary(modavgcor_fdpdes)$coefmat.subset,'~/DISENTANGLE/WP3/Arctic/AnalysisJan2018/modcoef_cond_gls_fdpdes.txt')


# RVI and MAC plots (Franklin) -------------------------------------------------------

#Figures

#Importance

#Diversity
isr<-read.table('~/DISENTANGLE/WP3/Arctic/AnalysisJan2018/importance_gls_sr.txt')
ipd<-read.table('~/DISENTANGLE/WP3/Arctic/AnalysisJan2018/importance_gls_pd.txt')
ifd<-read.table('~/DISENTANGLE/WP3/Arctic/AnalysisJan2018/importance_gls_fd.txt')
ifdpd<-read.table('~/DISENTANGLE/WP3/Arctic/AnalysisJan2018/importance_gls_fdpd_R1.txt')
#Give names
colnames(isr)<-'Species'
colnames(ipd)<-'Phylogenetic'
colnames(ifd)<-'Functional'
colnames(ifdpd)<-'Functional:Phylogenetic'

impdiv1<-cbind( isr[order(rownames(isr)),], ipd[order(rownames(ipd)),],ifd[order(rownames(ifd)),],ifdpd[order(rownames(ifdpd)),])
colnames(impdiv1)<-c('Species','Phylogenetic','Functional','Functional:Phylogenetic')
rownames(impdiv1)<-rownames(isr)[order(rownames(isr))]
impdiv1
impdivsort<-impdiv1[order(-impdiv1[,1]),]
rownames(impdivsort)<-(c('NDVI','Predator diversity (R)','Topographic heterogeneity','Habitat heterogeneity','Ice-free history','Winter minimum temperature','Zoogeographic region (F)'))
impdivsortx<-impdivsort[c(2,1,4,3,6,5,7),]
colsImp<-brewer.pal(5,'Blues')[2:5]
par(mar=c(5,12,1,1))
barplot(t(as.matrix(impdivsortx)), horiz=T,las=1,xlab='Relative variable importance',col=colsImp#,col=c('darkred','red','pink4',grey(0.8))
        ,beside=T,legend.text=T,args.legend=list(y=nrow(impdivsortx)-5,x=-0.05,title='Diversity',legend=rev(c('Species','Phylogenetic','Functional','Functional:Phylogenetic'))))


#Model averaged coefficients
#Diversity 
#Full coefficents
fullcoefsr<-read.table('~/DISENTANGLE/WP3/Arctic/AnalysisJan2018/modcoef_gls_sr.txt')
fullcoefpd<-read.table('~/DISENTANGLE/WP3/Arctic/AnalysisJan2018/modcoef_gls_pd.txt')
fullcoeffd<-read.table('~/DISENTANGLE/WP3/Arctic/AnalysisJan2018/modcoef_gls_fd.txt')
fullcoeffdpd<-read.table('~/DISENTANGLE/WP3/Arctic/AnalysisJan2018/modcoef_gls_fdpd_R1.txt')
#Reorder rows
fullcoefsr1<-fullcoefsr[order(rownames(fullcoefsr)),]
fullcoefpd1<-fullcoefpd[order(rownames(fullcoefpd)),]
fullcoeffd1<-fullcoeffd[order(rownames(fullcoeffd)),]
fullcoeffdpd1<-fullcoeffdpd[order(rownames(fullcoeffdpd)),]

macplot_est<-rbind(SR=fullcoefsr1[,1],PD=fullcoefpd1[,1],FD=fullcoeffd1[,1],FDPD=fullcoeffdpd1[,1])
colnames(macplot_est)<-rownames(fullcoefsr1)
macplot_est<-macplot_est[,c(6,4,7,3,8,9,5,2,3)]
colnames(macplot_est)<-c('Intercept','Predator diversity (R)','NDVI','Habitat heterogeneity','Topographic heterogeneity','Winter minimum temperature','Ice-free history','Region:Eur vs. Arc','Region: NA vs. Arc')
macplot_se<-rbind(SR=fullcoefsr1[,2],PD=fullcoefpd1[,2],FD=fullcoeffd1[,2],FDPD=fullcoeffdpd1[,2])
macplot_se<-macplot_se[,c(6,4,7,3,8,9,5,2,3)]
macplot_p<-rbind(SR=fullcoefsr1[,5],PD=fullcoefpd1[,5],FD=fullcoeffd1[,5],FDPD=fullcoeffdpd1[,5])
macplot_p<-macplot_p[,c(6,4,7,3,8,9,5,2,3)]
macplot_p[macplot_p>=0.05]<-1
macplot_p[macplot_p<0.05]<-16


tiff('PhylogeneticFunctionalAnalysisPicanteR/VarImpModAvgCoef.tif',width = 8,height=5,units='in',res=150,pointsize=8)
par(mfrow=c(1,2))
par(mar=c(5,13,1,1))
#Imp
#bI<-barplot(cbind(t(as.matrix(impdivsortx)),rep(NA,times=4)), horiz=T,las=1,xlab='Relative variable importance',col=colsImp#,col=c('darkred','red','pink4',grey(0.8))
#        ,beside=T,legend.text=T,args.legend=list(cex=0.9,'topr',title='Diversity',legend=rev(c('Species','Phylogenetic','Functional','Functional divergence'))))
bI<-barplot(cbind(t(as.matrix(impdivsortx)),rep(NA,times=4)), horiz=T,las=1,xlab='Relative variable importance',col=colsImp#,col=c('darkred','red','pink4',grey(0.8))
            ,beside=T,legend.text=T,args.legend=list(cex=0.9,'topr',title='Diversity',legend=rev(c('Species','Phylogenetic','Functional','Functional convergence'))),
            names.arg=c('Predator diversity (R)','Vegetation productivity \n (NDVI)','Habitat heterogeneity','Topographic heterogeneity','Climatic severity \n (winter minimum temperature)',
                        'Landscape history \n (time since glaciation)','Zoogeographic region (F)',''))
par(mar=c(5,2,1,1))
par(xpd=T)
bI1<-barplot(cbind(t(as.matrix(impdivsortx)),rep(NA,times=4)), horiz=T,beside=T,xlim=c(-0.58,0.5),col=F,border=F,
             xlab='Model averaged coefficients',las=1,names.arg=rep(NA,times=8))#,names.arg=colnames(macplot_est[,2:ncol(macplot_est)]))
#names.arg=c(rep(NA,times=6),colnames(macplot_est[,2:ncol(macplot_est)])[7:8]))
points(macplot_est[,2:ncol(macplot_est)],bI1,pch=macplot_p[,2:ncol(macplot_p)],col=colsImp,lwd=2,cex=1.5) #col=c('black','orange','blue','pink4')) 
arrows(macplot_est[,2:ncol(macplot_est)]+1.96*macplot_se[,2:ncol(macplot_est)],bI1,
       macplot_est[,2:ncol(macplot_est)]-1.96*macplot_se[,2:ncol(macplot_est)],bI1,
       angle=90,length=0.05,code=3,col=colsImp)#,col=c('black','orange','blue','pink4'))
#legend('topr',pch=c(1,16),col=colsImp[4],c('P>=0.05','P<0.05'),title='Significance',pt.cex=1.5,cex=0.9)
text(-0.47,colMeans(bI1),c(rep(NA,times=6),'Zoogeographic region \n (Eurasian vs. Arctico-Siberian)','Zoogeographic region \n (N. American vs Arctico-Siberian)'),cex=0.8)
axis(side=1)
title(xlab='Model averaged coefficients')
axis(side=2,pos=0,outer=F,lwd.ticks=NA,labels=F,lty=2)
mtext('(b)',side=3,adj=0)
dev.off()


#Effect size####
ipdes<-read.table('~/DISENTANGLE/WP3/Arctic/AnalysisJan2018/importance_gls_pdes.txt')
ifdes<-read.table('~/DISENTANGLE/WP3/Arctic/AnalysisJan2018/importance_gls_fdes.txt')
ifdpdes<-read.table('~/DISENTANGLE/WP3/Arctic/AnalysisJan2018/importance_gls_fdpdes.txt')

m1<-merge(t(ipdes),t(ifdes),all=T)
m2<-merge(m1,t(ifdpdes),all=T)

impes<-t(m2)
impes<-impes[order(-rowSums(impes)),]
impes
colnames(impes)<-c('FDPD','FD','PD')
par(mar=c(5,12,3,1))
barplot(t(as.matrix(impes)),col=c('darkred','red','pink4'), horiz=T,las=1,main='Relative variable importance'
        ,beside=T,legend.text=T,args.legend=list(y=nrow(impes)-4,x=-0.05,title='Effect size',legend=c('Phylogenetic','Functional','Functional:Phylogenetic')))


#Diversity
#Full coefficents
fullcoefsr<-read.table('~/DISENTANGLE/WP3/Arctic/AnalysisJan2018/modcoef_gls_sr.txt')
fullcoefpd<-read.table('~/DISENTANGLE/WP3/Arctic/AnalysisJan2018/modcoef_gls_pd.txt')
fullcoeffd<-read.table('~/DISENTANGLE/WP3/Arctic/AnalysisJan2018/modcoef_gls_fd.txt')
fullcoeffdpd<-read.table('~/DISENTANGLE/WP3/Arctic/AnalysisJan2018/modcoef_gls_fdpd.txt')
#Reorder rows
fullcoefsr<-fullcoefsr[c(1,4,3,5,2,7,6,8,9),]
fullcoefpd<-fullcoefpd[c(1,4,3,5,2,7,6,8,9),]
fullcoeffd<-fullcoeffd[c(1,4,3,5,2,7,6,8,9),]
fullcoeffdpd<-fullcoeffdpd[c(1,4,3,5,2,7,6,8,9),]

macplot_est<-rbind(SR=fullcoefsr[,1],PD=fullcoefpd[,1],FD=fullcoeffd[,1],FDPD=fullcoeffdpd[,1])
colnames(macplot_est)<-rownames(fullcoefsr)
colnames(macplot_est)<-c('Intercept','NDVI','Predator diversity (R)','Topographic heterogeneity','Habitat heterogeneity','Winter minimum temperature','Ice-free history','Region:Eur vs. Arc','Region: NA vs. Arc')
macplot_se<-rbind(SR=fullcoefsr[,2],PD=fullcoefpd[,2],FD=fullcoeffd[,2],FDPD=fullcoeffdpd[,2])
macplot_p<-rbind(SR=fullcoefsr[,5],PD=fullcoefpd[,5],FD=fullcoeffd[,5],FDPD=fullcoeffdpd[,5])
macplot_p[macplot_p>=0.05]<-1
macplot_p[macplot_p<0.05]<-16

par(mar=c(5,12,1,1))
b1<-barplot(macplot_est[,2:ncol(macplot_est)],horiz=T,col=F,border=F,las=1,xlim=c(-0.8,0.8)
            ,xlab='Model averaged coefficients',beside=T,space=c(1,10))
points(macplot_est[,2:ncol(macplot_est)],b1,pch=macplot_p[,2:ncol(macplot_p)],col=colsImp,lwd=2) #col=c('black','orange','blue','pink4')) 
arrows(macplot_est[,2:ncol(macplot_est)]+1.96*macplot_se[,2:ncol(macplot_est)],b1,
       macplot_est[,2:ncol(macplot_est)]-1.96*macplot_se[,2:ncol(macplot_est)],b1,
       angle=90,length=0.05,code=3,col=colsImp)#,col=c('black','orange','blue','pink4'))
abline(v=0,lty=2)
legend('topr',pch=16,col=rev(colsImp),legend=rev(c('Species','Phylogenetic','Functional','Functional:Phylogenetic')),title='Diversity',cex=0.7)
legend('topl',pch=c(1,16),col=colsImp[4],c('P>=0.05','P<0.05'),cex=0.7,title='Significance')

tiff('~/DISENTANGLE/WP3/Arctic/AnalysisJan2018/DivModAvgCoefFull.tif',res=100,units='in',width=6,height=6)
par(mar=c(5,12,1,1))
b1<-barplot(macplot_est[,2:ncol(macplot_est)],horiz=T,col=F,border=F,las=1,xlim=c(-0.7,0.9)
            ,xlab='Model averaged coefficients',beside=T,space=c(1,10))
points(macplot_est[,2:ncol(macplot_est)],b1,col=c('black','orange','blue','pink4'),pch=macplot_p[,2:ncol(macplot_p)]) 
arrows(macplot_est[,2:ncol(macplot_est)]+1.96*macplot_se[,2:ncol(macplot_est)],b1,
       macplot_est[,2:ncol(macplot_est)]-1.96*macplot_se[,2:ncol(macplot_est)],b1,
       angle=90,length=0.05,code=3,col=c('black','orange','blue','pink4'))
abline(v=0,lty=2)
legend('topr',pch=16,col=rev(c('black','orange','blue','pink4')),legend=rev(c('Species','Phylogenetic','Functional','Functional:Phylogenetic')),title='Diversity',cex=0.7)
legend('topl',pch=c(1,16),col=1,c('P>=0.05','P<0.05'),cex=0.7,title='Significance')
dev.off()


par(mfrow=c(1,2))
par(mar=c(5,12,1,1))
#Imp
barplot(t(as.matrix(impdiv)), horiz=T,las=1,xlab='Relative variable importance',col=colsImp#,col=c('darkred','red','pink4',grey(0.8))
        ,beside=T,legend.text=T,args.legend=list(y=nrow(impdiv)-5,x=-0.05,title='Diversity',legend=rev(c('Species','Phylogenetic','Functional','Functional:Phylogenetic'))))
#Model avg coef
b1<-barplot(macplot_est[,2:ncol(macplot_est)],horiz=T,col=F,border=F,las=1,xlim=c(-0.8,0.8)
            ,xlab='Model averaged coefficients',beside=T,space=c(1,10))
points(macplot_est[,2:ncol(macplot_est)],b1,pch=macplot_p[,2:ncol(macplot_p)],col=colsImp,lwd=2) #col=c('black','orange','blue','pink4')) 
arrows(macplot_est[,2:ncol(macplot_est)]+1.96*macplot_se[,2:ncol(macplot_est)],b1,
       macplot_est[,2:ncol(macplot_est)]-1.96*macplot_se[,2:ncol(macplot_est)],b1,
       angle=90,length=0.05,code=3,col=colsImp)#,col=c('black','orange','blue','pink4'))
abline(v=0,lty=2)
#legend('topr',pch=16,col=rev(colsImp),legend=rev(c('Species','Phylogenetic','Functional','Functional:Phylogenetic')),title='Diversity',cex=0.7)
legend('topl',pch=c(1,16),col=colsImp[4],c('P>=0.05','P<0.05'),cex=0.7,title='Significance')




#Diversity
#Conditional coefficents
condcoefsr<-read.table('~/DISENTANGLE/WP3/Arctic/AnalysisJan2018/modcoef_cond_gls_sr.txt')
condcoefpd<-read.table('~/DISENTANGLE/WP3/Arctic/AnalysisJan2018/modcoef_cond_gls_pd.txt')
condcoeffd<-read.table('~/DISENTANGLE/WP3/Arctic/AnalysisJan2018/modcoef_cond_gls_fd.txt')
condcoeffdpd<-read.table('~/DISENTANGLE/WP3/Arctic/AnalysisJan2018/modcoef_cond_gls_fdpd.txt')


macplot_est<-rbind(SR=condcoefsr[,1],PD=condcoefpd[,1],FD=condcoeffd[,1],FDPD=condcoeffdpd[,1])
colnames(macplot_est)<-rownames(condcoefsr)
macplot_se<-rbind(SR=condcoefsr[,2],PD=condcoefpd[,2],FD=condcoeffd[,2],FDPD=condcoeffdpd[,2])
macplot_p<-rbind(SR=condcoefsr[,5],PD=condcoefpd[,5],FD=condcoeffd[,5],FDPD=condcoeffdpd[,5])
macplot_p[macplot_p>=0.05]<-1
macplot_p[macplot_p<0.05]<-16

par(mar=c(5,12,1,1))
b1<-barplot(macplot_est[,2:ncol(macplot_est)],horiz=T,col=F,border=F,las=1,xlim=c(-0.8,0.8)
            ,xlab='Model averaged coefficients',beside=T,space=c(1,10))
points(macplot_est[,2:ncol(macplot_est)],b1,col=c('black','orange','blue','pink4'),pch=macplot_p[,2:ncol(macplot_p)]) 
arrows(macplot_est[,2:ncol(macplot_est)]+1.96*macplot_se[,2:ncol(macplot_est)],b1,
       macplot_est[,2:ncol(macplot_est)]-1.96*macplot_se[,2:ncol(macplot_est)],b1,
       angle=90,length=0.05,code=3,col=c('black','orange','blue','pink4'))
abline(v=0,lty=2)
legend('topr',pch=16,col=rev(c('black','orange','blue','pink4')),legend=rev(c('Species','Phylogenetic','Functional','Functional:Phylogenetic')),title='Diversity',cex=0.7)
legend('topl',pch=c(1,16),col=1,c('P>=0.05','P<0.05'),cex=0.7,title='Significance')


#Effect size
#Full coefficents
fullcoefpdes<-read.table('~/DISENTANGLE/WP3/Arctic/AnalysisJan2018/modcoef_gls_pdes.txt')
fullcoeffdes<-read.table('~/DISENTANGLE/WP3/Arctic/AnalysisJan2018/modcoef_gls_fdes.txt')
fullcoeffdpdes<-read.table('~/DISENTANGLE/WP3/Arctic/AnalysisJan2018/modcoef_gls_fdpdes.txt')

macplot_est<-rbind(PD=fullcoefpdes[,1],FD=fullcoeffdes[,1],FDPD=fullcoeffdpdes[,1])
colnames(macplot_est)<-rownames(fullcoefsr)
macplot_se<-rbind(PD=fullcoefpdes[,2],FD=fullcoeffdes[,2],FDPD=fullcoeffdpdes[,2])
macplot_p<-rbind(PD=fullcoefpdes[,5],FD=fullcoeffdes[,5],FDPD=fullcoeffdpdes[,5])
macplot_p[macplot_p>=0.05]<-1
macplot_p[macplot_p<0.05]<-16

par(mar=c(5,12,1,1))
b1<-barplot(macplot_est[,2:ncol(macplot_est)],horiz=T,col=F,border=F,las=1,xlim=c(-0.5,0.5)
            ,xlab='Model averaged coefficients',beside=T,space=c(1,10))
points(macplot_est[,2:ncol(macplot_est)],b1,col=c('orange','blue','pink4'),pch=macplot_p[,2:ncol(macplot_p)]) 
arrows(macplot_est[,2:ncol(macplot_est)]+1.96*macplot_se[,2:ncol(macplot_est)],b1,
       macplot_est[,2:ncol(macplot_est)]-1.96*macplot_se[,2:ncol(macplot_est)],b1,
       angle=90,length=0.05,code=3,col=c('orange','blue','pink4'))
abline(v=0,lty=2)
legend('topr',pch=16,col=rev(c('orange','blue','pink4')),legend=rev(c('Phylogenetic','Functional','Functional:Phylogenetic')),title='Effect size',cex=0.7)
legend('topl',pch=c(1,16),col=1,c('P>=0.05','P<0.05'),cex=0.7,title='Significance')

#Effect size
#Conditional coefficients
condcoefpdes<-read.table('~/DISENTANGLE/WP3/Arctic/AnalysisJan2018/modcoef_cond_gls_pdes.txt')
condcoeffdes<-read.table('~/DISENTANGLE/WP3/Arctic/AnalysisJan2018/modcoef_cond_gls_fdes.txt')
condcoeffdpdes<-read.table('~/DISENTANGLE/WP3/Arctic/AnalysisJan2018/modcoef_cond_gls_fdpdes.txt')

macplot_est<-rbind(PD=condcoefpdes[,1],FD=condcoeffdes[,1],FDPD=condcoeffdpdes[,1])
colnames(macplot_est)<-rownames(condcoefpdes)
macplot_se<-rbind(PD=condcoefpdes[,2],FD=condcoeffdes[,2],FDPD=condcoeffdpdes[,2])
macplot_p<-rbind(PD=condcoefpdes[,5],FD=condcoeffdes[,5],FDPD=condcoeffdpdes[,5])
macplot_p[macplot_p>=0.05]<-1
macplot_p[macplot_p<0.05]<-16

par(mar=c(5,12,1,1))
b1<-barplot(macplot_est[,2:ncol(macplot_est)],horiz=T,col=F,border=F,las=1,xlim=c(-0.5,0.5)
            ,xlab='Model averaged coefficients',beside=T,space=c(1,10))
points(macplot_est[,2:ncol(macplot_est)],b1,col=c('orange','blue','pink4'),pch=macplot_p[,2:ncol(macplot_p)]) 
arrows(macplot_est[,2:ncol(macplot_est)]+1.96*macplot_se[,2:ncol(macplot_est)],b1,
       macplot_est[,2:ncol(macplot_est)]-1.96*macplot_se[,2:ncol(macplot_est)],b1,
       angle=90,length=0.05,code=3,col=c('orange','blue','pink4'))
abline(v=0,lty=2)
legend('topr',pch=16,col=rev(c('orange','blue','pink4')),legend=rev(c('Phylogenetic','Functional','Functional:Phylogenetic')),title='Effect size',cex=0.7)
legend('topl',pch=c(1,16),col=1,c('P>=0.05','P<0.05'),cex=0.7,title='Significance')

# Clustering --------------------------------------------------------------


#PhyloS?rensen
require(ape)
#Range data
spplist<-list.files('RangeMaps',full.names=T)

herbstack<-stack(spplist)
herbstack

#Only all reindeer, and drop sheep
use1<-herbstack[[c(1:69,71:72,75:77)]]
names(use1)[3]<-'Rangifer_tarandus'
envvars1<-stack('EnvVars/Envvars1')
use1<-mask(use1,envvars1$CurrentIce,maskvalue=1)

sr_geo<-sum(use1,na.rm=T)
plot(sr_geo)

commdat<-getValues(use1)
#Replace NA with 0
commdat[is.na(commdat)]<-0

#FD
functree<-read.tree('~/DISENTANGLE/WP3/Arctic/AnalysisJan2018/Final trait classification/arcticherbivorefunctiontree.nwk')
plot(functree)
functree$tip.label[45]<-'Spermophilus_parryii' #Adjusting name to match phylogeny and range map

functree2<-functree
functree2$edge.length<-functree$edge.length/sum(functree$edge.length)

funcdata<-match.phylo.comm(functree2,commdat)

funcdiv<-pd(funcdata$comm,funcdata$phy,include.root=T)
#PD
phylogeny<-read.tree('~/DISENTANGLE/WP3/Arctic/AnalysisJan2018/Final phylogeny/rooted.nwk')
plot(phylogeny)
phylogeny$tip.label

phylogeny$tip.label[which(!phylogeny$tip.label%in%names(use1))]

#Standardise edge lengths
phylo2<-phylogeny
phylo2$edge.length<-phylogeny$edge.length/sum(phylogeny$edge.length)

#Use picante to trim community and phylogenetic data
phydata<-match.phylo.comm(phylo2,commdat)


#PhyloSor

phylosor1<-function (samp, tree) 
{
  if (is.null(tree$edge.length)) {
    stop("Tree has no branch lengths, cannot compute pd")
  }
  if (!is.rooted(tree)) {
    stop("Rooted phylogeny required for phylosor calculation")
  }
  samp <- as.matrix(samp)
  s <- nrow(samp)
  phylodist <- matrix(NA, s, s)
  rownames(phylodist) <- rownames(samp)
  colnames(phylodist) <- rownames(samp)
  samp_comb <- matrix(NA, s * (s - 1)/2, ncol(samp))
  colnames(samp_comb) <- colnames(samp)
  i <- 1
  for (l in 1:(s - 1)) {
    print(paste('l=',l))
    for (k in (l + 1):s) {
      samp_comb[i, ] <- samp[l, ] + samp[k, ]
      i <- i + 1
    }
  }
  pdsamp <- pd(samp, tree)
  pdsamp_comb <- pd(samp_comb, tree)
  i <- 1
  for (l in 1:(s - 1)) {
    pdl <- pdsamp[l, "PD"]
    for (k in (l + 1):s) {
      print(paste('k=',k))
      pdk <- pdsamp[k, "PD"]
      pdcomb <- pdsamp_comb[i, "PD"]
      pdsharedlk <- pdl + pdk - pdcomb
      phylodist[k, l] = 2 * pdsharedlk/(pdl + pdk)
      i <- i + 1
    }
  }
  return(as.dist(phylodist))
}

pd_sor<-phylosor1(phydata$comm[rowSums(phydata$comm)>1,],phydata$phy)#Only for sites with >1sppp
fd_sor<-phylosor1(funcdata$comm[rowSums(funcdata$comm)>1,],funcdata$phy)


#Convert phylogenetic similarity into distances
pd_sor_dist<-1-pd_sor
fd_sor_dist<-1-fd_sor

#Species clustering
require(vegan)
sppdata<-getValues(use1)
sppdata[is.na(sppdata)]<-0
sppdist<-vegdist(sppdata[rowSums(sppdata)>1,],method='bray',binary=T)

#write.csv(as.matrix(sppdist),'~/DISENTANGLE/WP3/Arctic/AnalysisJan2018/SpDistances.csv')
#write.csv(as.matrix(pd_sor_dist),'~/DISENTANGLE/WP3/Arctic/AnalysisJan2018/PdDistances.csv')
#write.csv(as.matrix(fd_sor_dist),'~/DISENTANGLE/WP3/Arctic/AnalysisJan2018/FdDistances.csv')


#To reimport distance matrices
sd2<-data.matrix(read.csv('~/DISENTANGLE/WP3/Arctic/AnalysisJan2018/SpDistances.csv'))
sd_dist<-as.dist(sd2[,2:ncol(sd2)])
pd2<-data.matrix(read.csv('~/DISENTANGLE/WP3/Arctic/AnalysisJan2018/PdDistances.csv'))
pd_dist<-as.dist(pd2[,2:ncol(pd2)])
fd2<-data.matrix(read.csv('~/DISENTANGLE/WP3/Arctic/AnalysisJan2018/FdDistances.csv'))
fd_dist<-as.dist(fd2[,2:ncol(fd2)])


#Cut cluster tree
spclusts<-hclust(sppdist)
#pdclusts<-hclust(pd_sor_dist)
#fdclusts<-hclust(fd_sor_dist)
pdclusts<-hclust(pd_dist)
fdclusts<-hclust(fd_dist)
plot(spclusts)
plot(pdclusts)
plot(fdclusts)

nclust<-c(9,4,12)
spgroups<-cutree(spclusts,nclust[1])
pdgroups<-cutree(pdclusts,nclust[2])
fdgroups<-cutree(fdclusts,nclust[3])

#Coloured dendrograms of site clusters
tiff('~/DISENTANGLE/WP3/Arctic/AnalysisJan2018/ColouredClusterDendrograms.tif',res=300,width=12,height=8,units='in')
par(mfrow=c(1,3))
d_sp<-as.dendrogram(spclusts)
darksdend<-brewer.pal(6,'Dark2')#[c(1,3,2,5,6,4)]#Same order here as for raster maps
spcol<-darksdend[spgroups]
spcol<-spcol[order.dendrogram(d_sp)]
spcol<-factor(spcol,unique(spcol))
plot(color_branches(d_sp,clusters=as.numeric(spcol),col=levels(spcol)),leaflab='none',main='Species Cluster',las=1)

d_pd<-as.dendrogram(pdclusts)
darksdend<-brewer.pal(6,'Dark2')#[c(4,1,3,2,5,6)]
pdcol<-darksdend[pdgroups]
pdcol<-pdcol[order.dendrogram(d_pd)]
pdcol<-factor(pdcol,unique(pdcol))
plot(color_branches(d_pd,clusters=as.numeric(pdcol),col=levels(pdcol)),leaflab='none',main='Phylogenetic Cluster',las=1)

d_fd<-as.dendrogram(fdclusts)
darksdend<-brewer.pal(6,'Dark2')#[c(4,1,5,3,4,2)]
fdcol<-darksdend[fdgroups]
fdcol<-fdcol[order.dendrogram(d_fd)]
fdcol<-factor(fdcol,unique(fdcol))
plot(color_branches(d_fd,clusters=as.numeric(fdcol),col=levels(fdcol)),leaflab='none',main='Functional Cluster',las=1)
dev.off()

#Set up vector to populate
r1<-raster(use1)
spprich<-rowSums(phydata$comm)
spprich[spprich<2]<-NA

spgroupvalues<-spprich
pdgroupvalues<-spprich
fdgroupvalues<-spprich
spgroupvalues[!is.na(spgroupvalues)]<-spgroups
pdgroupvalues[!is.na(pdgroupvalues)]<-pdgroups
fdgroupvalues[!is.na(fdgroupvalues)]<-fdgroups

sp_sorras<-setValues(r1,spgroupvalues)
pd_sorras<-setValues(r1,pdgroupvalues)
fd_sorras<-setValues(r1,fdgroupvalues)

levelplot(stack(sp_sorras,pd_sorras,fd_sorras),par.settings=viridisTheme,names.attr=c('SpeciesCluster','PhyloCluster','FunctionalCluster'),scales=(list(draw=F)))

writeRaster(sp_sorras,'~/DISENTANGLE/WP3/Arctic/AnalysisJan2018/sp_sorras')
writeRaster(pd_sorras,'~/DISENTANGLE/WP3/Arctic/AnalysisJan2018/pd_sorras')
writeRaster(fd_sorras,'~/DISENTANGLE/WP3/Arctic/AnalysisJan2018/fd_sorras')

#Ratify for factor plot
sp_sorras_rat<-ratify(sp_sorras)
pd_sorras_rat<-ratify(pd_sorras)
fd_sorras_rat<-ratify(fd_sorras)
sprat<-levels(sp_sorras_rat)[[1]]
pdrat<-levels(pd_sorras_rat)[[1]]
fdrat<-levels(fd_sorras_rat)[[1]]
sprat$cluster<-c('PA','Greenland','NA','Bylot','Sval','HighArc')
pdrat$cluster<-c('Sval','PA','Greenland','Bylot','Sval','HighArc')
fdrat$cluster<-c('PA','Greenland','NA','Bylot','Sval','HighArc')

levels(sp_sorras_rat) <- sprat
levels(pd_sorras_rat) <- pdrat
levels(fd_sorras_rat) <- fdrat

#Relevel for consequent cluster maps
#Same colour ordering as for dendrograms
spfac<-sp_sorras
spfac[spfac==1]<-60#PA
spfac[spfac==2]<-10#Greenland
spfac[spfac==3]<-30#NA
spfac[spfac==4]<-50#Bylot
spfac[spfac==5]<-40#Svalbard
spfac[spfac==6]<-20#Higharc

fdfac<-fd_sorras
fdfac[fdfac==1]<-40
fdfac[fdfac==2]<-10
fdfac[fdfac==3]<-50
fdfac[fdfac==4]<-30
fdfac[fdfac==5]<-60
fdfac[fdfac==6]<-20

pdfac<-pd_sorras
pdfac[pdfac==1]<-40
pdfac[pdfac==2]<-10
pdfac[pdfac==3]<-30
pdfac[pdfac==4]<-20
pdfac[pdfac==5]<-50
pdfac[pdfac==6]<-60

#mycol=rasterTheme(region=c('green','blue','purple','pink4','red','yellow'))

#Basic plottings
#Country outlines and north pole
#Lambert azimunthual equal area
laea<-'+proj=laea +lat_0=90 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0 '
charcount <- c('NO', 'SE', 'FI','RU','CA','US','IS','GL','SJ') 
allac2 <- do.call("bind", lapply(charcount, function(x)  raster::getData('GADM', country=x, level=0)))
allarc<-spTransform(allac2,CRS=crs(laea))
np<-SpatialPoints(cbind(0,0))


writeRaster(stack(spfac,pdfac,fdfac),'~/DISENTANGLE/WP3/Arctic/AnalysisJan2018/cluster6rasters',by.layer=T)

darks<-rasterTheme(region=brewer.pal(6,'Dark2'))
levelplot(stack(spfac,pdfac,fdfac),par.settings=darks,scales=list(draw=F),colorkey=F,names.attr=c('Species clusters','Phylogenetic clusters','Functional clusters'))
tiff('~/DISENTANGLE/WP3/Arctic/AnalysisJan2018/clustermap.tif',width = 8,height=3,units='in',res=300)
levelplot(stack(spfac,pdfac,fdfac),par.settings=darks,scales=list(draw=F),colorkey=F,names.attr=c('Species clusters','Phylogenetic clusters','Functional clusters'))+
  layer(sp.polygons(allarc,lwd=0.5))+
  layer(sp.points(np,col=1))
dev.off()


spfac_rat<-ratify(spfac)
pdfac_rat<-ratify(pdfac)
fdfac_rat<-ratify(fdfac)
spfacrat<-levels(spfac_rat)[[1]]
pdfacrat<-levels(pdfac_rat)[[1]]
fdfacrat<-levels(fdfac_rat)[[1]]
levels(spfac_rat) <- spfacrat
levels(pdfac_rat) <- pdfacrat
levels(fdfac_rat) <- fdfacrat
facstack<-stack(spfac_rat,pdfac_rat,pdfac_rat)


# Optimal number of clusters ----------------------------------------------


#Testing significance of clusters
require(NbClust)
indexchoice<-'cindex'#dunn
methodchoice<-'complete'
nbSpp<-NbClust(diss=sppdist,distance=NULL,method=methodchoice,min.nc=2,max.nc=12,index=indexchoice)
nbPD<-NbClust(diss=pd_sor_dist,distance=NULL,method=methodchoice,min.nc=2,max.nc=12,index=indexchoice)
nbFD<-NbClust(diss=fd_sor_dist,distance=NULL,method=methodchoice,min.nc=2,max.nc=12,index=indexchoice)

nbSpp
nbPD
nbFD
clusters<-c(nbSpp$Best.nc[1],nbPD$Best.nc[1],nbFD$Best.nc[1])

spclusts<-hclust(sppdist,method=methodchoice)
pdclusts<-hclust(pd_sor_dist,method=methodchoice)
fdclusts<-hclust(fd_sor_dist,method=methodchoice)
spgroups<-cutree(spclusts,nbSpp$Best.nc[1])
pdgroups<-cutree(pdclusts,nbPD$Best.nc[1])
fdgroups<-cutree(fdclusts,nbFD$Best.nc[1])

#Set up vector to populate
r1<-raster(use1)
spprich<-rowSums(phydata$comm)
spprich[spprich<2]<-NA

spgroupvalues<-spprich
pdgroupvalues<-spprich
fdgroupvalues<-spprich
spgroupvalues[!is.na(spgroupvalues)]<-spgroups
pdgroupvalues[!is.na(pdgroupvalues)]<-pdgroups
fdgroupvalues[!is.na(fdgroupvalues)]<-fdgroups

sp_sorras<-setValues(r1,spgroupvalues)
pd_sorras<-setValues(r1,pdgroupvalues)
fd_sorras<-setValues(r1,fdgroupvalues)

clusters<-rep(8,times=3)
clustpal<-rasterTheme(region=(brewer.pal(max(clusters),'Dark2')))

writeRaster(stack(sp_sorras,pd_sorras,fd_sorras),'~/DISENTANGLE/WP3/Arctic/AnalysisJan2018/cluster8rasters',by.layer=T)

tiff('~/DISENTANGLE/WP3/Arctic/AnalysisJan2018/clustermap8.tif',width = 8,height=3,units='in',res=300)
levelplot(stack(sp_sorras,pd_sorras,fd_sorras),par.settings=clustpal,colorkey=F,names.attr=c('SpeciesCluster','PhyloCluster','FunctionalCluster'),scales=(list(draw=F)))+
  layer(sp.polygons(allarc,lwd=0.5))+
  layer(sp.points(np,col=1))
dev.off()


# Regional cluster dendrograms --------------------------------------------


tiff('~/DISENTANGLE/WP3/Arctic/AnalysisJan2018/clusterdendros8.tif',width = 8,height=5,units='in',res=150)
par(mfrow=c(1,3))
d_sp<-as.dendrogram(spclusts)
darksdend<-(brewer.pal(max(clusters),'Dark2'))[1:clusters[1]]
spcol<-darksdend[spgroups]
spcol<-spcol[order.dendrogram(d_sp)]
spcol<-factor(spcol,unique(spcol))
plot(color_branches(d_sp,clusters=as.numeric(spcol),col=levels(spcol)),leaflab='none',main='Species Cluster',las=1)

d_pd<-as.dendrogram(pdclusts)
darksdend<-brewer.pal(max(clusters),'Dark2')
pdcol<-darksdend[pdgroups]
pdcol<-pdcol[order.dendrogram(d_pd)]
pdcol<-factor(pdcol,unique(pdcol))
plot(color_branches(d_pd,clusters=as.numeric(pdcol),col=levels(pdcol)),leaflab='none',main='Phylogenetic Cluster',las=1)

d_fd<-as.dendrogram(fdclusts)
darksdend<-brewer.pal(max(clusters),'Dark2')
fdcol<-darksdend[fdgroups]
fdcol<-fdcol[order.dendrogram(d_fd)]
fdcol<-factor(fdcol,unique(fdcol))
plot(color_branches(d_fd,clusters=as.numeric(fdcol),col=levels(fdcol)),leaflab='none',main='Functional Cluster',las=1)
dev.off()


# GLS of diversity ratios between zoo regions -----------------------------
divratdat<-read.csv('~/DISENTANGLE/WP3/Arctic/AnalysisJan2018/DiversitywithCoords.csv')
divratdat1<-divratdat[!is.na(divratdat$Regions) &!is.na(divratdat$Species.richness),]

glspdsr<- gls(Phylogenetic.diversity~Species.richness * as.factor(Regions),
              correlation=corExp(form=~x+y, nugget=T), method="ML",
              data=divratdat1)
glsfdsr<- gls(Functional.diversity~log(Species.richness) * as.factor(Regions),
              correlation=corExp(form=~x+y, nugget=T), method="ML",
              data=divratdat1)
glsfdpd<- gls(Functional.diversity~log(Phylogenetic.diversity) * as.factor(Regions),
              correlation=corExp(form=~x+y, nugget=T), method="ML",
              data=divratdat1)

glsfdispsr<- gls(Functional.Phylognetic.diversity~Species.richness * as.factor(Regions),
                 correlation=corExp(form=~x+y, nugget=T), method="ML",
                 data=divratdat1)

globmod_pdsr<-gls(Phylogenetic.diversity~Species.richness * as.factor(Regions),
                  correlation=corExp(form=~x+y,nugget=T,value=coef(glspdsr$modelStruct$corStruct,unconstrained=F),fixed=T), 
                  method="ML",data=divratdat1)
modset_pdsr<-dredge(globmod_pdsr)
modsel_pdsr<-model.sel(modset_pdsr)
modavg_pdsr<-model.avg(modsel_pdsr)

globmod_fdsr<-gls(Functional.diversity~log(Species.richness) * as.factor(Regions),
                  correlation=corExp(form=~x+y,nugget=T,value=coef(glsfdsr$modelStruct$corStruct,unconstrained=F),fixed=T), 
                  method="ML",data=divratdat1)
modset_fdsr<-dredge(globmod_fdsr)
modsel_fdsr<-model.sel(modset_fdsr)
modavg_fdsr<-model.avg(modsel_fdsr)

globmod_fdpd<-gls(Functional.diversity~log(Phylogenetic.diversity) * as.factor(Regions),
                  correlation=corExp(form=~x+y,nugget=T,value=coef(glsfdpd$modelStruct$corStruct,unconstrained=F),fixed=T), method="ML",data=divratdat1)
modset_fdpd<-dredge(globmod_fdpd)
modsel_fdpd<-model.sel(modset_fdpd)
modavg_fdpd<-model.avg(modsel_fdpd)



globmod_fdispsr<-gls(Functional.Phylognetic.diversity~Species.richness * as.factor(Regions),
                     correlation=corExp(form=~x+y,nugget=T,value=coef(glsfdispsr$modelStruct$corStruct,unconstrained=F),fixed=T), method="ML",data=divratdat1)
modset_fdispsr<-dredge(globmod_fdispsr)
modsel_fdispsr<-model.sel(modset_fdispsr)
modavg_fdispsr<-model.avg(modsel_fdispsr)

modavg_pdsr$importance
modavg_fdsr$importance
modavg_fdpd$importance
modavg_fdispsr$importance

write.table(summary(modavg_pdsr)$coefmat.full,'~/DISENTANGLE/WP3/Arctic/AnalysisJan2018/MAC_pdsrrat.csv')
write.table(summary(modavg_fdsr)$coefmat.full,'~/DISENTANGLE/WP3/Arctic/AnalysisJan2018/MAC_fdsrrat.csv')
write.table(summary(modavg_fdpd)$coefmat.full,'~/DISENTANGLE/WP3/Arctic/AnalysisJan2018/MAC_fdpdrat.csv')
write.table(summary(modavg_fdispsr)$coefmat.full,'~/DISENTANGLE/WP3/Arctic/AnalysisJan2018/MAC_fdispsrrat.csv')


macrat_mean<-cbind(PDSR=summary(modavg_pdsr)$coefmat.full[,1],FDSR=summary(modavg_fdsr)$coefmat.full[,1],FDPD=summary(modavg_fdpd)$coefmat.full[,1])
macrat_se<-cbind(PDSR=summary(modavg_pdsr)$coefmat.full[,2],FDSR=summary(modavg_fdsr)$coefmat.full[,2],FDPD=summary(modavg_fdpd)$coefmat.full[,2])

par(mfrow=c(1,2))
par(mar=c(3,10,3,1))
barplot(t(macrat_mean[2:6,]),horiz=T,beside=T,names.arg=c('Region: Eur vs. ArcSib','Region:NorAm vs. ArcSib','Species richness','Species richness:Region Eur vs ArcSib','Species richness: Region NorAm vs. ArcSib'),las=1)

#Predictions
lm1<-with(divratdat1,lm(Phylogenetic.diversity~Species.richness*as.factor(Regions),na.action=na.fail))
d1<-dredge(lm1)
modavg1<-model.avg(d1)

modavg_pdsr<-model.avg(modsel_pdsr,fit=T)
newdf<-expand.grid(Species.richness=seq(0,0.45,length.out = 100),Regions=c(3,6,12))
newdf$p1<-predict(modavg_pdsr,newdf)
newdf$pL<-predict(lm1,newdf)

with(divratdat1,plot(Species.richness,Phylogenetic.diversity,col=Regions))
with(newdf[newdf$Regions==3,],lines(Species.richness,p1,col=3))
with(newdf[newdf$Regions==6,],lines(Species.richness,p1,col=6))
with(newdf[newdf$Regions==12,],lines(Species.richness,p1,col=12))

with(divratdat1,plot(Species.richness,Phylogenetic.diversity,col=Regions))
with(newdf[newdf$Regions==3,],lines(Species.richness,pL,col=3))
with(newdf[newdf$Regions==6,],lines(Species.richness,pL,col=6))
with(newdf[newdf$Regions==12,],lines(Species.richness,pL,col=12))

#Figures

# GLS importance Figures --------------------------------------------------


# RVI and MAC plots (local) -----------------------------------------------
#Importance

#Diversity
isr<-read.table('S:/DISENTANGLE/WP3/ArcticHerbivoreFDPD/PhylogeneticFunctionalAnalysisPicanteR/importanceApr18/importance_gls_sr.txt')
ipd<-read.table('S:/DISENTANGLE/WP3/ArcticHerbivoreFDPD/PhylogeneticFunctionalAnalysisPicanteR/importanceApr18/importance_gls_pd.txt')
ifd<-read.table('S:/DISENTANGLE/WP3/ArcticHerbivoreFDPD/PhylogeneticFunctionalAnalysisPicanteR/importanceApr18/importance_gls_fd.txt')
ifdpd<-read.table('S:/DISENTANGLE/WP3/ArcticHerbivoreFDPD/PhylogeneticFunctionalAnalysisPicanteR/importanceApr18/importance_gls_fdpd.txt')
#Give names
colnames(isr)<-'Species'
colnames(ipd)<-'Phylogenetic'
colnames(ifd)<-'Functional'
colnames(ifdpd)<-'Functional:Phylogenetic'

impdiv1<-cbind( isr[order(rownames(isr)),], ipd[order(rownames(ipd)),],ifd[order(rownames(ifd)),],ifdpd[order(rownames(ifdpd)),])
colnames(impdiv1)<-c('Species','Phylogenetic','Functional','Functional:Phylogenetic')
rownames(impdiv1)<-rownames(isr)[order(rownames(isr))]
impdiv1
impdivsort<-impdiv1[order(-impdiv1[,1]),]
rownames(impdivsort)<-(c('NDVI','Predator diversity (R)','Topographic heterogeneity','Habitat heterogeneity','Ice-free history','Winter minimum temperature','Zoogeographic region (F)'))
impdivsortx<-impdivsort[c(2,1,4,3,6,5,7),]
colsImp<-brewer.pal(5,'Blues')[2:5]
par(mar=c(5,12,1,1))
barplot(t(as.matrix(impdivsortx)), horiz=T,las=1,xlab='Relative variable importance',col=colsImp#,col=c('darkred','red','pink4',grey(0.8))
        ,beside=T,legend.text=T,args.legend=list(y=nrow(impdivsortx)-5,x=-0.05,title='Diversity',legend=rev(c('Species','Phylogenetic','Functional','Functional:Phylogenetic'))))

#Model averaged coefficients
#Diversity 
#Full coefficents
fullcoefsr<-read.table('S:/DISENTANGLE/WP3/ArcticHerbivoreFDPD/PhylogeneticFunctionalAnalysisPicanteR/importanceApr18/modcoef_gls_sr.txt')
fullcoefpd<-read.table('S:/DISENTANGLE/WP3/ArcticHerbivoreFDPD/PhylogeneticFunctionalAnalysisPicanteR/importanceApr18/modcoef_gls_pd.txt')
fullcoeffd<-read.table('S:/DISENTANGLE/WP3/ArcticHerbivoreFDPD/PhylogeneticFunctionalAnalysisPicanteR/importanceApr18/modcoef_gls_fd.txt')
fullcoeffdpd<-read.table('S:/DISENTANGLE/WP3/ArcticHerbivoreFDPD/PhylogeneticFunctionalAnalysisPicanteR/importanceApr18/modcoef_gls_fdpd.txt')
#Reorder rows
fullcoefsr1<-fullcoefsr[order(rownames(fullcoefsr)),]
fullcoefpd1<-fullcoefpd[order(rownames(fullcoefpd)),]
fullcoeffd1<-fullcoeffd[order(rownames(fullcoeffd)),]
fullcoeffdpd1<-fullcoeffdpd[order(rownames(fullcoeffdpd)),]

macplot_est<-rbind(SR=fullcoefsr1[,1],PD=fullcoefpd1[,1],FD=fullcoeffd1[,1],FDPD=fullcoeffdpd1[,1])
colnames(macplot_est)<-rownames(fullcoefsr1)
macplot_est<-macplot_est[,c(1,5,7,2,8,9,6,2,3)]
colnames(macplot_est)<-c('Intercept','Predator diversity (R)','NDVI','Habitat heterogeneity','Topographic heterogeneity','Winter minimum temperature','Ice-free history','Region:Eur vs. Arc','Region: NA vs. Arc')
macplot_se<-rbind(SR=fullcoefsr1[,2],PD=fullcoefpd1[,2],FD=fullcoeffd1[,2],FDPD=fullcoeffdpd1[,2])
macplot_se<-macplot_se[,c(1,5,7,2,8,9,6,2,3)]
macplot_p<-rbind(SR=fullcoefsr1[,5],PD=fullcoefpd1[,5],FD=fullcoeffd1[,5],FDPD=fullcoeffdpd1[,5])
macplot_p<-macplot_p[,c(1,5,7,2,8,9,6,2,3)]
macplot_p[macplot_p>=0.05]<-1
macplot_p[macplot_p<0.05]<-16

par(mar=c(5,12,1,1))
b1<-barplot(macplot_est[,2:ncol(macplot_est)],horiz=T,col=F,border=F,las=1,xlim=c(-0.8,0.8)
            ,xlab='Model averaged coefficients',beside=T,space=c(1,10))
points(macplot_est[,2:ncol(macplot_est)],b1,pch=macplot_p[,2:ncol(macplot_p)],col=colsImp,lwd=2) #col=c('black','orange','blue','pink4')) 
arrows(macplot_est[,2:ncol(macplot_est)]+1.96*macplot_se[,2:ncol(macplot_est)],b1,
       macplot_est[,2:ncol(macplot_est)]-1.96*macplot_se[,2:ncol(macplot_est)],b1,
       angle=90,length=0.05,code=3,col=colsImp)#,col=c('black','orange','blue','pink4'))
abline(v=0,lty=2)
legend('topr',pch=16,col=rev(colsImp),legend=rev(c('Species','Phylogenetic','Functional','Functional:Phylogenetic')),title='Diversity',cex=0.7)
legend('topl',pch=c(1,16),col=colsImp[4],c('P>=0.05','P<0.05'),cex=0.7,title='Significance')


x11(12,6)
tiff('S:/DISENTANGLE/WP3/ArcticHerbivoreFDPD/PhylogeneticFunctionalAnalysisPicanteR/VarImpModAvgCoef_apr18.tif',width = 8,height=5,units='in',res=150,pointsize=8)
par(mfrow=c(1,2))
par(mar=c(5,13,1,1))
#Imp
#bI<-barplot(cbind(t(as.matrix(impdivsortx)),rep(NA,times=4)), horiz=T,las=1,xlab='Relative variable importance',col=colsImp#,col=c('darkred','red','pink4',grey(0.8))
#        ,beside=T,legend.text=T,args.legend=list(cex=0.9,'topr',title='Diversity',legend=rev(c('Species','Phylogenetic','Functional','Functional divergence'))))
bI<-barplot(cbind(t(as.matrix(impdivsortx)),rep(NA,times=4)), horiz=T,las=1,xlab='Relative variable importance',col=colsImp#,col=c('darkred','red','pink4',grey(0.8))
            ,beside=T,legend.text=T,args.legend=list(cex=0.9,'topr',title='Diversity',legend=rev(c('Species','Phylogenetic','Functional','Functional divergence'))),
            names.arg=c('Predator diversity (R)','Vegetation productivity \n (NDVI)','Habitat heterogeneity','Topographic heterogeneity','Climatic severity \n (winter minimum temperature)',
                        'Landscape history \n (time since glaciation)','Zoogeographic region (F)',''))
mtext('(a)',side=3,adj=0)
par(mar=c(5,2,1,1))
par(xpd=T)
bI1<-barplot(cbind(t(as.matrix(impdivsortx)),rep(NA,times=4)), horiz=T,beside=T,xlim=c(-0.58,0.5),col=F,border=F,
             xlab='Model averaged coefficients',las=1,names.arg=rep(NA,times=8))#,names.arg=colnames(macplot_est[,2:ncol(macplot_est)]))
             #names.arg=c(rep(NA,times=6),colnames(macplot_est[,2:ncol(macplot_est)])[7:8]))
points(macplot_est[,2:ncol(macplot_est)],bI1,pch=macplot_p[,2:ncol(macplot_p)],col=colsImp,lwd=2,cex=1.5) #col=c('black','orange','blue','pink4')) 
arrows(macplot_est[,2:ncol(macplot_est)]+1.96*macplot_se[,2:ncol(macplot_est)],bI1,
       macplot_est[,2:ncol(macplot_est)]-1.96*macplot_se[,2:ncol(macplot_est)],bI1,
       angle=90,length=0.05,code=3,col=colsImp)#,col=c('black','orange','blue','pink4'))
#legend('topr',pch=c(1,16),col=colsImp[4],c('P>=0.05','P<0.05'),title='Significance',pt.cex=1.5,cex=0.9)
text(-0.47,colMeans(bI1),c(rep(NA,times=6),'Zoogeographic region \n (Eurasian vs. Arctico-Siberian)','Zoogeographic region \n (N. American vs Arctico-Siberian)'),cex=0.8)
axis(side=1)
title(xlab='Model averaged coefficients')
axis(side=2,pos=0,outer=F,lwd.ticks=NA,labels=F,lty=2)
mtext('(b)',side=3,adj=0)
dev.off()


# Other indices -----------------------------------------------------------

#Phylo Sorensen


physor_pd<-phylosor(phydata$comm[sample(nrow(phydata$comm),100),],phydata$phy)
physor_pd[is.na(physor_pd)]<-0
a1<-hclust(physor_pd)
write.table(as.matrix(physor_pd),'S:/DISENTANGLE/WP3/ArcticHerbivoreFDPD/PhylogeneticFunctionalAnalysisPicanteR/PhySortest.txt')
clusts<-cutree(a1,k=5)
a2<-read.table('S:/DISENTANGLE/WP3/ArcticHerbivoreFDPD/PhylogeneticFunctionalAnalysisPicanteR/PhySortest.txt')

r1<-raster(use1)
pdpsv<-psv(phydata$comm,phydata$phy)
pdpdsras<-mask(setValues(r1,pdpsv$PSVs),sr_geo,maskvalue=0)

fdpsv<-psv(funcdata$comm,funcdata$phy)
fdpdsras<-mask(setValues(r1,fdpsv$PSVs),sr_geo,maskvalue=0)

psvstack<-stack(pdpdsras,fdpdsras)
names(psvstack)<-c('Phylogenetic species variability','Functional species variability')
levelplot(psvstack,margin=F,par.settings=YlOrRdTheme,scales=list(draw=FALSE))

levelplot(log(pdpdsras/fdpdsras+1),margin=F)

mpdpd<-mpd(phydata$comm,cophenetic(phydata$phy), abundance.weighted=F)
mpdras<-mask(setValues(r1,mpdpd),sr_geo,maskvalue=0)
levelplot(mpdras,margin=F)
mpdfd<-mpd(funcdata$comm,cophenetic(funcdata$phy), abundance.weighted=F)
mfdras<-mask(setValues(r1,mpdfd),sr_geo,maskvalue=0)
levelplot(mfdras,margin=F)
mpdstack<-stack((mpdras-cellStats(mpdras,stat='mean')),(mfdras-cellStats(mfdras,stat='mean')))


p<-levelplot(mpdstack,margin=F,names.attr=c('MPD','MFD'))+
  layer(sp.points(np,col=1))+
  layer(sp.polygons(allarc,col='grey',lwd=0.5))

diverge0(p,'RdBu')

#Regressions
envdata<-data.frame(getValues(extend(envvars,use1)))
require(brglm)
spr1<-sppregs(phydata$comm,envdata,phydata$phy,fam='binomial')

#Randomisation by richness
fd_es_rich<-ses.pd(funcdata$comm,funcdata$phy,null.model='richness',runs=1000)
pd_es_rich<-ses.pd(phydata$comm,phydata$phy,null.model='richness',runs=1000)

write.table(fd_es_rich,'S:/DISENTANGLE/WP3/ArcticHerbivoreFDPD/PhylogeneticFunctionalAnalysisPicanteR/fdes_rich.txt')
write.table(pd_es_rich,'S:/DISENTANGLE/WP3/ArcticHerbivoreFDPD/PhylogeneticFunctionalAnalysisPicanteR/pdes_rich.txt')


fdesrichras<-setValues(r1,fd_es_rich$pd.obs.p)
pdesrichras<-setValues(r1,pd_es_rich$pd.obs.p)
randrichstack<-stack(pdesrichras,fdesrichras)
levelplot(randrichstack)


#Phylostructure
ps1<-phylostruct(phydata$comm,phydata$phy,metric='psv')

# Cluster maps ------------------------------------------------------------
cluster6<-stack('S:/DISENTANGLE/WP3/ArcticHerbivoreFDPD/PhylogeneticFunctionalAnalysisPicanteR/cluster6rasters')
cluster6<-mask(cluster6,envvars1$CurrentIce,maskvalue=1)

cluster8<-stack('S:/DISENTANGLE/WP3/ArcticHerbivoreFDPD/PhylogeneticFunctionalAnalysisPicanteR/cluster8rasters_apr18')

darks<-rasterTheme(region=brewer.pal(8,'Dark2'))
levelplot(cluster8,par.settings=darks,scales=list(draw=F),colorkey=F,names.attr=c('Species clusters','Phylogenetic clusters','Functional clusters'))
tiff('S:/DISENTANGLE/WP3/ArcticHerbivoreFDPD/PhylogeneticFunctionalAnalysisPicanteR/clustermap.tif',width = 8,height=3,units='in',res=300)
levelplot(cluster8,par.settings=darks,scales=list(draw=F),colorkey=F,names.attr=c('Species clusters','Phylogenetic clusters','Functional clusters'))+
  layer(sp.polygons(allarc,lwd=0.5,col=grey(0.5)))+
  layer(sp.points(np,col=1))
dev.off()

levelplot(cluster6,par.settings=darks,scales=list(draw=F),colorkey=F,names.attr=c('Species clusters','Phylogenetic clusters','Functional clusters'))+
  layer(sp.polygons(allarc,lwd=0.5,col=grey(0.5)))+
  layer(sp.points(np,col=1))+
  layer(sp.polygons(rasterToPolygons(realms,dissolve=T),col='red',lwd=0.1))




# Mammals and birds -------------------------------------------------------
birdsdat<-use1[[c(5:11,15:19,21:26,35:37)]]#21 bird species
mammalsdat<-use1[[c(1:4,12:14,20,27:34,38:74)]]#53 mammal species
birdsr<-sum(birdsdat,na.rm=T)
mammalsr<-sum(mammalsdat,na.rm=T)

birdmam<-stack(birdsr,mammalsr)
names(birdmam)<-c('Birds','Mammals')
birdmam<-mask(birdmam,divstack2$Species.richness,maskvalue=NA,updatevalue=NA)
levelplot(birdmam,par.settings=YlOrRdTheme,scales=list(draw=FALSE))+
  layer(sp.points(np,col=1))#+
  layer(sp.polygons(allarc,col='grey',lwd=0.5))

birdsAndmammals<-birdsr
birdsAndmammals<-setValues(birdsAndmammals,rep(0,times=ncell(birdsr)))
birdsAndmammals[birdmam$Birds>0 & birdmam$Mammals>0]<-1
plot(birdsAndmammals)
cellStats(birdsAndmammals,'sum')     

#Proportion of cells with at least 1 bird and at least 1 mammal is >99%
cellStats(birdsAndmammals,'sum') / length(getValues(divstack2$Species.richness)[!is.na(getValues(divstack2$Species.richness))])
cellStats(birdsAndmammals,'sum') / length(getValues(divstack2$Species.richness)[getValues(divstack2$Species.richness)>0.02 & !is.na(getValues(divstack2$Species.richness))])

onlyonephy<-birdmam
onlyonephy[!is.na(onlyonephy)]<-0
onlyonephy$Birds[birdmam$Birds>=1& birdmam$Mammals==0]<-1
onlyonephy$Mammals[birdmam$Birds==0 & birdmam$Mammals>=1]<-1
plot(onlyonephy) #Some cells with only birds. No cells with only mammals
1-cellStats(onlyonephy$Birds,'sum')/(ncell(onlyonephy$Birds)-cellStats(onlyonephy$Birds,'countNA'))


#Bird PD
birdscom<-commdat[,c(5:11,15:19,21:26,35:37)]#21 bird species
birdphydata<-match.phylo.comm(phylo2,birdscom)
birdpd<-pd(birdphydata$com,birdphydata$phy)

pdbirdras<-raster(use1)
pdbirdras<-setValues(pdbirdras,birdpd$PD)
pdbirdras<-mask(pdbirdras,sr_geo,maskvalue=0)
plot(pdbirdras)
pdbird_p<-pdbirdras/sum(birdphydata$phy$edge.length)
summary(pdbird_p)

#Bird SR
srbirdras<-raster(use1)
srbirdras<-setValues(srbirdras,birdpd$SR)/21
srbirdras<-mask(srbirdras,sr_geo,maskvalue=0)
plot(srbirdras)

#Bird FD
birdfuncdata<-match.phylo.comm(functree2,birdscom)
birdfd<-pd(birdfuncdata$com,birdfuncdata$phy)
fdbirdras<-raster(use1)
fdbirdras<-setValues(fdbirdras,birdfd$PD)
fdbirdras<-mask(fdbirdras,sr_geo,maskvalue=0)
plot(fdbirdras)
fdbird_p<-fdbirdras/sum(birdfuncdata$phy$edge.length)
summary(fdbird_p)

#Bird FDPD
fdpdbird<-fdbird_p/pdbird_p
plot(fdpdbird)
funconbird<-1-fdpdbird
#Mammal PD
mamcom<-commdat[,c(1:4,12:14,20,27:34,38:74)]#53 mammal species
mamphydata<-match.phylo.comm(phylo2,mamcom)
mampd<-pd(mamphydata$com,mamphydata$phy)

pdmamras<-raster(use1)
pdmamras<-setValues(pdmamras,mampd$PD)
pdmamras<-mask(pdmamras,sr_geo,maskvalue=0)
plot(pdmamras)
pdmam_p<-pdmamras/sum(mamphydata$phy$edge.length)
summary(pdmam_p)

#Mammal SR
mamsr<-sum(mamscom)
srmamras<-raster(use1)
srmamras<-setValues(srmamras,mampd$SR)/53
srmamras<-mask(srmamras,sr_geo,maskvalue=0)
plot(srmamras)

#Mammal FD
mamfuncdata<-match.phylo.comm(functree2,mamcom)
mamfd<-pd(mamfuncdata$com,mamfuncdata$phy)
fdmamras<-raster(use1)
fdmamras<-setValues(fdmamras,mamfd$PD)
fdmamras<-mask(fdmamras,sr_geo,maskvalue=0)
plot(fdmamras)
fdmam_p<-fdmamras/sum(mamfuncdata$phy$edge.length)
summary(fdmam_p)

#Mammal FDPD
fdpdmam<-fdmam_p/pdmam_p
plot(fdpdmam)
funcconmam<-1-fdpdmam

#Bird and mammals diversity
birddiv<-stack(srbirdras,pdbird_p,fdbird_p,funconbird)
mamdiv<-stack(srmamras,pdmam_p,fdmam_p,funcconmam)
birddiv<-mask(birddiv,birddiv$layer.1,maskvalue=0,updatevalue=NA)#Masking cells with no birds
mamdiv<-mask(mamdiv,mamdiv$layer.1,maskvalue=0,updatevalue=NA)#Masking cells with no mammals
birdmamdiv<-stack(birddiv,mamdiv)
birdmamdiv<-mask(birdmamdiv,extend(envvars$CurrentIce,divstack2),maskvalue=1,updatevalue=NA)
names(birdmamdiv)<-c('Bird species richness','Bird phylogenetic diversity','Bird functional diversity','Bird functional convergence',
                     'Mammal species richness','Mammal phylogenetic diversity','Mammal functional diversity','Mammal functional convergence')
levelplot(birdmamdiv)

tiff('PhylogeneticFunctionalAnalysisPicanteR/BirdMamDiv.tif',width = 6,height=8,units='in',res=300)
p1 <- levelplot(birdmamdiv[[1]], at=my.at,par.settings=YlOrRdTheme(layout.heights=list(top.padding=0,bottom.padding=0)),scales=list(draw=FALSE), margin=FALSE,main=list('Bird species richness',cex=0.8))+
  layer(sp.points(np,col=1))+
  layer(sp.polygons(allarc,lwd=0.5,col=grey(0.5)))
p2 <- levelplot(birdmamdiv[[2]], at=my.at, par.settings=YlOrRdTheme(layout.heights=list(top.padding=0,bottom.padding=0)),scales=list(draw=FALSE), margin=FALSE,main=list('Bird phylogenetic diversity',cex=0.8))+
  layer(sp.points(np,col=1))+
  layer(sp.polygons(allarc,lwd=0.5,col=grey(0.5)))
p3 <- levelplot(birdmamdiv[[3]], at=my.at,par.settings=YlOrRdTheme(layout.heights=list(top.padding=0,bottom.padding=0)),scales=list(draw=FALSE), margin=FALSE,main=list('Bird functional diversity',cex=0.8))+
  layer(sp.points(np,col=1))+
  layer(sp.polygons(allarc,lwd=0.5,col=grey(0.5)))
p4 <- levelplot(birdmamdiv[[4]], par.settings=BTCTheme(region=rev(BTC(9)),layout.heights=list(top.padding=0,bottom.padding=0)),scales=list(draw=FALSE), margin=FALSE,main=list('Bird functional convergence',cex=0.8))+
  layer(sp.points(np,col=1))+
  layer(sp.polygons(allarc,lwd=0.5,col=grey(0.5)))
p5 <- levelplot(birdmamdiv[[5]], at=my.at,par.settings=YlOrRdTheme(layout.heights=list(top.padding=0,bottom.padding=0)),scales=list(draw=FALSE), margin=FALSE,main=list('Mammal species richness',cex=0.8))+
  layer(sp.points(np,col=1))+
  layer(sp.polygons(allarc,lwd=0.5,col=grey(0.5)))
p6 <- levelplot(birdmamdiv[[6]], at=my.at, par.settings=YlOrRdTheme(layout.heights=list(top.padding=0,bottom.padding=0)),scales=list(draw=FALSE), margin=FALSE,main=list('Mammal phylogenetic diversity',cex=0.8))+
  layer(sp.points(np,col=1))+
  layer(sp.polygons(allarc,lwd=0.5,col=grey(0.5)))
p7 <- levelplot(birdmamdiv[[7]], at=my.at,par.settings=YlOrRdTheme(layout.heights=list(top.padding=0,bottom.padding=0)),scales=list(draw=FALSE), margin=FALSE,main=list('Mammal functional diversity',cex=0.8))+
  layer(sp.points(np,col=1))+
  layer(sp.polygons(allarc,lwd=0.5,col=grey(0.5)))
p8 <- levelplot(birdmamdiv[[8]], par.settings=BTCTheme(region=rev(BTC(9)),layout.heights=list(top.padding=0,bottom.padding=0)),scales=list(draw=FALSE), margin=FALSE,main=list('Mammal functional convergence',cex=0.8))+
  layer(sp.points(np,col=1))+
  layer(sp.polygons(allarc,lwd=0.5,col=grey(0.5)))

print(p1, split=c(1, 1, 2, 4), more=TRUE)
print(p2, split=c(1, 2, 2, 4), more=TRUE)
print(p3, split=c(1, 3, 2, 4), more=TRUE)
print(p4, split=c(1, 4, 2, 4), more=TRUE)
print(p5, split=c(2, 1, 2, 4), more=TRUE)
print(p6, split=c(2, 2, 2, 4), more=TRUE)
print(p7, split=c(2, 3, 2, 4), more=TRUE)
print(p8, split=c(2, 4, 2, 4))
dev.off()

#Bird mammal and total richness raw
mamsr<-srmamras*53
birdsr<-srbirdras*21

speciesrich<-stack(srras,birdsr,mamsr)
speciesrich<-mask(speciesrich,srras,maskvalue=0,updatevalue=NA)
speciesrich<-mask(speciesrich,extend(envvars$CurrentIce,divstack2),maskvalue=1,updatevalue=NA)
names(speciesrich)<-c('Total species richness','Bird species richness','Mammal species richness')

tiff('PhylogeneticFunctionalAnalysisPicanteR/SpeciesRichness.tif',width = 6,height=2,units='in',res=200)
p1<-levelplot(speciesrich[[1]],par.settings=YlOrRdTheme(layout.heights=list(top.padding=0,bottom.padding=0)),scales=list(draw=FALSE), margin=FALSE,main=list('Total species richness',cex=0.8))+
  layer(sp.points(np,col=1))+
  layer(sp.polygons(allarc,lwd=0.5,col=grey(0.5)))
p2<-levelplot(speciesrich[[2]],par.settings=YlOrRdTheme(layout.heights=list(top.padding=0,bottom.padding=0)),scales=list(draw=FALSE), margin=FALSE,main=list('Bird species richness',cex=0.8))+
  layer(sp.points(np,col=1))+
  layer(sp.polygons(allarc,lwd=0.5,col=grey(0.5)))
p3<-levelplot(speciesrich[[3]],par.settings=YlOrRdTheme(layout.heights=list(top.padding=0,bottom.padding=0)),scales=list(draw=FALSE), margin=FALSE,main=list('Mammal species richness',cex=0.8))+
  layer(sp.points(np,col=1))+
  layer(sp.polygons(allarc,lwd=0.5,col=grey(0.5)))
print(p1, split=c(1, 1, 3, 1), more=TRUE)
print(p2, split=c(2, 1, 3, 1), more=TRUE)
print(p3, split=c(3, 1, 3, 1))
dev.off()

