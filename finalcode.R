#Load the libraries:
install.packages("spdep")

library(maptools) # loads sp library too
library(RColorBrewer) # creates nice color schemes
library(classInt) # finds class intervals for continuous variables
library(rgeos)
library(gpclib)
#library(rgdal)
library(gpclib)
#library(gstat)
library(fields)
library(spdep)
library(maps)
library(car)

#Read and fix the data:
CR.shp <- readShapePoly(file.choose(),proj4string=CRS("+proj=longlat"))
CR<-CR.shp[CR.shp@data$NAME2!="Isla",]
CR<-CR[CR@data$NAME2!="lag.arenal",]
CR@data$ID<-c(1:86)
CR<-CR[CR@data$ID!="45",]
dim(CR@data)
CR@data$ID<-c(1:85)

dat1<-read.csv(file.choose())
attach(dat1)

a<-merge(CR@data,dat1,by.x="DCODE",by.y="DCODE")
a[c(4,5,6,7),]
a[c(41,42,43,44),]
a[c(63,64,65,66),]
a[c(70,71,72,73,74,75),]

b<-(a[-(c(5,6,42,43,64,65,71,72,73,75)),])
b<-b[order(b$ID),]


CR@data<-b
attach(CR@data)


#####################
# Describe the data #
#####################

#All variables:
scatterplot.matrix(data.frame(cbind(var1,var2,CR@data[,-c(1:28)])))

#Response variables:

#Variable 1: GC rate of people resident on a county per each 10000 people.

var1<-(CR@data$res/CR@data$pop)*10000
var2<-((CR@data$born/CR@data$pop)*10000)

qqnorm(var1)
qqline(var1)

qqnorm(var2)
qqline(var2)

install.packages("truncgof")
library(truncgof)
ks.test(var1, "plnorm", list(meanlog=1.75, sdlog=1))

#Plot any variable:
plotvar <- var1
nclr <- 4
plotclr <- brewer.pal(nclr,"GnBu")
#plotclr <- plotclr[nclr:1] # reorder colors if appropriate
max.symbol.size=12
min.symbol.size=1
class <- classIntervals(plotvar, nclr, style="quantile")
colcode <- findColours(class, plotclr)
cutpts <- round(class$brks, digits=4)
symbol.size <- ((plotvar-min(plotvar))/(max(plotvar)-min(plotvar))*(max.symbol.size-min.symbol.size)+min.symbol.size)
# 4 plot the shapefile
plot(CR)
plot(CR, col=colcode, add=T)
# 5 -- add a title and legend
title("Costa Rica: Reported cases of Gastric Cancer \n per 10 000 population, 2005")
legend(x=-86, y=9, legend=leglabs(round(cutpts,4)), fill=colcode, bty="n",x.intersp = 1, y.intersp = 1)


par(mfrow=c(1,3))
plotvar <- IEV
nclr <- 4
plotclr <- brewer.pal(nclr,"GnBu")
#plotclr <- plotclr[nclr:1] # reorder colors if appropriate
max.symbol.size=12
min.symbol.size=1
class <- classIntervals(plotvar, nclr, style="quantile")
colcode <- findColours(class, plotclr)
cutpts <- round(class$brks, digits=4)
symbol.size <- ((plotvar-min(plotvar))/(max(plotvar)-min(plotvar))*(max.symbol.size-min.symbol.size)+min.symbol.size)
# 4 plot the shapefile
plot(CR)
plot(CR, col=colcode, add=T)
# 5 -- add a title and legend
title("Costa Rica: Life expentancy by county, 2005")
legend(x=-86, y=9, legend=leglabs(round(cutpts,4)), fill=colcode, bty="n",x.intersp = 1, y.intersp = 1)



#Centroids and adjancency:
plot(CR)
CR.cntr <- coordinates(CR)
points(CR.cntr, pch=16)
title("Costa Rica: Centroids per county")
legend(x=-86, y=9, legend=leglabs(round(cutpts,4)), fill=colcode, bty="n",x.intersp = 1, y.intersp = 1)


threshold<-20
Dist<-rdist.earth(CR.cntr,miles=T)
ADJbyDist<-ifelse(Dist<threshold,1,0)
nb2<-mat2listw(ADJbyDist)
summary(nb2)
plot(CR)
plot.listw(nb2,coords=CR.cntr,add=T,col=4,lwd=2)


threshold<-15
Dist<-rdist.earth(CR.cntr,miles=T)
ADJbyDist<-ifelse(Dist<threshold,1,0)
nb3<-mat2listw(ADJbyDist)
summary(nb3)
plot(CR)
plot.listw(nb3,coords=CR.cntr,add=T,col=3,lwd=2)
title("Costa Rica: Neighbor counties within 15 kilometers")


threshold<-10
Dist<-rdist.earth(CR.cntr,miles=T)
ADJbyDist<-ifelse(Dist<threshold,1,0)
nb4<-mat2listw(ADJbyDist)
summary(nb4)
plot(CR)
plot.listw(nb4,coords=CR.cntr,add=T,col=4,lwd=2)
title("Costa Rica: Neighbor counties within 10 kilometers")


map<-CR
nb.r<-poly2nb(map, queen=T)    #two areas are neighbors if share common edges with length >0
mat<-nb2mat(nb.r, style="B", zero.policy=T)        #mat is the 0/1 adjacency matrix
n.site<-dim(mat)[1]                 #n.site: number of areas
n.edge<-sum(mat)/2      #n.edge: number of unique pairs (ij and ji only count once).

nb1<-mat2listw(mat)
summary(nb1)

s<-coordinates(CR)
plot(CR)
plot.listw(nb1,coords=s,add=T,col=4,lwd=2)
title("Costa Rica: Neighbour counties by adjancency")


#Moran test for the data:

moran.test(var1,listw=nb1,randomisation=T)
moran.test(var1,listw=nb2,randomisation=T)
moran.test(var1,listw=nb3,randomisation=T)
moran.test(var1,listw=nb4,randomisation=T)

moran.test(var1,listw=nb1,randomisation=F)
moran.test(var1,listw=nb2,randomisation=F)
moran.test(var1,listw=nb3,randomisation=F)
moran.test(var1,listw=nb4,randomisation=F)

moran.mc(var1,listw=nb1,10000)
moran.mc(var1,listw=nb2,10000)
moran.mc(var1,listw=nb3,10000)
moran.mc(var1,listw=nb4,10000)

plot.spcor(sp.correlogram(nb.r, var2, order = 7, method = "I"), xlab = "Spatial lags", main = "Spatial correlogram: Autocorrelation CIs")

moran.mc(var2,listw=nb1,1000)

rho<-0.9999
sigma2<-1^2
M<-diag(rowSums(nb2))
Q<-(M-rho*W)/sigma2
COV<-solve(Q)
VAR<-diag(COV)
COR<-cov2cor(COV)

plot.nc(VAR,main="Variance",sig_digits=2)
plot.nc(COR[1,],main="Correlation with County 1",sig_digits=2,breaks=seq(0,1,.1))





#MOdels:

a1<-lm(var1~IEV+IBM+IC)
outlierTest(a1)
qqPlot(a1, main="QQ Plot")  
leveragePlots(a1) # leverage plots
# identify D values > 4/(n-k-1) 
cutoff <- 4/((nrow(CR@data)-length(a1$coefficients)-2)) 
plot(a1, which=4, cook.levels=cutoff)
# Influence Plot 
influencePlot(a1, id.method="identify", main="Influence Plot", sub="Circle size is proportial to Cook's Distance" )
#outlier: Leon Cortes.

# Normality of Residuals
# qq plot for studentized resid
par(mfrow=c(1,2))
qqPlot(a1, main="QQ Plot")
shapiro.test(a1$res)
# distribution of studentized residuals
library(MASS)
sresid <- studres(a1) 
hist(sresid, freq=FALSE,main="Distribution of Studentized Residuals")
xfit<-seq(min(sresid),max(sresid),length=40) 
yfit<-dnorm(xfit) 
lines(xfit, yfit)

# Evaluate homoscedasticity
# non-constant error variance test
ncvTest(a1)
# plot studentized residuals vs. fitted values 
spreadLevelPlot(a1, main="Spread-Level Plot")

vif(a1) # variance inflation factors 
sqrt(vif(a1)) > 2 # problem?

# Evaluate Nonlinearity
# component + residual plot 
crPlots(a1)
# Ceres plots 
ceresPlots(a1)

durbinWatsonTest(a1)

install.packages("gvlma")
library(gvlma)
gvmodel <- gvlma(a1) 
summary(gvmodel)


library(HH)
hov(var1~IEV+IBM+IC)
plot.hov(var1~IEV+IBM+IC)


#Plot any variable:
plotvar <- a1$res
nclr <- 4
plotclr <- brewer.pal(nclr,"GnBu")
#plotclr <- plotclr[nclr:1] # reorder colors if appropriate
max.symbol.size=12
min.symbol.size=1
class <- classIntervals(plotvar, nclr, style="quantile")
colcode <- findColours(class, plotclr)
cutpts <- round(class$brks, digits=4)
symbol.size <- ((plotvar-min(plotvar))/(max(plotvar)-min(plotvar))*(max.symbol.size-min.symbol.size)+min.symbol.size)
# 4 plot the shapefile
plot(CR)
plot(CR, col=colcode, add=T)
# 5 -- add a title and legend
title("Residuals for OLS model")
legend(x=-86, y=9, legend=leglabs(round(cutpts,4)), fill=colcode, bty="n",x.intersp = 1, y.intersp = 1)



moran.mc(var1,listw=nb1,10000)
moran.mc(var1,listw=nb2,10000)
moran.mc(var1,listw=nb3,10000)
moran.mc(var1,listw=nb4,10000)

moran.mc(a1$res,listw=nb1,10000)
moran.mc(a1$res,listw=nb2,10000)
moran.mc(a1$res,listw=nb3,10000)
moran.mc(a1$res,listw=nb4,10000)


summary(spautolm(var1~IEV+IBM+IC,listw=nb1,family="CAR"))
summary(spautolm(var1~IEV+IBM+IC,listw=nb2,family="CAR"))
summary(spautolm(var1~IEV+IBM+IC,listw=nb3,family="CAR"))
summary(spautolm(var1~IEV+IBM+IC,listw=nb4,family="CAR"))
#Best model: centroid 20 miles or adjacency.

a2<-spautolm(var1~IEV+IBM+IC,listw=nb2,family="CAR")

summary(spautolm(var1~IEV+IBM+IC,listw=nb1,family="SAR"))
summary(spautolm(var1~IEV+IBM+IC,listw=nb2,family="SAR"))
summary(spautolm(var1~IEV+IBM+IC,listw=nb3,family="SAR"))
summary(spautolm(var1~IEV+IBM+IC,listw=nb4,family="SAR"))
#Best model: centroid 20 miles.

a3<-spautolm(var1~IEV+IBM+IC,listw=nb2,family="SAR")

par(mfrow=c(1,3))
plotvar <- a1$res
nclr <- 4
plotclr <- brewer.pal(nclr,"GnBu")
#plotclr <- plotclr[nclr:1] # reorder colors if appropriate
max.symbol.size=12
min.symbol.size=1
class <- classIntervals(plotvar, nclr, style="quantile")
colcode <- findColours(class, plotclr)
cutpts <- round(class$brks, digits=4)
symbol.size <- ((plotvar-min(plotvar))/(max(plotvar)-min(plotvar))*(max.symbol.size-min.symbol.size)+min.symbol.size)
# 4 plot the shapefile
plot(CR)
plot(CR, col=colcode, add=T)
# 5 -- add a title and legend
title("Residuals for OLS model")
legend(x=-86, y=9, legend=leglabs(round(cutpts,4)), fill=colcode, bty="n",x.intersp = 1, y.intersp = 1)


#########################################
NEW MODEL INCLUDING BORN AS A PREDICTED VARIABLE
#########################################


a2<-lm(var1~var2+IEV+IBM+IC)
outlierTest(a2)
qqPlot(a2, main="QQ Plot")  
leveragePlots(a2) # leverage plots
# identify D values > 4/(n-k-1) 
cutoff <- 4/((nrow(CR@data)-length(a2$coefficients)-2)) 
plot(a2, which=4, cook.levels=cutoff)
# Influence Plot 
influencePlot(a2, id.method="identify", main="Influence Plot", sub="Circle size is proportial to Cook's Distance" )
#outlier: Leon Cortes.

# Normality of Residuals
# qq plot for studentized resid
par(mfrow=c(1,2))
qqPlot(a2, main="QQ Plot")
shapiro.test(a2$res)
# distribution of studentized residuals
library(MASS)
sresid <- studres(a2) 
hist(sresid, freq=FALSE,main="Distribution of Studentized Residuals")
xfit<-seq(min(sresid),max(sresid),length=40) 
yfit<-dnorm(xfit) 
lines(xfit, yfit)

# Evaluate homoscedasticity
# non-constant error variance test
ncvTest(a2)
# plot studentized residuals vs. fitted values 
spreadLevelPlot(a2, main="Spread-Level Plot")

vif(a2) # variance inflation factors 
sqrt(vif(a2)) > 2 # problem?

# Evaluate Nonlinearity
# component + residual plot 
crPlots(a2)
# Ceres plots 
ceresPlots(a2)

durbinWatsonTest(a2)

install.packages("gvlma")
library(gvlma)
gvmodel <- gvlma(a2) 
summary(gvmodel)


library(HH)
hov(var1~var2+IEV+IBM+IC)
plot.hov(var1~var2+IEV+IBM+IC)


#Plot any variable:
plotvar <- a2$res
nclr <- 4
plotclr <- brewer.pal(nclr,"GnBu")
#plotclr <- plotclr[nclr:1] # reorder colors if appropriate
max.symbol.size=12
min.symbol.size=1
class <- classIntervals(plotvar, nclr, style="quantile")
colcode <- findColours(class, plotclr)
cutpts <- round(class$brks, digits=4)
symbol.size <- ((plotvar-min(plotvar))/(max(plotvar)-min(plotvar))*(max.symbol.size-min.symbol.size)+min.symbol.size)
# 4 plot the shapefile
plot(CR)
plot(CR, col=colcode, add=T)
# 5 -- add a title and legend
title("Residuals for OLS model")
legend(x=-86, y=9, legend=leglabs(round(cutpts,4)), fill=colcode, bty="n",x.intersp = 1, y.intersp = 1)



moran.mc(var1,listw=nb1,10000)
moran.mc(var1,listw=nb2,10000)
moran.mc(var1,listw=nb3,10000)
moran.mc(var1,listw=nb4,10000)

moran.mc(a2$res,listw=nb1,10000)
moran.mc(a2$res,listw=nb2,10000)
moran.mc(a2$res,listw=nb3,10000)
moran.mc(a2$res,listw=nb4,10000)


summary(spautolm(var1~var2+IEV+IBM+IC,listw=nb1,family="CAR"))
summary(spautolm(var1~var2+IEV+IBM+IC,listw=nb2,family="CAR"))
summary(spautolm(var1~var2+IEV+IBM+IC,listw=nb3,family="CAR"))
summary(spautolm(var1~var2+IEV+IBM+IC,listw=nb4,family="CAR"))
#Best model: centroid 20 miles or adjacency.

a3<-spautolm(var1~var2+IEV+IBM+IC,listw=nb2,family="CAR")

summary(spautolm(var1~var2+IEV+IBM+IC,listw=nb1,family="SAR"))
summary(spautolm(var1~var2+IEV+IBM+IC,listw=nb2,family="SAR"))
summary(spautolm(var1~var2+IEV+IBM+IC,listw=nb3,family="SAR"))
summary(spautolm(var1~var2+IEV+IBM+IC,listw=nb4,family="SAR"))
#Best model: centroid 20 miles.

a4<-spautolm(var1~IEV+IBM+IC,listw=nb2,family="SAR")

par(mfrow=c(1,3))
plotvar <- a1$res
nclr <- 4
plotclr <- brewer.pal(nclr,"GnBu")
#plotclr <- plotclr[nclr:1] # reorder colors if appropriate
max.symbol.size=12
min.symbol.size=1
class <- classIntervals(plotvar, nclr, style="quantile")
colcode <- findColours(class, plotclr)
cutpts <- round(class$brks, digits=4)
symbol.size <- ((plotvar-min(plotvar))/(max(plotvar)-min(plotvar))*(max.symbol.size-min.symbol.size)+min.symbol.size)
# 4 plot the shapefile
plot(CR)
plot(CR, col=colcode, add=T)
# 5 -- add a title and legend
title("Residuals for OLS model")
legend(x=-86, y=9, legend=leglabs(round(cutpts,4)), fill=colcode, bty="n",x.intersp = 1, y.intersp = 1)






