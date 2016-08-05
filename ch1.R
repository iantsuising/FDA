library(fda)
library(ggplot2)
library(reshape)

source("multipolt.R")

#pdf("pic.in.ch1.pdf")


###### the grouth data
data(growth)
hgtm = growth$hgtm
hgtf = growth$hgtf                       
age = growth$age


knots  <- growth$age
norder <- 6
nbasis <- length(knots) + norder - 2
hgtbasis <- create.bspline.basis(c(1,18), nbasis, norder, knots)

Lfdobj <- 4
lambda <- 1e-1
growfdPar <- fdPar(hgtbasis, Lfdobj, lambda)

hgtffd <- smooth.basis(age, hgtf[,1:10], growfdPar)$fd

age<-growth$age

dat<-melt(growth$hgtf[,1:10])

ggplot(data=dat,
       aes(x=X1,y=value,colour=X2)) +
    geom_line()+
    geom_point()+
    xlab("Age(years)")+
    ylab("Height(cm)")+
    labs(colour="10 Girls")

D1hgtffd <- deriv.fd(hgtffd,1)
D2hgtffd <- deriv.fd(hgtffd,2)

plot(D2hgtffd,ylim=c(-4,2),xlab="Age (years)",ylab="Acceleration (cm/yr^2)")


###### refinery data

data(package="fda")

dat <- as.vector(t(nondurables))

plot(x=1:973,y=dat,type="l",xlab="Year", ylab="Nondurabel Goods Index")


p1<-ggplot(refinery,aes(x=Time,y=Tray47)) +
    geom_point(size=0.5)+
    ylab("Tray 47 Level")+
    xlab("")

p2 <-ggplot(refinery,aes(x=Time,y=Reflux)) +
    geom_point(size=0.5)+
    xlab("Time")+
    ylab("Reflux Flow")

multiplot(p1,p2,cols=1)

##### Canadian weather

summary(CanadianWeather)
summary(CanadianWeather$dailyAv)
summary(CanadianWeather$place)
CanadianWeather$place

colnames(CanadianWeather$dailyAv[,,1])

data <- CanadianWeather$dailyAv[,,1]

id <- cumsum(c(0,31,28,31,30,31,30,31,31,30,31,30,31))

monthdata <- data.frame(matrix(rep(NA,35*12),ncol=35))
colnames(monthdata) <- colnames(data)

for(i in 1:12){
    monthdata[i,] <- colMeans(data[(id[i]+1):id[i+1],])    
}

monthdata <- cbind(c(1:12),monthdata)

colnames(monthdata) <- c("Month",colnames(data))

Month <- monthdata$Month

data <- melt(monthdata[,c(1,13,24,30,36)],id="Month")

ggplot(data=data, aes(x=Month,y=value,colour=variable))+
    geom_smooth(se=FALSE,)+
    xlab("Month")+
    ylab("Mean temperature")+
    labs(colour="Locations") +
    geom_point(size=2,shape=2)

knots <- 1:12
norder <- 4
nbasis <- length(knots) + 4 - 2
mon.basis <- create.bspline.basis(c(1,12), nbasis, norder,knots)

Lfdobj <- 2
lambda <- 1e-1
growfdPar <- fdPar(hgtbasis, Lfdobj, lambda)


x <- 1:12
y <- as.matrix((monthdata[,c(13,24,30,36)]))

monthfd <- smooth.basis(x,y,growfdPar)$fd

plot(monthfd)

Dmonthfd <- deriv.fd(monthfd,1)
D2monthfd <- deriv.fd(monthfd,2)
D3monthfd <- deriv.fd(monthfd,3)

Ltemp <- D3monthfd+(pi/6)^2*Dmonthfd

plot(Ltemp,xlab="Month",ylab="L-temperature")

# Directary use monthly data


ave.mon <- melt(CanadianWeather$monthlyTemp[,c(12,23,29,35)])
monthdata <- CanadianWeather$monthlyTemp[,c(12,23,29,35)]

knots <- 1:12 -> Month
nbasis <- 9
mon.basis <- create.fourier.basis(c(1,12), nbasis)

x <- Month
y <- as.matrix(monthdata)

monthfd <- smooth.basis(x,y,mon.basis)$fd

plot(monthfd)

Dmonthfd <- deriv.fd(monthfd,1)
D2monthfd <- deriv.fd(monthfd,2)
D3monthfd <- deriv.fd(monthfd,3)

Ltemp <- D3monthfd+(pi/6)^2*Dmonthfd

plot(Ltemp,xlab="Month",ylab="L-temperature")


#plot(monthfd$coefs[,12],D2monthfd$coefs[,12])


pr.pre<-CanadianWeather$dailyAv[,"Pr. Rupert","Precipitation.mm"]

ggplot()+
    geom_smooth(aes(x=1:365,y=pr.pre),se=FALSE,method = "loess")+
    geom_point(aes(x=1:365,y=pr.pre))+
    xlab("Month")+
    ylab("Precipitation (mm)")+
    scale_x_continuous(breaks=cumsum(c(0,31,28,31,30,31,30,31,31,30,31,30)),
                      labels=c("J","F","M","A","M","J","J","A","S","O","N","D"))

smoothingSpline = smooth.spline(1:365, pr.pre, spar=1)
plot(1:365,pr.pre)
lines(smoothingSpline)

####gait data
data(package="fda")

time <- seq(0.025,0.975,0.05)


hipdata <- melt(gait[,,1])
kneedata <- melt(gait[,,2])

p.hip <- ggplot(data=hipdata,aes(x=X1,y=value,colour=X2)) +
    geom_line() +
    xlab("")+
    ylab("Hip angle (degree)") +
    theme(legend.position="none")

p.knee <- ggplot(data=kneedata,aes(x=X1,y=value,colour=X2)) +
    geom_line() +
    xlab("Time (proportion of gait cycle)")+
    ylab("Knee angle (degree)") +
    theme(legend.position="none")

multiplot(p.hip,p.knee)

inx <- c(1,4,8,12,16)

pp <- ggplot()+
    geom_point(aes(gait[,1,1],gait[,1,2]))+
    geom_path(aes(gait[,1,1],gait[,1,2]),arrow=arrow())

pp2 <- pp+
    geom_point(aes(rowMeans(gait[,,1]),rowMeans(gait[,,2])))+
    geom_path(aes(rowMeans(gait[,,1]),rowMeans(gait[,,2])),linetype=2,arrow=arrow()) +
    xlab("Hip angle (degrees)") +
    ylab("Knee angle (degrees)") +
    annotate("text",x=gait[inx,1,1],y=gait[inx,1,2],label=c("A","B","C","D","E"),size=8)+
    annotate("text",x=rowMeans(gait[,,1])[inx],y=rowMeans(gait[,,2])[inx],label=c("A","B","C","D","E"),size=8)


pp2

#### pinchraw data

pinchdata <- melt(pinchraw[,1:20])
pinchdata$X1 <- rep(pinchtime,20)
pinchdata$X2 <- as.factor(pinchdata$X2)
ggplot(data=pinchdata,aes(x=X1,y=value,colour=X2))+
    geom_line() +
    xlab("Seconds")+
    ylab("Force (N)") +
    theme(legend.position="none") +
    geom_hline(yintercept = 2,linetype=2)



nd <- log(nondurables[541:552])

knots  <- 1:12
norder <- 8
nbasis <- length(knots) + norder - 2
nondurablesbasis <- create.bspline.basis(c(1,12), nbasis, norder, knots)

Lfdobj <- 4
lambda <- 1e-1
growfdPar <- fdPar(nondurablesbasis, Lfdobj, lambda)

ndfd <- smooth.basis(1:12, nd, growfdPar)$fd


D1nd <- eval.fd(1:12, deriv.fd(ndfd,1))
D2nd <- eval.fd(1:12, deriv.fd(ndfd,2))

plot(D1nd,D2nd,xlab="Velocity",ylab="Acceleration")
text(cbind(D1nd,D2nd), labels=seq(12), pos=3)
xspline(D1nd,D2nd,shape = c(0,rep(-1,10),0),open=TRUE)

library(zoo)

nondurablesbasis <- create.bspline.basis(rangeval = c(1919,2000),
                                         nbasis=161,norder=8)
Lfdobjnondurabel <- int2Lfd(4)
logNondurSm <- smooth.basisPar(argvals=index(nondurables),
                               y=log10(coredata(nondurables)),fdobj = nondurablesbasis,
                               Lfdobj = Lfdobjnondurabel,lambda=1e-11)
phaseplanePlot(1964,logNondurSm$fd)


#dev.off()


