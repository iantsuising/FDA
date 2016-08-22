library(fda)
library(ggplot2)
library(zoo)
library(reshape)

source("multipolt.R")

#pdf("pic.in.ch2.pdf")

###### 1.15 in 2009
age <- growth$age


age.rng <- range(age)
agefine <- seq(age.rng[1], age.rng[2], length=501)

gr.basis = create.bspline.basis(norder=6, breaks=growth$age)


children = 1:10
ncasef   = length(children)

cvecf           = matrix(0, gr.basis$nbasis, ncasef)
dimnames(cvecf) = list(gr.basis$names,
              dimnames(growth$hgtf)[[2]][children])

gr.fd0  = fd(cvecf, gr.basis)

gr.Lfd    = 3


gr.lambda = 10^(-1.5)

#  Define the functional parameter object

gr.fdPar  = fdPar(gr.fd0, Lfdobj=gr.Lfd, lambda=gr.lambda)

#  Monotonically smooth the female data

hgtfmonfd   = with(growth, smooth.monotone(age, hgtf[,children], gr.fdPar))

hgtf.vel = predict(hgtfmonfd$yhatfd, agefine, 1)
hgtf.acc = predict(hgtfmonfd$yhatfd, agefine, 2)


plot(hgtf.vel, hgtf.acc, type='n', xlim=c(0, 12), ylim=c(-5, 2),
     xlab='Velocity (cm/yr)', ylab=expression(Acceleration (cm/yr^2)),
     las=1)
for(i in 1:10){
  lines(hgtf.vel[, i], hgtf.acc[, i])
  points(hgtf.vel[316, i], hgtf.acc[316, i])
}
abline(h=0, lty='dotted')




### pinch data

p1 <- ggplot()+
    geom_line(aes(x=pinchtime, y=rowMeans(pinch)))

p2 <- ggplot()+
    geom_line(aes(x=pinchtime, y=apply(pinch, 1, sd ) ))


matplot(pinchtime, pinch, type="l")

# experiment with lambda to get a plot that looks
# reasonably smooth but not too smooth
pinch.fd <- smooth.basisPar(pinchtime, pinch)
pinch.fd <- smooth.basisPar(pinchtime, pinch, lambda=1e-6)

plot(pinch.fd$fd)

plot(mean(pinch.fd$fd), xlab="Time (sec.)", ylab="Force mean (N)")
abline(v=.1, lty="dotted")

plot(sd.fd(pinch.fd$fd), ylim=c(0, 1), xlab="Time (sec.)",
     ylab="Force Std. Dev. (N)")
abline(v=.1, lty="dotted")

####gait
(gaittime <- as.numeric(dimnames(gait)[[1]]))

gaitbasis21 <- create.fourier.basis(nbasis=21)
# may be 11 or 15 basis

# A finite Fourier series must satisfy a harmonic
# differential equation, as described in the book.
# Give that a name.  
harmaccelLfd01 <- vec2Lfd(c(0, 0, (2*pi)^2, 0))

gaitfd <- smooth.basisPar(gaittime, gait,
       gaitbasis21, Lfdobj=harmaccelLfd01, lambda=1e-11)$fd

gaitfd$fdnames
# Change some of the names 
names(gaitfd$fdnames) = c("Normalized time", "Child", "Angle")
gaitfd$fdnames[[3]] = c("Hip", "Knee")

# Establish a finer grid so the plot appears smoother 
(gaitTimeSm <- seq(0, 1, length=41))
gaitCor <- cor.fd(gaitTimeSm, gaitfd)

gaitCorLbls <- (-5:5)/5

#jpeg("gait.jpeg",width = 800, height = 800)

#par(mfrow=c(2,2))

contour(gaitTimeSm, gaitTimeSm, gaitCor[,,1,1], labcex=1.2,
        levels=gaitCorLbls, xlab="Knee",
        ylab="Knee")
x01. <- seq(0, 1, .05)
text(x01., x01., 1)


contour(gaitTimeSm, gaitTimeSm, gaitCor[,,1,2], labcex=1.2,
        levels=gaitCorLbls, ylab="Knee",
        xlab="Hip")

contour(gaitTimeSm, gaitTimeSm, gaitCor[,,1,2], labcex=1.2,
        levels=gaitCorLbls, xlab="Knee",
        ylab="Hip")

contour(gaitTimeSm, gaitTimeSm, gaitCor[,,1,3], labcex=1.2,
        levels=gaitCorLbls, xlab="Hip",
        ylab="Hip")
text(x01., x01., 1)

#dev.off()

####Canadian


daybasis65 <- create.fourier.basis(c(0, 365), nbasis=65)
#  -----------  set up the harmonic acceleration operator  ----------
harmaccelLfd365 <- vec2Lfd(c(0,(2*pi/365)^2,0), c(0, 365))



TempSm6 <- smooth.basisPar(argvals=day.5, y=CanadianWeather$dailyAv[,,1],
     fdobj=daybasis65, Lfdobj=harmaccelLfd365, lambda=1e6)$fd

PrecipSm6 <- smooth.basisPar(argvals=day.5, y=CanadianWeather$dailyAv[,,3],
     fdobj=daybasis65, Lfdobj=harmaccelLfd365, lambda=1e6)$fd

lvl.1 <- (-10:10)/10

#jpeg("canadian.jpeg",width = 800, height = 800)

#par(mfrow=c(2,2))


contour(weeks, weeks, cor.fd(weeks, TempSm6), levels=lvl.1,
        xlab="Temp", 
        ylab="Temp", 
        axes=FALSE)
axisIntervals(1)
axisIntervals(2)
text(weeks,weeks,1)


contour(weeks, weeks, cor.fd(weeks, TempSm6, weeks, PrecipSm6),
        levels=lvl.1, xlab="Temp", 
        ylab="Prec", axes=FALS)
axisIntervals(1)
axisIntervals(2)

contour(weeks, weeks, cor.fd(weeks, PrecipSm6, weeks, TempSm6),
        levels=lvl.1, ylab="Temp", 
        xlab="Prec", axes=FALSE)
axisIntervals(1)
axisIntervals(2)

contour(weeks, weeks, cor.fd(weeks, PrecipSm6), levels=lvl.1,
        xlab="Prec", 
        ylab="Prec", 
        axes=FALSE)
axisIntervals(1)
axisIntervals(2)
text(weeks,weeks,1)

#dev.off()





### nondurable



plot(nondurables, log="y")
str(nondurables)

plot(log10(nondurables))
(log10.nondur.lm <- lm(log10(nondurables)~time(nondurables)))
abline(log10.nondur.lm, lty="dashed")

length(nondurables)

diff(range(time(nondurables)))

(nondurGrowth <- (nondurables[973]/nondurables[1]))


log(nondurables[973]/nondurables[1])/81


nondur1964.67 <- window(nondurables, 1964, 1967)

plot(log10(nondur1964.67), type="p", axes=FALSE, xlab="Year",
     ylab=expression(paste(log[10], " nondurable goods index")) )
axis(2)
axis(1, 1964:1967)
axis(1, seq(1964, 1967, by=0.5), labels=FALSE)


durtimefine <- seq(1964, 1967, length=181)

goodsbasis <- create.bspline.basis(rangeval=c(1919,2000),
                                   nbasis=979, norder=8)
LfdobjNonDur <- int2Lfd(4) 

logNondurSm <- smooth.basisPar(argvals=index(nondurables),
                y=log10(coredata(nondurables)), fdobj=goodsbasis,
                Lfdobj=LfdobjNonDur, lambda=1e-11)

logNondurSm1964.67 = eval.fd(durtimefine, logNondurSm$fd);
lines(durtimefine, logNondurSm1964.67)
abline(v=1965:1966, lty=2)

sin. <- expression(sin(2*pi*x))
D.sin <- D(sin., "x")
D2.sin <- D(D.sin, "x")

with(data.frame(x=seq(0, 1, length=46)),
     plot(eval(D.sin), eval(D2.sin), type="l",
          xlim=c(-10, 10), ylim=c(-50, 50), 
          xlab="Velocity", ylab="Acceleration") )
pi.2 <- (2*pi)

abline(h=0, lty="longdash")
pi.2.2 <- pi.2^2
lines(x=c(0,0), y=c(-pi.2.2, pi.2.2), lty="longdash")

text(c(0,0), c(-47, 47), rep("no kinetic, max potential", 2))
text(c(-8.5,8.5), c(0,0), rep("max kinetic\nno potential", 2))


goodsbasis <- create.bspline.basis(rangeval=c(1919,2000),
                                   nbasis=979, norder=8)
LfdobjNonDur <- int2Lfd(4) 

#logNondurSm <- smooth.basisPar(argvals=index(nondurables),
 #               y=log10(coredata(nondurables)), fdobj=goodsbasis,
  #              Lfdobj=LfdobjNonDur, lambda=1e-11)
phaseplanePlot(1964, logNondurSm$fd,xlim=c(-0.5,0.5),ylim=c(-7,10))

#dev.off()
