library(fda)


data(package='fda')


splinebasis <- create.bspline.basis(c(0,10), 13)

#pdf("pic.in.rmch3.pdf")
plot(splinebasis, xlab='t', ylab='Bspline basis functions B(t)', 
     las=1, lwd=2)

# Figure 3.2
 
basis2 <- create.bspline.basis(c(0,2*pi), 5, 2)
basis3 <- create.bspline.basis(c(0,2*pi), 6, 3)
basis4 <- create.bspline.basis(c(0,2*pi), 7, 4)

time <- seq(0, 2*pi, length=201)
sin.time <- sin(time)

sin2 <- Data2fd(time, sin.time, basis2)
sin3 <- Data2fd(time, sin.time, basis3)
sin4 <- Data2fd(time, sin.time, basis4)

sin2.time <- predict(sin2, time)
sin3.time <- predict(sin3, time)
sin4.time <- predict(sin4, time)

sinRng = range(sin2.time)
pi3    = ((1:3)*pi/2)

op = par(mfrow=c(3,2), mar=c(3,4,2,2)+.1)

plot(time, sin2.time, type='l', ylim=sinRng, xlab='', ylab='Order = 2',
     main='sine(t)' )
lines(time, sin.time, lty='dashed')
abline(v=pi3, lty='dotted')

Dsin2.time = predict(sin2, time, 1)
plot(time, Dsin2.time, type='l', ylim=sinRng, xlab='', ylab='',
     main='D sine(t)')
lines(time, cos(time), lty='dashed')
abline(v=pi3, lty='dotted')

plot(time, sin3.time, type='l', ylim=sinRng, xlab='', ylab='Order = 3')
lines(time, sin.time, lty='dashed')
abline(v=pi3, lty='dotted')

Dsin3.time = predict(sin3, time, 1)
plot(time, Dsin3.time, type='l', ylim=sinRng, xlab='', ylab='')
lines(time, cos(time), lty='dashed')
abline(v=pi3, lty='dotted')

plot(time, sin4.time, type='l', ylim=sinRng, xlab='t', ylab='Order = 4')
lines(time, sin.time, lty='dashed')
abline(v=pi3, lty='dotted')

Dsin4.time <- predict(sin4, time, 1)
plot(time, Dsin4.time, type='l', ylim=sinRng, xlab='t', ylab='')
lines(time, cos(time), lty='dashed')
abline(v=pi3, lty='dotted')

par(op)

#dev.off()

daybasis65 <- create.fourier.basis(c(0,365), 65)

coefmat <- matrix(0, 65, 35, dimnames=list(
     daybasis65$names, CanadianWeather$place) )
tempfd. <- fd(coefmat, daybasis65)


fdnames <- list("Age (years)", "Child", "Height (cm)")



Tempbasis <- create.fourier.basis(c(0,365), 65)
Tempfd <- smooth.basis(day.5,
          CanadianWeather$dailyAv[,,'Temperature.C'], Tempbasis)$fd
meanTempfd <- mean(Tempfd)
sumTempfd <- sum(Tempfd)

par(mfrow=c(1,1))
plot((meanTempfd-sumTempfd*(1/35)))

# round off error, as it should be.

#  plot the temperature for Resolute and add the Canadian mean

plot(Tempfd[35], lwd=1, ylim=c(-35,20))
lines(meanTempfd, lty=2)

#  evaluate the derivative of mean temperature and plot

DmeanTempVec <- eval.fd(day.5, meanTempfd, 1)
plot(day.5, DmeanTempVec, type='l')


harmaccelLfd <- vec2Lfd(c(0,c(2*pi/365)^2, 0), c(0, 365))
LmeanTempVec <- eval.fd(day.5, meanTempfd, harmaccelLfd)

par(mfrow=c(1,1))
plot(day.5, LmeanTempVec, type="l", cex=1.2,
     xlab="Day", ylab="Harmonic Acceleration")
abline(h=0)

#pdf("pic.in.fmch4.pdf")


dayOfYearShifted <- c(182:365, 1:181)

tempmat <- daily$tempav[dayOfYearShifted, ]
tempbasis <- create.fourier.basis(c(0,365),65)

temp.fd <- smooth.basis(day.5, tempmat, tempbasis)$fd

temp.fd$fdnames <- list("Day (July 2 to June 30)",
                      "Weather Station",
                      "Mean temperature (deg. C)")

plot(temp.fd, lwd=2, xlab='Day (July 1 to June 30)',
     ylab='Mean temperature (deg. C)')

#4.2

basis13 <- create.bspline.basis(c(0,10), 13)
tvec    <- seq(0,1,len=13)
sinecoef<- sin(2*pi*tvec)
sinefd  <- fd(sinecoef, basis13, list("t","","f(t)"))
op      <- par(cex=1.2)
plot(sinefd, lwd=2)
points(tvec*10, sinecoef, lwd=2)
par(op)

cat(MontrealTemp, file='MtlDaily.txt')

MtlDaily <- matrix(scan("MtlDaily.txt",0),34,365)
thawdata <- t(MtlDaily[,16:47])

daytime <- ((16:47)+0.5)
plot(daytime, apply(thawdata,1,mean), "b", lwd=2,
     xlab="Day", ylab="Temperature (deg C)", cex=1.2)

# 4.4

thawbasis <- create.bspline.basis(c(16,48),7)
thawbasismat <- eval.basis(thawbasis, daytime)

thawcoef <- solve(crossprod(thawbasismat),
    crossprod(thawbasismat,thawdata))
thawfd <- fd(thawcoef, thawbasis,
    list("Day", "Year", "Temperature (deg C)"))
plot(thawfd, lty=1, lwd=2, col=1)

# 4.5

plotfit.fd(thawdata[,1], daytime, thawfd[1],
           lty=1, lwd=2, main='')
##
## Section 4.4 The Linear Differential Operator or Lfd Class
##

#omega           <- 2*pi/365
#thawconst.basis <- create.constant.basis(thawbasis$rangeval)
#
#betalist      <- vector("list", 3)
#betalist[[1]] <- fd(0, thawconst.basis)
#betalist[[2]] <- fd(omega^2, thawconst.basis)
#betalist[[3]] <- fd(0, thawconst.basis)
#harmaccelLfd. <- Lfd(3, betalist)
#
#accelLfd <- int2Lfd(2)
#
#harmaccelLfd.thaw <- vec2Lfd(c(0,omega^2,0), thawbasis$rangeval)
#all.equal(harmaccelLfd.[-1], harmaccelLfd.thaw[-1])
#
#class(accelLfd)
#class(harmaccelLfd)
#
#Ltempmat  <- eval.fd(day.5, temp.fd, harmaccelLfd)
#
#D2tempfd <- deriv.fd(temp.fd, 2)
#Ltempfd  <- deriv.fd(temp.fd, harmaccelLfd)


#Bspl2 <- create.bspline.basis(nbasis=2, norder=1)
#Bspl3 <- create.bspline.basis(nbasis=3, norder=2)

#corrmat  <- array(1:6/6, dim=2:3)
#bBspl2.3 <- bifd(corrmat, Bspl2, Bspl3)

dev.off()
