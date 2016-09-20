library(fda)
library(zoo)
#library(ggplot2)

data(package="fda")

###Fig 4.1

Mbasis <- create.fourier.basis(rangeval = c(0,365),nbasis=109)

summer <- seq(150,250,1)
winter <- seq(335,435,1)
t <- seq(1,365,1)

data <- colMeans(MontrealTemp)

data <- as.numeric(c(data[150:365],data[1:149]))

summerTemp <- data[1:101]
winterTemp <- data[186:286]

MPar <- fdPar(fdobj = Mbasis)

Tempfd <- smooth.basis(t,data, MPar)$fd

Y <- eval.fd(1:101,Tempfd)

op <- par(mfcol = c(2,1))

plot(150:250,summerTemp,pch=19,cex=.5,xlab = "",ylab = "Temp. (deg C)",xaxt="n",yaxt = "n",ylim=c(15,25))
lines(150:250,Y,lwd=2)
axis(1, at = seq(150, 250, by = 25))
axis(2, at = seq(15, 25, by = 5))

plot(335:435,winterTemp,pch=19,cex=.5,xlab = "Days from January 1",ylab = "Temp. (deg C)",xaxt="n",yaxt = "n",ylim=c(-15,0))
lines(335:435,eval.fd(186:286,Tempfd),lwd=2)
axis(1, at = seq(335, 435, by = 25))
axis(2, at = seq(-15, 0, by = 5))

par(op)



###Fig 4.2
# NO Fels data


### Fig 4.7

Yhat <- eval.fd(t,Tempfd)
##sd(data[232:247]-Yhat[232:247])
sqrt(sum((data-Yhat)^2)/(365-109))

plot(15:30,data[232:247])
lines(15:30,Yhat[232:247])
lines(15:30,Yhat[232:247]+1.96*0.26,lty="dashed")
lines(15:30,Yhat[232:247]-1.96*0.26,lty="dashed")


var(Yhat-data)


########Tay basis
trydata <- colMeans(MontrealTemp)

tryvar <- rep(0,155)

for (i in 1:155) {
    k <- 2*i-1
    trybasis <- create.fourier.basis(rangeval = c(0,365),nbasis=k)
    tryfd <- smooth.basis(t,trydata, trybasis)$fd
    tryhat <- eval.fd(t,tryfd)
    tryvar[i] <- sum((tryhat-trydata)^2)/(365-k)
}

plot(seq(1,309,2),tryvar,ylim=c(0,0.61),type="l")






###Fig4.8

X <- (rowMeans(handwrit[,,1]))
Y <- (rowMeans(handwrit[,,2]))

plot(X,Y,)


t <- coredata(handwritTime/1000)

(delta.t <- t[3]-t[2])

x.jp1 <- X[-(1:2)]
x.jm1 <- X[-(1400:1401)]

x.D1 <- ((x.jp1-x.jm1)/(delta.t*2))*10

x.j <- X[-c(1,1401)]

x.D2 <- ((x.jp1+x.jm1-2*x.j)/(((delta.t)^2)))/10


y.jp1 <- Y[-(1:2)]
y.jm1 <- Y[-(1400:1401)]

y.D1 <- ((y.jp1-y.jm1)/(delta.t*2))*10

y.j <- Y[-c(1,1401)]

y.D2 <- ((y.jp1+y.jm1-2*y.j)/(((delta.t)^2)))/10

plot(ksmooth(t[-c(1,1401)], x.D2,"normal",bandwidth=.075),type = "l",ylim=c(-.75,.75))
lines(ksmooth(t[-c(1,1401)], y.D2,"normal",bandwidth=.075),lty = "dashed")
