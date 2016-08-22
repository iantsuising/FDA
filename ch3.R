library(fda)
library(ggplot2)
library(zoo)

pdf("pic.in.ch3.pdf")
###handwriting

dim(handwrit)

length(handwritTime)

X=1000*handwrit[,1,1]
Y=1000*handwrit[,1,2]
p3.1 <- ggplot()+
    geom_path(aes(x=X,y=Y))+
    labs(x="",y="")   
p3.1

p3.2.x <- ggplot()+
    geom_line(aes(x=coredata(handwritTime/1000),y=X))+
    labs(x="Time (Sec)", y="X coordinate")+
    annotate("text",x=c(0.4,1.23,2),y=c(39,39,39),label=c("f","d","a"),size=6)+
    geom_vline(xintercept=c(0.8,1.6),linetype="dashed")

p3.2.x


p3.2.y <- ggplot()+
    geom_line(aes(x=coredata(handwritTime/1000),y=Y))+
    labs(x="Time (Sec)", y="Y coordinate")+
    annotate("text",x=c(0.4,1.23,2),y=c(39,39,39),label=c("f","d","a"),size=6)+
    geom_vline(xintercept=c(0.8,1.6),linetype="dashed")

p3.2.y

##############3.3

X <- (rowMeans(handwrit[,,1]))


t <- coredata(handwritTime/1000)

(delta.t <- t[3]-t[2])

x.jp1 <- X[-(1:2)]
x.jm1 <- X[-(1400:1401)]

x.D1 <- ((x.jp1-x.jm1)/(delta.t*2))*10

x.j <- X[-c(1,1401)]

x.D2 <- ((x.jp1+x.jm1-2*x.j)/(((delta.t)^2)))/10

p3.3.1 <- ggplot()+
    geom_path(aes(x=t[-c(1,1401)],y=x.D1))+
    geom_hline(yintercept = 0, linetype="dashed")+
    labs(x="Time (Sce)", y= "X first difference/100")

p3.3.1

p3.3.2 <- ggplot()+
    geom_path(aes(x=t[-c(1,1401)],y=x.D2))+
    geom_hline(yintercept = 0, linetype="dashed")+
    labs(x="Time (Sce)", y= "X second difference/10000")


p3.3.2



##3.3 loop
#pdf("deletedme.3.3.pdf")
#for (i in 1:20){
#    X <- 1000*handwrit[,i,1]
#    t <- coredata(handwritTime/1000)
#
#   delta.t <- t[3]-t[2]
#
#    x.jp1 <- X[-(1:2)]
#    x.jm1 <- X[-(1400:1401)]
#
#    x.D1 <- ((x.jp1-x.jm1)/(delta.t*2))/100
#
#    plot(x=t[-c(1,1401)],y=x.D1,type="l")
#}
#dev.off()
##

X=1000*handwrit[,,1]
Y=1000*handwrit[,,2]

scriptBasis <- create.bspline.basis(rangeval = c(0,2.3), norder=6, breaks = t )

x.fd <- smooth.basisPar(argvals = t, y=X, fdobj = scriptBasis, Lfdobj = int2Lfd(4),lambda = 1e-11)$fd

y.fd <- smooth.basisPar(argvals = t, y=Y, fdobj = scriptBasis, Lfdobj = int2Lfd(4),lambda = 1e-11)$fd


x.fd.D1 <- deriv.fd(x.fd,1)

x.fd.D2 <- deriv.fd(x.fd,2)

y.fd.D1 <- deriv.fd(y.fd,1)

y.fd.D2 <- deriv.fd(y.fd,2)

plot(x.fd)
title("smoothed Xscript")

plot(y.fd)
title("smoothed Yscript")


plot(x.fd.D1)
title("smoothed Xscript fisrt dereiv")

plot(x.fd.D2)
title("smoothed Xscript second dereiv")


plot(y.fd.D1)
title("smoothed Yscript fisrt dereiv")

plot(y.fd.D2)
title("smoothed Yscript second dereiv")

dev.off()
