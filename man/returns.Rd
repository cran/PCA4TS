\name{returns}
\alias{returns}
\docType{data}
\title{
Daily returns of six stocks
}
\description{
The data contains the daily returns of the stocks
of Bank of America Corporation, Dell Inc., JPMorgan
Chase Co., FedEx Corporation, McDonald's Corp. and American
International Group in 2 January  2002 -- 10 July 2008. 
}
\usage{data("returns")}

\examples{
data(returns)

# Time series plots of the returns 

Y=returns
names(Y)=c("Y1","Y2","Y3","Y4","Y5","Y6")
location=c(1,253,505,757,1009,1260,1511)
year=c(2002,2003,2004,2005,2006,2007,2008)

par(mfrow=c(6,1),mar=c(2,2,0.5,1),cex=0.5, mgp = c(1.4, 0.6, 0))
plot(Y[,1],type='l',lty=1,xlab='',ylab='',col="blue",xaxt="n")
axis(1, at=location, labels=year, col.axis='brown')
plot(Y[,2],type='l',lty=1,xlab='',ylab='',col="blue",xaxt="n")
axis(1, at=location, labels=year, col.axis='brown')
plot(Y[,3],type='l',lty=1,xlab='',ylab='',col="blue",xaxt="n")
axis(1, at=location, labels=year, col.axis='brown')
plot(Y[,4],type='l',lty=1,xlab='',ylab='',col="blue",xaxt="n")
axis(1, at=location, labels=year, col.axis='brown')
plot(Y[,5],type='l',lty=1,xlab='',ylab='',col="blue",xaxt="n")
axis(1, at=location, labels=year, col.axis='brown')
plot(Y[,6],type='l',lty=1,xlab='',ylab='',col="blue",xaxt="n")
axis(1, at=location, labels=year, col.axis='brown')

}
\keyword{datasets}
