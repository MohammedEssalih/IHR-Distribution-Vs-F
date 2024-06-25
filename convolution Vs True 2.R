#install.packages("bayesmeta")
#library(bayesmeta)
## Not run: 
#  Skew-normal / logistic example:

dens1 <- function(x)
  # logistic distribution's density
{
  return(dlaplace(x, location=0, scale=0.1))
}

dens2 <- function(x) {
  0.0*dbeta(x,1,1)+0.2*dexp(x,1.5)+0.1*dlnorm(x,15,2.2)+0.7*dweibull(x,1,1)
}

cdf <- function(t) {
    integrate(dens2, 0, t)$value
}

conv <- convolve(dens1, dens2)

t1=seq(0.1,10,1)
Hazard_values1=-log(1-sapply(t1,cdf))
plot(t1,Hazard_values1,type="line")


t2=seq(-5,5,0.1)
Hazard_values2=-log(1-cumsum(conv$density(t2))/sum(conv$density(t2)))
plot(t2,Hazard_values2, type="line")
lines(t1,Hazard_values1,type="line",col=2)

#b=sqrt(var(0.4*dexp(t,1)+0.6*dbeta(t,1,1)))*0.4/sqrt(2)
#Hazard_values=-log(1-conv$cdf(t))
