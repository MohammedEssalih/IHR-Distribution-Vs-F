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
  0.7*dbeta(x,0.5,0.5)+0.3*dweibull(x,2.5,0.75)
}

cdf <- function(x) {
  0.5*pbeta(x,0.5,0.5)+0.5*pweibull(x,1.5,1)
}



convolve_cdf <- function(t_values) {
  sapply(t_values, function(t) {
    integrand <- function(x) {
      cdf(x) * dens1(t - x)
    }
    integrate(integrand, -Inf, Inf)$value
  })
}

#conv <- convolve(dens1, dens2)

t1=seq(0.1,4,0.01)
#convolve_cdf(t1)

#plot(t1,cdf(t1),type="line")
#lines(t1,convolve_cdf(t1),type = "line",col=2)

Hazard_values1=-log(1-cdf(t1))
Hazard_values2=-log(1-convolve_cdf(t1))
plot(t1,Hazard_values1,type="line",col=1)
lines(t1,Hazard_values2,type="line",col=2)