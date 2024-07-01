#install.packages("bayesmeta")
#library(bayesmeta)
## Not run: 
#  Skew-normal / logistic example:

lap <- function(x)
  # logistic distribution's density
{
  return(dlaplace(x, 0, 0.1))
}

Lap <- function(x)
  # logistic distribution's density
{
  return(plaplace(x, 0,0.1))
}

pdf <- function(x) {
  0.15*dexp(x,0.25)+0.3*dnormTrunc(x, mean = 0, sd = 1.25)+0.35*dweibull(x,1,1)+0.2*dunif(x,0,1)
}

cdf <- function(x) {
  0.15*pexp(x,0.25)+0.3*pnormTrunc(x, mean = 0, sd = 1.25)+0.35*pweibull(x,1,1)+0.2*punif(x,0,1)
}



convolve_cdf <- function(t_values) {
  sapply(t_values, function(t) {
    integrand <- function(x) {
      pdf(x) * Lap(t - x)
    }
    integrate(integrand, -100, 100,subdivisions = 1000)$value
  })
}
x=seq(0,20,0.01)
FnT <- data.frame(x=x,
                 y=cdf(x)
)

FnC <- data.frame(x=x,
                 y=convolve_cdf(x)
)

ST=splinefun(FnT$x, 1-FnT$y)
SC=splinefun(FnC$x, 1-FnC$y)
# Restrict the CDF to only positive interval
FnC$x=ifelse(FnC$x<0,0,FnC$x)
FnT$x=ifelse(FnT$x<0,0,FnT$x)
# Create the Qauntille estimator
QT=data.frame(x=FnT$y,y=FnT$x)
QC=data.frame(x=FnC$y,y=FnC$x)
# Remove duplicates 
#Q=Q[!duplicated(Q$x), ]
#Q=Q[!duplicated(Q$y), ]
# Create the Qauntile function
DC=splinefun(QC$x, QC$y)
DT=splinefun(QT$x, QT$y)
# Calculate the mean of the underlying data (X) (after deconvolution)
muC=integrate(DC, lower = 0, upper = 1,subdivisions=2000)$value
muT=integrate(DT, lower = 0, upper = 1,subdivisions=2000)$value
# Create the TTT data
Ttt=data.frame(
  x=QT$x,
  y=TTT(ST,QT$y)/muT,
  z=TTT(SC,QC$y)/muC
)
Ttt|>
  ggplot(aes(x=x))+
  geom_line(aes(y=y),color="blue")+
  geom_line(aes(y=z),color="red")
  