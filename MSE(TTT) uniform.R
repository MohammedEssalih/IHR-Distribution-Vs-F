#library("fdrtool")
#library(VGAM)
#library("decon")
#library(tidyverse)
n2=100
p1=1
p2=1
x2 <- rbeta(n2,p1,p2)
sig2=0.1
u2=ifelse(runif(n2) > 0.5, 1, -1) * rexp(n2, rate = 1/sig2)
w2 <- x2+u2
x=seq(0,1,0.01)
# estimate the bandwidth with the bootstrap method with re-sampling
bw2 <- bw.dnrd(w2,sig=sig2, error="laplacian")
# estimate the distribution function with measurement error
F2=DestimatePdf(x, w2, sig2, bw2)
F2 <- DeconCdf(w2,sig2,x,error='laplacian',bw=bw2)
# Create a smooth function
# Creat the Least concave majorant LCM
ll = gcmlcm(x,F2$y, type="lcm")
# Create a smooth function of the LCM
smooth_function <- splinefun(ll$x.knots, ll$y.knots,method = "natural")
# Create a data frame
M2=tibble(
  x = x, 
  y= F2$y,
  z=smooth_function(x),
  p=pbeta(x,p1,p2)
)
# Plote the data
ggplot(M2, aes(x = x)) +
  geom_line(aes(y = y), color = "red") +
  geom_line(aes(y=z),color="blue")+
  geom_line(aes(y=p),color="black")



###########################################################
#### Calculate the MSE for the Beta Distribution
###########################################################
m=50 # the number of simulations 
k=10 # the number of points support 
p1=1
p2=4
NSR=0.1 # the Noise to Signal Ration

varx=p1*p1^2 # the variance of Beta distribution 
sig2 =(varx/(2*NSR))^0.5
results=numeric(k-1)
for (j in 2:k-1) {
  p=j/k
  MSE1=0
  MSE2=0
  for (i in 1:m) {
    n2 <- 10
    x2 <- rbeta(n2,p1,p2)
    u2 <- rlaplace(n2, 0,sig2)
    w2 <- x2+u2
    
    # estimate the bandwidth with the bootstrap method with resampling
    bw2 <- bw.dboot2(w2,sig=sig2, error="normal")
    # estimate the distribution function with measurement error
    
    F2 <- DeconCdf(w2,sig2,error='normal',bw=bw2,from = 0)
    
    
    # Create a smooth function
    smooth_function1 <- splinefun(F2$x, F2$y)
    y1=smooth_function1(qbeta(p, p1, p2))
    #p1*(-log(1-p))^(1/p2)
    ll = gcmlcm(F2$x,F2$y, type="lcm")
    
    # Create a smooth function
    
    smooth_function2 <- splinefun(ll$x.knots, ll$y.knots,method = "natural")
    
    y2=smooth_function2(qbeta(p, p1, p2))
    
    MSE1=MSE1+(y1-p)^2
    MSE2=MSE2+(y2-p)^2
  }
  
  results[j] = MSE2/MSE1
}

M2=tibble(
  x = c(2:k-1)/k, 
  ratio = results
)

ggplot(M2, aes(x = x)) +
  geom_line(aes(y = ratio), color = "blue")  # Plot the second line

