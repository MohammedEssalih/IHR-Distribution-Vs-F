#library("decon")
#library("tidyverse")
# Plote and ensure
nsr=0.1
n=100
X <- rgamma(n,0.5,1)
varx=0.5
sig.x=sqrt(varx)
sig.e=nsr*sig.x
u <- rlaplace(n, 0,sig.e)
w <- X+u
##########
bw <- bw.dnrd(w,sig=sig.e, error="laplacian")
Fn <- DeconCdf(w,sig.e,error='laplacian',bw=bw)
#Modify the CDF
Fn$x=ifelse(Fn$x<0,0,Fn$x)
# Create and plot the survival functions
Sn=splinefun(Fn$x, 1-Fn$y)
S=splinefun(Fn$x, 1-pgamma(Fn$x,0.5,1))
# Create the Qauntille estimator
Qn=data.frame(x=Fn$y,y=Fn$x)
# Remove duplicates 
Qn=Qn[!duplicated(Qn$x), ]
# Create the Qauntille
Q=data.frame(x=Qn$x,y=qgamma(Qn$x,0.5,1))
#Q$y[Q$y== Inf] = 100
Q$y[Q$y== Inf] = min(Q$y[is.finite(Q$y)])
# Create the TTT function
TTT <- function(S,x) {
  d <- sapply(x, function(x) {
    integrate(S, lower = 0, upper = x)$value
  }) 
  d
}
# Create the TTT data
Ttt=data.frame(
  x=Qn$x,
  y=TTT(Sn,Qn$y),
  z=TTT(S,Q$y)
)
# Create the LCM of the TTT
ll = gcmlcm(sort(Ttt$x),Ttt$y, type="lcm")
# Create a smooth function
LCM_TTT <- splinefun(ll$x.knots, ll$y.knots,method = "natural")
Ttt=Ttt|>mutate(Lcm=LCM_TTT(Qn$x))
Ttt|>
  ggplot(aes(x=x))+
  geom_line(aes(y=y),color="red")+
  geom_line(aes(y=Lcm),color="blue")+
geom_line(aes(y=z),color="black")
#############################################################################
##Calcualte the Quantil
#############################################################################

Quantil=function(n,nsr){
  result1=numeric(m)
  for (i in 1:m) {
    sig.x=1
    sig.e=nsr*sig.x
    x <- rexp(n,1)
    u=ifelse(runif(n) > 0.5, 1, -1) * rexp(n, rate = 1/sig.e)
    w <- x+u
    # estimate the bandwidth with the bootstrap method with re-sampling
    bw <- bw.dnrd(w,sig=sig.e, error="laplacian")
    Fn <- DeconCdf(w,sig.e,error='laplacian',bw=bw)
    #Modify the CDF
    Fn$x=ifelse(Fn$x<0,0,Fn$x)
    # Create and the survival function and survival estimator
    Sn=splinefun(Fn$x, 1-Fn$y)
    # Create the Quantille data
    Qn=data.frame(x=Fn$y,y=Fn$x)
    # Remove duplicates 
    Qn=Qn[!duplicated(Qn$x), ]
    # Create the TTT estimator
    Ttt=data.frame(
      x=Qn$x,
      y=TTT(Sn,Qn$y)
    )
    # Create the LCM of the TTT 
    ll = gcmlcm(sort(Ttt$x),Ttt$y, type="lcm")
    # Create a smooth function of the LCM of TTT estimator
    LCM <- splinefun(ll$x.knots, ll$y.knots,method = "natural")
    # Add the LCM of the TTT estimator to the data 
    Ttt=Ttt|>mutate(Lcm=LCM(Qn$x))
    # The absolute deference between the estimator and its lcm
    result1[i] = max(abs(Ttt$Lcm-Ttt$y))
  }
  alpha=0.1
  Q=quantile(result1,1-alpha)
  return(Q)
}

#############################################################################"
# Calculate the Alternative
#############################################################################
Power=function(n,nsr){
  result2=numeric(m)
  for (i in 1:m) {
    X <- rbeta(n,0.5,1)
    varx=0.5/(1.5^2*2.5^2)
    sig.x=sqrt(varx)
    sig.e=nsr*sig.x
    e <- rlaplace(n, 0,sig.e)
    Y <- X+e
    ##########
    bw <- bw.dnrd(Y,sig=sig.e, error="laplacian")
    Fn <- DeconCdf(Y,sig.e,error='laplacian',bw=bw)
    #Modify the CDF
    Fn$x=ifelse(Fn$x<0,0,Fn$x)
    # Create and plot the survival functions
    Sn=splinefun(Fn$x, 1-Fn$y)
    # Create the Qauntille estimator
    Qn=data.frame(x=Fn$y,y=Fn$y)
    # Remove duplicates 
    Qn=Qn[!duplicated(Qn$x), ]
    # Create the TTT function
    Ttt=data.frame(
      x=Qn$x,
      y=TTT(Sn,Qn$x.1)
    )
    # Create the LCM of the TTT estimator
    ll = gcmlcm(sort(Ttt$x),Ttt$y, type="lcm")
    # Create a smooth function of the LCM of the TTT estimator
    LCM_TTT <- splinefun(ll$x.knots, ll$y.knots,method = "natural")
    Ttt=Ttt|>mutate(Lcm=LCM_TTT(Qn$x))
    # The absolute deference between the estimator and its lcm
    result2[i] = max(abs(Ttt$Lcm-Ttt$y))
  }
  mean(result2>Quantil(n,nsr))
}
Power(100,0.1)
