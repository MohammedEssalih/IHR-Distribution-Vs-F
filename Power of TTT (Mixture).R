#library(fdrtool)
#library("VGAM")
#library("decon")
#library("tidyverse")
#library(cobs)

m=100
nsr=0.8
n=100
########################
p1=0.5
p2=1
##############
X <- rexp(n,1)
sig.e=nsr
e <- rlaplace(n, 0,sig.e)
Y<- X+e
##########
bw <- bw.dboot2(Y,sig=sig.e, error="laplacian")
Fn <- DeconCdf(Y,sig.e,error='laplacian',bw=bw)
# Create the survival function
Sn=splinefun(Fn$x, 1-Fn$y)
# Restrict the CDF to only positive interval
Fn$x=ifelse(Fn$x<0,0,Fn$x)
# Create the Qauntille estimator
Qn=data.frame(x=Fn$y,y=Fn$x)
# Remove duplicates 
Qn=Qn[!duplicated(Qn$x), ]
Qn=Qn[!duplicated(Qn$y), ]
# Create the Qauntile function
D=splinefun(Qn$x, Qn$y)
# Calculate the mean of the underlying data (X) (after deconvolution)
mu=integrate(D, lower = 0, upper = 1,subdivisions=2000)$value
# Create the TTT data
Ttt=data.frame(
  x=Qn$x,
  y=TTT(Sn,Qn$y)/mu
)
# Create the LCM of the TTT
ll = gcmlcm(sort(Ttt$x),Ttt$y, type="lcm")
linear_interp <- approx(ll$x.knots, ll$y.knots, xout = seq(0, 1, by = 0.001))
LCM_TTT <- splinefun(linear_interp$x, linear_interp$y,method = "natural")
#LCM_TTT <- splinefunH(ll$x.knots, ll$y.knots,m=ll$slope.knots)
Ttt=Ttt|>mutate(Lcm=LCM_TTT(Qn$x))
Ttt|>
  ggplot(aes(x=x))+
  geom_line(aes(y=y),color="red")+
  geom_line(aes(y=Lcm),color="blue")+
  xlab(max(Ttt$Lcm-Ttt$y))
#############################################################################"
# Calculate the critical value from the null (Least favorable case)
#############################################################################
R1=numeric(m)
for (i in 1:m) {
  X <- rexp(n,1)
  sig.e=nsr
  e <- rlaplace(n, 0,sig.e)
  Y<- X+e
  ##########
  bw <- bw.dboot2(Y,sig=sig.e, error="laplacian")
  Fn <- DeconCdf(Y,sig.e,error='laplacian',bw=bw)
  # Create the survival function
  Sn=splinefun(Fn$x, 1-Fn$y)
  # Restrict the CDF to only positive interval
  Fn$x=ifelse(Fn$x<0,0,Fn$x)
  # Create the Qauntille estimator
  Qn=data.frame(x=Fn$y,y=Fn$x)
  # Remove duplicates 
  Qn=Qn[!duplicated(Qn$x), ]
  Qn=Qn[!duplicated(Qn$y), ]
  # Create the Qauntile function
  D=splinefun(Qn$x, Qn$y)
  # Calculate the mean of the underlying data (X) (after deconvolution)
  mu=integrate(D, lower = 0, upper = 1,subdivisions=2000)$value
  # Create the TTT data
  Ttt=data.frame(
    x=Qn$x,
    y=TTT(Sn,Qn$y)/mu
  )
  # Create the LCM of the TTT
  ll = gcmlcm(sort(Ttt$x),Ttt$y, type="lcm")
  linear_interp <- approx(ll$x.knots, ll$y.knots, xout = seq(0, 1, by = 0.001))
  LCM_TTT <- splinefun(linear_interp$x, linear_interp$y,method = "natural")
  #LCM_TTT <- splinefunH(ll$x.knots, ll$y.knots,m=ll$slope.knots)
  Ttt=Ttt|>mutate(Lcm=LCM_TTT(Qn$x))
  R1[i] = max(Ttt$Lcm-Ttt$y)
}

alpha=0.1
q=quantile(R1,1-alpha)
###############################################################
n=1000
X <- 0.6*rbeta(n,0.5,0.5)+0.4*rweibull(n,2.5,0.75)
varx=var(X)
sig.x=sqrt(varx)
sig.e=0.005
e <- rlaplace(n, 0,sig.e)
Y<- X+e
##########
bw <- bw.dboot2(Y,sig=sig.e, error="laplacian")
Fn <- DeconCdf(Y,sig.e,error='laplacian',bw=bw)
# Create the survival function
Sn=splinefun(Fn$x, 1-Fn$y)
#Restrict the CDF to the positive interval
Fn$x=ifelse(Fn$x<0,0,Fn$x)
# Create the Qauntille estimator
Qn=data.frame(x=Fn$y,y=Fn$x)
# Remove duplicates 
Qn=Qn[!duplicated(Qn$x), ]
Qn=Qn[!duplicated(Qn$y), ]
D=splinefun(Qn$x,Qn$y)
# Calculate the mean of the underlying variable (After deconvolution)
mu=integrate(D,0,1,subdivisions = 2000)$value
# Create the TTT data
Ttt=data.frame(
  x=Qn$x,
  y=TTT(Sn,Qn$y)/mu
)
# Create the LCM of the TTT
ll = gcmlcm(sort(Ttt$x),Ttt$y, type="lcm")
linear_interp <- approx(ll$x.knots, ll$y.knots, xout = seq(0, 1, by = 0.001))
LCM_TTT <- splinefun(linear_interp$x, linear_interp$y,method = "natural")
Ttt=Ttt|>mutate(Lcm=LCM_TTT(Qn$x))
Ttt|>
  ggplot(aes(x=x))+
  geom_line(aes(y=y),color="red")+
  geom_line(aes(y=Lcm),color="blue")+
  xlab(max(Ttt$Lcm-Ttt$y))
#############################################################################"
# Calculate the Alternative
#############################################################################
#Power=function(n,nsr){

R2=numeric(m)
for (i in 1:m) {
  X <- 0.6*rbeta(n,0.5,0.5)+0.4*rweibull(n,2.5,0.75)
  varx=var(X)
  sig.x=sqrt(varx)
  sig.e=nsr*sig.x
  e <- rlaplace(n, 0,sig.e)
  Y<- X+e
  ##########
  bw <- bw.dboot2(Y,sig=sig.e, error="laplacian")
  Fn <- DeconCdf(Y,sig.e,error='laplacian',bw=bw)
  # Create the survival function
  Sn=splinefun(Fn$x, 1-Fn$y)
  #Restrict the CDF to the positive interval
  Fn$x=ifelse(Fn$x<0,0,Fn$x)
  # Create the Qauntille estimator
  Qn=data.frame(x=Fn$y,y=Fn$x)
  # Remove duplicates 
  Qn=Qn[!duplicated(Qn$x), ]
  Qn=Qn[!duplicated(Qn$y), ]
  D=splinefun(Qn$x,Qn$y)
  # Calculate the mean of the underlying variable (After deconvolution)
  mu=integrate(D,0,1,subdivisions = 2000)$value
  # Create the TTT data
  Ttt=data.frame(
    x=Qn$x,
    y=TTT(Sn,Qn$y)/mu
  )
  # Create the LCM of the TTT
  ll = gcmlcm(sort(Ttt$x),Ttt$y, type="lcm")
  linear_interp <- approx(ll$x.knots, ll$y.knots, xout = seq(0, 1, by = 0.001))
  LCM_TTT <- splinefun(linear_interp$x, linear_interp$y,method = "natural")
  Ttt=Ttt|>mutate(Lcm=LCM_TTT(Qn$x))
  R2[i] = max(Ttt$Lcm-Ttt$y)
}
mean(R2>q)


