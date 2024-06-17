# Load necessary libraries
#library(splines)
#library(stats)
#library(LaplacesDemon)
#library(decon)
#library(numDeriv)

# Set parameters
nsr=0.5
n=50
X <- rexp(n,1)
varx=1
sig.x=sqrt(varx)
sig.e=nsr*sig.x
e <- rlaplace(n, 0,sig.e)
Y<- X+e
##########
bw <- bw.dboot2(Y,sig=sig.e, error="laplacian")
Fn <- DeconCdf(Y,sig.e,error='laplacian',bw=bw)
# Create the survival functions
Sn=splinefun(Fn$x, 1-Fn$y)
S <- splinefun(Fn$x, 1 - pexp(Fn$x,1))
#Modify the CDF
Fn$x=ifelse(Fn$x<0,0,Fn$x)
# Create the Qauntille estimator
Qn=data.frame(x=Fn$y,y=Fn$x)
# Remove duplicates 
Qn=Qn[!duplicated(Qn$x), ]

Q <- qexp(Qn$x,1)
Q[is.infinite(Q)] <- max(Q[is.finite(Q)])

# Create the TTT data

Ttt=data.frame(
  x=Qn$x,
  y=TTT(Sn,Qn$y),
  z=TTT(S,Q)
)
# remove Duplicates in order to apply LCM
Ttt=Ttt[!duplicated(Ttt$x), ]
# Create the LCM of the TTT
ll <- gcmlcm(sort(Ttt$x), Ttt$y, type = "lcm")
linear_interp <- approx(ll$x.knots, ll$y.knots,xout = seq(0,1,0.001))
LCM_TTT <- splinefun(linear_interp$x, linear_interp$y, method = "natural")
# mutate the data
Ttt=mutate(Ttt,t=LCM_TTT(Ttt$x))

Ttt|>
  ggplot(aes(x=x))+
  geom_line(aes(y=y),color="red")+
  geom_line(aes(y=z),color="black")+
  geom_line(aes(y=t),color="blue")

###########################################
############ HR Vs DHR ####################
###########################################
# The derivative of the LCM operator
DRV_LCM=splinefun(Ttt$x,LCM_TTT(Ttt$x,deriv=1),method = "natural")
# Define the Deconvolution Hazard rate estimator
DHR=splinefun(Fn$x,1/DRV_LCM(Fn$y),method = "natural")
Data=data.frame(
  x=Fn$x,
  y=1,
  z=DHR(Fn$x)
)
Data|>ggplot(aes(x=x))+
  geom_line(aes(y=y),color="red")+
  geom_line(aes(y=z),color="black")+
  xlim(c(0,5))+
  ylim(c(0,5))
###########################################
########### CHR Vs DCHR ###################
###########################################
# Define the Cumulative operator
Cum <- function(S, x) {
  sapply(x, function(x) {
    integrate(S, 0, x)$value
  })
}
Data=data.frame(
  x=Fn$x,
  y=Fn$x,
  z=Cum(DHR,Fn$x)
)

Data|>ggplot(aes(x=x))+
  geom_line(aes(y=y),color="red")+
  geom_line(aes(y=z),color="black")+
  xlim(c(0,5))+
  ylim(c(0,5))
###########################################
########### IHR F Vs F ####################
###########################################

S1 <- function(x) { x*0+1 }
# Example vector x
x <- seq(0, 10, 0.1)
Data=data.frame(
  x=x,
  y=E(Cum(S1, x)),
  z=E(Cum(DHR, x)),
  t=DeconCdf(Y,sig.e,x,error='laplacian',bw=bw)$y
)
Data|>ggplot(aes(x=x))+
  geom_line(aes(y=y),color="black")+
  geom_line(aes(y=z),color="red")+
  geom_line(aes(y=t),color="blue")
max(abs(Data$y-Data$z)) 
max(abs(Data$y-Data$t))

