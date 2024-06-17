#library("decon")
#library("tidyverse")
# Plote and ensure
m=500
nsr=0.5
n=100
X <- rexp(n,1)
varx=1
sig.x=sqrt(varx)
sig.e=nsr*sig.x
e <- rlaplace(n, 0,sig.e)
Y<- X+e
##########
bw <- bw.dboot2(Y,sig=sig.e, error="laplacian")
Fn <- DeconCdf(Y,sig.e,error='laplacian',bw=bw)
# Create the survival function
Sn=splinefun(Fn$x, 1-Fn$y)
#Modify the CDF
Fn$x=ifelse(Fn$x<0,0,Fn$x)
# Create the Qauntille estimator
Qn=data.frame(x=Fn$y,y=Fn$x)
# Remove duplicates 
Qn=Qn[!duplicated(Qn$x), ]
# Create the TTT data
Ttt=data.frame(
  x=Qn$x,
  y=TTT(Sn,Qn$y)
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
R1=numeric(m)
for (i in 1:m) {
  X <- rexp(n,1)
  varx=1
  sig.x=sqrt(varx)
  sig.e=nsr*sig.x
  e <- rlaplace(n, 0,sig.e)
  Y<- X+e
  ##########
  bw <- bw.dboot2(Y,sig=sig.e, error="laplacian")
  Fn <- DeconCdf(Y,sig.e,error='laplacian',bw=bw)
  #Modify the CDF
  Sn=splinefun(Fn$x, 1-Fn$y)
  Fn$x=ifelse(Fn$x<0,0,Fn$x)
  # Create and plot the survival functions
  
  # Create the Qauntille estimator
  Qn=data.frame(x=Fn$y,y=Fn$x)
  # Remove duplicates 
  Qn=Qn[!duplicated(Qn$x), ]
  # Create the TTT data
  Ttt=data.frame(
    x=Qn$x,
    y=TTT(Sn,Qn$y)
  )
  # Create the LCM of the TTT
  ll = gcmlcm(sort(Ttt$x),Ttt$y, type="lcm")
  linear_interp <- approx(ll$x.knots, ll$y.knots, xout = seq(0, 1, by = 0.001))
  LCM_TTT <- splinefun(linear_interp$x, linear_interp$y,method = "natural")
  Ttt=Ttt|>mutate(Lcm=LCM_TTT(Qn$x))
  
  R1[i] = max(Ttt$Lcm-Ttt$y)
}
alpha=0.1
q=quantile(R1,1-alpha)
###############################################################
x <- rlnorm(n,0.5,1)
varx=4.67078
#varx=0.5/(1.5^2*2.5)
sig.x=sqrt(varx)
sig.e=nsr*sig.x
u <- rlaplace(n, 0,sig.e)
w <- x+u
##########
bw <- bw.dboot2(w,sig=sig.e, error="laplacian")
Fn <- DeconCdf(w,sig.e,error='laplacian',bw=bw)
#Modify the CDF
Sn=splinefun(Fn$x, 1-Fn$y)
Fn$x=ifelse(Fn$x<0,0,Fn$x)
# Create the Qauntille estimator
Qn=data.frame(x=Fn$y,y=Fn$x)
# Remove duplicates 
Qn=Qn[!duplicated(Qn$x), ]
# Create the TTT data
Ttt=data.frame(
  x=Qn$x,
  y=TTT(Sn,Qn$y)
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

result2=numeric(m)
for (i in 1:m) {
  x <- rlnorm(n,0.5,1)
  varx=4.67078
  #varx=0.5/(1.5^2*2.5)
  sig.x=sqrt(varx)
  sig.e=nsr*sig.x
  u <- rlaplace(n, 0,sig.e)
  w <- x+u
  ##########
  bw <- bw.dboot2(w,sig=sig.e, error="laplacian")
  Fn <- DeconCdf(w,sig.e,error='laplacian',bw=bw)
  #Modify the CDF
  Sn=splinefun(Fn$x, 1-Fn$y)
  Fn$x=ifelse(Fn$x<0,0,Fn$x)
  # Create the Qauntille estimator
  Qn=data.frame(x=Fn$y,y=Fn$x)
  # Remove duplicates 
  Qn=Qn[!duplicated(Qn$x), ]
  # Create the TTT data
  Ttt=data.frame(
    x=Qn$x,
    y=TTT(Sn,Qn$y)
  )
  # Create the LCM of the TTT
  ll = gcmlcm(sort(Ttt$x),Ttt$y, type="lcm")
  linear_interp <- approx(ll$x.knots, ll$y.knots, xout = seq(0, 1, by = 0.001))
  LCM_TTT <- splinefun(linear_interp$x, linear_interp$y,method = "natural")
  Ttt=Ttt|>mutate(Lcm=LCM_TTT(Qn$x))
  result2[i] = max(Ttt$Lcm-Ttt$y)
}
mean(result2>=q)
