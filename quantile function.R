m=100
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

TTT <- function(S,x) {
  d <- sapply(x, function(x) {
    integrate(S, lower = 0, upper = x)$value
  }) 
  d
}

Quantil(1000,0.1)
