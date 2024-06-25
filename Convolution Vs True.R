pdf <- function(x, weights, params) {
  # Unpack the parameters
  beta_params <- params$beta
  lognormal_params <- params$lognormal
  weibull_params <- params$weibull
  
  # Calculate individual densities
  beta_density <- dbeta(x, beta_params$shape1, beta_params$shape2)
  lognormal_density <- dlnorm(x, lognormal_params$meanlog, lognormal_params$sdlog)
  weibull_density <- dweibull(x, weibull_params$shape, weibull_params$scale)
  
  # Combine densities with weights
  mixture_density <- weights[1] * beta_density + 
    weights[2] * lognormal_density + 
    weights[3] * weibull_density
  
  return(mixture_density)
}

# Define the mixture CDF function
cdf <- function(t, weights, params) {
  sapply(t, function(t_val) {
    integrate(function(x) pdf(x, weights, params), 0, t_val)$value
  })
}


Hazard <- function(t, weights, params) {
  -log(1-cdf(t, weights, params))
}

# Example usage
weights <- c(0.5, 0.5, 0)  # Ensure the weights sum to 1
# Define parameters for each distribution
params <- list(
  beta = list(shape1 = 0.5, shape2 = 1),
  lognormal = list(meanlog = 0, sdlog = 1),
  weibull = list(shape = 1, scale = 1)
)

# Define the range of x values
x1 <- seq(0.1, 20, by = 0.01)

# plot the mix_pdf function
plot(x1,pdf(x1, weights, params),type="line",col=1)
# plot the mix_cdf function
plot(x1,cdf(x1, weights, params),type="line",col=1)

data=data.frame(x=x1,
                y=Hazard(x1, weights, params))
dat<- na.omit(data)
plot(data$x,data$y,type="line",col=1)
##################################################
######## Define the convolution function #########
##################################################
# Define the Laplace distribution function
lap <- function(x, mu, b) {
  return(dlaplace(x,mu,b))
}


convolve_pdf <- function(t, weights, params, mu, b) {
  integrand <- function(x) {
    pdf(x, weights, params) * lap(t - x, mu, b)
  }
  integrate(integrand, -Inf, Inf, subdivisions = 1000)$value
}
x2 <- seq(-10, 20, by = 0.1)
mu <- 0
b <- 0.1
convolve_pdf=splinefun(x2, sapply(x2, convolve_pdf, weights, params, mu, b))

convolve_cdf <- function(t) {
  integrand <- function(x) {
    convolve_pdf(x)
  }
  integrate(integrand, -Inf, t, subdivisions = 1000)$value
}

plot(x2,sapply(x2, convolve_pdf),type="line",col=2)
lines(x1,pdf(x1, weights, params),type="line",col=1)


plot(x2,sapply(x2, convolve_cdf),type="line",col=2)
lines(x1,cdf(x1, weights, params),type="line",col=1)

data=data.frame(x=x1,
                y=Hazard(x1, weights, params))
data<- na.omit(data)
plot(data$x,data$y,type="line",col=1)
lines(x2, -log(1-sapply(x2, convolve_cdf)),type="line",col=2)

