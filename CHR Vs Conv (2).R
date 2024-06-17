# Define the function to compute the log-normal density
PDF <- function(x, a, b ) {
  density_values <- dbeta(x, a, b)
  return(density_values)
}

Lap <- function(x, a, b) {
  density_values <- dlaplace(x, a, b)
  return(density_values)
}
# Define the convolution function
convolve_pdf <- function(t, a, b,al) {
  integrand <- function(x) {
    PDF(x, a, b) * Lap(t - x, 0, al)
  }
  integrate(integrand, -Inf, Inf)$value
}

# Example: Parameters for the Beta and Laplace distributions
a <- 0.5  # Beta parameter
b <- 1  # Beta parameter
al=0.5
# Compute the convolution for a sequence of values
s <- seq(-5, 5, length.out = 100)

Data1=data.frame(x=s,y=PDF(s,a,b),z=sapply(s, convolve_pdf, a = a, b = b,al=al))
ggplot(Data1, aes(x = x)) +
  geom_line(aes(y = y), color = "red")+
  geom_line(aes(y=z),color="blue")

###############################################
######### CDF
###############################################
# The True CDF

CDF <- function(x, a, b ) {
  density_values <- pbeta(x, a, b)
  return(density_values)
}

# Define the Convolution CDF (Integrating the PDF)
convolve_CDF <- function(t, a, b,al) {
  integrate(function(x) sapply(x, convolve_pdf, a=a, b=b,al=al), -Inf, t)$value
}
s <- seq(-5, 5, length.out = 100)

Data1=data.frame(x=s,y=CDF(s,a,b),z=sapply(s, convolve_CDF, a = a, b = b,al=al))
ggplot(Data1, aes(x = x)) +
  geom_line(aes(y = y), color = "red")+
  geom_line(aes(y=z),color="blue")
#################################################
######## HAZARD
#################################################
# Define the Hazard rate function
Hazard=function(x, a, b) {
  -log(1 - CDF(x,a,b))
}
# Define the convoluted Hazard
convolve_Hazard=function(x, a, b,al) {
  -log(1 - sapply(x, convolve_cdf, a = a, b= b,al))
}
# Example: Parameters for the distributions
s <- seq(-1.5, 1.5, length.out = 100)

Data3=data.frame(x=s,y=Hazard(s,a,b))
Data3|>
  ggplot(aes(x = x)) +
  geom_line(aes(y = y), color = "black",size=1)+
  labs(
    title = "The Cumulative Hazard function",
    x = "X-axis",
    y = "Cumulative Hazard function"
  ) +
  theme_minimal()+
xlim(0,1.05)

Data4=data.frame(x=s,y=sapply(s, convolve_Hazard, a = a, b= b,al=al))
Data4|>
  ggplot(aes(x = x)) +
  geom_line(aes(y = y), color = "black",size=1)+
  labs(
    title = "The Cumulative Hazard function",
    x = "X-axis",
    y = "Cumulative Hazard function"
  ) +
  theme_minimal()+
  xlim(-5,5)
