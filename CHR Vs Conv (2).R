# Define the function to compute the log-normal density
PDF <- function(x, a, b ) {
  density_values <- dlnorm(x, a, b)
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
  integrate(integrand, 0, +Inf,subdivisions = 1000)$value
}
# The True CDF
CDF <- function(x, a, b ) {
  density_values <- plnorm(x, a, b)
  return(density_values)
}
# Define the Convolution CDF (Integrating the PDF)
convolve_CDF <- function(t, a, b,al) {
  integrate(function(x) sapply(x, convolve_pdf, a=a, b=b,al=al), -Inf, t)$value
}
# Define the Hazard rate function
Hazard=function(x, a, b) {
  -log(1 - CDF(x,a,b))
}
# Define the convoluted Hazard
convolve_Hazard=function(x, a, b,al) {
  -log(1 - sapply(x, convolve_CDF, a = a, b= b,al))
}
#################################################
#### Examples ###################################
#################################################
a <- 0  # Beta parameter
b <- 1 # Beta parameter
al=0.5
# Compute the convolution for a sequence of values
s <- seq(0, 10, length.out = 100)

Data=data.frame(x=s,y=PDF(s,a,b),z=sapply(s, convolve_pdf, a = a, b = b,al=al))
ggplot(Data, aes(x = x)) +
  geom_line(aes(y = y), color = "red")+
  geom_line(aes(y=z),color="blue")


Data=data.frame(x=s,y=CDF(s,a,b),z=sapply(s, convolve_CDF, a = a, b = b,al=al))
ggplot(Data1, aes(x = x)) +
  geom_line(aes(y = y), color = "red")+
  geom_line(aes(y=z),color="blue")

# Example: Parameters for the distributions
#s <- seq(0, 1, length.out = 100)

NSR1=0.2
NSR2=0.4
NSR3=0.5

a <- 0  # Beta parameter
b <- 1  # Beta parameter

sig.x=sqrt((a*b)/((a+b)^2*(a+b+1)))
sig.e=sqrt(2)*0.1
sig.e/sig.x
al1=NSR1*sig.x/sqrt(2)
al2=NSR2*sig.x/sqrt(2)
al3=NSR3*sig.x/sqrt(2)

Data3=data.frame(x=s,y=Hazard(s,a,b))
Data3|>
  ggplot(aes(x = x)) +
  geom_line(aes(y = y), color = "black",size=1)+
  labs(
    title = "The Cumulative Hazard function",
    x = "X-axis",
    y = "Cumulative Hazard function"
  ) 


#+
#  theme_minimal()+
#xlim(0,1.05)

Data4=data.frame(
  x=s,
  y=Hazard(s,a,b),
  z=sapply(s, convolve_Hazard, a = a, b= b,al=al1),
  t=sapply(s, convolve_Hazard, a = a, b= b,al=al2),
  k=sapply(s, convolve_Hazard, a = a, b= b,al=al3)
  )

Data4|>
  ggplot(aes(x = x)) +
  geom_line(aes(y = y), color = "black",size=0.5)+
  geom_line(aes(y = z), color = "red",size=0.5)+
  geom_line(aes(y = t), color = "green",size=0.5)+
  geom_line(aes(y = k), color = "blue",size=0.5)

+

  
  labs(
    title = "The Cumulative Hazard function",
    x = "X-axis",
    y = "Cumulative Hazard function"
  ) +
  theme_minimal()+
  xlim(0,1)+
  ylim(0,5)
