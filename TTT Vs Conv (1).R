# Define the Beta PDF that handles vector inputs using sapply
PDF <- function(x, a, b) {
  sapply(x, function(xi) {
    if (xi >= 0 && xi <= 1) {
      (xi^(a - 1) * (1 - xi)^(b - 1)) / beta(a, b)
    } else {
      0
    }
  })
}


CDF <- function(t, a, b) {
  integrate(function(x) PDF(x, a, b), 0, t)$value
}

Hazard=function(x, a, b) {
  -log(1 - CDF(x,a,b))
}
# Example: Parameters for the Beta distribution
a <- 0.5
b <- 1
s <- seq(-10, 10, length.out = 100)
pdf_values <- sapply(s, PDF, a = a, b = b)
cdf_values <- sapply(s, CDF, a = a, b = b)
hazard_values <- sapply(s, Hazard, a = a, b = b)
# Plot the PDF
plot(t_values, pdf_values, type = "l", main = "CDF of Beta Distribution", xlab = "t", ylab = "CDF")
# Plot the CDF
plot(t_values, cdf_values, type = "l", main = "CDF of Beta Distribution", xlab = "t", ylab = "CDF")
# Plot the Hazard function
plot(s, hazard_values, type = "l", main = "Cumulative Hazard rate of Beta Distribution", xlab = "t", ylab = "Hazard")

# Laplace
Lap <- function(x,a,b) {
  1/(2*b) * exp(-abs(x-a)/b)
}
# Define the convoluted PDF
conv_pdf <- function(t, a, b) {
  integrand <- function(x) {
    Lap(t-x,0,1) * PDF(x, a, b)
  }
  integrate(integrand, -50, 50)$value
}
# Define the convoluted CDF
conv_cdf <- function(t, a, b) {
  integrate(function(x) sapply(x, conv_pdf, a=a, b=b), -50, t)$value
}
# Define the convoluted Hazard
conv_Hazard=function(x, a, b) {
  -log(1 - conv_cdf(x,a,b))
}
# Example: Parameters for the distributions

conv_pdf_values <- sapply(s, conv_pdf, a = a, b= b)
conv_cdf_values <- sapply(s, conv_cdf, a = a, b= b)
conv_hazard_values <- sapply(s, conv_Hazard, a = a, b= b)
plot(s,conv_pdf_values)

Data1=data.frame(x=s,y=PDF(s,a,b),z=sapply(s, conv_pdf, a = a, b= b))
ggplot(Data1, aes(x = x)) +
  geom_line(aes(y = y), color = "red")+
  geom_line(aes(y=z),color="blue")


Data2=data.frame(x=s,y=CDF(s,a,b),z=sapply(s, conv_cdf, a = a, b= b))
ggplot(Data2, aes(x = x)) +
  geom_line(aes(y = y), color = "red")+
  geom_line(aes(y=z),color="blue")

Data3=data.frame(x=s,y=Hazard(s,a,b),z=sapply(s, conv_Hazard, a = a, b= b))
ggplot(Data3, aes(x = x)) +
  geom_line(aes(y = y), color = "black")+
  geom_line(aes(y=z),color="red")
