# Exponentiated Weibull PDF and CDF
PDF <- function(x, a, b) {
  (1/beta(a,b))*x^(a-1)*(1-b)^(b-1)*(x>=0)
}

CDF <- function(x, a, b) {
  (1 - exp(-x^a))^b*(x>=0)
}
Hazard=function(x, a, b) {
  -log(1 - CDF(x,a,b))
}
# Laplace PDF and CDF
Lap <- function(x,a,b) {
  1/(2*b) * exp(-abs(x-a)/b)
}

s=seq(0,6,0.01)
a=0.5
b=1
plot(s,PDF(s,a,b),type="line")
plot(s,CDF(s,a,b),type="line")
plot(s,Hazard(s,a,b),type="line")


# Define the convolution function for the PDFs
conv_pdf <- function(t, a, b) {
  integrand <- function(x) {
    Lap(t-x,0,1) * PDF(x, a, b)
  }
  integrate(integrand, -50, 50)$value
}
# Define the convolved CDF
conv_cdf <- function(t, a, b) {
  integrate(function(x) sapply(x, conv_pdf, a=a, b=b), -50, t)$value
}

conv_Hazard=function(x, a, b) {
  -log(1 - conv_cdf(x,a,b))
}
# Example: Parameters for the distributions

conv_pdf_values <- sapply(s, conv_pdf, a = a, b= b)
conv_cdf_values <- sapply(s, conv_cdf, a = a, b= b)
conv_hazard_values <- sapply(s, conv_Hazard, a = a, b= b)

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