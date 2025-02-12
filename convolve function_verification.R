# Load necessary libraries
install.packages("bayesmeta")
library(bayesmeta)
library(VGAM)
library(ggplot2)
library(statmod)  # For dlaplace function

# Define the function to compute the log-normal density
log_normal_pdf <- function(x, meanlog = meanlog, sdlog = sdlog) {
  dlnorm(x, meanlog, sdlog)
}
log_normal_pdf(1)
# Define the function to compute the Laplace density
laplace_pdf <- function(x, location = location, scale = scale) {
  dlaplace(x, location, scale)
}

# Define the convolution function for PDFs
convolve_pdf <- Vectorize(function(t, meanlog=meanlog, sdlog=sdlog, scale=scale) {
  integrand <- function(x) {
    log_normal_pdf(x, meanlog, sdlog) * laplace_pdf(t - x, 0, scale)
  }
  integrate(integrand, -Inf, Inf, subdivisions = 1000)$value
})

# Set parameters
meanlog <- 0
sdlog <- 1
scale <- 0.5
location <- 0

# Generate a sequence of values for evaluation
y <- seq(0, 5, length = 1000)

# Compute the convolution result
z <- convolve_pdf(y, meanlog = meanlog, sdlog = sdlog, scale = scale)

# Plot the convolution result
plot(y, z, type = "l", lwd = 2, col = "blue",
     main = "Density of the Convolution of Log-Normal and Laplace",
     xlab = "x", ylab = "Density")

# Verify with bayesmeta::convolve function
conv <- convolve(log_normal_pdf, laplace_pdf)

# Add the bayesmeta::convolve result to the plot
lines(y, conv$density(y), type = "l", lwd = 2, col = "red")

# Add a legend
legend("topright", legend = c("Custom Convolution", "bayesmeta::convolve"),
       col = c("blue", "red"), lwd = 2)

### The plots are identical, so the function convolve from the package {bayesmeta} is correct