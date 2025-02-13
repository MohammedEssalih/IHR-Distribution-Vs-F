log_normal_pdf <- function(x, meanlog = 0, sdlog =1) {
  dlnorm(x, meanlog, sdlog)
}
log_normal_pdf(1)
# Define the function to compute the Laplace density
laplace_pdf <- function(x, location = 0, scale = 0.1) {
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
scale <- 0.1
location <- 0

# Generate a sequence of values for evaluation
y <- seq(0, 10, length = 1000)

# Compute the convolution result
z1 <- convolve_pdf(y, meanlog = meanlog, sdlog = sdlog, scale = scale)

# Plot the convolution result
plot(y, z, type = "l", lwd = 2, col = "blue",
     main = "Density of the Convolution of Log-Normal and Laplace",
     xlab = "x", ylab = "Density")

# Verify with bayesmeta::convolve function
conv <- convolve(log_normal_pdf, laplace_pdf)

# Add the bayesmeta::convolve result to the plot
z2=conv$density(y)
lines(y, z2, type = "l", lwd = 2, col = "red")

# Add a legend
legend("topright", legend = c("Custom Convolution", "bayesmeta::convolve"),
       col = c("blue", "red"), lwd = 2)

# Define the Convolution CDF (Integrating the PDF)
convolve_cdf <- Vectorize(function(t, meanlog, sdlog, scale) {
  integrate(function(x) sapply(x, convolve_pdf, meanlog = meanlog, sdlog = sdlog, scale = scale), 0, t)$value
})
d1=convolve_cdf(y, meanlog = meanlog, sdlog = sdlog, scale = scale)
d2=conv$cdf(y)

plot(y, d2, type = "l", lwd = 2, col = "blue",
     main = "CDF of the Convolution of Log-Normal and Laplace",
     xlab = "x", ylab = "Cumulative Probability")
lines(y, d1, type = "l", lwd = 2, col = "red")

plot(y, -log(1-d2), type = "l", lwd = 2, col = "blue",
     main = "CDF of the Convolution of Log-Normal and Laplace",
     xlab = "x", ylab = "Cumulative Probability")
lines(y, -log(1 - log_normal_cdf(y, meanlog, sdlog)), type = "l", lwd = 2, col = "red")

# Define the Hazard rate function
hazard_rate <- function(x, meanlog, sdlog) {
  -log(1 - log_normal_cdf(x, meanlog, sdlog))
}


# Define the convoluted Hazard rate function
convolve_hazard_rate <- function(x, meanlog, sdlog, scale) {
  -log(1 - sapply(x, convolve_cdf, meanlog = meanlog, sdlog = sdlog, scale = scale))
}

# Example parameters
meanlog <- 0
sdlog <- 1
scale <- 0.5

# Compute the convolution for a sequence of values
s <- seq(0, 10, length.out = 100)

# Create a data frame for PDF comparison
pdf_data <- data.frame(
  x = s,
  log_normal_pdf = log_normal_pdf(s, meanlog, sdlog),
  convoluted_pdf = sapply(s, convolve_pdf, meanlog = meanlog, sdlog = sdlog, scale = scale)
)

# Plot the PDFs
ggplot(pdf_data, aes(x = x)) +
  geom_line(aes(y = log_normal_pdf), color = "red") +
  geom_line(aes(y = convoluted_pdf), color = "blue") +
  labs(
    title = "Log-Normal PDF vs Convoluted PDF",
    x = "X-axis",
    y = "Density"
  ) +
  theme_minimal()

# Create a data frame for CDF comparison
cdf_data <- data.frame(
  x = s,
  log_normal_cdf = log_normal_cdf(s, meanlog, sdlog),
  convoluted_cdf = sapply(s, convolve_cdf, meanlog = meanlog, sdlog = sdlog, scale = scale)
)

# Plot the CDFs
ggplot(cdf_data, aes(x = x)) +
  geom_line(aes(y = log_normal_cdf), color = "red") +
  geom_line(aes(y = convoluted_cdf), color = "blue") +
  labs(
    title = "Log-Normal CDF vs Convoluted CDF",
    x = "X-axis",
    y = "Cumulative Probability"
  ) +
  theme_minimal()

# Example: Parameters for the distributions
nsr1 <- 0.2
nsr2 <- 0.4
nsr3 <- 0.5

meanlog <- 0.7
sdlog <- 1

sig_x <- sqrt((meanlog * sdlog) / ((meanlog + sdlog)^2 * (meanlog + sdlog + 1)))
sig_e <- sqrt(2)
scale1 <- nsr1 * sig_x / sqrt(2)
scale2 <- nsr2 * sig_x / sqrt(2)
scale3 <- nsr3 * sig_x / sqrt(2)

# Create a data frame for Hazard rate comparison
hazard_data <- data.frame(
  x = s,
  log_normal_hazard = hazard_rate(s, meanlog, sdlog),
  convoluted_hazard1 = sapply(s, convolve_hazard_rate, meanlog = meanlog, sdlog = sdlog, scale = scale1),
  convoluted_hazard2 = sapply(s, convolve_hazard_rate, meanlog = meanlog, sdlog = sdlog, scale = scale2),
  convoluted_hazard3 = sapply(s, convolve_hazard_rate, meanlog = meanlog, sdlog = sdlog, scale = scale3)
)

# Plot the Hazard rates
ggplot(hazard_data, aes(x = x)) +
  geom_line(aes(y = log_normal_hazard), color = "black", size = 0.5) +
  geom_line(aes(y = convoluted_hazard1), color = "red", size = 0.5) +
  geom_line(aes(y = convoluted_hazard2), color = "green", size = 0.5) +
  geom_line(aes(y = convoluted_hazard3), color = "blue", size = 0.5) +
  labs(
    title = "Cumulative Hazard Function",
    x = "X-axis",
    y = "Cumulative Hazard Function"
  ) +
  theme_minimal() +
  xlim(0, 1) +
  ylim(0, 5)

