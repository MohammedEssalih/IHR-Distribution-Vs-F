#############################################
#### Hazard Function & CDF for Pareto Distribution with NSR
#############################################

# Load required packages
library(ggplot2)
library(LaplacesDemon)  # For dlaplace function

# Define the CDF of the Pareto distribution
cdf_pareto <- function(x, alpha, x_m) {
  ppareto(x, scale = x_m, shape = alpha)
}

# Define the hazard function of the Pareto distribution
hazard_pareto <- function(x, alpha, x_m) {
  ifelse(x >= x_m, alpha / x, 0)
}

# Define parameters
x_m <- 1  # Scale parameter (minimum value)
alpha <- 2  # Shape parameter

# 1. Mean and Variance of Pareto Distribution
mu <- if (alpha > 1) x_m * alpha / (alpha - 1) else NA
var_pareto <- if (alpha > 2) (x_m^2 * alpha) / ((alpha - 1)^2 * (alpha - 2)) else NA

# 2. Noise-to-signal ratios (NSRs)
nsr_values <- c(0.1, 0.2, 0.5)
sigma_x <- sqrt(var_pareto)
b_values <- nsr_values * sigma_x / sqrt(2)  # Scale parameter for Laplace error

# Define Laplace PDF
pdf2 <- function(x, b) dlaplace(x, location = 0, scale = b)

convolve_cdf <- Vectorize(function(t, b) {
  integrate(function(x) pdf2(t-x, b) * cdf_pareto(x, alpha, x_m), 1, 10, subdivisions = 1000)$value
})

# Generate evaluation points
y <- seq(0.5, 5, length.out = 1000)
d <- cdf_pareto(y, alpha, x_m)

# Compute convoluted CDFs for different NSRs
cdf_list <- lapply(b_values, function(b) convolve_cdf(y, b))
names(cdf_list) <- c("10%", "20%", "50%")

# Compute hazard functions
H <- -log(1 - d)

# Compute hazard functions for contaminated distributions
CH_values <- lapply(cdf_list, function(d) -log(1 - d))

# Remove NA values for consistency
valid_indices <- complete.cases(H, CH_values[[1]], CH_values[[2]], CH_values[[3]])
y <- y[valid_indices]
H <- H[valid_indices]
CH_values <- lapply(CH_values, function(CH) CH[valid_indices])

# Create DataFrame for plotting
Data2 <- data.frame(
  x = rep(y, times = 4),
  y = c(H, CH_values[[1]], CH_values[[2]], CH_values[[3]]),
  group = rep(c("0%", "10%", "20%", "50%"), each = length(y))
)

# Plot hazard functions
plot_hazard <- ggplot(Data2, aes(x = x, y = y, linetype = group)) +
  geom_line() +
  labs(y = "Hazard function", x = "x", linetype = "NSR") +
  scale_linetype_manual(values = c("0%" = "solid", "10%" = "dashed", "20%" = "twodash", "50%" = "dotted")) +
  theme_minimal() +
  ylim(0, 5) +
  xlim(0.5, 5)

# Plot CDFs
Data3 <- data.frame(
  x = rep(y, times = 4),
  y = c(d, cdf_list[[1]], cdf_list[[2]], cdf_list[[3]]),
  group = rep(c("0%", "10%", "20%", "50%"), each = length(y))
)

plot_cdf <- ggplot(Data3, aes(x = x, y = y, linetype = group)) +
  geom_line() +
  labs(y = "CDF", x = "x", linetype = "NSR") +
  scale_linetype_manual(values = c("0%" = "solid", "10%" = "dashed", "20%" = "twodash", "50%" = "dotted")) +
  theme_minimal() +
  ylim(0, 1) +
  xlim(0.5, 5)

# Display plots
print(plot_hazard)
print(plot_cdf)
