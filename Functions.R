# ---------------------------------------------------------
# Title: Total Time on Test (TTT) and Generalized TTT (GTTT)
# Author: Benjrada
# Description: Functions for computing classical and generalized TTT
# ---------------------------------------------------------

# Classical TTT function (for Exponential Distribution)
TTT <- function(S, x) {
  sapply(x, function(x) integrate(S, lower = 0, upper = x)$value)
}

# Generalized TTT function
GTTT <- function(G, CDF, x) {
  if (G == "U") {
    return(x)  # Uniform case: TTT(x) = x
  }
  
  # Define the function to integrate based on G
  integral_func <- switch(G,
                          "EXP" = function(s) (1 - CDF(s)),  # Exponential case
                          "LL"  = function(s) (1 - CDF(s))^2,  # Log-logistic case
                          stop("Unsupported distribution. Use 'U', 'EXP', or 'LL'.")
  )
  
  sapply(x, function(t) integrate(integral_func, lower = 0, upper = t)$value)
}

# ---------------------------------------------------------
# Example Usage
# ---------------------------------------------------------

# Define the CDF of the Log-Logistic(1,1) distribution
CDF_loglogistic <- function(x) pllogis(x, 1, 1)

# Define the CDF of the Exponential(1) distribution
CDF_exponential <- function(x) pexp(x, 1)

# Generate x values
x_vals <- seq(0.01, 0.99, length.out = 100)

# Compute quantiles
y_vals_loglogistic <- qllogis(x_vals, 1, 1)
y_vals_exponential <- qexp(x_vals, 1)

# Compute GTTT for different cases and plot results
plot_GTTT <- function(G, CDF, y_vals, title) {
  result <- GTTT(G = G, CDF = CDF, x = y_vals)
  plot(x_vals, result, type = "l", col = "blue", lwd = 2,
       main = title, xlab = "x", ylab = "GTTT(x)")
}

# Plot for Uniform Case
plot_GTTT(G = "U", CDF = CDF_loglogistic, y_vals = y_vals_loglogistic, 
          title = "GTTT for Uniform Case")

# Plot for Log-Logistic Case
plot_GTTT(G = "LL", CDF = CDF_loglogistic, y_vals = y_vals_loglogistic, 
          title = "GTTT for Log-Logistic Case")

# Plot for Exponential Case
plot_GTTT(G = "EXP", CDF = CDF_exponential, y_vals = y_vals_exponential, 
          title = "GTTT for Exponential Case")
