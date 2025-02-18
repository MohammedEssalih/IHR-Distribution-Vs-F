# ---------------------------------------------------------
# Title: Total Time on Test (TTT) and Generalized TTT (GTTT)
# Author: Ben Jrada Mohammed Es-salih
# Description: Functions for computing classical and generalized TTT
# ---------------------------------------------------------

# Classical TTT function (for Exponential Distribution)
TTT <- function(S, x) {
  sapply(x, function(t) integrate(S, lower = 0, upper = t)$value)
}

# Generalized TTT function
Gttt <- function(G, Q, PDF, x) {
  if (G == "U") {
    return(dunif(qunif(x)) / PDF(Q(x)))  
  }
  
  if (G == "EXP") {
    return(dexp(qexp(x)) / PDF(Q(x)))  
  }
  
  if (G == "LL") {
    return(dllogis(qllogis(x, 1, 1), 1, 1) / PDF(Q(x)))  
  }
  
  stop("Unsupported distribution. Use 'U', 'EXP', or 'LL'.")
}

# Quantile and PDF functions for Uniform distribution
Q_unif <- function(x) qunif(x)
PDF_unif <- function(x) dunif(x)

# Example usage of Gttt function
Gttt("U", Q_unif, PDF_unif, x = seq(0.01, 0.99, length.out = 100))

# ---------------------------------------------------------
# Generalized TTT function
# ---------------------------------------------------------

GTTT <- function(G, CDF, x) {
  if (G == "U") {
    return(x)  # Uniform case: TTT(x) = x
  }
  
  integral_func <- switch(G,
                          "EXP" = function(s) (1 - CDF(s)),  
                          "LL"  = function(s) (1 - CDF(s))^2,  
                          stop("Unsupported distribution. Use 'U', 'EXP', or 'LL'.")
  )
  
  sapply(x, function(t) integrate(integral_func, lower = 0, upper = t)$value)
}

# ---------------------------------------------------------
# Example Usage
# ---------------------------------------------------------

# Define CDF functions
CDF_loglogistic <- function(x) pllogis(x, 1, 1)
CDF_exponential <- function(x) pexp(x, 2)

# Generate x values
x_vals <- seq(0.01, 0.99, length.out = 100)

# Compute quantiles
y_vals_loglogistic <- qllogis(x_vals, 1, 1)
y_vals_exponential <- qexp(x_vals, 2)

# Function to plot GTTT results
plot_GTTT <- function(G, CDF, y_vals, title) {
  result <- GTTT(G = G, CDF = CDF, x = y_vals)
  plot(x_vals, result, type = "l", col = "blue", lwd = 2,
       main = title, xlab = "x", ylab = "GTTT(x)")
}

# Plot for different distributions
plot_GTTT("U", CDF_loglogistic, y_vals_loglogistic, "GTTT for Uniform Case")
plot_GTTT("LL", CDF_loglogistic, y_vals_loglogistic, "GTTT for Log-Logistic Case")
plot_GTTT("EXP", CDF_exponential, y_vals_exponential, "GTTT for Exponential Case")

# ---------------------------------------------------------
# Generalized Hazard Function (GHF) and Generalized Hazard Rate (GHR)
# ---------------------------------------------------------

# Generalized Hazard Function
GHF <- function(G, CDF, x) {
  if (G == "U") {
    return(CDF(x))  
  }
  
  if (G == "EXP") {
    return(qexp(CDF(x)))  
  }
  
  if (G == "LL") {
    return(qllogis(CDF(x), 1, 1))  
  }
  
  stop("Unsupported distribution. Use 'U', 'EXP', or 'LL'.")
}

# Generalized Hazard Rate (GHR)
GHR <- function(G, CDF, PDF, x) {
  if (G == "U") {
    return(PDF(x))  
  }
  
  if (G == "EXP") {
    return(PDF(x) / dexp(qexp(CDF(x))))  
  }
  
  if (G == "LL") {
    return(PDF(x) / dllogis(qllogis(CDF(x), 1, 1), 1, 1))  
  }
  
  stop("Unsupported distribution. Use 'U', 'EXP', or 'LL'.")
}

# Define CDF and PDF functions for Log-Logistic(1,1)
CDF_loglogistic <- function(x) pllogis(x, 1, 1)
PDF_loglogistic <- function(x) dllogis(x, 1, 1)

# Compute and plot the difference GHF("LL", CDF, x) - x
x_vals <- seq(0.01, 5, length.out = 100)
ghf_vals <- sapply(x_vals, function(x) GHF("LL", CDF_loglogistic, x))

plot(x_vals, ghf_vals, type = "l", col = "red", lwd = 2,
     main = "GHF - x for Log-Logistic Case", xlab = "x", ylab = "GHF(x) - x")
grid()
