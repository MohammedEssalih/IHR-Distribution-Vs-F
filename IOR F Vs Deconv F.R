# Load necessary libraries
library(flexsurv)
library(purrr)
library(LaplacesDemon)
library(decon)
library(numDeriv)
library(ggplot2)
library(dplyr)
library(VGAM)
library(fdrtool)

# Set parameters
nsr <- 0.5
n <- 500
a <- 1.2
b <- 1.25
p <- seq(0, 0.98, length.out = 100)  # Probability grid
x <- seq(0, 100, 0.1)

# Variance and standard deviation calculations
term <- (pi / b) / sin(pi / b)
sigma_x <- sqrt(a^2 * (term^2 - 1))
sigma_e <- nsr * sigma_x / sqrt(2)

# Generate data
X <- rllogis(n, a, b)
e <- rlaplace(n, 0, sigma_e)
Y <- X + e

# Bandwidth estimation and deconvolution CDF/PDF
bw <- bw.dboot2(Y, sig = sigma_e, error = "laplacian")
Fn <- DeconCdf(Y, sigma_e,x, error = "laplacian", bw = bw)
fn <- DeconPdf(Y, sigma_e,x, error = "laplacian", bw = bw)

# Ensure non-negative CDF values
Fn$x <- pmax(Fn$x, 0)

# Define survival and density functions
CDF_n <- splinefun(Fn$x, Fn$y)
PDF_n <- splinefun(fn$x, fn$y)
CDF <- function(x) pllogis(x, a, b)
PDF <- function(x) dllogis(x, a, b)
Q <- function(x) qllogis(x, a, b)
Qn <- splinefun(Fn$y, Fn$x, method = "natural")

# --- Plot CDF and PDF Comparisons ---
ggplot() +
  geom_line(aes(x = Fn$x, y = Fn$y, color = "Deconvolution CDF")) +
  geom_line(aes(x = Fn$x, y = CDF(Fn$x), color = "True CDF")) +
  scale_color_manual(values = c("Deconvolution CDF" = "black", "True CDF" = "red")) +
  labs(y = "CDF", x = "x", color = "Legend") +
  xlim(0,100)+
  theme_minimal()

ggplot() +
  geom_line(aes(x = fn$x, y = fn$y, color = "Deconvolution PDF")) +
  geom_line(aes(x = fn$x, y = PDF(fn$x), color = "True PDF")) +
  scale_color_manual(values = c("Deconvolution PDF" = "black", "True PDF" = "red")) +
  labs(y = "PDF", x = "x", color = "Legend") +
  theme_minimal()

# --- Create GTTT Data ---
GTtt <- data.frame(
  x = p,
  y = pmax(GTTT("LL", CDF_n, Qn, p),0),
  z = GTTT("LL", CDF, Q, p)
)

ggplot(GTtt, aes(x = x)) +
  geom_line(aes(y = y, color = "Estimated TTT")) +
  geom_line(aes(y = z, color = "True TTT")) +
  ylim(0, 3) +
  theme_minimal()

# --- Compute Least Concave Majorant (LCM) ---
GTtt <- GTtt %>% distinct(x, .keep_all = TRUE) %>% arrange(x)
ll <- gcmlcm(GTtt$x, GTtt$y, type = "lcm")
linear_interp <- approx(ll$x.knots, ll$y.knots, xout = seq(0, 1, 0.001))
LCM_TTT <- splinefun(linear_interp$x, linear_interp$y, method = "natural")

GTtt <- GTtt %>% mutate(t = LCM_TTT(x))

# --- Plot GTTT Functions with LCM ---
ggplot(GTtt, aes(x = x)) +
  geom_line(aes(y = y, color = "Estimated TTT")) +
  geom_line(aes(y = z, color = "True TTT")) +
  geom_line(aes(y = t, color = "LCM of TTT")) +
  ylim(0, 3) +
  scale_color_manual(values = c("Estimated TTT" = "red", "True TTT" = "black", "LCM of TTT" = "blue")) +
  labs(y = "TTT Function", x = "x", color = "Legend") +
  theme_minimal()

# --- Gttt Estimation ---
DRV_LCM <- splinefun(p, LCM_TTT(p, deriv = 1), method = "natural")
Data=data.frame(
  x=p,
  y=Gttt("LL", Q, PDF, p),
  z=Gttt("LL", Qn, PDF_n, p),
  t=DRV_LCM(p)
)

Data|>
  ggplot(aes(x=x))+
  geom_line(aes(y=y,color="True Gttt"))+
  geom_line(aes(y=z,color="Deconvolution Gttt"))+
  geom_line(aes(y=t,color="Constrained Gttt"))+
  scale_color_manual(values=c("True Gttt"="black","Deconvolution Gttt"="red","Constrained Gttt"="green"))+
  labs(y="Gttt",x="p",color="Legend")+
  xlim(0,0.9)+
  ylim(0,3)+
  theme_minimal()

# --- Define the Deconvolution Hazard Rate Estimator ---
DRV <- splinefun(p, 1 / DRV_LCM(p), method = "natural")


Data=data.frame(
  x=Fn$x,
  y=GHR("LL", CDF,PDF, Fn$x),
  z=GHR("LL", CDF_n,PDF_n, Fn$x),
  t=DRV(Fn$y)
  
)
Data|>
  ggplot(aes(x=x))+
  geom_line(aes(y=y,color="True GHR"))+
  geom_line(aes(y=z,color="Deconvolution GHR"))+
  geom_line(aes(y=t,color="Constrained GHR"))+
  scale_color_manual(values=c("True GHR"="black","Deconvolution GHR"="red","Constrained GHR"="green"))+
  labs(y="GHR",x="x",color="Legend")+
  xlim(0,60)+
  ylim(0,15)+
  theme_minimal()




#################################################
# --- Compute and compare the Hazard function ---
#################################################
# --- Cumulative Hazard Rate Estimation ---
Cum <- function(S, x) {
  sapply(x, function(xi) integrate(S, lower = 0, upper = xi)$value)
}
Examp1 <- splinefun(Fn$x, DRV(Fn$y), method = "natural")

Data_CHR_DCHR <- data.frame(
  x = Fn$x,
  y = GHF("LL", CDF, Fn$x),
  z = GHF("LL", CDF_n, Fn$x),
  t = Cum(Examp1, Fn$x)
)

ggplot(Data_CHR_DCHR, aes(x = x)) +
  geom_line(aes(y = y, color = "True CHR")) +
  geom_line(aes(y = z, color = "Deconvolution CHR")) +
  geom_line(aes(y = t, color = "Constrained")) +
  scale_color_manual(values = c("True CHR" = "red", "Deconvolution CHR" = "black", "Constrained" = "blue")) +
  labs(y = "Cumulative HR", x = "x", color = "Legend") +
  xlim(0, 25) + ylim(0, 100) +
  theme_minimal()

# --- Integrated Hazard Rate (IHR) ---
Cumulative_DHR <- Cum(Examp1, Fn$x)

Data_IHR <- data.frame(
  x = Fn$x,
  y = pllogis(Fn$x, a, b),
  z = pllogis(Cumulative_DHR, a, b),
  t = Fn$y
)

ggplot(Data_IHR, aes(x = x)) +
  geom_line(aes(y = y, color = "True IHR")) +
  geom_line(aes(y = z, color = "Deconvolution IHR")) +
  geom_line(aes(y = t, color = "Deconvolution CDF")) +
  scale_color_manual(values = c("True IHR" = "black", "Deconvolution IHR" = "red", "Deconvolution CDF" = "blue")) +
  labs(y = "Integrated HR", x = "x", color = "Legend") +
  theme_minimal()



##############################################
##############################################
##############################################

# Set parameters
nsr <- 0.5
a=1.5
b=1
var_x <- b^2 * (gamma(1 + 2/a) - (gamma(1 + 1/a))^2)
sigma_x <- sqrt(var_x)
sigma_e <- nsr * sigma_x / sqrt(2)


# Generate data
n <- 100
X <- rweibull(n,a,b)
e <- rlaplace(n, 0, sigma_e)
Y <- X + e

# Bandwidth estimation and deconvolution CDF
bw <- bw.dboot2(Y, sig = sigma_e, error = "laplacian")
Fn <- DeconCdf(Y, sigma_e, error = "laplacian", bw = bw)

# Modify the CDF to ensure non-negative values
Fn$x <- pmax(Fn$x, 0)

# Create survival functions
Sn <- splinefun(Fn$x, 1 - Fn$y)
S <- splinefun(Fn$x, 1 - pweibull(Fn$x, a,b))

# Create quantile estimator
Qn <- data.frame(x = Fn$y, y = Fn$x) %>%
  distinct(x, .keep_all = TRUE)  # Remove duplicates

Q <- qweibull(Qn$x,a,b)
Q[is.infinite(Q)] <- max(Q[is.finite(Q)])  # Handle infinite values

# Create TTT data
Ttt <- data.frame(
  x = Qn$x,
  y = TTT(Sn, Qn$y),
  z = TTT(S, Q)
)

# Remove duplicates before applying LCM
Ttt <- Ttt %>% distinct(x, .keep_all = TRUE)

# Compute the LCM (Least Concave Majorant)
ll <- gcmlcm(sort(Ttt$x), Ttt$y, type = "lcm")
linear_interp <- approx(ll$x.knots, ll$y.knots, xout = seq(0, 1, 0.001))

# Create LCM function and apply it to TTT
LCM_TTT <- splinefun(linear_interp$x, linear_interp$y, method = "natural")

###########################################
############ HR Vs DHR ####################
###########################################
# Compute the derivative of the LCM operator
DRV_LCM <- splinefun(Ttt$x, LCM_TTT(Ttt$x, deriv = 1), method = "natural")

# Define the Deconvolution Hazard Rate estimator
DHR <- splinefun(Fn$x, 1 / DRV_LCM(Fn$y), method = "natural")

###########################################
########### CHR Vs DCHR ###################
###########################################

# Define cumulative function (integral approximation)
Cum <- function(S, x) {
  sapply(x, function(xi) integrate(S, lower = 0, upper = xi)$value)
}

###########################################
########### IHR F Vs F ####################
###########################################

Cumulative_DHR <- Cum(DHR, x)

z = pexp(Cumulative_DHR)  # Deconvolved IHR
t = DeconCdf(Y, sigma_e, x, error = "laplacian", bw = bw)$y  # Estimated CDF


plot(x, z, type = "l", col = "red", xlab = "x", ylab = "IHR")
lines(x, t, col = "blue")

###########################################################
#### Calculate the MSE for the Beta Distribution
###########################################################
library(dplyr)
library(ggplot2)

# Set parameters
m <- 500  # Number of simulations
k <- 10  # Number of support points
NSR <- 0.1  # Noise-to-Signal Ratio
n <- 500  # Sample size

nsr <- 0.1
a=1.25
b=1
var_x <- b^2 * (gamma(1 + 2/a) - (gamma(1 + 1/a))^2)
sigma_x <- sqrt(var_x)
sigma_e <- nsr * sigma_x / sqrt(2)

# Define a function for cumulative integration
Cum <- function(S, x) {
  sapply(x, function(xi) integrate(S, lower = 0, upper = xi)$value)
}

# Simulation loop
results <- numeric(k - 1)
for (j in 2:(k - 1)) {
  p <- j / k
  MSE1 <- 0
  MSE2 <- 0
  
  for (i in 1:m) {
    # Generate data
    X <- rweibull(n, shape = a, scale = b)
    e <- rlaplace(n, 0, sigma_e)
    Y <- X + e
    
    # Bandwidth estimation and deconvolution
    bw <- bw.dboot2(Y, sig = sigma_e, error = "laplacian")
    Fn <- DeconCdf(Y, sigma_e, error = "laplacian", bw = bw)
    Fn$x <- pmax(Fn$x, 0)  # Ensure non-negative values
    
    # Survival functions
    Sn <- splinefun(Fn$x, 1 - Fn$y)
    S <- splinefun(Fn$x, 1 - pweibull(Fn$x, shape = a, scale = b))
    
    # Quantile estimator
    Qn <- distinct(data.frame(x = Fn$y, y = Fn$x), x, .keep_all = TRUE)
    Q <- qweibull(Qn$x, shape = a, scale = b)
    Q[is.infinite(Q)] <- max(Q[is.finite(Q)])  # Handle infinities
    
    # TTT data
    Ttt <- data.frame(
      x = Qn$x,
      y = TTT(Sn, Qn$y),
      z = TTT(S, Q)
    ) %>% distinct(x, .keep_all = TRUE)  # Remove duplicates
    
    # Least Concave Majorant (LCM)
    ll <- gcmlcm(sort(Ttt$x), Ttt$y, type = "lcm")
    linear_interp <- approx(ll$x.knots, ll$y.knots, xout = seq(0, 1, 0.001))
    LCM_TTT <- splinefun(linear_interp$x, linear_interp$y, method = "natural")
    
    # Hazard rate estimation
    DRV_LCM <- splinefun(Ttt$x, LCM_TTT(Ttt$x, deriv = 1), method = "natural")
    DHR <- splinefun(Fn$x, 1 / DRV_LCM(Fn$y), method = "natural")
    
    # Cumulative hazard rate
    Cumulative_DHR <- Cum(DHR, Fn$x)
    
    # Create smooth functions for estimation
    smooth_function1 <- splinefun(Fn$x, Fn$y, method = "natural")
    smooth_function2 <- splinefun(Fn$x, pexp(Cumulative_DHR), method = "natural")
    
    # Evaluate at quantiles
    y1 <- smooth_function1(qweibull(p, shape = a, scale = b))
    y2 <- smooth_function2(qweibull(p, shape = a, scale = b))
    
    # Compute MSE
    MSE1 <- MSE1 + (y1 - p)^2
    MSE2 <- MSE2 + (y2 - p)^2
  }
  
  # Store results
  results[j] <- MSE2 / MSE1
}

# Store and visualize results
M2 <- tibble(x = (1:(k - 1)) / k, ratio = results)

ggplot(M2, aes(x = x, y = ratio)) +
  geom_line(color = "blue") +
  labs(x = "Quantile", y = "MSE Ratio", title = "MSE Ratio across Quantiles")

plot(x,pbeta(x,p1,p2),type="l",col="black")

