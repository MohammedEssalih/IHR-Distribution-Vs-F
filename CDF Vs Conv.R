
# Set the number of samples
n <- 100000
t=seq(-1,1,0.01)
# Generate samples from Beta(0.5, 1) distribution
beta_samples <- pbeta(t, 0.5, 1)
# Generate samples from Laplace(0, 0.5) distribution
laplace_samples <- plaplace(t, location = 0, scale = 0.5)

# Sum the samples to get the convolution
convolution_samples <- beta_samples + laplace_samples

# Estimate the ECDF of the convolution
ecdf_convolution <- ecdf(convolution_samples)

# Create a data frame for ggplot
Data <- data.frame(
  x = t,
  y = laplace_samples
)

# Plot the estimated ECDF using ggplot2
Data|>
ggplot(aes(x = x, y = y)) +
  geom_line(color = "black", size = 2) +
  labs(
    title = "The CDF of the error",
    x = "X-axis",
    y = "Cumulative Probability"
  ) +
  theme_minimal()

Data|>
  ggplot(aes(x=x,y=y))+
  geom_line(color = "black", size = 2) +
  
  # Add a Dash of Style: Customize the aesthetics for a visually appealing touch
  theme_minimal() +
  labs(
    title = "The CDF of the error",
    x = "X-axis",
    y = "Cumulative Probability"
  )

###################################33
n1=1000
x <- rbeta(n1,p1,p2)
u <- rlaplace(n1, 0,sig2)
w <- x+u
y=dbeta(t,0.5,1)
t=seq(-1,5,0.01)
density(x)

plot(density(x))
lines(t,y)
# Theoretical CDF function
theoretical_cdf <- function(y) {
  y
}

# Values for the theoretical CDF
x <- seq(0, 1, length.out = 1000)


Data=data.frame(x=x,y=theoretical_cdf(x))

Data|>
  ggplot(aes(x=x,y=y))+
  geom_line(color = "black", size = 2) +
  
  # Add a Dash of Style: Customize the aesthetics for a visually appealing touch
  theme_minimal() +
  labs(
    title = "The CDF of the error",
    x = "X-axis",
    y = "Cumulative Probability"
  )


