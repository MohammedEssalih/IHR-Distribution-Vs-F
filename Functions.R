TTT <- function(S,x) {
  d <- sapply(x, function(x) {
    integrate(S, lower = 0, upper = x)$value
  }) 
  d
}
