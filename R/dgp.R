
std <- 0.3
mid <- 0.8

gamma_shape <- function(x) {
  return(0.5 + x ^ 2)
}

m <- function(x) {
  return(1 + 2 * (x - mid))
}

p <- function(x) {
  return(exp(-10 * (x - mid) ^ 2))
}

mu0 <- function(x) {
  return(5 * exp(15 * (x - 0.5)) / (1 + exp(15 * (x - 0.5))) - 4 * x)
}

# Generation function
r <- function(x) {
  n <- length(x)
  z <- rbinom(n, 1, p(x))
  return(z * rnorm(n, m(x), std) + 
           (1 - z) * log(rgamma(n, gamma_shape(x), 1)) + mu0(x))
}

f1 <- function(y, x) {
  return(dnorm(y, mu0(x) + m(x), std))
}

f2 <- function(y, x) {
  return(dgamma(exp(y - mu0(x)), gamma_shape(x), 1) * exp(y - mu0(x)))
}

F1 <- function(y, x) {
  return(pnorm(y, mu0(x) + m(x), std))
}

F2 <- function(y, x) {
  return(pgamma(exp(y - mu0(x)), gamma_shape(x), 1))
}

pp <- function(y, x) {
  p <- p(x)
  return(p * F1(y, x) + (1 - p) * F2(y, x))
}

# Density function
d <- function(y, x) {
  p <- p(x) 
  return(p * f1(y, x) + (1 - p) * f2(y, x))
}

## Functions for DR-BART -- SBART-DS comparison. 
mu_quad <- function(x, a) {
  return(a * (x - 0.5) ^ 2)
}

r_quad <- function(x, a) {
  n <- length(x)
  z <- rbinom(n, 1, p(x))
  z * rnorm(n, m(x), std) +
    (1 - z) * log(rgamma(n, gamma_shape(x), 1)) + mu_quad(x, a)
}

f1_quad <- function(y, x, a) {
  return(dnorm(y, mu_quad(x, a) + m(x), std))
}

f2_quad <- function(y, x, a) {
  dgamma(exp(y - mu_quad(x, a)), gamma_shape(x), 1) * exp(y - mu_quad(x, a)) 
}


d_quad <- function(y, x, a) {
  p <- p(x)
  dens <- p * f1_quad(y, x, a) + (1 - p) * f2_quad(y, x, a)
  return(dens)
}
##

gg <- function(x) {
  x_grid[findInterval(uu, pp(x_grid, x))]
}