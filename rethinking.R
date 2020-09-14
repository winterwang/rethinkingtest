install.packages(c("devtools","mvtnorm","loo","coda"),dependencies=TRUE)
library(devtools)
install_github("rmcelreath/rethinking",ref="Experimental")



##%######################################################%##
#                                                          #
####                    R code  0.3                     ####
#                                                          #
##%######################################################%##

(log(0.01^200))
200*log(0.01)


##%######################################################%##
#                                                          #
####                     R code 2.2                     ####
#                                                          #
##%######################################################%##

dbinom(6, size = 9, prob = 0.5)


##%######################################################%##
#                                                          #
####                     Rcode  2.3                     ####
#                                                          #
##%######################################################%##


# define grid 

p_grid <- seq(from = 0, to = 1, length.out = 20)

# define prior 

prior <- rep(1, 20)
prior <- ifelse(p_grid < 0.5, 0, 1)
prior <- exp(-5 * abs( p_grid - 0.5))


# compute the likelihood at each value in grid

likelihood <- dbinom(6, size = 9, prob = p_grid)

# compute the product of likelihood and prior 

unstd.posterior <- likelihood*prior

# standardize the posterior, so it sums to 1

posterior <- unstd.posterior / sum(unstd.posterior)


##%######################################################%##
#                                                          #
####                     R code 2.4                     ####
#                                                          #
##%######################################################%##


plot(p_grid, posterior, type = "b", 
     xlab = "probability of water", 
     ylab = "posterior probability")

mtext("20 points")

##%######################################################%##
#                                                          #
####                     R code 2.6                     ####
#                                                          #
##%######################################################%##

library(rethinking)

globe.qa <- quap(
  alist(
    W ~ dbinom(W + L, p) , # binomial likelihood
    p ~ dunif(0, 1)        # uniform prior
  ),
  data = list(W = 6, 
              L = 3)
)

# Display summary of quadratic approximation

precis(globe.qa)


##%######################################################%##
#                                                          #
####                     R code 2.7                     ####
#                                                          #
##%######################################################%##


# analytical calculation 

W <-  6
L <-  3 

curve(dbeta(x, W + 1, L + 1), from = 0, to = 1)

# quadratic approximation
curve(dnorm( x, 0.67, 0.16), lty = 2, add = TRUE)



##%######################################################%##
#                                                          #
####                    R code  2.8                     ####
#                                                          #
##%######################################################%##

n_samples <- 1000

p <- rep(NA, n_samples)

p[1] <- 0.5
W <- 6 
L <- 3

for (i in 2:n_samples) {
  p_new <- rnorm(1, p[i - 1], 0.1)
  if (p_new < 0) p_new <- abs(p_new)
  if (p_new > 1) p_new <- 2 - p_new
  q0 <- dbinom(W, W+L, p[i - 1])
  q1 <- dbinom(W, W+L, p_new)
  p[i] <- ifelse( runif(1) < q1/q0, p_new, p[i-1])
}

##%######################################################%##
#                                                          #
####                     Rcode 2.9                      ####
#                                                          #
##%######################################################%##
dens(p, xlim=c(0,1))
curve(dbeta(x, W + 1, L + 1), lty = 2, add = TRUE)


##%######################################################%##
#                                                          #
####                    R code  3.2                     ####
#                                                          #
##%######################################################%##




p_grid <- seq(from = 0, to = 1, length.out = 1000)

prob_p <- rep(1, 1000)

prob_data <- dbinom(6, size = 9, prob = p_grid)

posterior <- prob_data * prob_p

posterior <- posterior / sum(posterior)

##%######################################################%##
#                                                          #
####                     R code 3.3                     ####
#                                                          #
##%######################################################%##

samples <- sample(p_grid, prob = posterior, size = 10000, replace = TRUE)

plot(samples)

dens(samples)


##%######################################################%##
#                                                          #
####                     R Code 3.6                     ####
#                                                          #
##%######################################################%##

# add up posterior probability where p < 0.05

sum(posterior[ p_grid < 0.5])

sum(samples < 0.5) / 10000

sum( samples > 0.5 & samples < 0.75) / 10000

quantile( samples, 0.8)

quantile( samples, c(0.1, 0.9))


##%######################################################%##
#                                                          #
####                    R code 3.11                     ####
#                                                          #
##%######################################################%##


p_grid <- seq(from = 0, to = 1, length.out = 1000)
prior <- rep(1, 1000)
likelihood <- dbinom(3, size = 3, prob = p_grid)
likelihood <- dbinom(6, size = 9, prob = p_grid)

posterior <- likelihood * prior

posterior <- posterior / sum(posterior)

samples <- sample( p_grid, size = 10000, replace = TRUE, prob = posterior)

plot(samples)
dens(samples)

PI(samples , prob = 0.5)
HPDI(samples, prob = 0.5)
PI(samples , prob = 0.8)
HPDI(samples, prob = 0.8)
PI(samples , prob = 0.95)
HPDI(samples, prob = 0.95)


p_grid[ which.max(posterior)]

chainmode( samples, adj = 0.01)


##%######################################################%##
#                                                          #
####                    R code 3.17                     ####
#                                                          #
##%######################################################%##


sum( posterior*abs( 0.5 - p_grid))

loss <- sapply( p_grid, function(d) sum(posterior * abs(d - p_grid)))

p_grid[ which.min(loss)]

median(samples)



##%######################################################%##
#                                                          #
####                    R code 3.20                     ####
#                                                          #
##%######################################################%##


dbinom(0 : 2, size = 2, prob = 0.7)


rbinom(1, size = 2, prob = 0.7)

rbinom(10, size = 2, prob = 0.7)


dummy_w <- rbinom(100000, size = 2, prob = 0.7)

table(dummy_w) / 100000

dummy_w <- rbinom(100000, size = 9, prob = 0.7)

simplehist(dummy_w, xlab = "dummary water count")



##%######################################################%##
#                                                          #
####                    R code 3.25                     ####
#                                                          #
##%######################################################%##

w <- rbinom(10000, size = 9, prob = 0.6)


simplehist(w)


w <- rbinom(10000, size = 9, prob = samples)

simplehist(w)

##%######################################################%##
#                                                          #
####                     R code 4.1                     ####
#                                                          #
##%######################################################%##

pos <- replicate( 1000, sum( runif(16, -1, 1)))
hist(pos)
plot(density(pos))


##%######################################################%##
#                                                          #
####                     R code 4.2                     ####
#                                                          #
##%######################################################%##



prod( 1 + runif(12, 0, 0.1))

growth <- replicate( 10000, prod( 1 + runif(12, 0, 0.1)))

dens( growth, norm.comp = TRUE)

##%######################################################%##
#                                                          #
####                        4.4                         ####
#                                                          #
##%######################################################%##


big <- replicate( 10000, prod( 1 + runif(12, 0, 0.5)))

small <- replicate( 10000, prod( 1 + runif(12, 0, 0.01)))

dens( big, norm.comp = TRUE)


dens( small, norm.comp = TRUE)


##%######################################################%##
#                                                          #
####                        4.5                         ####
#                                                          #
##%######################################################%##

log.big <- replicate( 10000, log(prod( 1 + runif(12, 0, 0.5))))

dens( log.big, norm.comp = TRUE)

dnorm(0, 0, 0.1)



##%######################################################%##
#                                                          #
####                        4.7                         ####
#                                                          #
##%######################################################%##


data("Howell1")

d <- Howell1
str(d)

precis(d)

dens(d[ d$age >= 18, ])



