# optim example of a log liklihood of normal distribution and a cooler sphere example that shows 3d "probability space"

# logliklihood of normal
x =rnorm(10000) # sample 100 numbers from normal distribution mean =0 std =1

# log liklihood of a normal ("Gaussian distribution") is goven by
# L(X|μ,σ2)  =−(n2)ln(2πσ2) − 1/2σ2 ∑ni(xi−μ)2

llik = function(x,par){ # par is mean and std estimates
  m=par[1]
  s=par[2]
  n=length(x)# log of the normal likelihood
  # -n/2 * log(2*pi*s^2) + (-1/(2*s^2)) * sum((x-m)^2)
  ll = -(n/2)*(log(2*pi*s^2)) + (-1/(2*s^2)) *sum((x-m)^2)# return the negative to maximize rather than minimize
  return(-ll)
}
hist(x, freq=FALSE,col='tan')
lines(density(x),col='red',lwd=2)

plot(seq(-5,5,.1),-1*sapply(seq(-5,5,.1),FUN=llik,par=c(0,1)),type='l',ylab='log liklihood',xlab='Xi')
plot(seq(-5,5,.1),-1*sapply(seq(-5,5,.1),FUN=llik,par=c(-.5,1)),type='l',ylab='log liklihood',xlab='Xi')
plot(seq(-5,5,.1),-1*sapply(seq(-5,5,.1),FUN=llik,par=c(.5,1)),type='l',ylab='log liklihood',xlab='Xi')

# optim allows us to find the par that maximises the liklihood for the data
# we can choose and initial values
optim(par=c(.5,.5), llik, x=x) # maximise function
optim(par=c(.5,.5), llik, x=x,control=list(fnscale=-1)) # this tells optim to find params that minimise func 
# we see that we almost get a mean =0 and std = 1
mean(x)
sd(x)
##optim example https://statacumen.com/teach/SC1/SC1_17_optim.pdf
f.name <- "Sphere function"
f.sphere <- function(x){# make x a matrix so this function works for plotting and for optimizing
  x <- matrix(x, ncol=2)# calculate the function value for each row of x
  f.x <- apply(x^2, 1, sum)# return function value
  return(f.x)
}

x1 <- seq(-10, 10, length = 101)
x2 <- seq(-10, 10, length = 101)
X <- as.matrix(expand.grid(x1, x2))
colnames(X) <- c("x1", "x2")# evaluate function
y <- f.sphere(X)# put X and y values in a data.frame for plotting
df <- data.frame(X, y)# plot the function
library(lattice)# use the lattice package
wireframe(y ~ x1 * x2,
          data = df
          , main = f.name
          , shade = TRUE, scales = list(arrows = FALSE), screen = list(z = -50, x = -70))

optim(c(1,1), f.sphere, method = "Nelder-Mead")
optim(c(1,1), f.sphere, method = "SANN")
