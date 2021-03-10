# shameless copy job possibly from matt stephens ()
# This example we are sampling from an exponential distribution it would be easy to convert this so we are sampling from some state space we would then 
# use the acceptance from the notes aij = min[1,(πj * qji/πi * qij)] where qij is from the proposal matrix

# 
EXP_dist = function(x){
  if(x<0){
    return(0)}
  else {
    return( exp(-x))
  }
}

x = rep(0,5000) # chain
x[1] = 3 # random init
for(i in 2:5000){
  currentx = x[i-1] # current value
  proposedx = currentx + rnorm(1,mean=0,sd=1) # add a random number from a normal distribution
  A = EXP_dist(proposedx)/EXP_dist(currentx) # This is our alpha value
  if(runif(1)<A){
    x[i] = proposedx       # accept move with probabily min(1,A)
  } else {
    x[i] = currentx        # otherwise "reject" move, and stay where we are
  }
}
x1 = x

#lets repeat 3 times.
for(i in 2:5000){
  currentx = x[i-1]
  proposedx = currentx + rnorm(1,mean=0,sd=1)
  A = EXP_dist(proposedx)/EXP_dist(currentx) 
  if(runif(1)<A){
    x[i] = proposedx       
  } else {
    x[i] = currentx 
  }
}
x2 <- x

for(i in 2:5000){
  currentx = x[i-1]
  proposedx = currentx + rnorm(1,mean=0,sd=1)
  A = EXP_dist(proposedx)/EXP_dist(currentx) 
  if(runif(1)<A){
    x[i] = proposedx       # accept move with probabily min(1,A)
  } else {
    x[i] = currentx        # otherwise "reject" move, and stay where we are
  }
}
x3 <- x
# plot our results
plot(x1,type="l")
lines(x2, col=2)
lines(x3, col=3)

## sample allele frequency

prior = function(p){
  if((p<0) || (p>1)){  # || here means "or"
    return(0)}
  else{
    return(1)}
}

likelihood = function(p, nAA, nAa, naa){
  return(p^(2*nAA) * (2*p*(1-p))^nAa * (1-p)^(2*naa))
}

psampler = function(nAA, nAa, naa, niter, pstartval, pproposalsd){
  p = rep(0,niter)
  p[1] = pstartval
  for(i in 2:niter){
    currentp = p[i-1]
    newp = currentp + rnorm(1,0,pproposalsd)
    A = prior(newp)*likelihood(newp,nAA,nAa,naa)/(prior(currentp) * likelihood(currentp,nAA,nAa,naa))
    if(runif(1)<A){
      p[i] = newp       # accept move with probabily min(1,A)
    } else {
      p[i] = currentp        # otherwise "reject" move, and stay where we are
    }
  }
  return(p)
}


z=psampler(50,21,29,10000,0.5,0.01)
hist(z , probability = T)

x=seq(0,1,length=1000)
hist(z,prob=T)
lines(x,dbeta(x,122, 80))
