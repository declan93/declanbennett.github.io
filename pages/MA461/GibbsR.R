library(expm)
gibbs_func = function(nsteps = 10000) {
  
  
  #P(T = T1 | D = g): .4
  #P(T = T1 | D = r): .6
  #P(D = g | T = T1): .5
  #P(D = g | T = T2): .6
  ####
  ##
  
  #Initialize
  Tw = "T2"
  D = "g"
  
  count = cbind(rep(0,4))
  count1 = cbind(rep(0,16))
  rownames(count) = c("T1g","T2g","T1r","T2r")
  
  x <- c("T1g","T2g","T1r","T2r")
  d1 <- combn(x,2)
  
  d2 <- expand.grid(x,x)
  nmes <- c()
  for (i in 1:length(d2[,1])){
    nmes <- c(nmes,paste(d2[i,1],d2[i,2],sep=""))
  }
  rownames(count1) <-nmes
  pair_s <- c() 
  
  for(i in 1:nsteps) {
    # update the S variable
    r = runif(1)
    if((D == "g" & r < 0.4) | (D == "r" & r < .6)) {
      Tw = "T1"
    }
    else {
      Tw = "T2"
    }
    
    state1 = paste(T,D,sep="")
    count[state1,1] = count[state1,1] + 1
    
    # update the H variable
    r = runif(1)
    if((Tw == "T2" & r < .6) | (Tw == "T1" & r < .5)) {
      D <- "g"
    }
    else {
      D <- "r"
    }
    
    state2 = paste(T,D,sep="")
    count[state2,1] = count[state2,1] + 1
    # get chain 
    pair_s <- append(pair_s, state1)
    pair_s <- append(pair_s, state2)
  }
  
  long_term <- count/(2*nsteps)
  # count transitions in chain.
for (j in 2:length(pair_s) -1){
  trans <- paste(pair_s[j],pair_s[j+1],sep="")
  count1[trans,1] = count1[trans,1] + 1
}
  tpm <- matrix(count1, 4) 
  # you can work out the transition probabilities from the conditionals directly. ie, T1rT2r -> T1rT2g = 1/n * P(T2g | T1r) and so on...
  # I think this is confucing so I'm just going to work it out empirically from the sampled chain
  tpm1 <- tpm
  for (i in 1:4){
    tpm1[i,] <- tpm[i,]/sum(tpm[i,])
  }
  
  colnames(tpm1) <- c("rr","rg","gr","gg")
  rownames(tpm1) <- c("rr","rg","gr","gg")
  return(list(long_term=long_term, tpm=tpm1))
}

X <-gibbs_func()

# calculate limiting dist and stationary dist 
tpm <- X$tpm
lim_d <- X$long_term
lim_d[,1]%*%tpm

eig <- eigen(t(tpm))
eig
#eigenvectors
e_vec <- eig$vectors

# eigenvalue
e_val <- eig$values
e_val

# rescale to get probabilities where lambda = 1 (col 1) 
pi_vec_eig <- e_vec[,1]/sum(e_vec[,1]) 


pi_vec_eig
lim_d
# as nsteps is increased the limiting distribution and staionaruy distribution should converge. if our conditional probabilities work (I've chosen random probabilities)

