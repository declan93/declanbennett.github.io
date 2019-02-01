library(expm)
# Download R script
download.file("https://declan93.github.io/pages/MA461/MA461_tut1_stat_dist.R",destfile = "MA461_tut1_stat_dist.R")


# generate TPM
TPM <- matrix(c(.29,.21,.05,.45,
                .58,.06,.22,.14,
                .10,.28,.37,.25,
                .34,.19,.11,.36),4, byrow = T)
TPM
# stationary distribution can be found from eigenvectors with eigenvalue = 1 as it is a left eigenvector of the TPM
# get eigenvector corresponding to eigenvalue = 1 of transposed TPM 
eig <- eigen(t(TPM))
eig
#eigenvectors
e_vec <- eig$vectors

# eigenvalue
e_val <- eig$values
e_val

# rescale to get probabilities where lambda = 1 (col 1) 
pi_vec_eig <- e_vec[,1]/sum(e_vec[,1]) 
pi_vec_eig

# When our chain is ergodic (irreducible and aperiodic)
# The limiting distribution is equal to stationary distribution
# E.g limiting distribution P()**N = pi where N -> inf (** is raise to power)

pi_vec <- TPM%^%100

# Go Compare
pi_vec[1,]
pi_vec_eig

# Demonstrate that for unique pi
#  pi = pi * TPM
pi_vec[1,]%*%TPM


# Recall chapman kolmogorov eq. P(n) = P**n i.e probability transition in n steps corresponds to P**n
# Rename rows cols
colnames(TPM) <- c("i","j","k","l")
rownames(TPM) <- c("i","j","k","l")
TPM

# Probability of going from i to j where nsteps = 1 
TPM["i","j"]

# Probability of going from i to j where nsteps = 2
(TPM%^%2)["i","j"]

# and so on ..
 
# Demonstrate that the stationary distribution does NOT always == limiting distribution
P <- matrix(c(0,1,1,0),2)
P

# no convergence on a limiting distribution
P%^%10000010
P%^%10000010

## but P does have a stationary distribution
P_eig <- eigen(t(P))
P_eig
#eigenvectors
P_e_vec <- P_eig$vectors

# eigenvalue
P_e_val <- P_eig$values
e_val

# rescale for probabilities = lambda is col 1 
P_vec_eig <- P_e_vec[,1]/sum(P_e_vec[,1]) 
P_vec_eig

###########################################################
## Symmetric matrix always satisfies detailed balance
## doubly stochastic (rows and cols sum to 1) has a special quality that pi is uniform 1/N (N = no. states)
TPM2 <- matrix(c(.2, .1, .4, .3,
                .1, .2, .3, .4,
                .4, .3, .2, .1,
                .3, .4, .1, .2),4, byrow = T)
TPM2
t(TPM2)
TPM2 == t(TPM2)

# Eigenvalue decomp
eig2 <- eigen(t(TPM2))
eig2
#eigenvectors
e_vec2 <- eig2$vectors
# eigenvalue
e_val2 <- eig2$values

# rescale for probabilities = lambda is col 1 
pi_vec_eig2 <- e_vec2[,1]/sum(e_vec2[,1]) 

# limiting distribution
pi_vec2 <- TPM2%^%100

# Compare
pi_vec_eig2
pi_vec2[1,]

##############################################
## Example of 3x3 stochastic that is nonsymmetric (non uniform stationary)
TMP3 <- matrix(c(.3,.4,.5,.3,.4,.3,.4,.2,.2),3) 
TMP3
t(TMP3)
TMP3%^%100



##############################################
# Reversibility. Needs to satisfy detailed balance. pi_i Pij = pi_j Pji
# Q is TPM reversible
TPM
pi_vec_eig[1]*TPM["i","j"]
pi_vec_eig[2]*TPM["j","i"]

# How about TPM2? Symmetric matrices are always reversible
pi_vec_eig2[1]*TPM2[1,2] == pi_vec_eig2[2]*TPM2[2,1]

# aside Pij =  pi_j/ pi_i * Pji
pi_vec_eig2[2]/pi_vec_eig2[1] * TPM2[2,1]

