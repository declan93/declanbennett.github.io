library(expm)

# generate TPM
TPM <- matrix(c(.29,.21,.05,.45,
                .58,.06,.22,.14,
                .10,.28,.37,.25,
                .34,.19,.11,.36),4,byrow = T)

# stadionary distribution can be found from eigenvectors with eigenvalue = 1 as it is a left eigenvector of the TPM
# get eigenvector corresponding to eigenvalue = 1 of transposed TPM 
eig <- eigen(t(TPM2))

#eigenvectors
e_vec <- eig$vectors

# eigenvalue
e_val <- eig$values
e_val

# rescale for probabilities = lambda is col 1 
pi_vec_eig <- e_vec[,1]/sum(e_vec[,1]) 


# When our chain is ergodic (irreducible and aperiodic)
# The limiting distribution is equal to stationary distribution
# E.g limiting distribution P()**N = pi where N -> inf (** is raise to power)

pi_vec <- TPM%^%100

# go compare
pi_vec[1,]
pi_vec_eig

# Proof that for unique pi
#  pi = pi * TPM
pi_vec[1,]
pi_vec[1,]%*%TPM
