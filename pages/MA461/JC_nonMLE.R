library(expm) # only library required

## sequences must be of same length
seq1 <- "CAGCAGAGAGTCTCCCATGATAACTGCGTT"
seq2 <- "CAGCATACTGGAATCCGTGCTAACCGCATT"

Get_dist<- function(seq1,seq2){
  if (nchar(seq1) != nchar(seq2)){
    stop("Error: Difference is sequence length")
    geterrmessage()}
  
  seq_l1 <- as.list(strsplit(seq1, ""))[[1]]
  seq_l2 <- as.list(strsplit(seq2, ""))[[1]]
  sequences <- mapply(c,seq_l1,seq_l2) # make matrix 
  change <- c("AG", "GA", "CT", "TC","AT", "AC", "CA", "CG", "GT", "GC", "TG", "TA")
  No_change <- c("CC","AA","GG","TT")
  
  alpha<- 0 # Number of changes
  i <- 0
  # loop over sequence matrix and counting changes this part can be changed to calculate num of Ti and num of Tv. 
  # i.e split change into transitions (p) and transversions (q)
  
  for (i in 1:length(sequences[1,])){
    if (paste(sequences[,i],collapse = "") %in% change){
      alpha <- alpha +1
    }else if (paste(sequences[,i],collapse = "") %in% No_change){}
    i<-i+1
  }
  # make relative to sequence length
  p <- (alpha/(i-1)) # % fraction of sites differing
  
  # Generator matrix if you want you can set the off diagonals to be what ever you like 
  # JC69 all are equal so will not matter when you set the diagonals and scale
  # kimura 2 param uses different Ti and Tv #!! k != Ti/Tv !!#
  G = matrix(c((-3),1,1,1,1,(-3),1,1,1,1,(-3),1,1,1,1,(-3)),4)
  
  # bells and whistles
  Nucleotides <- c("A","C","G","T")
  rownames(G) <- Nucleotides
  colnames(G) <- Nucleotides
  
  # Beta scaling factor for G wikipedia "models of DNA evolution" scaling branch lengths. time is now "expected mutations per site"
  # Scaling is easy as we assume equal nucleotide frequencies pi = .25
  # Our scaling equates to 1 / (-1 * .25*-12)
  Beta<-(1/(-1*(0.25*sum(G[1,1],G[2,2],G[3,3],G[4,4]))))
  G1<- Beta*G # G1 is the scaled generator matrix diagonals ==  -1 
  
  # calculate genetic distance using JC69 parameter distance equation
  # You could remove this step and pass as an argument and optimise to find d 
  d <- -1*(3/4)*(log(1-((4/3)*(p))))
  
  # create transition probability matrix for time equal to d 
  TPM <- expm(d*G1)
  scaleG =d*G1
  # Log-Likelihood of the sequences
  accum <- 0
  # now to start off probability of sequence
  for (i in length(sequences[1,])){
    accum = accum + log(TPM[seq_l1[1],seq_l2[i]])
  }
  
  return(list(Transition_Probability_Matrix_time_d = TPM, Associated_Genetic_Distance = d,
              Likelihood_of_seq=exp(accum),difference_Percentage= p, Differnece_sites=alpha,
              Seq_length=i-1, Generator= G,scaled_Generator = G1, beta_scaling_factor= Beta, Before_exponential_d_G1=scaleG))
}

Get_dist(seq1,seq2)

hld = Get_dist(seq1,seq2)
expm(hld$d_G1)
