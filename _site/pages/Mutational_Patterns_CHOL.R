#setwd("/home/declan/teaching/Cancer_genomics/")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("BSgenome")
#BiocManager::install("MutationalPatterns")
library("MutationalPatterns")
library("BSgenome")
#BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
library("NMF")
library("curl")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
available.genomes()
ref_genome<- "BSgenome.Hsapiens.UCSC.hg19"
BiocManager::install(ref_genome)
library(ref_genome, character.only = TRUE)

### Get data 
# url to data
url_dat = "https://declan93.github.io/pages/CHOL_vcfs.tar.gz"
download.file(url_dat,destfile="tmp.tar.gz")
untar("tmp.tar.gz",list=TRUE)  ## check contents
untar("tmp.tar.gz")
setwd("./CHOL_vcfs/")
getwd()
# load data 
vcf_files <- list.files("./", pattern = ".vcf") # list vcf files
sample_names <- unlist(lapply(strsplit(unlist(vcf_files),".vcf"), "[[", 1)) # list sample names
tissue <- c(rep("CHOL",length(sample_names)))

vcfs <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome) # load vcfs as grange objects             
muts = mutations_from_vcf(vcfs[[1]]) # look at first sample
head(muts,12)
types = mut_type(vcfs[[1]])
head(types,12)
context = mut_context(vcfs[[1]], ref_genome)
head(context)
type_context = type_context(vcfs[[1]], ref_genome)
lapply(type_context,head,12)

type_occurrences <- mut_type_occurrences(vcfs, ref_genome = ref_genome) # counts table
type_occurrences

#Plots # Rstudio drop
plot_spectrum(type_occurrences)
plot_spectrum(type_occurrences,CT=TRUE)
plot_spectrum(type_occurrences,CT=TRUE,legend=F)
plot_spectrum(type_occurrences, by = tissue, CT = TRUE, legend = TRUE)

mut_mat <- mut_matrix(vcf_list = vcfs, ref_genome = ref_genome) # context matrix for NMF
head(mut_mat)

plot_96_profile(mut_mat[,c(1,7)]) # look at spectrum for samples 1 & 7

estimate <- nmf(mut_mat, rank=2:5, method="brunet", nrun=10, seed=123456)
plot(estimate)

res <- extract_signatures(mut_mat, rank = 3, nrun = 10) # function to call NMF return list of vectors
head(res$signatures,12) # r x k
head(res$contribution,12)# k x n
colnames(res$signatures) <- c("Sig_A","Sig_B","Sig_C")
rownames(res$contribution) <- c("Sig_A","Sig_B","Sig_C")
plot_96_profile(res$signatures, condensed = TRUE)

plot_contribution(res$contribution, res$signature, mode = "relative")
plot_contribution(res$contribution, res$signature, mode = "absolute")

plot_contribution_heatmap(res$contribution,sig_order = c("Signature A", "Signature B", "Signature C"))

plot_compare_profiles(mut_mat[,1], res$reconstructed[,1],profile_names = c("Original", "Reconstructed"),condensed = TRUE)

### using curl for known sigs. having and issue with documentation
sp_url <- paste("http://cancer.sanger.ac.uk/cancergenome/assets/","signatures_probabilities.txt", sep = "")
cancer_signatures = read.table(sp_url, sep = "\t", header = TRUE)
cancer_signatures <-read.table(curl(sp_url),sep="\t",h=T)


new_order = match(row.names(mut_mat), cancer_signatures$Somatic.Mutation.Type)
cancer_signatures = cancer_signatures[as.vector(new_order),]
row.names(cancer_signatures) = cancer_signatures$Somatic.Mutation.Type
cancer_signatures = as.matrix(cancer_signatures[,4:33])

plot_96_profile(cancer_signatures[,1:2], condensed = TRUE, ymax = 0.3)

hclust_cosmic = cluster_signatures(cancer_signatures, method = "average")
cosmic_order = colnames(cancer_signatures)[hclust_cosmic$order]

plot(hclust_cosmic)


# similarity
cos_sim_samples_signatures = cos_sim_matrix(mut_mat, cancer_signatures)
cos_sim_sig_signatures = cos_sim_matrix(res$signatures, cancer_signatures)
plot_cosine_heatmap(cos_sim_samples_signatures, col_order = cosmic_order,cluster_rows = F)
plot_cosine_heatmap(cos_sim_samples_signatures, col_order = cosmic_order,cluster_rows = TRUE)
plot_cosine_heatmap(cos_sim_sig_signatures, col_order = cosmic_order,cluster_rows = F,plot_values = T)
fit_res <- fit_to_signatures(mut_mat, cancer_signatures)
select <- which(rowSums(fit_res$contribution) > 10)

plot_cosine_heatmap(cos_sim_samples_signatures, col_order = cosmic_order,cluster_rows = TRUE)
fit_res <- fit_to_signatures(mut_mat, cancer_signatures)
select <- which(rowSums(fit_res$contribution) > 10)

plot_contribution(fit_res$contribution[select,],cancer_signatures[,select], coord_flip = FALSE,mode = "absolute")

plot_contribution_heatmap(fit_res$contribution,cluster_samples = TRUE,method = "complete")

# Strand Bias

genes_hg19 <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
mut_mat_s <- mut_matrix_stranded(vcfs, ref_genome, genes_hg19)
strand_counts <- strand_occurrences(mut_mat_s, by=tissue)
strand_bias <- strand_bias_test(strand_counts)

plot_strand(strand_counts, mode = "relative")
plot_strand_bias(strand_bias)

