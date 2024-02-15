args <- commandArgs(TRUE)
q <- as.numeric(args[1])
outfile <- as.character(args[2])

library(snpR); library(ggplot2); library(data.table);

## ------------------------------------------------------------------------------------------------------------------------------
# Generate some example data:
## define parameters
seed <- 231266
threads <- 25
allele_frequencies <- seq(.001, .1, by = .001) # allele frequencies
n1 <- 5000 # p1 gene copy number
n2 <- 500 # p2 gene copy number
k <- seq(10, 100, by = 10) # numbers to rarefact to
mcmcs <- seq(1000, 10000, by = 1000) # number of mcmc draws
max_missing <- 30 # maximum percentage of missing data in a locus, from a uniform dist


## derived
set.seed(seed)
nal <- length(allele_frequencies)
miss_percent <- floor(runif(nal, 0, max_missing))/100


## ------------------------------------------------------------------------------------------------------------------------------
# generate genotypes
g1 <- matrix(0, nrow = nal, ncol = n1) # init
g2 <- matrix(0, nrow = nal, ncol = n2)

# add missing data
ms1 <- lapply(n1*miss_percent, function(x) sample(n1, x, FALSE))  # which are missing
ms2 <- lapply(n2*miss_percent, function(x) sample(n2, x, FALSE))

# determine minor allele placements
ma1 <- ma2 <- vector("list", length(ms1))
for(i in 1:nal){
  ma1[[i]] <- sample(n1*(1-miss_percent)[i], (n1*(1-miss_percent)*allele_frequencies)[i], FALSE)
  ma2[[i]] <- sample(n2*(1-miss_percent)[i], (n2*(1-miss_percent)*allele_frequencies)[i], FALSE)

}

# populate
## function to replace minors
populate <- function(g, missing, minors){
  if(length(missing) > 0){
    g[missing] <- NA
    g[-missing][minors] <- 1
  }
  else{
    g[minors] <- 1
  }
  return(g)
}
## do for every locus
for(i in 1:nal){ 
  g1[i,] <- populate(g1[i,], ms1[[i]], ma1[[i]])
  g2[i,] <- populate(g2[i,], ms2[[i]], ma2[[i]])
}

## combine into individuals
g1 <- g1[,1:(ncol(g1)/2)] + g1[,((ncol(g1)/2) + 1):ncol(g1)]
g2 <- g2[,1:(ncol(g2)/2)] + g2[,((ncol(g2)/2) + 1):ncol(g2)]



## ------------------------------------------------------------------------------------------------------------------------------
counts <- function(g){
  return(data.frame(majhom = rowSums(g == 0, na.rm = TRUE), 
                    het = rowSums(g == 1, na.rm = TRUE), 
                    minhom = rowSums(g == 2, na.rm = TRUE)))
}
rarefact <- function(counts, k){
  if(length(k) > 1){
    warning("length(k) > 1, using first element only.\n")
    k <- k[1]
  }
  return(1 - ((choose(counts[,1], k) + choose(counts[,3], k))/choose(rowSums(counts), k)))
}
varseg <- function(p) p*(1-p)

c1 <- counts(g1) # geno counts
c2 <- counts(g2)

# rarefact
pseg1 <- lapply(k, function(tk) rarefact(c1, tk))
pseg2 <- lapply(k, function(tk) rarefact(c2, tk))

vseg1 <- lapply(pseg1, varseg)
vseg2 <- lapply(pseg2, varseg)

p1 <- lapply(pseg1, function(p) sum(p)/nal)
p2 <- lapply(pseg2, function(p) sum(p)/nal)


## ------------------------------------------------------------------------------------------------------------------------------
set.seed(seed + q) # reset seed for cluster

rarefact_mcmc <- function(g, k){
  if(length(k) > 1){
    warning("length(k) > 1, using first element only.\n")
    k <- k[1]
  }
  
  miss <- is.na(g)
  if(any(miss)){
    g <- g[-which(miss)]
  }
  
  samp <- sample(g, k, FALSE)
  seg <- ifelse(any(samp == 1), 1, # got hets
                ifelse(any(samp == 0) & any(samp == 2), 1, 0)) # got both homs
  return(seg)
}

# draw for each k/mcmc count combination
combs <- expand.grid(k, mcmcs)
colnames(combs) <- c("k", "iters")

tmcmcs <- combs$iters[q]
tk <- combs$k[q]

seg1 <- seg2 <- matrix(0, nal, tmcmcs)
for(i in 1:tmcmcs){
  for(j in 1:nal){
    seg1[j,i] <- rarefact_mcmc(g1[j,], k = tk)
    seg2[j,i] <- rarefact_mcmc(g2[j,], k = tk)
  }
}

mpseg1 <- binom::binom.confint(rowSums(seg1), rep(tmcmcs, nrow(seg1)), 
                               method = "agresti-coull")[,c("mean", "lower", "upper")]
mpseg2 <- binom::binom.confint(rowSums(seg2), rep(tmcmcs, nrow(seg1)),
                               method = "agresti-coull")[,c("mean", "lower", "upper")]

colnames(mpseg1) <- colnames(mpseg2) <- c("mean", "lower", "upper")

res <- list(seg1 = seg1, seg2 = seg2, mpseg1 = mpseg1, mpseg2 = mpseg2, k = tk, iters = tmcmcs)

saveRDS(res, paste0(outfile, "_", q, ".RDS"))
