library(ggplot2); library(snpR)

d <- GeneArchEst::process_ms("../genomic_prediction_accuracy/data/theta4k_1000_10_rho40k.txt", 10000000)

keeps <- sample(nrow(d$meta), 10000, FALSE)
d$x <- d$x[keeps,]
d$meta <- d$meta[keeps,]
genotypes <- list(d$x[,1:250])

meta <- d$meta
effects <- rep(0, nrow(meta))

effects <- cbind(effects, effects, effects)
inbreeding <- c(10, -10, 0)

target_size <- 50
set.seed(123343)
tno <- simulate_populations(genotypes[rep(1, 3)], meta = meta, effects = effects, h = .5, gens = 10, chr_length = 10000000,
                            growth_function = function(n) BL_growth(n, 2), # use negative binomial instead
                            var_theta = 0,
                            K_thin_post_surv = target_size,
                            thin = FALSE,
                            thin_fixed = TRUE, 
                            sampling_point = "offspring",
                            inbreeding = inbreeding,
                            starting_surv_opt = c(0, 0, 0),
                            survival_function = function(phenotypes, ...) rep(1, length(phenotypes)),
                            selection_shift_function = function(opt, ...) opt,
                            verbose = TRUE)


genos <- cbind(tno$genotypes[[1]], tno$genotypes[[2]], tno$genotypes[[3]])
genos <- GeneArchEst::convert_2_to_1_column(genos)
genos <- t(genos)
d <- import.snpR.data(as.data.frame(genos), snp.meta = tno$meta, 
                      sample.meta = data.frame(pop = c(rep("Inbred", ncol(tno$genotypes[[1]])/2),
                                                       rep("Outbred", ncol(tno$genotypes[[2]])/2),
                                                       rep("Random", ncol(tno$genotypes[[3]])/2))))



set.seed (3218434)
d <- filter_snps(d, non_poly = TRUE)
d <- calc_fis(d, "pop")
get.snpR.stats(d, "pop", "fis")$weighted.means


n <- 25

d <- calc_seg_sites(d, "pop", g = n)
d <- calc_tajimas_d(d, "pop", global = TRUE)
res <- get.snpR.stats(d, "pop", c("tsd", "seg_sites"))$weighted.means
sfs1 <- calc_sfs(d, "pop", "Inbred", projection = n*2)
sfs2 <- calc_sfs(d, "pop", "Outbred", projection = n*2)
sfs3 <- calc_sfs(d, "pop", "Random", projection = n*2)

one_sub <- function(d, n){
  d_sub <- d[,c(sample(which(sample.meta(d)$pop == "Inbred"), n), 
                sample(which(sample.meta(d)$pop == "Outbred"), n), 
                sample(which(sample.meta(d)$pop == "Random"), n))]
  d_sub <- filter_snps(d_sub, non_poly = TRUE, verbose = FALSE)
  d_sub <- calc_seg_sites(d_sub, "pop")
  res_sub <- get.snpR.stats(d_sub, "pop", c("tsd", "seg_sites"))$weighted.means
  
  return(res_sub)
}


boots <- 1000
sub_res <- vector("list", boots)

cl <- parallel::makePSOCKcluster(5)
doParallel::registerDoParallel(cl)

pboot <- split(1:boots, rep(1:5, length.out = length(1:boots), each = ceiling(length(1:boots)/5)))


sub_res <- foreach(q = 1:5, .packages = "snpR") %dopar% {
  set.seed(3218214 + q)
  
  ti <- pboot[[q]]
  tres <- vector("list", length(ti))
  for(i in 1:length(ti)){
    tres[[i]] <- one_sub(d, n)
    tres[[i]]$iter <- i
  }
  tres <- rbindlist(tres)
  tres
}
parallel::stopCluster(cl)


sub_res <- data.table::rbindlist(sub_res)

colours <- khroma::color("batlow")
colors <- colours(4, range=c(0.2,1))

p1 <- ggplot(sub_res, aes(x = subfacet, y = seg_sites)) +
  geom_violin(draw_quantiles = c(.25, .5, .75)) +
  geom_point(data = data.frame(subfacet = c("Inbred", "Outbred", "Random"),
                               method = rep(c("E(N[S])", "Projection"), each = 3),
                               seg_sites = c(res$seg_sites,
                                             sum(sfs1[-1]),
                                             sum(sfs2[-1]),
                                             sum(sfs3[-1]))),
             aes(color = method), size = 2) +
  scale_color_manual(values = colors[c(1,4)]) +
  theme_bw() +
  ylab("S") +
  xlab("Demography")

sub_res <- as.data.table(sub_res)
rand_only <- sub_res[subfacet == "Random",]
sub_res[,ref :=  rand_only$seg_sites[match(sub_res$iter, rand_only$iter)], ]
sub_res[,prop := seg_sites/ref]

res <- as.data.table(res)

K <- n*2 #average sample size for ws.theta, as in hohenlohe 2010. Alternative would make this into a function, then use lapply on the vector of K values
#if(is.nan(K) == TRUE){browser()}
a1 <- sum(1/seq(1:(K-1))) #get a1
res[,ws.theta_adj := seg_sites/a1]


p2 <- ggplot(sub_res[subfacet != "Random",], aes(x = subfacet, y = prop)) +
  geom_violin(draw_quantiles = c(.25, .5, .75)) +
  geom_point(data = data.table(subfacet = c("Inbred", "Outbred", "Random"),
                               method = rep(c("E(N[S])", "Watterson's Theta", "E(N[S])-adjusted Watterson's Theta", "Projection"), each = 3),
                               prop = c(res$seg_sites/res$seg_sites[res$subfacet == "Random"],
                                        res$global_ws.theta/res$global_ws.theta[res$subfacet == "Random"],
                                        res$ws.theta_adj/res$ws.theta_adj[res$subfacet == "Random"],
                                        c(sum(sfs1[-1]),
                                          sum(sfs2[-1]),
                                          sum(sfs3[-1]))/sum(sfs3[-1])))[subfacet != "Random",],
             aes(color = method), position = position_dodge(width = .1), size = 2) +
  scale_color_manual(values = colors) +
  theme_bw() +
  ylab("Proportion S, relative to Neutral Demography") +
  xlab("Demography")


cowplot::plot_grid(p1, p2, nrow = 1, align = "hv", axis = "brlt")

