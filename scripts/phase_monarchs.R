library(snpR); library(ggplot2)
d <- readRDS("../ConGen/microConGen_RW/Data/monarch_nomaf.RDS")

chrs <- table(snp.meta(d)$group)
bad_chrs <- names(chrs)[which(chrs == 1)]
bad_snps <- which(snp.meta(d)$group %in% bad_chrs)
d <- import.snpR.data(genotypes(d)[-bad_snps,], snp.meta(d)[-bad_snps,], sample.meta(d))

format_snps(d, "vcf", chr = "group", outfile = "data/monarchs.vcf")
cmd <- "java -Xmx4g -jar ~/bin/beagle.22Jul22.46e.jar gt=data/monarchs.vcf out=data/monarchs_phased"
system(cmd)
system("gunzip data/monarchs_phased.vcf.gz")


dv <- data.table::fread("data/monarchs_phased.vcf")
pops <- colnames(dv)[-c(1:9)]
pops <- substr(pops, 1, 3)


meta <- dv[,1:2]
dv <- dv[,-c(1:9), with = FALSE]


genotypes <- list(nam = dv[,which(pops == "NAM"), with = FALSE],
                  qld = dv[,which(pops == "QLD"), with = FALSE])
genotypes <- lapply(genotypes, function(x){
  x1 <- matrix(as.numeric(substr(as.matrix(x), 1, 1)), nrow = nrow(x))
  x2 <- matrix(as.numeric(substr(as.matrix(x), 3, 3)), nrow = nrow(x))
  return(cbind(x1, x2)[,order(c(1:ncol(x1), 1:ncol(x2)))])
})

effects <- rep(0, nrow(dv))
effects <- cbind(effects, effects)


target_size <- 50
migration <- matrix(.5, 2, 2)
chr_lengths <- tapply(meta$POS, meta$`#CHROM`, max)
colnames(meta) <- c("chr", "position")
set.seed(123343)
tno <- simulate_populations(genotypes, meta = meta, effects = effects, h = .5, gens = 1, chr_length = chr_lengths,
                            growth_function = function(n) BL_growth(n, 2),
                            var_theta = 0,
                            K_thin_post_surv = target_size,
                            thin = FALSE,
                            thin_fixed = FALSE, 
                            sampling_point = "offspring",
                            starting_surv_opt = c(0, 0),
                            migration = migration,
                            survival_function = function(phenotypes, ...) rep(1, length(phenotypes)),
                            selection_shift_function = function(opt, ...) opt,
                            verbose = TRUE)

tno2 <- simulate_populations(tno$genotypes[1], meta = tno$meta, effects = tno$effects[,1, drop = FALSE],
                             h = .5, gens = 5, chr_length = chr_lengths,
                             growth_function = function(n) BL_growth(n, 2),
                             var_theta = 0,
                             K_thin_post_surv = target_size,
                             thin = FALSE,
                             thin_fixed = FALSE, 
                             sampling_point = "offspring",
                             starting_surv_opt = c(0),
                             survival_function = function(phenotypes, ...) rep(1, length(phenotypes)),
                             selection_shift_function = function(opt, ...) opt,
                             verbose = TRUE)


genos <- cbind(genotypes[[1]], genotypes[[2]], tno$genotypes[[1]], tno$genotypes[[2]], tno2$genotypes[[1]])
genos <- GeneArchEst::convert_2_to_1_column(genos)
genos <- t(genos)
d <- import.snpR.data(as.data.frame(genos), snp.meta = tno$meta, 
                      sample.meta = data.frame(pop = c(rep("NAM", ncol(genotypes[[1]])/2),
                                                       rep("QLD", ncol(genotypes[[2]])/2), 
                                                       rep("OutbredA", ncol(tno$genotypes[[1]])/2),
                                                       rep("OutbredB", ncol(tno$genotypes[[2]])/2),
                                                       rep("Random", ncol(tno2$genotypes[[1]])/2))))




d <- d[pop = c("NAM", "OutbredB", "Random")]
d <- filter_snps(d, non_poly = TRUE)
set.seed (3218434)
d <- d[sample(nrow(d), 10000, FALSE),]
gc();

d <- calc_fis(d, "pop")
get.snpR.stats(d, "pop", "fis")$weighted.means




n <- 25

d <- calc_seg_sites(d, "pop", g = n)
d <- calc_tajimas_d(d, "pop", global = TRUE)
d <- calc_fis(d, "pop")
res <- get.snpR.stats(d, "pop", c("tsd", "seg_sites", "fis"))$weighted.means
sfs1 <- calc_sfs(d, "pop", "NAM", projection = n*2)
sfs2 <- calc_sfs(d, "pop", "OutbredB", projection = n*2)
sfs3 <- calc_sfs(d, "pop", "Random", projection = n*2)

one_sub <- function(d, n){
  d_sub <- d[,c(sample(which(sample.meta(d)$pop == "NAM"), n), 
                sample(which(sample.meta(d)$pop == "OutbredB"), n), 
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
sub_res[,fis := res$weighted_mean_fis[match(sub_res$subfacet, res$subfacet)]]
saveRDS(sub_res, "data/sub_res_theta_proj_comp.RDS")
saveRDS(res, "data/res_theta_proj_comp.RDS")


colours <- khroma::color("batlow")
colors <- colours(4, range=c(0.2,1))

pd <- data.table(subfacet = c("NAM", "OutbredB", "Random"),
                 `"E(" * "N"["S"] * ")"` = res$seg_sites,
                 `"Projection"` = c(sum(sfs1[-1]),
                                sum(sfs2[-1]),
                                sum(sfs3[-1])),
                 fis = res$weighted_mean_fis,
                 `"Watterson" * "'" * "s " * theta` = res$global_ws.theta)



K <- n*2 #average sample size for ws.theta, as in hohenlohe 2010. Alternative would make this into a function, then use lapply on the vector of K values
#if(is.nan(K) == TRUE){browser()}
a1 <- sum(1/seq(1:(K-1))) #get a1
pd[, `"Watterson" * "'" * "s " * theta["adjusted"]` := res$seg_sites/a1]
pd[, `"Watterson" * "'" * "s " * theta["adjusted"]` :=`"Watterson" * "'" * "s " * theta["adjusted"]`/mean(`"Watterson" * "'" * "s " * theta["adjusted"]`)]
pd <- melt(pd, id.vars = c("subfacet", "fis"))
pd[,value := value/mean(value), by = variable]


sub_res <- as.data.table(sub_res)
sub_res[,value := seg_sites/mean(seg_sites)]

p1 <- ggplot(sub_res, aes(x = fis, y = log10(value), group = fis)) +
  geom_violin(draw_quantiles = c(.25, .5, .75)) +
  geom_point(data = pd,
             aes(color = variable), size = 2, position = position_dodge2(width = .005)) +
  scale_color_manual(values = colors, labels = scales::parse_format()) +
  geom_vline(xintercept = 0, linetype = "dashed") + 
  theme_bw() %+replace%
  theme(panel.border = element_blank(),
        axis.line = element_line(colour = "black", 
                     linewidth = rel(1))) +
  ggbreak::scale_x_break(c(0.005, 0.135), ticklabels = c(.135, .14, .145)) +
  scale_x_continuous(breaks = c(-.015, -.010, -.005, 0)) +
  ylab(expression(theta * ' or S/' * mu[theta * ' or S' ])) +
  xlab(expression(F[italic(IS)])) +
  guides(color = guide_legend(title = "Method"))
p1
