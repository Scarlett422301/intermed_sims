# parameters
ns <- c(2000, 4000)
seed <- 1:1000
parm <- expand.grid(n = ns, seed = seed)

save_dir <- "/Users/scarlett/Desktop/nonparametric/submission2" 
n_out <- 157

rslt <- matrix(NA, nrow = nrow(parm), ncol = n_out)
for(i in 1:nrow(parm)){
  tmp <- tryCatch({
    load(paste0(save_dir, "/all_fits/jialufit_", i, ".RData"))
    out
  }, error = function(e){
    rep(NA, n_out)
  })
  rslt[i, ] <- tmp
}

out <- data.frame(rslt)
colnames(out) <- c("n", "seed",
                   "total_plugin_est", "direct_plugin_est", "indirectM1_plugin_est", "indirectM2_plugin_est", "covarM1M2_plugin_est",
                   "total_parametric_plugin_est", "direct_parametric_plugin_est", "indirectM1_parametric_plugin_est", "indirectM2_parametric_plugin_est", "covarM1M2_parametric_plugin_est",
                   "total_tmle_est", "total_tmle_cil", "total_tmle_ciu",
                   "direct_tmle_est","direct_tmle_cil","direct_tmle_ciu",
                   "indirectM1_tmle_est","indirectM1_tmle_cil","indirectM1_tmle_ciu",
                   "indirectM2_tmle_est","indirectM2_tmle_cil","indirectM2_tmle_ciu",
                   "covarM1M2_tmle_est","covarM1M2_tmle_cil","covarM1M2_tmle_ciu",
                   "total_aiptw_est", "total_aiptw_cil", "total_aiptw_ciu",
                   "direct_aiptw_est","direct_aiptw_cil","direct_aiptw_ciu",
                   "indirectM1_aiptw_est","indirectM1_aiptw_cil","indirectM1_aiptw_ciu",
                   "indirectM2_aiptw_est","indirectM2_aiptw_cil","indirectM2_aiptw_ciu",
                   "covarM1M2_aiptw_est","covarM1M2_aiptw_cil","covarM1M2_aiptw_ciu",
                   "total_tmle_cover", "direct_tmle_cover", "indirectM1_tmle_cover", "indirectM2_tmle_cover", "covarM1M2_tmle_cover",
                   "total_aiptw_cover", "direct_aiptw_cover", "indirectM1_aiptw_cover", "indirectM2_aiptw_cover", "covarM1M2_aiptw_cover",
                   "total_parametric_tmle_est", "total_parametric_tmle_cil", "total_parametric_tmle_ciu",
                   "direct_parametric_tmle_est","direct_parametric_tmle_cil","direct_parametric_tmle_ciu",
                   "indirectM1_parametric_tmle_est","indirectM1_parametric_tmle_cil","indirectM1_parametric_tmle_ciu",
                   "indirectM2_parametric_tmle_est","indirectM2_parametric_tmle_cil","indirectM2_parametric_tmle_ciu",
                   "covarM1M2_parametric_tmle_est","covarM1M2_parametric_tmle_cil","covarM1M2_parametric_tmle_ciu",
                   "total_parametric_aiptw_est", "total_parametric_aiptw_cil", "total_parametric_aiptw_ciu",
                   "direct_parametric_aiptw_est","direct_parametric_aiptw_cil","direct_parametric_aiptw_ciu",
                   "indirectM1_parametric_aiptw_est","indirectM1_parametric_aiptw_cil","indirectM1_parametric_aiptw_ciu",
                   "indirectM2_parametric_aiptw_est","indirectM2_parametric_aiptw_cil","indirectM2_parametric_aiptw_ciu",
                   "covarM1M2_parametric_aiptw_est","covarM1M2_parametric_aiptw_cil","covarM1M2_parametric_aiptw_ciu",
                   "total_parametric_tmle_cover", "direct_parametric_tmle_cover", "indirectM1_parametric_tmle_cover", "indirectM2_parametric_tmle_cover", "covarM1M2_parametric_tmle_cover",
                   "total_parametric_aiptw_cover", "direct_parametric_aiptw_cover", "indirectM1_parametric_aiptw_cover", "indirectM2_parametric_aiptw_cover", "covarM1M2_parametric_aiptw_cover",
                   "total_gcomp_est0", "direct_gcomp_est0", "indirectM1_gcomp_est0", "indirectM2_gcomp_est0", "covarM1M2_gcomp_est0",
                   "total_gcomp_cil0", "direct_gcomp_cil0", "indirectM1_gcomp_cil0", "indirectM2_gcomp_cil0", "covarM1M2_gcomp_cil0",
                   "total_gcomp_ciu0", "direct_gcomp_ciu0", "indirectM1_gcomp_ciu0", "indirectM2_gcomp_ciu0", "covarM1M2_gcomp_ciu0",
                   "total_gcomp_cover0", "direct_gcomp_cover0", "indirectM1_gcomp_cover0", "indirectM2_gcomp_cover0", "covarM1M2_gcomp_cover0",
                   "total_gcomp_est1", "direct_gcomp_est1", "indirectM1_gcomp_est1", "indirectM2_gcomp_est1", "covarM1M2_gcomp_est1",
                   "total_gcomp_cil1", "direct_gcomp_cil1", "indirectM1_gcomp_cil1", "indirectM2_gcomp_cil1", "covarM1M2_gcomp_cil1",
                   "total_gcomp_ciu1", "direct_gcomp_ciu1", "indirectM1_gcomp_ciu1", "indirectM2_gcomp_ciu1", "covarM1M2_gcomp_ciu1",
                   "total_gcomp_cover1", "direct_gcomp_cover1", "indirectM1_gcomp_cover1", "indirectM2_gcomp_cover1", "covarM1M2_gcomp_cover1",
                   "total_gcomp_est2", "direct_gcomp_est2", "indirectM1_gcomp_est2", "indirectM2_gcomp_est2", "covarM1M2_gcomp_est2",
                   "total_gcomp_cil2", "direct_gcomp_cil2", "indirectM1_gcomp_cil2", "indirectM2_gcomp_cil2", "covarM1M2_gcomp_cil2",
                   "total_gcomp_ciu2", "direct_gcomp_ciu2", "indirectM1_gcomp_ciu2", "indirectM2_gcomp_ciu2", "covarM1M2_gcomp_ciu2",
                   "total_gcomp_cover2", "direct_gcomp_cover2", "indirectM1_gcomp_cover2", "indirectM2_gcomp_cover2", "covarM1M2_gcomp_cover2",
                   "truth_total", "truth_direct", "truth_indirectM1", "truth_indirectM2", "truth_covarM1M2"
) 

save(out, file = "/Users/scarlett/Desktop/nonparametric/submission2/all_out_intermed_new1.RData")

get_bias <- function(out, which_eff = "total"){
  by(out, out$n, function(x){
    tmle <- mean(x[ , grep(paste0(which_eff, "_tmle_est"), colnames(x))] - x[,paste0("truth_", which_eff)], na.rm = TRUE)
    aiptw <- mean(x[ , grep(paste0(which_eff, "_aiptw_est"), colnames(x))] - x[,paste0("truth_", which_eff)], na.rm = TRUE)
    tmle2 <- mean(x[ , grep(paste0(which_eff, "_parametric_tmle_est"), colnames(x))] - x[,paste0("truth_", which_eff)], na.rm = TRUE)
    aiptw2 <- mean(x[ , grep(paste0(which_eff, "_parametric_aiptw_est"), colnames(x))] - x[,paste0("truth_", which_eff)], na.rm = TRUE)
    plugin <- mean(x[ , grep(paste0(which_eff, "_plugin_est"), colnames(x))] - x[,paste0("truth_", which_eff)], na.rm = TRUE)
    plugin2 <- mean(x[ , grep(paste0(which_eff, "_parametric_plugin_est"), colnames(x))] - x[,paste0("truth_", which_eff)], na.rm = TRUE)
    gcomp0 <- mean(x[ , grep(paste0(which_eff, "_gcomp_est0"), colnames(x))] - x[,paste0("truth_", which_eff)], na.rm = TRUE)
    gcomp1 <- mean(x[ , grep(paste0(which_eff, "_gcomp_est1"), colnames(x))] - x[,paste0("truth_", which_eff)], na.rm = TRUE)
    gcomp2 <- mean(x[ , grep(paste0(which_eff, "_gcomp_est2"), colnames(x))] - x[,paste0("truth_", which_eff)], na.rm = TRUE)
    return(c(tmle, aiptw, tmle2, aiptw2, plugin, plugin2, gcomp0, gcomp1, gcomp2))
  })
}

get_rootnbias <- function(out, which_eff = "total"){
  by(out, out$n, function(x){
    tmle <- mean(x$n[1]^(1/2)*(x[ , grep(paste0(which_eff, "_tmle_est"), colnames(x))] - x[,paste0("truth_", which_eff)]), na.rm = TRUE)
    aiptw <- mean(x$n[1]^(1/2)*(x[ , grep(paste0(which_eff, "_aiptw_est"), colnames(x))] - x[,paste0("truth_", which_eff)]), na.rm = TRUE)
    tmle2 <- mean(x$n[1]^(1/2)*(x[ , grep(paste0(which_eff, "_parametric_tmle_est"), colnames(x))] - x[,paste0("truth_", which_eff)]), na.rm = TRUE)
    aiptw2 <- mean(x$n[1]^(1/2)*(x[ , grep(paste0(which_eff, "_parametric_aiptw_est"), colnames(x))] - x[,paste0("truth_", which_eff)]), na.rm = TRUE)
    plugin <- mean(x$n[1]^(1/2)*(x[ , grep(paste0(which_eff, "_plugin_est"), colnames(x))] - x[,paste0("truth_", which_eff)]), na.rm = TRUE)
    plugin2 <- mean(x$n[1]^(1/2)*(x[ , grep(paste0(which_eff, "_parametric_plugin_est"), colnames(x))] - x[,paste0("truth_", which_eff)]), na.rm = TRUE)
    gcomp0 <- mean(x$n[1]^(1/2)*(x[ , grep(paste0(which_eff, "_gcomp_est0"), colnames(x))] - x[,paste0("truth_", which_eff)]), na.rm = TRUE)
    gcomp1 <- mean(x$n[1]^(1/2)*(x[ , grep(paste0(which_eff, "_gcomp_est1"), colnames(x))] - x[,paste0("truth_", which_eff)]), na.rm = TRUE)
    gcomp2 <- mean(x$n[1]^(1/2)*(x[ , grep(paste0(which_eff, "_gcomp_est2"), colnames(x))] - x[,paste0("truth_", which_eff)]), na.rm = TRUE)
    return(c(tmle, aiptw, tmle2, aiptw2, plugin, plugin2, gcomp0, gcomp1, gcomp2))
  })
}

get_mse <- function(out, which_eff = "total"){
  by(out, out$n, function(x){
    tmle <- mean((x[ , grep(paste0(which_eff, "_tmle_est"), colnames(x))] - x[,paste0("truth_", which_eff)])^2, na.rm = TRUE)
    aiptw <- mean((x[ , grep(paste0(which_eff, "_aiptw_est"), colnames(x))] - x[,paste0("truth_", which_eff)])^2, na.rm = TRUE)
    tmle2 <- mean((x[ , grep(paste0(which_eff, "_parametric_tmle_est"), colnames(x))] - x[,paste0("truth_", which_eff)])^2, na.rm = TRUE)
    aiptw2 <- mean((x[ , grep(paste0(which_eff, "_parametric_aiptw_est"), colnames(x))] - x[,paste0("truth_", which_eff)])^2, na.rm = TRUE)
    plugin <- mean((x[ , grep(paste0(which_eff, "_plugin_est"), colnames(x))] - x[,paste0("truth_", which_eff)])^2, na.rm = TRUE)
    plugin2 <- mean((x[ , grep(paste0(which_eff, "_parametric_plugin_est"), colnames(x))] - x[,paste0("truth_", which_eff)])^2, na.rm = TRUE)
    gcomp0 <- mean((x[ , grep(paste0(which_eff, "_gcomp_est0"), colnames(x))] - x[,paste0("truth_", which_eff)])^2, na.rm = TRUE)
    gcomp1 <- mean((x[ , grep(paste0(which_eff, "_gcomp_est1"), colnames(x))] - x[,paste0("truth_", which_eff)])^2, na.rm = TRUE)
    gcomp2 <- mean((x[ , grep(paste0(which_eff, "_gcomp_est2"), colnames(x))] - x[,paste0("truth_", which_eff)])^2, na.rm = TRUE)
    return(c(tmle, aiptw, tmle2, aiptw2, plugin, plugin2, gcomp0, gcomp1, gcomp2))
  })
}

get_nmse <- function(out, which_eff = "total"){
  by(out, out$n, function(x){
    tmle <- x$n[1]*mean((x[ , grep(paste0(which_eff, "_tmle_est"), colnames(x))] - x[,paste0("truth_", which_eff)])^2, na.rm = TRUE)
    aiptw <- x$n[1]*mean((x[ , grep(paste0(which_eff, "_aiptw_est"), colnames(x))] - x[,paste0("truth_", which_eff)])^2, na.rm = TRUE)
    tmle2 <- x$n[1]*mean((x[ , grep(paste0(which_eff, "_parametric_tmle_est"), colnames(x))] - x[,paste0("truth_", which_eff)])^2, na.rm = TRUE)
    aiptw2 <- x$n[1]*mean((x[ , grep(paste0(which_eff, "_parametric_aiptw_est"), colnames(x))] - x[,paste0("truth_", which_eff)])^2, na.rm = TRUE)
    plugin <- x$n[1]*mean((x[ , grep(paste0(which_eff, "_plugin_est"), colnames(x))] - x[,paste0("truth_", which_eff)])^2, na.rm = TRUE)
    plugin2 <- x$n[1]*mean((x[ , grep(paste0(which_eff, "_parametric_plugin_est"), colnames(x))] - x[,paste0("truth_", which_eff)])^2, na.rm = TRUE)
    gcomp0 <- x$n[1]*mean((x[ , grep(paste0(which_eff, "_gcomp_est0"), colnames(x))] - x[,paste0("truth_", which_eff)])^2, na.rm = TRUE)
    gcomp1 <- x$n[1]*mean((x[ , grep(paste0(which_eff, "_gcomp_est1"), colnames(x))] - x[,paste0("truth_", which_eff)])^2, na.rm = TRUE)
    gcomp2 <- x$n[1]*mean((x[ , grep(paste0(which_eff, "_gcomp_est2"), colnames(x))] - x[,paste0("truth_", which_eff)])^2, na.rm = TRUE)
    return(c(tmle, aiptw, tmle2, aiptw2, plugin, plugin2, gcomp0, gcomp1, gcomp2))
  })
}
get_sd <- function(out, which_eff = "total"){
  by(out, out$n, function(x){
    tmle <- sd(x[ , grep(paste0(which_eff, "_tmle_est"), colnames(x))], na.rm = TRUE)
    aiptw <- sd(x[ , grep(paste0(which_eff, "_aiptw_est"), colnames(x))] , na.rm = TRUE)
    tmle2 <- sd(x[ , grep(paste0(which_eff, "_parametric_tmle_est"), colnames(x))], na.rm = TRUE)
    aiptw2 <- sd(x[ , grep(paste0(which_eff, "_parametric_aiptw_est"), colnames(x))] , na.rm = TRUE)
    plugin <- sd(x[ , grep(paste0(which_eff, "_plugin_est"), colnames(x))], na.rm = TRUE)
    plugin2 <- sd(x[ , grep(paste0(which_eff, "_parametric_plugin_est"), colnames(x))] , na.rm = TRUE)
    gcomp0 <- sd(x[ , grep(paste0(which_eff, "_gcomp_est0"), colnames(x))] , na.rm = TRUE)
    gcomp1 <- sd(x[ , grep(paste0(which_eff, "_gcomp_est1"), colnames(x))] , na.rm = TRUE)
    gcomp2 <- sd(x[ , grep(paste0(which_eff, "_gcomp_est2"), colnames(x))] , na.rm = TRUE)
    return(c(tmle, aiptw, tmle2, aiptw2, plugin, plugin2, gcomp0, gcomp1, gcomp2))
  })
}
get_rootnsd <- function(out, which_eff = "total"){
  by(out, out$n, function(x){
    tmle <- x$n[1]^(1/2)*sd(x[ , grep(paste0(which_eff, "_tmle_est"), colnames(x))], na.rm = TRUE)
    aiptw <- x$n[1]^(1/2)*sd(x[ , grep(paste0(which_eff, "_aiptw_est"), colnames(x))] , na.rm = TRUE)
    tmle2 <- x$n[1]^(1/2)*sd(x[ , grep(paste0(which_eff, "_parametric_tmle_est"), colnames(x))], na.rm = TRUE)
    aiptw2 <- x$n[1]^(1/2)*sd(x[ , grep(paste0(which_eff, "_parametric_aiptw_est"), colnames(x))] , na.rm = TRUE)
    plugin <- x$n[1]^(1/2)*sd(x[ , grep(paste0(which_eff, "_plugin_est"), colnames(x))], na.rm = TRUE)
    plugin2 <- x$n[1]^(1/2)*sd(x[ , grep(paste0(which_eff, "_parametric_plugin_est"), colnames(x))] , na.rm = TRUE)
    gcomp0 <- x$n[1]^(1/2)*sd(x[ , grep(paste0(which_eff, "_gcomp_est0"), colnames(x))] , na.rm = TRUE)
    gcomp1 <- x$n[1]^(1/2)*sd(x[ , grep(paste0(which_eff, "_gcomp_est1"), colnames(x))] , na.rm = TRUE)
    gcomp2 <- x$n[1]^(1/2)*sd(x[ , grep(paste0(which_eff, "_gcomp_est2"), colnames(x))] , na.rm = TRUE)
    return(c(tmle, aiptw, tmle2, aiptw2, plugin, plugin2, gcomp0, gcomp1, gcomp2))
  })
}
get_se_est <- function(out, which_eff = "total"){
  by(out, out$n, function(x){
    tmle <- mean(x[ , grep(paste0(which_eff, "_tmle_ciu"), colnames(x))] - x[ , grep(paste0(which_eff, "_tmle_cil"), colnames(x))]) / (2 * 1.96)
    aiptw <- mean(x[ , grep(paste0(which_eff, "_aiptw_ciu"), colnames(x))] - x[ , grep(paste0(which_eff, "_aiptw_cil"), colnames(x))]) / (2 * 1.96)
    tmle2 <- mean(x[ , grep(paste0(which_eff, "_parametric_tmle_ciu"), colnames(x))] - x[ , grep(paste0(which_eff, "_parametric_tmle_cil"), colnames(x))]) / (2 * 1.96)
    aiptw2 <- mean(x[ , grep(paste0(which_eff, "_parametric_aiptw_ciu"), colnames(x))] - x[ , grep(paste0(which_eff, "_parametric_aiptw_cil"), colnames(x))]) / (2 * 1.96)
    # plugin <- mean(x[ , grep(paste0(which_eff, "_plugin_ciu"), colnames(x))] - x[ , grep(paste0(which_eff, "_plugin_cil"), colnames(x))]) / (2 * 1.96)
    # plugin2 <- mean(x[ , grep(paste0(which_eff, "_parametric_plugin_ciu"), colnames(x))] - x[ , grep(paste0(which_eff, "_parametric_plugin_cil"), colnames(x))]) / (2 * 1.96)
    gcomp0 <- mean(x[ , grep(paste0(which_eff, "_gcomp_ciu0"), colnames(x))] - x[ , grep(paste0(which_eff, "_gcomp_cil0"), colnames(x))]) / (2 * 1.96)
    gcomp1 <- mean(x[ , grep(paste0(which_eff, "_gcomp_ciu1"), colnames(x))] - x[ , grep(paste0(which_eff, "_gcomp_cil1"), colnames(x))]) / (2 * 1.96)
    gcomp2 <- mean(x[ , grep(paste0(which_eff, "_gcomp_ciu2"), colnames(x))] - x[ , grep(paste0(which_eff, "_gcomp_cil2"), colnames(x))]) / (2 * 1.96)
    # return(c(tmle, aiptw, tmle2, aiptw2, plugin, plugin2, gcomp0, gcomp1, gcomp2))
    return(c(tmle, aiptw, tmle2, aiptw2, gcomp0, gcomp1, gcomp2))
  })
}
get_cover <- function(out, which_eff = "total"){
  by(out, out$n, function(x){
    tmle <- mean(x[ , grep(paste0(which_eff, "_tmle_cover"), colnames(x))], na.rm = TRUE)
    aiptw <- mean(x[ , grep(paste0(which_eff, "_aiptw_cover"), colnames(x))], na.rm = TRUE)
    # tmle2 <- mean(x[ , grep(paste0(which_eff, "_parametric_tmle_cover"), colnames(x))], na.rm = TRUE)
    # aiptw2 <- mean(x[ , grep(paste0(which_eff, "_parametric_aiptw_cover"), colnames(x))], na.rm = TRUE)
    gcomp0 <- mean(x[ , grep(paste0(which_eff, "_gcomp_cover0"), colnames(x))], na.rm = TRUE)
    gcomp1 <- mean(x[ , grep(paste0(which_eff, "_gcomp_cover1"), colnames(x))], na.rm = TRUE)
    gcomp2 <- mean(x[ , grep(paste0(which_eff, "_gcomp_cover2"), colnames(x))], na.rm = TRUE)
    return(c(tmle, aiptw, gcomp0, gcomp1, gcomp2))
  })
}

get_density_orc <- function(out, which_eff = "total"){
  by(out, out$n, function(x){
    tmle <- density((x[ , grep(paste0(which_eff, "_tmle_est"), colnames(x))] - x[,paste0("truth_", which_eff)]) / sd(x[ , grep(paste0(which_eff, "_tmle_est"), colnames(x))]), na.rm = TRUE)
    aiptw <- density((x[ , grep(paste0(which_eff, "_aiptw_est"), colnames(x))] - x[,paste0("truth_", which_eff)]) / sd(x[ , grep(paste0(which_eff, "_aiptw_est"), colnames(x))]), na.rm = TRUE)
    tmle2 <- density((x[ , grep(paste0(which_eff, "_parametric_tmle_est"), colnames(x))] - x[,paste0("truth_", which_eff)]) / sd(x[ , grep(paste0(which_eff, "_parametric_tmle_est"), colnames(x))]), na.rm = TRUE)
    aiptw2 <- density((x[ , grep(paste0(which_eff, "_parametric_aiptw_est"), colnames(x))] - x[,paste0("truth_", which_eff)]) / sd(x[ , grep(paste0(which_eff, "_parametric_aiptw_est"), colnames(x))]), na.rm = TRUE)
    return(list(tmle, aiptw, tmle2, aiptw2))
  })
}  
get_density_est <- function(out, which_eff = "total"){
  by(out, out$n, function(x){
    tmle_se <- (x[ , grep(paste0(which_eff, "_tmle_ciu"), colnames(x))] - x[ , grep(paste0(which_eff, "_tmle_cil"), colnames(x))]) / (2 * 1.96)
    aiptw_se <- (x[ , grep(paste0(which_eff, "_aiptw_ciu"), colnames(x))] - x[ , grep(paste0(which_eff, "_aiptw_cil"), colnames(x))]) / (2 * 1.96)
    tmle <- density((x[ , grep(paste0(which_eff, "_tmle_est"), colnames(x))] - x[,paste0("truth_", which_eff)]) / tmle_se, na.rm = TRUE)
    aiptw <- density((x[ , grep(paste0(which_eff, "_aiptw_est"), colnames(x))] - x[,paste0("truth_", which_eff)]) / aiptw_se, na.rm = TRUE)
    tmle_se2 <- (x[ , grep(paste0(which_eff, "_parametric_tmle_ciu"), colnames(x))] - x[ , grep(paste0(which_eff, "_parametric_tmle_cil"), colnames(x))]) / (2 * 1.96)
    aiptw_se2 <- (x[ , grep(paste0(which_eff, "_parametric_aiptw_ciu"), colnames(x))] - x[ , grep(paste0(which_eff, "_parametric_aiptw_cil"), colnames(x))]) / (2 * 1.96)
    tmle2 <- density((x[ , grep(paste0(which_eff, "_parametric_tmle_est"), colnames(x))] - x[,paste0("truth_", which_eff)]) / tmle_se2, na.rm = TRUE)
    aiptw2 <- density((x[ , grep(paste0(which_eff, "_parametric_aiptw_est"), colnames(x))] - x[,paste0("truth_", which_eff)]) / aiptw_se2, na.rm = TRUE)
    return(list(tmle, aiptw, tmle2, aiptw2))
  })
}

# check oracle CI coverage
check_oracle_ci <- function(out, which_eff = "total"){
  by(out, out$n, function(x){
    oracle_se_tmle <- sd(x[, paste0(which_eff, "_tmle_est")], na.rm = TRUE)
    oracle_se_aiptw <- sd(x[, paste0(which_eff, "_aiptw_est")], na.rm = TRUE)

    tmle_cil <- x[, paste0(which_eff, "_tmle_est")] - 1.96 * oracle_se_tmle
    tmle_ciu <- x[, paste0(which_eff, "_tmle_est")] + 1.96 * oracle_se_tmle
    aiptw_cil <- x[, paste0(which_eff, "_aiptw_est")] - 1.96 * oracle_se_aiptw
    aiptw_ciu <- x[, paste0(which_eff, "_aiptw_est")] + 1.96 * oracle_se_aiptw
    tmle_cover <- tmle_cil < x[, paste0("truth_", which_eff)] & tmle_ciu > x[, paste0("truth_", which_eff)]
    aiptw_cover <- aiptw_cil < x[, paste0("truth_", which_eff)] & aiptw_ciu > x[, paste0("truth_", which_eff)]

    oracle_se_gcomp0 <- sd(x[, paste0(which_eff, "_gcomp_est0")], na.rm = TRUE)
    oracle_se_gcomp1 <- sd(x[, paste0(which_eff, "_gcomp_est1")], na.rm = TRUE)
    oracle_se_gcomp2 <- sd(x[, paste0(which_eff, "_gcomp_est2")], na.rm = TRUE)

    gcomp0_cil <- x[, paste0(which_eff, "_gcomp_est0")] - 1.96 * oracle_se_gcomp0
    gcomp0_ciu <- x[, paste0(which_eff, "_gcomp_est0")] + 1.96 * oracle_se_gcomp0
    gcomp1_cil <- x[, paste0(which_eff, "_gcomp_est1")] - 1.96 * oracle_se_gcomp1
    gcomp1_ciu <- x[, paste0(which_eff, "_gcomp_est1")] + 1.96 * oracle_se_gcomp1
    gcomp2_cil <- x[, paste0(which_eff, "_gcomp_est2")] - 1.96 * oracle_se_gcomp2
    gcomp2_ciu <- x[, paste0(which_eff, "_gcomp_est2")] + 1.96 * oracle_se_gcomp2

    gcomp0_cover <- gcomp0_cil < x[, paste0("truth_", which_eff)] & gcomp0_ciu > x[, paste0("truth_", which_eff)]
    gcomp1_cover <- gcomp1_cil < x[, paste0("truth_", which_eff)] & gcomp1_ciu > x[, paste0("truth_", which_eff)]
    gcomp2_cover <- gcomp2_cil < x[, paste0("truth_", which_eff)] & gcomp2_ciu > x[, paste0("truth_", which_eff)]

    return(c(mean(tmle_cover, na.rm = TRUE), mean(aiptw_cover, na.rm = TRUE),
             mean(gcomp0_cover, na.rm = TRUE), mean(gcomp1_cover, na.rm = TRUE), mean(gcomp2_cover, na.rm = TRUE)))
  })
}

format_out <- function(out, summary_fn,
                       all_eff = c("total", "direct", "indirectM1", 
                                   "indirectM2", "covarM1M2")){
  if(!grepl("get_density", summary_fn)){    
    all_out <- sapply(all_eff, function(x){
      suppressWarnings(a <- data.frame(Reduce(rbind, 
                                              do.call(summary_fn, 
                                                      args = list(out = out, 
                                                                  which_eff = x))),
                                       n = c(2000, 4000),
                                       stringsAsFactors = FALSE))
      # 2 means parametric
      num_cols <- ncol(a)
      if(ncol(a) == 5){
          colnames(a) <- c("tmle", "os", "tmle2", "os2", "n")
      }else if(ncol(a) == 6){
          colnames(a) <- c("tmle", "os", "gcomp0", "gcomp1", "gcomp2", "n")
      }else if(ncol(a) == 8){
          colnames(a) <- c("tmle", "os", "tmle2", "os2", "gcomp0", "gcomp1", "gcomp2", "n")
      }else{
          colnames(a) <- c("tmle", "os", "tmle2", "os2", "plugin", "plugin2", "gcomp0", "gcomp1", "gcomp2", "n")
      }
      row.names(a) <- NULL
      return(a)
    }, simplify = FALSE)
  }else{
    all_out <- sapply(all_eff, function(x){
      do.call(summary_fn, args = list(out = out, which_eff = x))
    })
    return(all_out)
  }
}

load("/Users/scarlett/Desktop/nonparametric/submission2/all_out_intermed_new1.RData")

all_bias <- format_out(out, "get_bias")
all_rootnbias <- format_out(out, "get_rootnbias")
all_mse <- format_out(out, "get_mse")
all_nmse <- format_out(out, "get_nmse")
# !!! could get true asymptotic variance by simulating large data set
# !!! using true versions of nuisance estimators and evaluating variance of
# !!! EIF. 
all_sd <- format_out(out, "get_sd")
all_rootnsd <- format_out(out, "get_rootnsd")
all_cover <- format_out(out, "get_cover")
all_cover_orc <- format_out(out, "check_oracle_ci")
all_se <- format_out(out, "get_se_est")

all_density_orc <- format_out(out, "get_density_orc")
all_density_est <- format_out(out, "get_density_est")


blank_plot <- function(...){
  plot(1e-10, 1e-10, pch = "", bty = "n", xaxt = "n", ...)
}

plot_one_est_row <- function(which_eff = "total",
                             all_bias, all_mse, all_sd, 
                             all_cover, all_cover_orc, all_se,
                             bias_ylim = NULL){
  # six panel plot, two rows
  # top row = bias, sd, MSE
  # bottom row = samp dist. TMLE, samp. dist AIPTW, coverage ()
  # bias plot
  bias_range <- range(c(all_bias[[which_eff]]$tmle, all_bias[[which_eff]]$os, all_bias[[which_eff]]$gcomp0,
                        all_bias[[which_eff]]$gcomp1, all_bias[[which_eff]]$gcomp2))
  max_bias_val <- max(abs(bias_range))
  
  if(is.null(bias_ylim)){
    bias_ylim <- c(-1.05 * max_bias_val, 1.05 * max_bias_val)
  }
  blank_plot(xlim = c(0.5, 2.5), ylim = bias_ylim,
             xlab = "n", ylab = "Bias")
  axis(side = 1, at = 1:2, labels = all_bias[[which_eff]]$n)
  abline(h = 0, lty = 3)
  points(y = all_bias[[which_eff]]$tmle, x = 1:2, type = "b")
  points(y = all_bias[[which_eff]]$os, x = 1:2, type = "b", pch = 2)
  points(y = all_bias[[which_eff]]$gcomp0, x = 1:2, type = "b", pch = 3, col = "red")
  points(y = all_bias[[which_eff]]$gcomp1, x = 1:2, type = "b", pch = 4, col = "red")
  points(y = all_bias[[which_eff]]$gcomp2, x = 1:2, type = "b", pch = 5, col = "red")
  
  # sd plot
  # sd_range <- range(c(all_sd[[which_eff]]$tmle, all_sd[[which_eff]]$os, all_sd[[which_eff]]$gcomp0,
  #                    all_sd[[which_eff]]$gcomp1, all_sd[[which_eff]]$gcomp2))
  # max_sd_val <- max(abs(sd_range))
  blank_plot(xlim = c(0.5, 2.5), ylim = c(0, 0.13),
             xlab = "n", ylab = "Standard deviation")
  axis(side = 1, at = 1:2, labels = all_sd[[which_eff]]$n)
  # abline(h = 0, lty = 3)
  points(y = all_sd[[which_eff]]$tmle, x = 1:2, type = "b")
  points(y = all_sd[[which_eff]]$os, x = 1:2, type = "b", pch = 2)
  points(y = all_sd[[which_eff]]$gcomp0, x = 1:2, type = "b", pch = 3, col = "red")
  points(y = all_sd[[which_eff]]$gcomp1, x = 1:2, type = "b", pch = 4, col = "red")
  points(y = all_sd[[which_eff]]$gcomp2, x = 1:2, type = "b", pch = 5, col = "red")
  
  # mean squared error plot
  mse_range <- range(c(all_mse[[which_eff]]$tmle, all_mse[[which_eff]]$os, all_mse[[which_eff]]$gcomp0,
                       all_mse[[which_eff]]$gcomp1, all_mse[[which_eff]]$gcomp2))
  max_mse_val <- max(abs(mse_range))
  # min_mse_val <- min(abs(mse_range))
  blank_plot(xlim = c(0.5, 2.5), ylim = c(0, max_mse_val * 1.1), 
             xlab = "n", ylab = "Mean squared error")
  axis(side = 1, at = 1:2, labels = all_mse[[which_eff]]$n)
  points(y = all_mse[[which_eff]]$tmle, x = 1:2, type = "b")
  points(y = all_mse[[which_eff]]$os, x = 1:2, type = "b", pch = 2)
  points(y = all_mse[[which_eff]]$gcomp0, x = 1:2, type = "b", pch = 3, col = "red")
  points(y = all_mse[[which_eff]]$gcomp1, x = 1:2, type = "b", pch = 4, col = "red")
  points(y = all_mse[[which_eff]]$gcomp2, x = 1:2, type = "b", pch = 5, col = "red")
}

# plot results for point estimates
pdf("/Users/scarlett/Desktop/nonparametric/submission2/estimation_new1.pdf", height = 7, width = 7)
layout(matrix(c(1, 1, 1, 2:13), nrow =5, byrow = TRUE),
       heights = c(0.25, 1, 1, 1, 1))
par(oma = c(0, 2.1, 0, 0), mar = c(2.9, 2.9, 0.1, 0.1), mgp = c(1.8, 0.5, 0))
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend(x = 0.17, y = 0.5, ncol = 5, pch = c(2, 1, 3, 4, 5), title = "Estimator",
       legend = c("One-step", "TMLE", "Gcomp0", "Gcomp1", "Gcomp2"), xpd = NA,
       col = c("black", "black", "red", "red", "red"))
at_loc <- seq(0, 1, length = 6)[-c(1,6)] + c(0.052, 0.0125, -0.022, -0.055)[4:1]
plot_one_est_row(which_eff = "direct", all_bias = all_bias, all_mse = all_mse, all_sd = all_sd,
                 bias_ylim = c(-0.15, 0.1))
mtext(side = 2, outer = TRUE, line = 0, expression(psi[A]), at = at_loc[4])
plot_one_est_row(which_eff = "indirectM1", all_bias = all_bias, all_mse = all_mse, all_sd = all_sd,
                 bias_ylim = c(-0.15, 0.1))
mtext(side = 2, outer = TRUE, line = 0, expression(psi[M[1]]), at = at_loc[3])
plot_one_est_row(which_eff = "indirectM2", all_bias = all_bias, all_mse = all_mse, all_sd = all_sd,
                 bias_ylim = c(-0.15, 0.1))
mtext(side = 2, outer = TRUE, line = 0, expression(psi[M[2]]), at = at_loc[2])
plot_one_est_row(which_eff = "covarM1M2", all_bias = all_bias, all_mse = all_mse, all_sd = all_sd,
                 bias_ylim = c(-0.15, 0.1))
mtext(side = 2, outer = TRUE, line = 0, expression(psi[M[1]*","*M[2]]), at = at_loc[1])
dev.off()

## generate table for coverage -----------------------------

# coverage result
all_cover <- format_out(out, "get_cover")
all_cover_orc <- format_out(out, "check_oracle_ci")
# extract estimated coverage (n = 2000, n = 4000)
est_cover = Reduce(rbind, all_cover)[-c(1, 2),]
est_cover$type = factor(rep(c("direct", "indirectM1", "indirectM2", "covarM1M2"), each = 2), 
                        levels = c("direct", "indirectM1", "indirectM2", "covarM1M2"))
est_cover_2000 = est_cover[est_cover$n == 2000, ]
est_cover_4000 = est_cover[est_cover$n == 4000, ]
# extract oracle coverage (n = 2000, n = 4000)
oracle_cover = Reduce(rbind, all_cover_orc)[-c(1, 2),]
oracle_cover$type = factor(rep(c("direct", "indirectM1", "indirectM2", "covarM1M2"), each = 2), 
                           levels = c("direct", "indirectM1", "indirectM2", "covarM1M2"))
oracle_cover_2000 = oracle_cover[oracle_cover$n == 2000, ]
oracle_cover_4000 = oracle_cover[oracle_cover$n == 4000, ]
# row bind n = 2000 coverages
cover_2000 = rbind(est_cover_2000, oracle_cover_2000)
cover_2000 = cover_2000[order(cover_2000$type),]
# row bind n = 4000 coverages
cover_4000 = rbind(est_cover_4000, oracle_cover_4000)
cover_4000 = cover_4000[order(cover_4000$type),]
# column bind n = 2000 and n = 4000 coverages
mat = data.frame(type = rep(c("Estimated", "Oracle"),4))
mat = cbind(mat, cover_2000[, -c(6, 7)], cover_4000[, -c(6, 7)])
# change os and tml order
mat = mat[c(1,3,2,4:6,8,7,9:11)]

# print the table
# column for multirow describing scenario
n_sample_sizes = 2
add_multirow_scen <- function(scenario){
  c(paste0("\\multirow{", n_sample_sizes, "}{3em}{", scenario, "}"), rep("", n_sample_sizes - 1))
}

versions = c("$\\psi_A$", "$\\psi_{M_1}$", "$\\psi_{M_2}$", "$\\psi_{M_1, M_2}$")
version_rows = c(sapply(versions, add_multirow_scen))

# row names
addtorow = list()
addtorow$pos = list(0, 0)
addtorow$command = c("& & \\multicolumn{5}{c}{n = 2000} & \\multicolumn{5}{c}{n = 4000} \\\\ \n",
                     "& & One-step & TMLE & Gcomp0 & Gcomp1 & Gcomp2 & One-step & TMLE & Gcomp0 & Gcomp1 & Gcomp2 \\\\ \n")

df = data.frame(Estimator = version_rows, mat)

library(xtable)
print(xtable(df, digits=c(0,0,0,3,3,3,3,3,3,3,3,3,3)), include.rownames = F, 
      hline.after = c(-1, 0, seq(n_sample_sizes, nrow(df), by = n_sample_sizes)),
      include.colnames = F,
      sanitize.colnames.function = identity,
      sanitize.text.function = identity,
      add.to.row = addtorow)


