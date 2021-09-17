###########################################################################
# Model diagnostics
###########################################################################

# load(path/to/model/output/model_output.RData)
# Here our model output name is cc_052021


## Convergence and distribution check ----
###############################################

# checking Rhat values, if there ara any over 1.1
hist(cc_052021$summary[,8])
length(which(cc_052021$summary[,8]>1.1))

# visual traceplot checking
par(mfrow=c(3,3))
traceplot(cc_052021, parameters = c("alphapsi" , "alphagamma" , "betagamma" , "alphaeps" , "alphap" , "betap"))

par(mfrow = c(1,1))
whiskerplot(cc_052021, parameters = c("alphapsi" , "alphagamma" , "betagamma" , "alphaeps" , "alphap" , "betap"))
whiskerplot(cc_052021, parameters = c("n.occ"))
# whiskerplot(cc_052021, parameters = c("p"))
# whiskerplot(cc_052021, parameters = c("z"))

# visual parameter distribution checking

par(mfrow = c(3,3))

hist(cc_052021$sims.list$alphapsi, col="gray", main="")
abline(v=cc_052021$mean$alphapsi, col="red", lwd=2)

hist(cc_052021$sims.list$alphagamma, col="gray", main="")
abline(v=cc_052021$mean$alphagamma, col="red", lwd=2)

hist(cc_052021$sims.list$betagamma, col="gray", main="")
abline(v=cc_052021$mean$betagamma, col="red", lwd=2)

hist(cc_052021$sims.list$alphaeps, col="gray", main="")
abline(v=cc_052021$mean$alphaeps, col="red", lwd=2)

hist(cc_052021$sims.list$alphap, col="gray", main="")
abline(v=cc_052021$mean$alphap, col="red", lwd=2)

hist(cc_052021$sims.list$betap[,1], col="gray", main="")
abline(v=cc_052021$mean$betap[1], col="red", lwd=2)

hist(cc_052021$sims.list$betap[,2], col="gray", main="")
abline(v=cc_052021$mean$betap[2], col="red", lwd=2)

##  calculating WAIC ----
#########################

niter <- cc_052021$mcmc.info$n.samples # 67500
nobs <- nrow(dat) # 2826

# Probability of obtaining y detections after J surveys when the probability of success is p*z
ppd.save <- matrix(0, nrow = niter, ncol = nobs)
lppd.save <- matrix(0, nrow = niter, ncol = nobs)
for (i in 1:niter){
 for (t in 1:nobs){
   ppd.save[i,t] <- dbinom(dat$y[t],1, cc_052021$sims.list$p[i,t]*cc_052021$sims.list$z[i, dat$site[t], dat$pocc[t]])
   lppd.save[i,t] <- dbinom(dat$y[t],1, cc_052021$sims.list$p[i,t]*cc_052021$sims.list$z[i, dat$site[t], dat$pocc[t]], log = TRUE)
 }
}

# WAIC (Gelman recommends to use WAIC.2)
# computed log pointwise preditive density (Gelman et al 2014)
ppd_log <- log(apply(ppd.save, 2, mean)) 
ppd_ll <- -2 * sum(ppd_log)
# adjustments for WAIC
pD1 <- 2 * sum(ppd_log - apply(lppd.save, 2, mean))
pD2 <- sum(apply(lppd.save, 2, var))
WAIC1 <- ppd_ll + 2 * pD1
WAIC2 <- ppd_ll + 2 * pD2
WAIC2 # 2986.098
