###########################################################################
# Predictions and visualization
###########################################################################

# load(path/to/model/output/model_output.RData)
# Here our model output name is cc_052021

roedeer <- png::readPNG("phylopic_png/roedeer.png") # need to download a png file first

str(cc_052021$sims.list)

## Colonization probability ~ season ----
#########################################

# length.out X seasons X mcmc list
col <- matrix(NA, 2, cc_052021$mcmc.info$n.samples)

col[1,] <- plogis(cc_052021$sims.list$alphagamma)       # summer = 0 so betagamma will be gone

col[2,] <- plogis(cc_052021$sims.list$alphagamma 
                  + cc_052021$sims.list$betagamma*1)    # winter = 1
str(col)

# then we take the mean of the mcmc list
pm.col <- apply(col, 1, mean)
str(pm.col)

# then calculate the credible intervals
CRI.col <- apply(col, 1, function(x) quantile(x, c(0.025, 0.975)))
str(CRI.col)

col.prob <- data.frame(season=c("Winter to Summer", "Summer to Winter"))

col.prob$pred <- c(pm.col[2], pm.col[1])
col.prob$pred <- round(col.prob$pred,3)

col.prob$lower <- c(CRI.col[1,2], CRI.col[1,1])

col.prob$upper <- c(CRI.col[2,2], CRI.col[2,1])

col.prob$species <- "Roe Deer"

head(col.prob)

ggplot(col.prob)+
  geom_errorbar(aes(season,ymin = lower, ymax = upper, col=season), size = 1.5, width = 0.2) +
  geom_point(aes(season, pred, col=season), size=10) +
  #scale_color_hue(direction = -1) +
  # geom_text(data=col.prob[1:2,],aes(1:2,pred,label=pred),hjust=0.5, vjust=-3, size=6, fontface="bold")+
  add_phylopic(roedeer, x=1.5, y=0.50, ysize=1.5,alpha = 0.3, color = "black") +
  ylim(c(0,1)) + 
  theme_minimal() + 
  theme(panel.grid.minor = element_blank())+
  ylab(expression(paste("Colonization Probabilitiy", " ", (gamma)))) +
  xlab('Season Transition') +
  theme(axis.text = element_text(size=18),
        axis.title = element_text(size=20),
        legend.position = "none")


## Desertion probability ~ 1 ----
#################################
# since there are no effects on this parameter, we plot with seasons to show there is no seasonality

# length.out X seasons X mcmc list
ext <- matrix(NA, 2, cc_052021$mcmc.info$n.samples)

ext[1,] <- plogis(cc_052021$sims.list$alphaeps)   # summer = 0

ext[2,] <- plogis(cc_052021$sims.list$alphaeps
                 + cc_052021$sims.list$betaeps*1))   # winter = 1

str(ext)

# then we take the mean of the mcmc list
pm.ext <- apply(ext, 1, mean)
str(pm.ext)

# then calculate the credible intervals
CRI.ext <- apply(ext, 1, function(x) quantile(x, c(0.025, 0.975)))
str(CRI.ext)

ext.prob <- data.frame(season=c("Winter to Summer", "Summer to Winter"))

ext.prob$pred <- c(pm.ext[2], pm.ext[1])
ext.prob$pred <- round(ext.prob$pred,3)

ext.prob$lower <- c(CRI.ext[1,2], CRI.ext[1,1])

ext.prob$upper <- c(CRI.ext[2,2], CRI.ext[2,1])

ext.prob$species <- "Roe Deer"

head(ext.prob)

ggplot(ext.prob)+
  geom_errorbar(aes(season,ymin = lower, ymax = upper, col=season), size = 1.5, width = 0.2) +
  geom_point(aes(season, pred, col=season), size=8) +
  #scale_color_hue(direction = -1) +
  # geom_text(data=col.prob[1:2,],aes(1:2,pred,label=pred),hjust=0.5, vjust=-3, size=6, fontface="bold")+
  add_phylopic(roedeer, x=1.5, y=0.50, ysize=1.5,alpha = 0.3, color = "black") +
  ylim(c(0,1)) + 
  theme_minimal() + 
  theme(panel.grid.minor = element_blank())+
  ylab(expression(paste("Desertion Probabilitiy", " ", (epsilon)))) +
  xlab('Season Transition') +
  theme(axis.text = element_text(size=18),
        axis.title = element_text(size=20),
        legend.position = "none")

## Detection probability ~ season + distpop ----
################################################
# (variable_data$DistPop) is in meters, we standardized it by dividing it with 10000

distpop_d <- seq(min(variable_data$DistPop)/10000,max(variable_data$DistPop)/10000, length.out = 1000)

# season X length.out.distpop X mcmc list
det <- array(NA, dim = c(1000,2,cc_052021$mcmc.info$n.samples))

for (i in 1:1000){
  det[i,1,] <- plogis(cc_052021$sims.list$alphap                      # summer = 0 so betap1 will be gone
                      + cc_052021$sims.list$betap[,2]*distpop_d[i])
  
  det[i,2,] <- plogis(cc_052021$sims.list$alphap 
                      + cc_052021$sims.list$betap[,1]*1               # winter = 1
                      + cc_052021$sims.list$betap[,2]*distpop_d[i])
}
str(det)

# then we take the mean of the mcmc list
pm.det <- apply(det, c(1,2), mean)
str(pm.det)

# then calculate the credible intervals
CRI.det <- apply(det, c(1,2), function(x) quantile(x, c(0.025, 0.975)))
str(CRI.det)

par(mfrow=c(1,1))
plot(distpop_d, pm.det[,1], ylim=c(0,1), col="red")
segments(distpop_d, CRI.det[1,,1], distpop_d, CRI.det[2,,1], lwd=0.1, col='red')
points(distpop_d, pm.det[,2], ylim=c(0,1), col="blue")
segments(distpop_d, CRI.det[1,,2], distpop_d, CRI.det[2,,2], lwd=0.1, col='blue')

# ggplot

# converting distpop into km here
det.prob <- data.frame(distpop=rep(distpop_d*10,2), season=rep(c("Winter to Summer", "Summer to Winter"), each = 1000, length.out = 2000))

det.prob$pred[1:1000] <- pm.det[,2]
det.prob$pred[1001:2000] <- pm.det[,1]

det.prob$lower[1:1000] <- CRI.det[1,,2]
det.prob$lower[1001:2000] <- CRI.det[1,,1]

det.prob$upper[1:1000] <- CRI.det[2,,2]
det.prob$upper[1001:2000] <- CRI.det[2,,1]

head(det.prob)

ggplot(det.prob)+
  geom_ribbon(aes(x= distpop, ymin = lower, ymax = upper, fill=season),alpha=0.5) +
  #scale_fill_brewer(type = "qual",palette=6,direction = -1)+
  geom_line(aes(distpop, pred, colour=season), lwd=1, linetype=1)+
  #scale_colour_brewer(type = "qual",palette=6,direction = -1)+
  add_phylopic(roedeer, x=35, y=0.50, ysize=40,alpha = 0.3, color = "black") +
  ylim(c(0,1)) + 
  theme_minimal() + 
  theme(panel.grid.minor = element_blank())+
  ylab(expression(paste("Detection Probabilitiy", " ", (p)))) +
  xlab('Distance to populated places (km)') +
  labs(fill = "Season Transition") + 
  labs(col = "Season Transition") + 
  theme(axis.text = element_text(size=18),
        axis.title = element_text(size=20),
        legend.text = element_text(size=18),
        legend.title = element_text(size=20),
        legend.position = c(0.15, 0.85),
        legend.direction = "vertical")

## Proportion of used sites ----
#########################################

plot(cc_052021$mean$n.prop)

prop_occu <- data.frame(season=1:22)
prop_occu$pred <- cc_052021$mean$n.prop
prop_occu$lower <- cc_052021$q2.5$n.prop
prop_occu$upper <- cc_052021$q97.5$n.prop
prop_occu$species <- "Roe Deer"

head(prop_occu)
# write.csv(fin_occu, "cc_052021_finoccu.csv")

prop_occu <- read.csv("/Users/dilsaddagtekin/Dropbox/PhD/North_Anatolia/NA_singlesp_2021/cloud/ccapreolus/cc_052021_propoccu.csv")[,-1]

ggplot(prop_occu) + 
  geom_ribbon(aes(x=1:22, ymin=lower, ymax=upper), fill="#cc4c02") +
  geom_point(aes(1:22, pred)) + 
  geom_line(aes(1:22, pred))+ 
  #add_phylopic(wildboar,x=11.4, y=120, ysize=60,alpha = 0.3, color = "black") +
  scale_x_continuous(breaks = 1:22) + theme_minimal() + 
  ylim(0,1) +
  theme(panel.grid.minor = element_blank())+
  ylab("Proportion of used sites") +
  xlab('Seasons') +
  theme(axis.text = element_text(size=18),
        axis.title = element_text(size=20),
        legend.position = "none")
