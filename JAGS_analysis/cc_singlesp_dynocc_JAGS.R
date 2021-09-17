############################################################################
# This script is for the re-analysis of large-mammals in North Anatolia.
# Capreolus capreolus: single-species dynamic occupancy models.
# Data is from 11/2007 to 11/2018. Data prep file: data_prep_052021.R

# We decided to do seasons as:
# Active/Summer season: May-June-July-Aug-Sept-Oct
# Inactive/Winter season: Dec-Jan-Feb-March
# we dropped April and November data to allow for transition periods
# so data needs to be equal from 2007 December to 2018 October: 24096 - 24226
# This makes 11 Winter and 11 Summer. From 2008 Winter until 2018 Summer.

# FINAL MODEL AFTER NESTED MODEL SELECTION BASED ON WAIC AND WAICw
# from previous analysis.

# Date: 03/05/2021
# Author: Dilsad Dagtekin
###########################################################################

###########################################################################
# 1. House keeping and loading libraries and data
###########################################################################

## 1.1. House keeping ----
##########################

rm(list = ls())

## 1.2. Loading libraries ----
##############################

load.libraries <- function(){
  library(dplyr)
  library(jagsUI)
  library(ggplot2)
  library(ggthemes)
  library(rphylopic)
}

load.libraries()

## 1.3. Loading data ----
#########################

setwd("/Users/dilsaddagtekin/Dropbox/PhD/North_Anatolia/NA_singlesp_2021")

variable_data <- read.csv("variable_data.csv")[,-1]
head(variable_data)

stations <- read.csv("stations.csv")[,-1]
# remove TUL08 and TUL12 stations
stations <- stations[-c(165, 169)]
length(stations) # 171

# sessions: 12/2007 - 11/2018
aprils <- ((2008:2018)*12) + 4
novembers <- ((2008:2018)*12) + 11

sessions <- 24096:24227
sessions <- sessions[-which(sessions %in% c(aprils,novembers))]
length(sessions) # 110

cc_dh <- read.csv("cc_dh_052021.csv")[,-1]
rownames(cc_dh) <- stations
colnames(cc_dh) <- sessions
head(cc_dh)

dat <- read.csv("cc_dat_052021.csv")[,-1]

## 1.4. Prepare data for JAGS ----
##################################

# we do not have any records on TUL08 and TUL12 stations, so we decided to discard them here 
variable_data <- variable_data[-c(165,169),]

# arrange everything that needs to be numeric
str(variable_data)

variable_data$Ruggedness <- as.numeric(variable_data$Ruggedness)
variable_data$Elevation <- as.numeric(variable_data$Elevation)

# We standardized the *elevation*, *populatiton density*, *distance to populated places* and *distance to nearest road* by making them as close to 0 as possible. 

habitat       <- variable_data$Habitat
aspect        <- variable_data$Aspect_Cat2
elevation     <- variable_data$Elevation/1000
slope         <- variable_data$Slope
rugg          <- variable_data$Ruggedness
popden        <- variable_data$PopDen/10
roadden       <- variable_data$RoadDen
distpop       <- variable_data$DistPop/10000
distroad      <- variable_data$DistRoad/10000

# prepare variables for season
# we start with Dec 2007 so winter and end with Oct 2018 so summer
# 11 winters and 11 summers: 22 seasons, 11 years

season <- matrix(c("W", "S"), 171,22, byrow=T)     # all seasons


# season as numeric -- 2 level dummy variable: 1/0
season_n <- matrix(NA, nrow(season), ncol(season))
for (i in 1:nrow(season)){
  for (j in 1:ncol(season)) {
    season_n[i,j] <- ifelse(season[i,j]=="W", 1,0)
  }
}
str(season_n)

# habitat -- more than 2 levels
habitat_n <- as.integer(as.factor(habitat)) ; table(habitat_n)
habitat_1 <- NA; habitat_2 <- NA; habitat_3 <- NA; habitat_4 <- NA

for (i in 1:171){
  habitat_1[i] <- ifelse(habitat_n[i]==1,1,0)
  habitat_2[i] <- ifelse(habitat_n[i]==2,1,0)
  habitat_3[i] <- ifelse(habitat_n[i]==3,1,0)
  habitat_4[i] <- ifelse(habitat_n[i]==4,1,0)
}
habitat_n
habitat_4

# aspect -- 2 level dummy variable: 1/0
aspect_n <- NA
for (i in 1:length(aspect)){
  aspect_n[i] <- ifelse(aspect[i]=="N",1,0)
}
aspect_n

###########################################################################
# 2. JAGS analysis
###########################################################################

# dh: detection history
# nsite = nrow(cc_dh)             # 171
# nprimary= (ncol(cc_dh)/10)*2    # 22
# nsecondary = ncol(cc_dh)        # 110

# Bundle data
str(jags.data <- list(y              = dat$y, 
                      nsite          = nrow(cc_dh), 
                      nprimary       = (ncol(cc_dh)/10)*2,
                      nsecondary     = length(unique(dat$socc)),
                      nobs           = nrow(dat),
                      pocc           = dat$pocc,
                      site           = dat$site,
                      season.p       = dat$season,
                      habitat.p      = dat$habitat, 
                      elevation.p    = dat$elevation, 
                      slope.p        = dat$slope, 
                      aspect.p       = dat$aspect, 
                      rugg.p         = dat$rugg,
                      popden.p       = dat$popden, 
                      roadden.p      = dat$roadden,
                      distpop.p      = dat$distpop, 
                      distroad.p     = dat$distroad,
                      season         = season_n,
                      habitat        = habitat_n,
                      elevation      = elevation,
                      slope          = slope,
                      aspect         = aspect_n,
                      rugg           = rugg,
                      popden         = popden,
                      roadden        = roadden,
                      distpop        = distpop,
                      distroad       = distroad))

# Specify model in BUGS language
cat(file = "cc_052021.txt"," 
model {

# Priors

# psi ~ 1
alphapsi <- logit(psi.mean)                                 # constant, no beta value
psi.mean ~ dunif(0, 1)

# gamma ~ season
alphagamma <- logit(gamma.mean)                             # intercept
gamma.mean ~ dunif(0, 1)
betagamma ~ dnorm(0,0.01)                                   # betagamma season effect
                                                        
# eps ~ 1
alphaeps <- logit(eps.mean)                                 # intercept
eps.mean ~ dunif(0, 1)

# p ~ season + distpop
alphap <- logit(p.mean)                                     # intercept
p.mean ~ dunif(0, 1)
for(b in 1:2){                                              # we have 2 betas:
  betap[b] ~  dnorm(0,0.01)                                 # betap[1] season effect
}                                                           # betap[2] distpop effect
  
# Ecological submodel: intial occupancy as derived parameter
for (i in 1:nsite){
  z[i,1] ~ dbern(psi1[i])
  logit(psi1[i]) <- alphapsi 
}

# State transitions: colonization and extinction
# for col and ext
for (i in 1:nsite){
  for (s in 2:nprimary){
    logit(gamma[i,s-1]) <- alphagamma 
    + betagamma*season[i,s-1]
    
    logit(eps[i,s-1]) <- alphaeps
  }
}

# for z: true occupancy
for(i in 1:nsite){
  for (t in 2:nprimary){
    
    z[i,t] ~ dbern((z[i,t-1]*(1-eps[i,t-1])) + ((1-z[i,t-1])*gamma[i,t-1]))
  }
}

# Derived parameters
# Compute population and sample occupancy

n.occ[1] <- sum(z[1:nsite,1])  # Number of occupied sites in sample
n.prop[1] <- n.occ[1]/nsite

for (t in 2:nprimary){
  n.occ[t] <- sum(z[1:nsite,t])
  n.prop[t] <- n.occ[t]/nsite
}

# Observation model 
for (i in 1:nobs){
  logit(p[i]) <- alphap 
  + betap[1]*season.p[site[i]] 
  + betap[2]*distpop.p[site[i]]

    y[i] ~ dbern(z[site[i],pocc[i]]*p[i])  # y: observed occupancy
}

}")

# Initial values
inits<-function(){list(z=matrix(1,171,22))}

# Parameters monitored
params <- c("alphapsi", "alphagamma", "betagamma", "alphaeps","alphap","betap","n.occ","n.prop","p","z")

# MCMC settings 
# na <- 50000
# nb <- 50000
# nt <- 20
# ni <- 500000
# nc <- 3

na <- 50
nb <- 50
nt <- 20
ni <- 500
nc <- 3

parallel:::setDefaultClusterOptions(setup_strategy = "sequential") # this is needed for parallel run

# call JAGS
cc_052021 <- jags(jags.data, inits, params, "cc_052021.txt", n.chains = nc, n.adapt = na, n.iter = ni, n.burnin = nb, n.thin = nt, parallel = TRUE)

print(cc_052021, 3)
str(cc_052021$sims.list)

