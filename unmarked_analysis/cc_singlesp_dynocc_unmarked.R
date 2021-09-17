# Preliminary analysis with unmarked to determine candidate models
# Multi-season dynamic occupancy model for Capreolus capreolus 2007-2018 data
library(plyr)
library(unmarked)

setwd("")

# load data --------------------------------------

# Capreolus capreolus detection history
cc_dh <- read.csv("cc_dh_unmarked.csv")[,-1]

# variable data
variable_data <- read.csv("variable_data.csv")
variable_data <- variable_data[,-1]
variable_data$Ruggedness <- as.numeric(variable_data$Ruggedness)
str(variable_data)

# get variables for the model

season <- cbind(matrix('W',173,1), matrix(c('S','W'), 173,22, byrow=T))

habitat <- variable_data$Habitat
elevation <- variable_data$Elevation
slope <- variable_data$Slope
aspect <- variable_data$Aspect_Cat2
rugg <- variable_data$Ruggedness

# models ------------------------------------------------------------------
ccapreolus <- unmarkedMultFrame(y = cc_dh,
                                siteCovs = data.frame(habitat=habitat, elevation=elevation, slope=slope, 
                                                      aspect=aspect, rugg = rugg),
                                yearlySiteCovs = list(season=season),numPrimary = 23)

str(ccapreolus)
summary(ccapreolus)

# null model
cc_dm0 <- colext(psiformula = ~ 1,
                 gammaformula = ~ 1,
                 epsilonformula = ~ 1,
                 pformula = ~ 1,
                 data = ccapreolus)


# detection probability models --------------------------------------------

cc_dm1 <- colext(psiformula = ~ habitat,
                 gammaformula = ~ season*habitat,
                 epsilonformula = ~ season*habitat,
                 pformula = ~ season,
                 data = ccapreolus)

cc_dm2 <- colext(psiformula = ~ habitat,
                 gammaformula = ~ season*habitat,
                 epsilonformula = ~ season*habitat,
                 pformula = ~ habitat,
                 data = ccapreolus)

cc_dm3 <- colext(psiformula = ~ habitat,
                 gammaformula = ~ season*habitat,
                 epsilonformula = ~ season*habitat,
                 pformula = ~ elevation,
                 data = ccapreolus)

cc_dm4 <- colext(psiformula = ~ habitat,
                 gammaformula = ~ season*habitat,
                 epsilonformula = ~ season*habitat,
                 pformula = ~ slope,
                 data = ccapreolus)

cc_dm5 <- colext(psiformula = ~ habitat,
                 gammaformula = ~ season*habitat,
                 epsilonformula = ~ season*habitat,
                 pformula = ~ aspect,
                 data = ccapreolus)

cc_dm6 <- colext(psiformula = ~ habitat,
                 gammaformula = ~ season*habitat,
                 epsilonformula = ~ season*habitat,
                 pformula = ~ rugg,
                 data = ccapreolus)

cc_dm7 <- colext(psiformula = ~ habitat,
                 gammaformula = ~ season*habitat,
                 epsilonformula = ~ season*habitat,
                 pformula = ~ season + habitat,
                 data = ccapreolus)

cc_dm8 <- colext(psiformula = ~ habitat,
                 gammaformula = ~ season*habitat,
                 epsilonformula = ~ season*habitat,
                 pformula = ~ season + elevation,
                 data = ccapreolus)

cc_dm9 <- colext(psiformula = ~ habitat,
                 gammaformula = ~ season*habitat,
                 epsilonformula = ~ season*habitat,
                 pformula = ~ season + slope,
                 data = ccapreolus)

cc_dm10 <- colext(psiformula = ~ habitat,
                  gammaformula = ~ season*habitat,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ season + aspect,
                  data = ccapreolus)

cc_dm11 <- colext(psiformula = ~ habitat,
                  gammaformula = ~ season*habitat,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ season + rugg,
                  data = ccapreolus)

cc_dm12 <- colext(psiformula = ~ habitat,
                  gammaformula = ~ season*habitat,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ habitat + elevation,
                  data = ccapreolus)

cc_dm13 <- colext(psiformula = ~ habitat,
                  gammaformula = ~ season*habitat,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ habitat + slope,
                  data = ccapreolus)

cc_dm14 <- colext(psiformula = ~ habitat,
                  gammaformula = ~ season*habitat,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ habitat + aspect,
                  data = ccapreolus)

cc_dm15 <- colext(psiformula = ~ habitat,
                  gammaformula = ~ season*habitat,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ habitat + rugg,
                  data = ccapreolus)

cc_dm16 <- colext(psiformula = ~ habitat,
                  gammaformula = ~ season*habitat,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ elevation + slope,
                  data = ccapreolus)

cc_dm17 <- colext(psiformula = ~ habitat,
                  gammaformula = ~ season*habitat,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ elevation + aspect,
                  data = ccapreolus)

cc_dm18 <- colext(psiformula = ~ habitat,
                  gammaformula = ~ season*habitat,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ elevation + rugg,
                  data = ccapreolus)

cc_dm19 <- colext(psiformula = ~ habitat,
                  gammaformula = ~ season*habitat,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ slope + aspect,
                  data = ccapreolus)

cc_dm20 <- colext(psiformula = ~ habitat,
                  gammaformula = ~ season*habitat,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ slope + rugg,
                  data = ccapreolus)

cc_dm21 <- colext(psiformula = ~ habitat,
                  gammaformula = ~ season*habitat,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ aspect + rugg,
                  data = ccapreolus)

cc_dm22 <- colext(psiformula = ~ habitat,
                  gammaformula = ~ season*habitat,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ season + habitat + season:habitat,
                  data = ccapreolus)

cc_dm23 <- colext(psiformula = ~ habitat,
                  gammaformula = ~ season*habitat,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ season + elevation + season:elevation,
                  data = ccapreolus)

cc_dm24 <- colext(psiformula = ~ habitat,
                  gammaformula = ~ season*habitat,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ season + slope + season:slope,
                  data = ccapreolus)

cc_dm25 <- colext(psiformula = ~ habitat,
                  gammaformula = ~ season*habitat,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ season + aspect + season:aspect,
                  data = ccapreolus)

cc_dm26 <- colext(psiformula = ~ habitat,
                  gammaformula = ~ season*habitat,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ season + rugg + season:rugg,
                  data = ccapreolus)

cc_dm27 <- colext(psiformula = ~ habitat,
                  gammaformula = ~ season*habitat,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ habitat + elevation + habitat:elevation,
                  data = ccapreolus)

cc_dm28 <- colext(psiformula = ~ habitat,
                  gammaformula = ~ season*habitat,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ habitat + slope + habitat:slope,
                  data = ccapreolus)

cc_dm29 <- colext(psiformula = ~ habitat,
                  gammaformula = ~ season*habitat,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ habitat + aspect + habitat:aspect,
                  data = ccapreolus)

cc_dm30 <- colext(psiformula = ~ habitat,
                  gammaformula = ~ season*habitat,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ habitat + rugg + habitat:rugg,
                  data = ccapreolus)

cc_dm31 <- colext(psiformula = ~ habitat,
                  gammaformula = ~ season*habitat,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ elevation + slope + elevation:slope,
                  data = ccapreolus)

cc_dm32 <- colext(psiformula = ~ habitat,
                  gammaformula = ~ season*habitat,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ elevation + aspect + elevation:aspect,
                  data = ccapreolus)

cc_dm33 <- colext(psiformula = ~ habitat,
                  gammaformula = ~ season*habitat,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ elevation + rugg + elevation:rugg,
                  data = ccapreolus)

cc_dm34 <- colext(psiformula = ~ habitat,
                  gammaformula = ~ season*habitat,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ slope + aspect + slope:aspect,
                  data = ccapreolus)

cc_dm35 <- colext(psiformula = ~ habitat,
                  gammaformula = ~ season*habitat,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ slope + rugg + slope:rugg,
                  data = ccapreolus)

cc_dm36 <- colext(psiformula = ~ habitat,
                  gammaformula = ~ season*habitat,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ aspect + rugg + aspect:rugg,
                  data = ccapreolus)

cc_dm37 <- colext(psiformula = ~ habitat,
                  gammaformula = ~ season*habitat,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ 1,
                  data = ccapreolus)

models <- fitList('cc_dm1'=cc_dm1,'cc_dm2'=cc_dm2,'cc_dm3'=cc_dm3,'cc_dm4'=cc_dm4,'cc_dm5'=cc_dm5,
                  'cc_dm6'=cc_dm6,'cc_dm7'=cc_dm7,'cc_dm8'=cc_dm8,'cc_dm9'=cc_dm9,'cc_dm10'=cc_dm10,
                  'cc_dm11'=cc_dm11,'cc_dm12'=cc_dm12,'cc_dm13'=cc_dm13,'cc_dm14'=cc_dm14,'cc_dm15'=cc_dm15,
                  'cc_dm16'=cc_dm16,'cc_dm17'=cc_dm17,'cc_dm18'=cc_dm18,'cc_dm19'=cc_dm19,'cc_dm20'=cc_dm20,
                  'cc_dm21'=cc_dm21,'cc_dm22'=cc_dm22,'cc_dm23'=cc_dm23,'cc_dm24'=cc_dm24,'cc_dm25'=cc_dm25,
                  'cc_dm26'=cc_dm26,'cc_dm27'=cc_dm27,'cc_dm28'=cc_dm28,'cc_dm29'=cc_dm29,'cc_dm30'=cc_dm30,
                  'cc_dm31'=cc_dm31,'cc_dm32'=cc_dm32,'cc_dm33'=cc_dm33,'cc_dm34'=cc_dm34,'cc_dm35'=cc_dm35,
                  'cc_dm36'=cc_dm36,'cc_dm37'=cc_dm37)

modSel(models)
# cc_dm8

# initial occupancy probabilitty models -----------------------------------

cc_dm38 <- colext(psiformula = ~ habitat,
                  gammaformula = ~ season*habitat,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ season + elevation,
                  data = ccapreolus)

cc_dm39 <- colext(psiformula = ~ elevation,
                  gammaformula = ~ season*habitat,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ season + elevation,
                  data = ccapreolus)

cc_dm40 <- colext(psiformula = ~ slope,
                  gammaformula = ~ season*habitat,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ season + elevation,
                  data = ccapreolus)

cc_dm41 <- colext(psiformula = ~ aspect,
                  gammaformula = ~ season*habitat,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ season + elevation,
                  data = ccapreolus)

cc_dm42 <- colext(psiformula = ~ rugg,
                  gammaformula = ~ season*habitat,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ season + elevation,
                  data = ccapreolus)

cc_dm43 <- colext(psiformula = ~ habitat + elevation,
                  gammaformula = ~ season*habitat,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ season + elevation,
                  data = ccapreolus)

cc_dm44 <- colext(psiformula = ~ habitat + slope,
                  gammaformula = ~ season*habitat,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ season + elevation,
                  data = ccapreolus)

cc_dm45 <- colext(psiformula = ~ habitat + aspect,
                  gammaformula = ~ season*habitat,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ season + elevation,
                  data = ccapreolus)

cc_dm46 <- colext(psiformula = ~ habitat + rugg,
                  gammaformula = ~ season*habitat,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ season + elevation,
                  data = ccapreolus)

cc_dm47 <- colext(psiformula = ~ elevation + slope,
                  gammaformula = ~ season*habitat,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ season + elevation,
                  data = ccapreolus)

cc_dm48 <- colext(psiformula = ~ elevation + aspect,
                  gammaformula = ~ season*habitat,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ season + elevation,
                  data = ccapreolus)

cc_dm49 <- colext(psiformula = ~ elevation + rugg,
                  gammaformula = ~ season*habitat,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ season + elevation,
                  data = ccapreolus)

cc_dm50 <- colext(psiformula = ~ slope + aspect,
                  gammaformula = ~ season*habitat,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ season + elevation,
                  data = ccapreolus)

cc_dm51 <- colext(psiformula = ~ slope + rugg,
                  gammaformula = ~ season*habitat,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ season + elevation,
                  data = ccapreolus)

cc_dm52 <- colext(psiformula = ~ aspect + rugg,
                  gammaformula = ~ season*habitat,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ season + elevation,
                  data = ccapreolus)

cc_dm53 <- colext(psiformula = ~ habitat + elevation + habitat:elevation,
                  gammaformula = ~ season*habitat,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ season + elevation,
                  data = ccapreolus)

cc_dm54 <- colext(psiformula = ~ habitat + slope + habitat:slope,
                  gammaformula = ~ season*habitat,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ season + elevation,
                  data = ccapreolus)

cc_dm55 <- colext(psiformula = ~ habitat + aspect + habitat:aspect,
                  gammaformula = ~ season*habitat,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ season + elevation,
                  data = ccapreolus)

cc_dm56 <- colext(psiformula = ~ habitat + rugg + habitat:rugg,
                  gammaformula = ~ season*habitat,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ season + elevation,
                  data = ccapreolus)

cc_dm57 <- colext(psiformula = ~ elevation + slope + elevation:slope,
                  gammaformula = ~ season*habitat,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ season + elevation,
                  data = ccapreolus)

cc_dm58 <- colext(psiformula = ~ elevation + aspect + elevation:aspect,
                  gammaformula = ~ season*habitat,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ season + elevation,
                  data = ccapreolus)

cc_dm59 <- colext(psiformula = ~ elevation + rugg + elevation:rugg,
                  gammaformula = ~ season*habitat,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ season + elevation,
                  data = ccapreolus)

cc_dm60 <- colext(psiformula = ~ slope + aspect + slope:aspect,
                  gammaformula = ~ season*habitat,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ season + elevation,
                  data = ccapreolus)

cc_dm61 <- colext(psiformula = ~ slope + rugg + slope:rugg,
                  gammaformula = ~ season*habitat,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ season + elevation,
                  data = ccapreolus)

cc_dm62 <- colext(psiformula = ~ aspect + rugg + aspect:rugg,
                  gammaformula = ~ season*habitat,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ season + elevation,
                  data = ccapreolus)

cc_dm63 <- colext(psiformula = ~ 1,
                  gammaformula = ~ season*habitat,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ season + elevation,
                  data = ccapreolus)

models <- fitList('cc_dm38'=cc_dm38,'cc_dm39'=cc_dm39,'cc_dm40'=cc_dm40,'cc_dm41'=cc_dm41,'cc_dm42'=cc_dm42,
                  'cc_dm43'=cc_dm43,'cc_dm44'=cc_dm44,'cc_dm45'=cc_dm45,'cc_dm46'=cc_dm46,'cc_dm47'=cc_dm47,
                  'cc_dm48'=cc_dm48,'cc_dm49'=cc_dm49,'cc_dm50'=cc_dm50,'cc_dm51'=cc_dm51,'cc_dm52'=cc_dm52,
                  'cc_dm53'=cc_dm53,'cc_dm54'=cc_dm54,'cc_dm55'=cc_dm55,'cc_dm56'=cc_dm56,'cc_dm57'=cc_dm57,
                  'cc_dm58'=cc_dm58,'cc_dm59'=cc_dm59,'cc_dm60'=cc_dm60,'cc_dm61'=cc_dm61,'cc_dm62'=cc_dm62,
                  'cc_dm63'=cc_dm63)

modSel(models)
# cc_dm63

# colonization probability models -----------------------------------------

cc_dm64 <- colext(psiformula = ~ 1,
                  gammaformula = ~ season,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ season + elevation,
                  data = ccapreolus)

cc_dm65 <- colext(psiformula = ~ 1,
                  gammaformula = ~ habitat,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ season + elevation,
                  data = ccapreolus)

cc_dm66 <- colext(psiformula = ~ 1,
                  gammaformula = ~ elevation,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ season + elevation,
                  data = ccapreolus)

cc_dm67 <- colext(psiformula = ~ 1,
                  gammaformula = ~ slope,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ season + elevation,
                  data = ccapreolus)

cc_dm68 <- colext(psiformula = ~ 1,
                  gammaformula = ~ aspect,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ season + elevation,
                  data = ccapreolus)

cc_dm69 <- colext(psiformula = ~ 1,
                  gammaformula = ~ rugg,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ season + elevation,
                  data = ccapreolus)

cc_dm70 <- colext(psiformula = ~ 1,
                  gammaformula = ~ season + habitat,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ season + elevation,
                  data = ccapreolus)

cc_dm71 <- colext(psiformula = ~ 1,
                  gammaformula = ~ season + elevation,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ season + elevation,
                  data = ccapreolus)

cc_dm72 <- colext(psiformula = ~ 1,
                  gammaformula = ~ season + slope,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ season + elevation,
                  data = ccapreolus)

cc_dm73 <- colext(psiformula = ~ 1,
                  gammaformula = ~ season + aspect,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ season + elevation,
                  data = ccapreolus)

cc_dm74 <- colext(psiformula = ~ 1,
                  gammaformula = ~ season + rugg,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ season + elevation,
                  data = ccapreolus)

cc_dm75 <- colext(psiformula = ~ 1,
                  gammaformula = ~ habitat + elevation,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ season + elevation,
                  data = ccapreolus)

cc_dm76 <- colext(psiformula = ~ 1,
                  gammaformula = ~ habitat + slope,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ season + elevation,
                  data = ccapreolus)

cc_dm77 <- colext(psiformula = ~ 1,
                  gammaformula = ~ habitat + aspect,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ season + elevation,
                  data = ccapreolus)

cc_dm78 <- colext(psiformula = ~ 1,
                  gammaformula = ~ habitat + rugg,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ season + elevation,
                  data = ccapreolus)

cc_dm79 <- colext(psiformula = ~ 1,
                  gammaformula = ~ elevation + slope,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ season + elevation,
                  data = ccapreolus)

cc_dm80 <- colext(psiformula = ~ 1,
                  gammaformula = ~ elevation + aspect,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ season + elevation,
                  data = ccapreolus)

cc_dm81 <- colext(psiformula = ~ 1,
                  gammaformula = ~ elevation + rugg,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ season + elevation,
                  data = ccapreolus)

cc_dm82 <- colext(psiformula = ~ 1,
                  gammaformula = ~ slope + aspect,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ season + elevation,
                  data = ccapreolus)

cc_dm83 <- colext(psiformula = ~ 1,
                  gammaformula = ~ slope + rugg,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ season + elevation,
                  data = ccapreolus)

cc_dm84 <- colext(psiformula = ~ 1,
                  gammaformula = ~ aspect + rugg,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ season + elevation,
                  data = ccapreolus)

cc_dm85 <- colext(psiformula = ~ 1,
                  gammaformula = ~ season + habitat + season:habitat,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ season + elevation,
                  data = ccapreolus)

cc_dm86 <- colext(psiformula = ~ 1,
                  gammaformula = ~ season + elevation + season:elevation,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ season + elevation,
                  data = ccapreolus)

cc_dm87 <- colext(psiformula = ~ 1,
                  gammaformula = ~ season + slope + season:slope,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ season + elevation,
                  data = ccapreolus)

cc_dm88 <- colext(psiformula = ~ 1,
                  gammaformula = ~ season + aspect + season:aspect,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ season + elevation,
                  data = ccapreolus)

cc_dm89 <- colext(psiformula = ~ 1,
                  gammaformula = ~ season + rugg + season:rugg,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ season + elevation,
                  data = ccapreolus)

cc_dm90 <- colext(psiformula = ~ 1,
                  gammaformula = ~ habitat + elevation + habitat:elevation,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ season + elevation,
                  data = ccapreolus)

cc_dm91 <- colext(psiformula = ~ 1,
                  gammaformula = ~ habitat + slope + habitat:slope,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ season + elevation,
                  data = ccapreolus)

cc_dm92 <- colext(psiformula = ~ 1,
                  gammaformula = ~ habitat + aspect + habitat:aspect,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ season + elevation,
                  data = ccapreolus)

cc_dm93 <- colext(psiformula = ~ 1,
                  gammaformula = ~ habitat + rugg + habitat:rugg,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ season + elevation,
                  data = ccapreolus)

cc_dm94 <- colext(psiformula = ~ 1,
                  gammaformula = ~ elevation + slope + elevation:slope,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ season + elevation,
                  data = ccapreolus)

cc_dm95 <- colext(psiformula = ~ 1,
                  gammaformula = ~ elevation + aspect + elevation:aspect,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ season + elevation,
                  data = ccapreolus)

cc_dm96 <- colext(psiformula = ~ 1,
                  gammaformula = ~ elevation + rugg + elevation:rugg,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ season + elevation,
                  data = ccapreolus)

cc_dm97 <- colext(psiformula = ~ 1,
                  gammaformula = ~ slope + aspect + slope:aspect,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ season + elevation,
                  data = ccapreolus)

cc_dm98 <- colext(psiformula = ~ 1,
                  gammaformula = ~ slope + rugg + slope:rugg,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ season + elevation,
                  data = ccapreolus)

cc_dm99 <- colext(psiformula = ~ 1,
                  gammaformula = ~ aspect + rugg + aspect:rugg,
                  epsilonformula = ~ season*habitat,
                  pformula = ~ season + elevation,
                  data = ccapreolus)

cc_dm100 <- colext(psiformula = ~ 1,
                   gammaformula = ~ 1,
                   epsilonformula = ~ season*habitat,
                   pformula = ~ season + elevation,
                   data = ccapreolus)

models <- fitList('cc_dm64'=cc_dm64,'cc_dm65'=cc_dm65,'cc_dm66'=cc_dm66,'cc_dm67'=cc_dm67,'cc_dm68'=cc_dm68,
                  'cc_dm69'=cc_dm69,'cc_dm70'=cc_dm70,'cc_dm71'=cc_dm71,'cc_dm72'=cc_dm72,'cc_dm73'=cc_dm73,
                  'cc_dm74'=cc_dm74,'cc_dm75'=cc_dm75,'cc_dm76'=cc_dm76,'cc_dm77'=cc_dm77,'cc_dm78'=cc_dm78,
                  'cc_dm79'=cc_dm79,'cc_dm80'=cc_dm80,'cc_dm81'=cc_dm81,'cc_dm82'=cc_dm82,'cc_dm83'=cc_dm83,
                  'cc_dm84'=cc_dm84,'cc_dm85'=cc_dm85,'cc_dm86'=cc_dm86,'cc_dm87'=cc_dm87,'cc_dm88'=cc_dm88,
                  'cc_dm89'=cc_dm89, 'cc_dm90'=cc_dm90,'cc_dm91'=cc_dm91,'cc_dm92'=cc_dm92,'cc_dm93'=cc_dm93,
                  'cc_dm94'=cc_dm94,'cc_dm95'=cc_dm95,'cc_dm96'=cc_dm96,'cc_dm97'=cc_dm97,'cc_dm98'=cc_dm98,
                  'cc_dm99'=cc_dm99, 'cc_dm100'=cc_dm100)

modSel(models)
# cc_dm68

# extinction probability models -------------------------------------------

cc_dm101 <- colext(psiformula = ~ 1,
                   gammaformula = ~ aspect,
                   epsilonformula = ~ season,
                   pformula = ~ season + elevation,
                   data = ccapreolus)

cc_dm102 <- colext(psiformula = ~ 1,
                   gammaformula = ~ aspect,
                   epsilonformula = ~ habitat,
                   pformula = ~ season + elevation,
                   data = ccapreolus)

cc_dm103 <- colext(psiformula = ~ 1,
                   gammaformula = ~ aspect,
                   epsilonformula = ~ elevation,
                   pformula = ~ season + elevation,
                   data = ccapreolus)

cc_dm104 <- colext(psiformula = ~ 1,
                   gammaformula = ~ aspect,
                   epsilonformula = ~ slope,
                   pformula = ~ season + elevation,
                   data = ccapreolus)

cc_dm105 <- colext(psiformula = ~ 1,
                   gammaformula = ~ aspect,
                   epsilonformula = ~ aspect,
                   pformula = ~ season + elevation,
                   data = ccapreolus)

cc_dm106 <- colext(psiformula = ~ 1,
                   gammaformula = ~ aspect,
                   epsilonformula = ~ rugg,
                   pformula = ~ season + elevation,
                   data = ccapreolus)

cc_dm107 <- colext(psiformula = ~ 1,
                   gammaformula = ~ aspect,
                   epsilonformula = ~ season + habitat,
                   pformula = ~ season + elevation,
                   data = ccapreolus)

cc_dm108 <- colext(psiformula = ~ 1,
                   gammaformula = ~ aspect,
                   epsilonformula = ~ season + elevation,
                   pformula = ~ season + elevation,
                   data = ccapreolus)

cc_dm109 <- colext(psiformula = ~ 1,
                   gammaformula = ~ aspect,
                   epsilonformula = ~ season + slope,
                   pformula = ~ season + elevation,
                   data = ccapreolus)

cc_dm110 <- colext(psiformula = ~ 1,
                   gammaformula = ~ aspect,
                   epsilonformula = ~ season + aspect,
                   pformula = ~ season + elevation,
                   data = ccapreolus)

cc_dm111 <- colext(psiformula = ~ 1,
                   gammaformula = ~ aspect,
                   epsilonformula = ~ season + rugg,
                   pformula = ~ season + elevation,
                   data = ccapreolus)

cc_dm112 <- colext(psiformula = ~ 1,
                   gammaformula = ~ aspect,
                   epsilonformula = ~ habitat + elevation,
                   pformula = ~ season + elevation,
                   data = ccapreolus)

cc_dm113 <- colext(psiformula = ~ 1,
                   gammaformula = ~ aspect,
                   epsilonformula = ~ habitat + slope,
                   pformula = ~ season + elevation,
                   data = ccapreolus)

cc_dm114 <- colext(psiformula = ~ 1,
                   gammaformula = ~ aspect,
                   epsilonformula = ~ habitat + aspect,
                   pformula = ~ season + elevation,
                   data = ccapreolus)

cc_dm115 <- colext(psiformula = ~ 1,
                   gammaformula = ~ aspect,
                   epsilonformula = ~ habitat + rugg,
                   pformula = ~ season + elevation,
                   data = ccapreolus)

cc_dm116 <- colext(psiformula = ~ 1,
                   gammaformula = ~ aspect,
                   epsilonformula = ~ elevation + slope,
                   pformula = ~ season + elevation,
                   data = ccapreolus)

cc_dm117 <- colext(psiformula = ~ 1,
                   gammaformula = ~ aspect,
                   epsilonformula = ~ elevation + aspect,
                   pformula = ~ season + elevation,
                   data = ccapreolus)

cc_dm118 <- colext(psiformula = ~ 1,
                   gammaformula = ~ aspect,
                   epsilonformula = ~ elevation + rugg,
                   pformula = ~ season + elevation,
                   data = ccapreolus)

cc_dm119 <- colext(psiformula = ~ 1,
                   gammaformula = ~ aspect,
                   epsilonformula = ~ slope + aspect,
                   pformula = ~ season + elevation,
                   data = ccapreolus)

cc_dm120 <- colext(psiformula = ~ 1,
                   gammaformula = ~ aspect,
                   epsilonformula = ~ slope + rugg,
                   pformula = ~ season + elevation,
                   data = ccapreolus)

cc_dm121 <- colext(psiformula = ~ 1,
                   gammaformula = ~ aspect,
                   epsilonformula = ~ aspect + rugg,
                   pformula = ~ season + elevation,
                   data = ccapreolus)

cc_dm122 <- colext(psiformula = ~ 1,
                   gammaformula = ~ aspect,
                   epsilonformula = ~ season + habitat + season:habitat,
                   pformula = ~ season + elevation,
                   data = ccapreolus)

cc_dm123 <- colext(psiformula = ~ 1,
                   gammaformula = ~ aspect,
                   epsilonformula = ~ season + elevation + season:elevation,
                   pformula = ~ season + elevation,
                   data = ccapreolus)

cc_dm124 <- colext(psiformula = ~ 1,
                   gammaformula = ~ aspect,
                   epsilonformula = ~ season + slope + season:slope,
                   pformula = ~ season + elevation,
                   data = ccapreolus)

cc_dm125 <- colext(psiformula = ~ 1,
                   gammaformula = ~ aspect,
                   epsilonformula = ~ season + aspect + season:aspect,
                   pformula = ~ season + elevation,
                   data = ccapreolus)

cc_dm126 <- colext(psiformula = ~ 1,
                   gammaformula = ~ aspect,
                   epsilonformula = ~ season + rugg + season:rugg,
                   pformula = ~ season + elevation,
                   data = ccapreolus)

cc_dm127 <- colext(psiformula = ~ 1,
                   gammaformula = ~ aspect,
                   epsilonformula = ~ habitat + elevation + habitat:elevation,
                   pformula = ~ season + elevation,
                   data = ccapreolus)

cc_dm128 <- colext(psiformula = ~ 1,
                   gammaformula = ~ aspect,
                   epsilonformula = ~ habitat + slope + habitat:slope,
                   pformula = ~ season + elevation,
                   data = ccapreolus)

cc_dm129 <- colext(psiformula = ~ 1,
                   gammaformula = ~ aspect,
                   epsilonformula = ~ habitat + aspect + habitat:aspect,
                   pformula = ~ season + elevation,
                   data = ccapreolus)

cc_dm130 <- colext(psiformula = ~ 1,
                   gammaformula = ~ aspect,
                   epsilonformula = ~ habitat + rugg + habitat:rugg,
                   pformula = ~ season + elevation,
                   data = ccapreolus)

cc_dm131 <- colext(psiformula = ~ 1,
                   gammaformula = ~ aspect,
                   epsilonformula = ~ elevation + slope + elevation:slope,
                   pformula = ~ season + elevation,
                   data = ccapreolus)

cc_dm132 <- colext(psiformula = ~ 1,
                   gammaformula = ~ aspect,
                   epsilonformula = ~ elevation + aspect + elevation:aspect,
                   pformula = ~ season + elevation,
                   data = ccapreolus)

cc_dm133 <- colext(psiformula = ~ 1,
                   gammaformula = ~ aspect,
                   epsilonformula = ~ elevation + rugg + elevation:rugg,
                   pformula = ~ season + elevation,
                   data = ccapreolus)

cc_dm134 <- colext(psiformula = ~ 1,
                   gammaformula = ~ aspect,
                   epsilonformula = ~ slope + aspect + slope:aspect,
                   pformula = ~ season + elevation,
                   data = ccapreolus)

cc_dm135 <- colext(psiformula = ~ 1,
                   gammaformula = ~ aspect,
                   epsilonformula = ~ slope + rugg + slope:rugg,
                   pformula = ~ season + elevation,
                   data = ccapreolus)

cc_dm136 <- colext(psiformula = ~ 1,
                   gammaformula = ~ aspect,
                   epsilonformula = ~ aspect + rugg + aspect:rugg,
                   pformula = ~ season + elevation,
                   data = ccapreolus)

cc_dm137 <- colext(psiformula = ~ 1,
                   gammaformula = ~ aspect,
                   epsilonformula = ~ 1,
                   pformula = ~ season + elevation,
                   data = ccapreolus)

models <- fitList('cc_dm101'=cc_dm101,'cc_dm102'=cc_dm102,'cc_dm103'=cc_dm103,'cc_dm104'=cc_dm104,'cc_dm105'=cc_dm105,
                  'cc_dm106'=cc_dm106,'cc_dm107'=cc_dm107,'cc_dm108'=cc_dm108,'cc_dm109'=cc_dm109,'cc_dm110'=cc_dm110,
                  'cc_dm111'=cc_dm111,'cc_dm112'=cc_dm112,'cc_dm113'=cc_dm113,'cc_dm114'=cc_dm114,'cc_dm115'=cc_dm115,
                  'cc_dm116'=cc_dm116,'cc_dm117'=cc_dm117,'cc_dm118'=cc_dm118,'cc_dm119'=cc_dm119,'cc_dm120'=cc_dm120,
                  'cc_dm121'=cc_dm121,'cc_dm122'=cc_dm122,'cc_dm123'=cc_dm123,'cc_dm124'=cc_dm124,'cc_dm125'=cc_dm125,
                  'cc_dm126'=cc_dm126,'cc_dm127'=cc_dm127,'cc_dm128'=cc_dm128,'cc_dm129'=cc_dm129,'cc_dm130'=cc_dm130,
                  'cc_dm131'=cc_dm131,'cc_dm132'=cc_dm132,'cc_dm133'=cc_dm133,'cc_dm134'=cc_dm134,'cc_dm135'=cc_dm135,
                  'cc_dm136'=cc_dm136,'cc_dm137'=cc_dm137)

modSel(models)
# cc_dm111

# estimate predictions ----------------------------------------------------

# best model according to AIC was 
# colext(psiformula = ~ 1,
#        gammaformula = ~ aspect,
#        epsilonformula = ~ season + rugg,
#        pformula = ~ season + elevation,
#        data = ccapreolus)

# best parsimonous model is:
# colext(psiformula = ~ 1,
#        gammaformula = ~ aspect,
#        epsilonformula = ~ season + rugg,
#        pformula = ~ season + elevation,
#        data = ccapreolus)

cc_dm <- colext(psiformula = ~ 1,
                gammaformula = ~ aspect,
                epsilonformula = ~ season + rugg,
                pformula = ~ season + elevation,
                data = ccapreolus)

summary(cc_dm)

# predictions
library(ggplot2)

# predict detection probability: season + elevation
nd_det <- expand.grid(season=c('S', 'W'), elevation = seq(741, 2134))

dat_det = cbind(nd_det,predict(cc_dm, type='det', newdata=nd_det))

ggplot(dat_det,aes(elevation, Predicted, fill=season))+
  #geom_point(position = position_dodge(width = 0.5))+
  geom_line(position = position_dodge(width = 0.5))+
  geom_ribbon(aes(ymin = lower, ymax = upper),alpha=0.5,position = position_dodge(width = 0.5))+
  #geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, position = position_dodge(width = 0.5)) +
  ylim(c(0,1)) +
  ylab('C.capreolus - Detection')

# predict colonization probability: aspect
nd_col <- expand.grid(aspect = c('N','S'))

dat_col = cbind(nd_col,predict(cc_dm, type='col', newdata=nd_col))

ggplot(dat_col,aes(aspect, Predicted, col=aspect))+
  geom_point(position = position_dodge(width = 0.5))+
  #geom_ribbon(aes(ymin = lower, ymax = upper),alpha=0.5,position = position_dodge(width = 0.5), fill='blue') +
  geom_line(position = position_dodge(width = 0.5))+
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, position = position_dodge(width = 0.5)) +
  ylim(c(0,1)) +
  ylab('C.capreolus - Colonization')

# predict extinction probability: season + elevation
nd_ext <- expand.grid(season=c('S', 'W'), rugg = seq(2, 7))

dat_ext = cbind(nd_ext,predict(cc_dm, type='ext', newdata=nd_ext))

ggplot(dat_ext,aes(rugg, Predicted, col=season))+
  geom_point(position = position_dodge(width = 0.5))+
  geom_line(position = position_dodge(width = 0.5))+
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, position = position_dodge(width = 0.5)) +
  ylim(c(0,1)) +
  ylab('C.capreolus - Extinction')

