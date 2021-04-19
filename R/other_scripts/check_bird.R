# Compare radiative transfer simulations with bird model
library(dplyr)
source("R/BIO_bird.R")

########################################################
# Simulation from radiative transfer code
########################################################
tab <- read.table("C:/Users/devrede/Desktop/TMP/Students/Kitty_Kam/coart_simulation_Kitty.txt",skip = 12,sep="")

names(tab) <- c("wv","sza","z","dif_down","dir_down","Ed","Eu","Eu/Ed")
tab$wv = tab$wv*1000 # from micrometers to nanometers
# We need to convert irradiance from W m-2 microm-1 to W m-2 nm-1, so we divide by 1000
#tab = tab %>% mutate(dif_down = dif_down,
#                     dir_down = dir_down,
#                     Ed = Ed,
#                     Eu = Eu
#)

zerom = tab %>% filter (z == -1)
zerop = tab %>% filter (z == 0)


# test at 45 degrees
zend = 45. 
zen = zend * pi /180.
Edp45 = zerop %>% filter(sza == zend)

#######################################################
# Brid model 
#######################################################
source("R/BIO_bird.R")
source("R/BIO_irradiance_corrections.R")

lambda = seq(0.4,0.7,0.005)
irr.bird = bird(ZEN = zen, lambda)
irr.seacor = I0_corr_season(I0_direct = irr.bird$DIRECT,I0_diffuse = irr.bird$DIFFUSE,Day = 180)
irr.clcor = I0_corr_cloud_BIO(I0_direct = irr.seacor$I0_direct,I0_diffuse = irr.seacor$I0_diffuse,ZEN = zen,Cloud = 0)

#######################################################
# Gregg and Carder model (From Hydrolight)
#######################################################
source("R/Gcmod.r")




# correct bird for season
plot(Edp45$wv,Edp45$dir_down,type="l",ylim = c(0,1600))
lines(Edp45$wv,Edp45$dif_down,lty=2)
lines(lambda*1000,irr.bird$DIRECT,col=2)
lines(lambda*1000,irr.bird$DIFFUSE,col=2,lty=2)
lines(Edp45$wv,Edp45$Ed,lwd=2)
lines(lambda*1000,irr.bird$DIFFUSE+irr.bird$DIRECT,lwd=2,col=2)
lines(lambda*1000,irr.seacor$I0_direct+irr.seacor$I0_diffuse,lwd=2,col=3)
lines(lambda*1000,irr.clcor$I0_direct+irr.clcor$I0_diffuse,lwd=2,col=4)
