library(dplyr)
library(oce)
library(lubridate)
library(stringr)
library(insol)

source("R/BIO_zenith.R")
source("R/BIO_bird.R")
source("R/BIO_irradiance_corrections.R")
source("R/Gcmod.r")
# There are various ways to compute the sun zenith angle, let's check that first:
daypxl = 120
yearpxl = 2019
hrpxl = 14
latipxl = 58
datepxl = doyday(yearpxl,daypxl) + hours(hrpxl)
julianpxl = JD(doyday(yearpxl,daypxl))
datepxl = str_replace(as.character(datepxl),"AST","UTC") # Here I put the date to U
#Here we take the mean longitude for the area of interest
#it would have to be changed 
zendR = 90 - sunAngle(datepxl,longitude = 0, latitude = latipxl)$altitude

## BIO way:
res = daylat(daypxl,latipxl)
zenbio = zenith_bio(hrpxl,0,delta = res$delta, phi = res$phi)
zendbio = zenbio * 180./ pi

## Hydrolight, manu recoded:
zendhm = sunang(iday = daypxl, hr = hrpxl, rad = rad, xlon = 0, ylat = latipxl)$sunz

print(c(zendR,zendbio, zendhm))

# CONCLUSION: BOTH MODEL PROVIDE SAME ANGLES.

# simulation were made using an online radiative coupled ocean-atmosphere
# radiative transfer code from NASA:
# https://satcorps.larc.nasa.gov/jin/coart.html

########################################################
# Simulation from radiative transfer code
########################################################
tab <- read.table("C:/Users/devrede/Desktop/TMP/Students/Kitty_Kam/coart_simulation_check.txt",skip = 12,sep="")

names(tab) <- c("wv","sza","z","dif_down","dir_down","Ed","Eu","Eu/Ed")
tab$wv = tab$wv*1000 # from micrometers to nanometers
# We need to convert irradiance from W m-2 microm-1 to W m-2 nm-1, so we divide by 1000
#tab = tab %>% mutate(dif_down = dif_down,
#                     dir_down = dir_down,
#                     Ed = Ed,
#                     Eu = Eu
#)

zerom = tab %>% filter (z == -0.1)
zerop = tab %>% filter (z == 0)
zerot = tab %>% filter (z == 100.0)

# test at 45 degrees
zend = 45. 
zen = zend * pi /180.
Edp45 = zerop %>% filter(sza == zend)

#######################################################
# Bird model 
#######################################################
source("R/BIO_bird.R")
source("R/BIO_irradiance_corrections.R")

lambda = seq(0.4,0.7,0.005)
irr.bird = bird(ZEN = zend, lambda)
irr.seacor = I0_corr_season(I0_direct = irr.bird$DIRECT,I0_diffuse = irr.bird$DIFFUSE,Day = 180)
irr.clcor = I0_corr_cloud_BIO(I0_direct = irr.seacor$I0_direct,I0_diffuse = irr.seacor$I0_diffuse,ZEN = zen,Cloud = 0)

#######################################################
# Gregg and Carder model (From Hydrolight)
#######################################################
source("R/Gcmod.r")

rh = 0.8  # humidity default value
am = 1 # aerosol model (1-10), 1 is default value (marine model)
wv = 1.5 # precipitable water
wsm = 4.0 # Average wind speed over 24 hours, default value (m/s)
ws = 5.0 # instant wind speed (m/s)
vis = 25  # visibility in km; enter 15 for default
ro3 = 100 # ozone in Dobson Units, or enter -99 for climatological ozone
# values (computed from the day and Earth position)
iblw = 0 # = 0 to compute above-surface values
# = 1 to compute below-surface values
p0 = 29.92 # reference atmospheric pressure at sea level
pres = 29.9 # Real time atmospheric pressure

rad = 180/pi
pi2 <- 2*pi

sco3 = comp_sco3(ylat = latipxl, xlon = 0,jday = daypxl,ro3 = ro3)
Fo = solar_irr(Fobar,daypxl)

irr.gc = atmodd(iblw=iblw,rad = rad, lam = lam, theta = zend, oza = oza, ag = ag, 
       aw = aw, sco3 = sco3, p=p0, wv = wv, rh = rh, am = am, wsm = wsm,
       ws = ws, vis = vis, Fo = Fo)

irrm.gc = atmodd(iblw=1,rad = rad, lam = lam, theta = 45., oza = oza, ag = ag, 
                aw = aw, sco3 = sco3, p=p0, wv = wv, rh = rh, am = am, wsm = wsm,
                ws = ws, vis = vis, Fo = Fo)
#####################################################################

plot(Edp45$wv,Edp45$Ed,type="l",ylim = c(400,1500),lwd=2,
     main = "Irradiance above the sea surface",
     xlab = "Wavelength (nm)",
     ylab = "Irradiance (mW.m-2.?)")
#lines(lambda*1000,irr.bird$DIFFUSE+irr.bird$DIRECT,lwd=2,col=2)
lines(lambda*1000,irr.clcor$I0_direct+irr.clcor$I0_diffuse,lwd=2,col="red")
lines(lam,irr.gc$Ed,col="Green",lwd=2)
legend("bottomright",c("COART","Bird","G&C"),pch=NA,lwd=2,lty=1,bty="n",
       col = c("black","red","Green"))
#lines(lam[101:401],irrm.gc$Ed[101:401]*1000,col="blue")
lines(zerop$wv,zerop$Ed,col="green")

#################################################################

