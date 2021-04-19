
#*******************************************************************************
# Load libraries, custom source code, and define variables

library(insol)
library(lubridate)
library(oce)
library(stringr)
library(fields)

# GCmod.r is the Gregg and Carder model, that is used
# to retrieve the spectral irradiance just under the 
# sea surface: E_d(0-, lambda) in Einstiein.m-2
# We first use the visibility (in km) as a parameter  to 
# scale above-water Satellite PAR to the above-water modeled
# PAR, once the visibility is found, we derive E_d(0-,lambda)
source("R/Gcmod.r")
source("R/Watt_to_quanta.R")

# here We download clear sky irradiance, it will be used to
# derive Ed(lambda,0-) as we sill scale the visibility to match
# model and satellite PAR, the variable is PARday
load("data/LUT_PAR_lati_day_vis.Rdata")
# These are the steps for the 
latibin = seq(35,80,0.25) # x axis for PARday 
daybin = 1:365 # y axis for PARday 
visibility = 10^seq(-1.6,1.8,,30) # z axis for PARday 

# see above for year, day and latitude
# Here we convert PAR from Einstein.m-2.d-1 back to mW.cm-2.microm-2
# see seawifs tech report vol 22 page 49 
SatPAR = 24/1.193 

# We define lambda for the visible spectrum (i.e. PAR)
lambda = 400:700
# Atmospheric parameters
rh = 0.8  # humidity default value
am = 1 # aerosol model (1-10), 1 is default value (marine model)
wv = 1.5 # precipitable water
wsm = 0.0 # Average wind speed over 24 hours, default value (m/s)
ws = 0.0 # instant wind speed (m/s)
vis = 40  # visibility in km; enter 15 for default
ro3 = 100 # ozone in Dobson Units, or enter -99 for climatological ozone
# values (computed from the day and Earth position)
iblw = 1 # = 0 to compute above-surface values
# = 1 to compute below-surface values
p0 = 29.92 # reference atmospheric pressure at sea level
pres = 29.9 # Real time atmospheric pressure

rad = 180/pi
pi2 <- 2*pi

nw = 1.34 # refractive index of seawater



#*******************************************************************************
# FUNCTION

PAR_resolved <- function(latipxl, daypxl, yearpxl) {
    
    # computer sun angle for day
    xhr = 1:23
    datepxl = doyday(yearpxl,daypxl) + hours(xhr)
    julianpxl = JD(doyday(yearpxl,daypxl))
    datepxl = str_replace(as.character(datepxl),"AST","UTC") # Here I put the date to U
    
    # Here we compute the sun zenith angle in degree, which is 90 - elevation
    zendR = 90 - sunAngle(datepxl,longitude = 0, latitude = latipxl)$altitude # Function from "oce" package
    
    # We compute the ozone content, it can be derived from climatology
    # or use a real time value
    sco3 = comp_sco3(ylat = latipxl, xlon = 0,jday = daypxl,ro3 = ro3)
    Fo = solar_irr(Fobar,daypxl) # Extraterrestiral irradiance in W.m-2.nm-1
    
    # Find the visibility *pixvil) that corresponds to  SatPAR for latitude 
    # and day of year:
    latidx = round(latipxl-35)/0.25 # we need to remove 35 deg because the look-up-table starts
    # at 35 deg.
    PARvis = PARday[latidx,daypxl,]
    pxlvis = approx(PARvis,visibility,xout = SatPAR)$y
    
    # Ed(0+,lambda)
    Ed = matrix(0,23,301) # total in mW.m-2.nm-1
    
    #Ed(0-, lambda)
    Edirw = matrix(0,23,301) # direct 
    Edifw = matrix(0,23,301) # diffus
    Edw = matrix(0,23,301) # total in mW.m-2.nm-1
    
    #Eq(0+)
    Eqdir = matrix(0,23,301) # direct
    Eqdif = matrix(0,23,301) # diffus
    Eqd = matrix(0,23,301) # total in Einstein.m-2
    
    # Eq(0-)
    Eqdirw = matrix(0,23,301) # direct
    Eqdifw = matrix(0,23,301) # diffus
    Eqdw = matrix(0,23,301) # total in Einstein.m-2
    
    zendw = rep(23,0)
    
    for (i in 1:23)
    {
        zendw[i] = asin( 1/nw * sin(zendR[i] / rad))
        #    if (zendR[i] >= 90){PARhr[i] = 0}
        if (zendR[i] < 90.)
        {
            
            # above water spectral irradiance        
            irr.gc = atmodd(iblw=0,rad = rad, lam = lam, theta = zendR[i], oza = oza, ag = ag, 
                            aw = aw, sco3 = sco3, p=p0, wv = wv, rh = rh, am = am, wsm = wsm,
                            ws = ws, vis = pxlvis, Fo = Fo)
            Ed[i,] = irr.gc$Ed[51:351]
            
            Eqdir[i,] = Watt_perm_2_to_microMol_per_m2_per_s(Ed = irr.gc$Edir[51:351],
                                                             lambda = lam[51:351])
            Eqdif[i,] = Watt_perm_2_to_microMol_per_m2_per_s(Ed = irr.gc$Edif[51:351],
                                                             lambda = lam[51:351])
            # Under water spectral irradiance.
            irrm.gc = atmodd(iblw=1,rad = rad, lam = lam, theta = zendR[i], oza = oza, ag = ag, 
                             aw = aw, sco3 = sco3, p=p0, wv = wv, rh = rh, am = am, wsm = wsm,
                             ws = ws, vis = pxlvis, Fo = Fo)
            Edirw[i,] = irrm.gc$Edir[51:351] #w.m-2.nm-1
            Edifw[i,] = irrm.gc$Edif[51:351] #w.m-2.nm-1
            
            #We divide the direct irradiance by the in-water sun zenith angle to convert from 
            # planar to scalar
            Eqdirw[i,] = Watt_perm_2_to_microMol_per_m2_per_s(Ed = irrm.gc$Edir[51:351]/cos(zendw[i]),
                                                              lambda = lam[51:351])
            #We divide the diffuse irradiance by 0.833 to convert from 
            # planar to scalar, as in Platt and Sathyendranath 1997, page 2624 last paragraph
            Eqdifw[i,] = Watt_perm_2_to_microMol_per_m2_per_s(Ed = irrm.gc$Edif[51:351]/0.83,
                                                              lambda = lam[51:351])
            
        }
    }
    
    Eqd = Eqdir + Eqdif
    
    # Now we use Ed and SatPar to scale Edw
    # Note that extraterrestrial solar irradiance in Gregg and Carder is in w.cm-2.microm and we want
    # it in mW.m-2.nm, so we *1000 (for W to mW), we /100000 (for cm-2 to m-2) and we /1000 (for microm to nm)
    
    ScalingPARfactor = SatPAR / (sum(Ed)*36)
    
    Eqd = Eqd *SatPAR /(sum(Ed)*3600/1000)
    #Here we have the underwater light field in Ein.m-2 for a given hour.
    Eqdirw = Eqdirw * SatPAR / (sum(Ed)*3600/1000)
    Eqdifw = Eqdifw * SatPAR / (sum(Ed)*3600/1000)
    Eqdw = Eqdirw + Eqdifw
    
    return(list(Eqdifw=Eqdifw, Eqdirw=Eqdirw, Eqdw=Eqdw, zendR=zendR, zendw=zendw, i=i))
    
}
