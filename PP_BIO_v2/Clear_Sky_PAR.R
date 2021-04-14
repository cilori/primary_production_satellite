# Ratio of direct/diffuse irradiance as a function of visibility 
# for 3 sun angle: 20, 40 and 60:
library(fields)
library(insol)
library(lubridate)
library(oce)
library(stringr)

source("Gcmod.r")

# Atmospheric parameters
rh = 0.8  # humidity default value
am = 1 # aerosol model (1-10), 1 is default value (marine model)
wv = 4.5 # precipitable water
wsm = 4.0 # Average wind speed over 24 hours, default value (m/s)
ws = 5.0 # instant wind speed (m/s)
vis = 30.  # visibility in km; enter 15 for default
ro3 = 100 # ozone in Dobson Units, or enter -99 for climatological ozone
# values (computed from the day and Earth position)
iblw = 0 # = 0 to compute above-surface values
# = 1 to compute below-surface values
p0 = 29.92 # reference atmospheric pressure at sea level
pres = 29.9 # Real time atmospheric pressure

rad = 180/pi
pi2 <- 2*pi

sco3 = comp_sco3(ylat = latipxl, xlon = 0,jday = daypxl, ro3 = ro3)
Fo = solar_irr(Fobar,daypxl)

visibility = seq(0.5,50,0.5)
nbv = length(visibility)
ratioEd = matrix(NaN,nbv,301)
Ed550 = rep(NaN,nbv)
for (i in 1:nbv)
{
    vis = visibility[i]
    irr.gc = atmodd(iblw=0,rad = rad, lam = lam, theta = 10., oza = oza, ag = ag, 
                 aw = aw, sco3 = sco3, p=p0, wv = wv, rh = rh, am = am, wsm = wsm,
                 ws = ws, vis = vis, Fo = Fo)
    ratioEd[i,] = irr.gc$Edir[51:351]/irr.gc$Edif[51:351]
    Ed550[i] = irr.gc$Ed[201]
}

image.plot(visibility,400:700,ratioEd,main = "E_dir/E_diff")

plot(visibility,Ed550)


# Let's get the theoretical clear sky PAR for a given latitude:

latibin = seq(35,80,0.25)
daybin = 1:365
visbin = 10^seq(-1.6,1.8,,30)
visibility = visbin
PARday = array(NaN,c(length(latibin),length(daybin),length(visbin)))
xhr = 1:23
yearpxl = 2018

for (k in 1:length(latibin)) # start of loop on latitude
{
    latipxl = latibin[k]
    print(k)
    for (j in 1:length(daybin)) # start of loon on day of year
    {
    daypxl = daybin[j]
    datepxl = doyday(yearpxl,daypxl) + hours(xhr)
    julianpxl = JD(doyday(yearpxl,daypxl))
    datepxl = str_replace(as.character(datepxl),"AST","UTC") # Here I put the date to U

    #Here we compute the sun zenith angle in degree, which is 90 - elevation
    zendR = 90 - sunAngle(datepxl,longitude = 0, latitude = latipxl)$altitude # Function from "oce" package
    for(l in 1:length(visbin))
    { # start loop on visibility
        vis = visbin[l]
        for (i in 1:23)
        {
            if (zendR[i] >= 90){PARhr[i] = 0}
            if (zendR[i] < 90.)
        {
        PARhr[i] =  sum(atmodd(iblw=0,rad = rad, lam = lam, theta = zendR[i], oza = oza, ag = ag, 
                aw = aw, sco3 = sco3, p=p0, wv = wv, rh = rh, am = am, wsm = wsm,
                ws = ws, vis = vis, Fo = Fo)$Ed[51:351])
        }
        } # end of loop on time of day
        res.int = approxfun(xhr,PARhr)
        PARday[k,j,l] = integrate(res.int,lower = 1, upper = 23)$value/100000
    } # end of loop on visibility
    } # end of loop on day of year
} # end of loop on latitudes

image.plot(daybin,latibin,t(PARday[,,30]))


