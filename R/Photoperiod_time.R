# TESTING - COMPARING TO THE HFX STANFIELD INTL AIRPORT SUN RISE/SET TIMES AT
# ENVIRONMENT CANADA WEATHER WEBSITE
# https://weather.gc.ca/city/pages/ns-19_metric_e.html

# > Photoperiod_time(latit=44.886021,longit=-63.509424,y=2019,m=2,day=14,tzone=-4)
# $riset
# [1] 7.322821 = 7:19:37AM (Environment Canada: sunrise = 7:18)
# 
# $settm
# [1] 17.61813 = 5:37:09PM (Environment Canada: sunset = 5:40)
# 
# $noont
# [1] 12.47047
# 
# $daylen
# [1] 10.29531

##################################################################################

# tzone = hours relative to GMT (i.e. local time = GMT + tzone)
Photoperiod_time <- function(latit=latit,longit=longit,y=y,m=m,day=day,tzone=NA) {
    
    h <- 12
    degs <- 180.0/pi
    rads <- pi/180.0
    SunDia <- 0.53            # Sunradius degrees
    AirRefr <- 34.0/60.0      # Atmospheric refraction degrees
    
    d <- days_to_J2000(y, m, day, h)
    
    # Use sun_lon to find the ecliptic longitude of the Sun
    
    lambda <- sun_lon(d,rads)$DD
    L <- sun_lon(d,rads)$L
    
    # Obliquity of the ecliptic
    
    obliq <- 23.439 * rads - .0000004 * rads * d
    
    # Find the RA and DEC of the Sun
    
    alpha <- atan2(cos(obliq) * sin(lambda), cos(lambda))
    delta <- asin(sin(obliq) * sin(lambda))
    
    # Find the Equation of Time in minutes
    # Correction suggested by David Smith
    
    LL <- L - alpha;
    if (L < pi) LL <- LL + (2*pi)
    equation <- 1440.0 * (1.0 - LL / (2*pi))
    
    ha <- hour_angle(latit,delta,rads,SunDia,AirRefr)
    
    if (is.na(ha)) { 
        daylen <- 0
        riset <- 0
        settm <- 0
        noont <- 0
    } 
    
    if (!is.na(ha)) {
        
        hb <- hour_angle_twilight(latit,delta,rads)
        twx <- hb - ha;   # length of twilight in radians
        twx <- 12.0*twx/pi;      # length of twilight in degrees
        
        # Conversion of angle to hours and minutes
        daylen <- degs * ha / 7.5;
        if (daylen < 0.0001) { daylen = 0.0 }
        
        riset <- 12.0 - 12.0 * ha/pi + tzone - longit/15.0 + equation/60.0
        settm <- 12.0 + 12.0 * ha/pi + tzone - longit/15.0 + equation/60.0
        noont <- riset + 12.0 * ha/pi
        altmax <- 90.0 + delta * degs - latit
        
        # Correction suggested by David Smith to express as degrees from the N horizon        
        if (delta * degs > latit ) { altmax <- 90.0 + latit - delta * degs }    
        twam <- riset - twx      # morning twilight begin
        twpm <- settm + twx      # evening twilight end
        if (riset > 24.0) { riset == 24.0 }
        if (settm > 24.0) { settm == 24.0 }
        
    }
    
    list (riset=riset,settm=settm,noont=noont,daylen=daylen)
    
}
