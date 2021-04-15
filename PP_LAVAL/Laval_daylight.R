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

# Get the days to J2000
# h is UT in decimal hours
# FNday only works between 1901 to 2099 - see Meeus chapter 7
FNday <- function(y, m, day, h) {
    y <- as.integer(y)
    m <- as.integer(m)
    day <- as.integer(day)
    luku <- as.integer(- 7 * (y + (m + 9)/12)/4 + 275*m/9 + day)
    # Typecasting needed for TClite on PC DOS at least, to avoid product overflow
    luku <- luku + as.integer(y*367)
    luku <- luku - 730531.5 + h/24.0
    return (luku)
}

# The function below returns an angle in the range 0 to 2*pi
FNrange <- function(x) {
    b <- x / (2*pi)
    a <- (2*pi) * (b - as.integer(b))
    if (a < 0) a <- (2*pi) + a
    return (a)
}

# Calculating the hour angle
f0 <- function(lat, declin, rads, SunDia, AirRefr) {
    
    # Correction: different sign at S HS
    dfo <- rads*(0.5*SunDia + AirRefr)
    
    if (lat < 0.0) dfo <- -dfo
    fo <- tan(declin + dfo) * tan(lat*rads)
    
    if (fo > 0.99999) fo <- 1.0; # to avoid overflow
    fo <- asin(fo) + pi/2.0
    return (fo)
    
}

# Calculating the hour angle for twilight times
f1 <- function(lat,declin, rads) {
    
    # Correction: different sign at S HS
    df1 <- rads * 6.0; if (lat < 0.0) df1 = -df1;
    fi <- tan(declin + df1) * tan(lat*rads);
    
    if (fi > 0.99999) fi <- 1.0 # to avoid overflow
    fi <- asin(fi) + pi/2.0
    return (fi)
    
}

# Find the ecliptic longitude of the Sun 
FNsun <- function (d, rads) {
    
    # mean longitude of the Sun
    L <- FNrange(280.461 * rads + .9856474 * rads * d)
    
    # mean anomaly of the Sun
    g <- FNrange(357.528 * rads + .9856003 * rads * d)
    
    DD <- FNrange(L + 1.915 * rads * sin(g) + .02 * rads * sin(2 * g))
    
    # Ecliptic longitude of the Sun and the mean longitude of the Sun
    list (DD=DD, L=L)
    
}

# tzone = hours relative to GMT (i.e. local time = GMT + tzone)
Photoperiod_time <- function(latit=latit,longit=longit,y=y,m=m,day=day,tzone=NA) {
    
    h <- 12
    degs <- 180.0/pi
    rads <- pi/180.0
    SunDia <- 0.53            # Sunradius degrees
    AirRefr <- 34.0/60.0      # Atmospheric refraction degrees
    
    d <- FNday(y, m, day, h)
    
    # Use FNsun to find the ecliptic longitude of the Sun
    
    lambda <- FNsun(d,rads)$DD
    L <- FNsun(d,rads)$L
    
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
    
    ha <- f0(latit,delta,rads,SunDia,AirRefr)
    
    if (is.na(ha)) { 
        daylen <- 0
        riset <- 0
        settm <- 0
        noont <- 0
    } 
    
    if (!is.na(ha)) {
        
        hb <- f1(latit,delta,rads)
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

daylat <- function(day,lat) {
    
    # Jan1 = Day 0
    
    # From Fortran zenith angle subroutine:
    # delta  = declination
    #          (error in delta equation is less than 3',
    #           maximum error is 35")
    # lat    = latitude
    # phi    = latitude in radians
    
    theta = (2 * pi * day)/365
    delta = (0.006918 - 0.399912 * cos(theta) +  0.070257 * sin(theta)
             - 0.006758 * cos(2*theta) + 0.000907 * sin(2*theta)
             - 0.002697 * cos(3*theta) + 0.001480 * sin(3*theta))
    phi = lat * (pi/180)
    phidel = -tan(phi) * tan(delta)
    
    return(list(theta=theta,delta=delta,phi=phi,phidel=phidel))
    
}

