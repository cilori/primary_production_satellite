# Calculating the hour angle
hour_angle <- function(lat, declin, rads, SunDia, AirRefr) {
    
    # Correction: different sign at S HS
    dfo <- rads*(0.5*SunDia + AirRefr)
    
    if (lat < 0.0) dfo <- -dfo
    fo <- tan(declin + dfo) * tan(lat*rads)
    
    if (fo > 0.99999) fo <- 1.0; # to avoid overflow
    fo <- asin(fo) + pi/2.0
    return (fo)
    
}

# Calculating the hour angle for twilight times
hour_angle_twilight <- function(lat, declin, rads) {
    
    # Correction: different sign at S HS
    df1 <- rads * 6.0; if (lat < 0.0) df1 = -df1;
    fi <- tan(declin + df1) * tan(lat*rads);
    
    if (fi > 0.99999) fi <- 1.0 # to avoid overflow
    fi <- asin(fi) + pi/2.0
    return (fi)
    
}
