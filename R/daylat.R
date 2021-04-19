# Calculate angles based on day and latitude, to be used in zenith angle calculation.
daylat <- function(day,lat) {
    
    # Jan1 = Day 0
    
    # delta  = declination (error in delta equation is less than 3', maximum error is 35")
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
