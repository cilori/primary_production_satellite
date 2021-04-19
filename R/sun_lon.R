# Find the ecliptic longitude of the Sun 
sun_lon <- function (d, rads) {
    
    # mean longitude of the Sun
    L <- pos_angle(280.461 * rads + .9856474 * rads * d)
    
    # mean anomaly of the Sun
    g <- pos_angle(357.528 * rads + .9856003 * rads * d)
    
    DD <- pos_angle(L + 1.915 * rads * sin(g) + .02 * rads * sin(2 * g))
    
    # Ecliptic longitude of the Sun and the mean longitude of the Sun
    list (DD=DD, L=L)
    
}
