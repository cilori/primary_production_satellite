Isurface_laval <- function(LUT,thetas,solz) {
    
    # For the selected solar zenith angle, get surface irradiance for all lambda
    # and solz (solar zenith angles), given the look-up table provided by Laval.
    # UNITS: uE/m^2/s
    # SOURCE: ATBD Belanger et al 2011
    
    # NOTE:
    #   If I0 are all NA, the solar zenith angle might be too high (LUT columns
    #   for theta=75-90 are empty, so anything > 70 will be NA). To "fix" this
    #   temporarily, use rule=2 in the approx function (for values outside the
    #   interval of interpolation, the value at the closest extreme is returned,
    #   instead of NA)
    
    I0 <- apply(LUT,2,function(y) approx(thetas,y,xout=solz,method="linear",rule=2)$y)
    
    return(I0)
    
}
