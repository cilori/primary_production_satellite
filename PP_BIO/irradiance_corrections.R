# CORRECT FOR SEASONAL VARIATION IN SOLAR ENERGY
# Input: surface direct and diffuse irradiance (each a vector of values across
# wavelengths), and the Day, a numeric ZERO-INDEXED day of year (0-364).
I0_corr_season <- function(I0_direct,I0_diffuse,Day) {
    
    # DayNo   = DAY NUMBER
    DayNo = c(0,3,31,42,59,78,90,93,120,133,151,170,181,183,206,212,243,265,273,
              277,304,306,334,355,365)
    
    # SolIr   = EXTRATERRESTRIAL SOLAR IRRADIANCE
    # Values are from Thekaekara 1977 (ref in JGR 1988 - s & p)
    SolIr = c(1399.,1399.,1393.,1389.,1378.,1364.,1355.,1353.,1332.,1324.,1316.,1310.,
              1309.,1309.,1312.,1313.,1329.,1344.,1350.,1353.,1347.,1375.,1392.,1398.,1399.)
    
    # Find nearest day index that is greater or equal to satellite Day.
    D = which(DayNo >= Day)[1]
    
    # Set value of solar irradiance according to SolIr data array
    if (D==1) {
        sol = SolIr[1]
    } else {
        sol = SolIr[D-1]-(SolIr[D-1]-SolIr[D]) * ((Day-DayNo[D-1])/(DayNo[D]-DayNo[D-1]))
    }
    
    # Adjust direct and diffuse components by fractional correction
    # Note: 1353. is max value (100%)
    I0_direct_corrected = I0_direct * sol/1353
    I0_diffuse_corrected = I0_diffuse * sol/1353
    
    return(list("I0_direct"=I0_direct_corrected,"I0_diffuse"=I0_diffuse_corrected))
    
}


# CLOUD COVER CORRECTION (BIO METHOD)
# Source: Platt 1991, equation 6 (taken from Paltridge and Platt, 1976?)
# Input: surface direct and diffuse irradiance (each a vector of values across
# wavelengths), a numeric value for ZEN (the solar zenith angle in air), and
# Cloud, a numeric value between 0 and 1 representing fraction of cloud cover.
I0_corr_cloud_BIO <- function(I0_direct,I0_diffuse,ZEN,Cloud) {
    
    # Components of clear sky irradiance.
    Idir_clear = sum(I0_direct) * cos(ZEN)
    Idif_clear = sum(I0_diffuse)
    
    # Direct surface irradiance, accounting for cloud cover
    Rd = 1 - Cloud
    Idir = Idir_clear * Rd
    
    # Atmospheric albedo
    Albedo = 0.28/(1 + 6.43 * cos(ZEN))
    
    # Ratio of total surface irradiance to clear sky total surface irradiance
    R = (1-0.5*Cloud)*(0.82-Albedo*(1-Cloud)) / (0.82 - Albedo)
    
    # Diffuse surface irradiance, accounting for cloud cover
    Itotal_clear = Idir_clear + Idif_clear
    Idif = Itotal_clear * R - Idir
    Rs = Idif/Idif_clear
    
    # Adjust direct and diffuse components of irradiance across spectrum
    I0_direct_corrected = I0_direct * Rd
    I0_diffuse_corrected = I0_diffuse * Rs
    
    return(list("I0_direct"=I0_direct_corrected,"I0_diffuse"=I0_diffuse_corrected))
    
}


# CALCULATE REFLECTION LOSSES AT AIR-SEA INTERFACE
# Sathyendranath et al 1998, "The Model" (page 9270)
# Input: surface direct and diffuse irradiance (each a vector of values across
# wavelengths), and a numeric value for ZEN (the solar zenith angle in air)
I0_corr_refl_loss <- function(I0_direct,I0_diffuse,ZEN) {
    
    # Solar zenith angle in water
    ZENW = asin(sin(ZEN)/1.333)
    
    # Reflection - Fresnel's equation
    # Found in Kirk 1983, Light and Photosynthesis in Aquatic Ecosystems, equation 2.12 in section 2.5 (page 38)
    ref = 0.5*(sin(ZEN-ZENW))^2/(sin(ZEN+ZENW))^2 + 0.5*(tan(ZEN-ZENW))^2/(tan(ZEN+ZENW))^2
    
    # Adjust surface irradiance across spectrum for reflection losses
    I0_direct_corrected = I0_direct * (1 - ref)
    I0_diffuse_corrected = I0_diffuse * 0.945
    
    return(list("I0_direct"=I0_direct_corrected,"I0_diffuse"=I0_diffuse_corrected, "ZENW"=ZENW))
    
}
