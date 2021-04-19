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
