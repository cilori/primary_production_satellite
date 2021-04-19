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
