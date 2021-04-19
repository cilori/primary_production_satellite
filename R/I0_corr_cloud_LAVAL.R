# CLOUD COVER CORRECTION (LAVAL METHOD)
# Srikanth et al 2018, equation 4: Adjust I0 based on total cloud cover,
# using the I0 based on a clear sky (taucl=0)
# NOTE: IF USING O3=300 AND TAUCL=5, THIS MUST BE CHANGED******
# ALSO NOTE: This is for TOTAL surface irradiance, not separated into direct
# and diffuse components.
I0_corr_cloud_LAVAL <- function(I0,TCC) {
    return((1 - 0.29*(TCC + TCC*TCC)) * I0)
}
