# SURFACE IRRADIANCE (direct and diffuse components at different wavelengths)

bird <- function(ZEN,wl) {
    
    # ZEN   = solar zenith angle, radians
    # wl    = wavelengths where direct and diffuse should be computed, micrometres
    # 
    # Compute direct and diffuse components of spectral irradiance at sea level,
    # given solar zenith angle in air (< 80 degrees). The direct component is
    # calculated normal to the sun's direction.
    # UNITS: watts * metre^-2 * micrometre^-1
    # SOURCE: Bird 1984, Bird and Riordan 1986
    # 
    # AIRMAS = AIRMASS
    # C1     = C, INTERPOLATED FOR THE ZENITH ANGLE OF INTEREST
    # C2     = C1, INTERPOLATED FOR THE 24 LAMBDAS.
    # DIF    = DIFFUSE SPECTRAL IRRADIANCE AT 24 WAVELENGTHS
    # DIFUSE = DIFFUSE SPECTRAL IRRADIANCE AT 61 WAVELENGTHS
    # DIR    = DIRECT  SPECTRAL IRRADIANCE AT 24 WAVELENGTHS
    # DIRECT = DIRECT  SPECTRAL IRRADIANCE AT 61 WAVELENGTHS
    # EM0    = M(0), THE AIRMASS EXPRESSION FOR OZONE
    # ROS    = RO (S), THE ALBEDO OF AIR.
    # TA     = AEROSOL TRANSMITTANCE
    # TO     = OZONE TRANSMITTANCE
    # TR     = RAYLEIGH TRANSMITTANCE
    # TW     = WATER VAPOUR TRANSMITTANCE
    # ZEND   = ZENITH ANGLE IN DEGREES
    
    #***************************************************************************
    
    # Value for alpha2 corrected March 1993, changed again Feb 1998
    ALPHA1 <- 1.0274
    BETA1 <- 0.13238
    ALPHA2 <- 1.206
    BETA2 <- 0.116981
    
    # LAMBDA = WAVELENGTHS AT WHICH TRANSMITTANCE IS CALCULATED 
    LAMBDA <- c(seq(400,550,by=10),570,593,610,630,656,667.6,690,710)
    LAMPTH <- LAMBDA/1000
    
    # EXTER = EXTRA TERRESTRIAL SPECTRAL IRRADIANCE
    EXTER <- c(1479.1,1701.3,1740.4,1587.2,1837.0,2005.0,2043.0,1987.0,2027.0,
               1896.0,1909.0,1927.0,1831.0,1891.0,1898.0,1892.0,1840.0,1768.0,
               1728.0,1658.0,1524.0,1531.0,1420.0,1399.0)
    
    # AV = WATER VAPOUR ABSORPTION COEFFICIENT
    # Feb 1998 - corrected value at wavelength 710 from 0.05 to 0.0125
    AV <- c(rep(0,17),0.075,0,0,0,0,0.016,0.0125)
    
    # AO = OZONE ABSORPTION COEFFICIENT 
    # Feb 1998 - corrected value at wavelengths 550, 610, 630, 6676 from 0.095 to 0.085,
    # from 0.132 to 0.120, from 0.120 to 0.090, and from 0.060 to 0.051 respectively  
    AO <- c(0,0,0,0,0,0.003,0.006,0.009,0.014,0.021,0.030,0.040,0.048,0.063,
            0.075,0.085,0.120,0.119,0.120,0.090,0.065,0.051,0.028,0.018)
    
    # AU = ABSORPTION COEFFICIENT AND GASEOUS AMOUNT
    AU <- c(rep(0,22),0.15,0)
    
    
    #***************************************************************************
    # DIRECT NORMAL IRRADIANCE FOR 24 WAVELENGTHS
    
    # Compute direct irradiance transmittance components using AIRMAS=1.9,
    # then use them to compute ROS (air albedo).
    
    # Bird and Riordan, 1986: M0 = (1 + ho/6370) / sqrt( (cos(ZEN)^2) + 2*ho/6370),
    # where ho=22km, approx max height of ozone
    # Below: TO = exp(-AO * 0.344 * EM0)
    EM0 <- 35. / sqrt(1224. * (cos(ZEN))^2 + 1.)
    
    AIRMAS <- 1.9
    W <- 2 # precipitable water vapor in a vertical path (cm), value from George's script
    
    # RAYLEIGH SCATTERING
    TR <- exp( -AIRMAS / (LAMPTH^4 * ( 115.6406 -1.335/LAMPTH^2 )))
    # WATER VAPOUR ABSORPTION
    # Feb 1998 - corrected TW equation as per Bird and Riordan 1986
    TW <- exp( (-0.2385 * AV * W *  AIRMAS) / ((1. + 20.07 * AV * W * AIRMAS)^0.45))
    # OZONE
    # From Campbell and Aarup 1989, assume O3=0.344cm (corrected Feb 1998),
    # where O3 = ozone amount in atmosphere
    TO <- exp (-AO*0.344*EM0)
    # AEROSOL SCATTERING AND ABSORPTION (TURBIDITY ASSUMED TO BE 0.27)
    TA <- rep(NA,24)
    TA[1:10] <- exp(-BETA1 * LAMPTH[1:10]^(-ALPHA1) * AIRMAS)
    TA[11:24] <- exp(-BETA2 * LAMPTH[11:24]^(-ALPHA2) * AIRMAS)
    # UNIFORMLY MIXED GAS ABSORPTION
    # Note: AU = 0 for all but lambda=690nm, so TU = 1 for those values
    # Feb 1998 - Corrected Leckner's value from 118.93 to 118.3 as described in
    # Bird and Riordan, 1986
    # Also, AU(690) changed from .30 to .15
    TU <- exp(-1.41 * AU * AIRMAS / ((1. + 118.3 * AU * AIRMAS)^0.45) )
    # ALBEDO OF AIR
    ROS <- TO * TW * (TA * (1. - TR) *0.5 + TR * (1. - TA) * 0.22 * 0.928) * TU
    
    
    # Get ZEND (zenith angle in degrees), and recompute AIRMAS and all direct
    # irradiance transmittance components that use the AIRMAS variable.
    # ***NOTE: do NOT recompute ROS, just use the airmas=1.9 ROS after this
    ZEND = ZEN * 180/pi
    AIRMAS <- 1 / (cos(ZEN) + 0.15*(93.885 - ZEND)^(-1.253))
    
    if (!is.finite(AIRMAS)) {return(list(DIRECT=NA,DIFFUSE=NA))}
    
    if (AIRMAS < 1) {AIRMAS <- 1}
    
    TR <- exp(-AIRMAS / (LAMPTH^4 * ( 115.6406 -1.335/(LAMPTH^2) )))
    TW <- exp( (-0.2385 * AV * W *  AIRMAS) / ((1. + 20.07 * AV * W * AIRMAS)^0.45))
    TA <- rep(NA,24)
    TA[1:10] <- exp(-BETA1 * LAMPTH[1:10]^(-ALPHA1) * AIRMAS)
    TA[11:24] <- exp(-BETA2 * LAMPTH[11:24]^(-ALPHA2) * AIRMAS)
    TU <- exp (-1.41 * AU * AIRMAS / ((1. + 118.3 * AU * AIRMAS)^0.45) )
    # TOTAL DIRECT IRRADIANCE
    DIR <- EXTER * TR * TA * TW * TO * TU
    
    
    # #***************************************************************************
    # # CORRECTIONS - Nov 2020, no longer using this
    # 
    # # C = CORRECTION FACTRS FOR DIFFUSE IRRADIANCE, FOR 5 WAVELENGTHS AND 7 ZENITH ANGLES
    # C <- c(1.11,1.04,1.15,1.12,1.32,1.13,1.05,1.00,0.96,1.12,1.18,1.09,1.00,0.96,
    #        1.07,1.24,1.11,0.99,0.94,1.02,1.46,1.24,1.06,0.99,1.10,1.70,1.34,1.07,
    #        0.96,0.90,2.61,1.72,1.22,1.04,0.80)
    # # 5 rows (wavelengths), 7 columns (zenith angles)
    # # Fortran writes the vector above to a 5x7 matrix by column, like R
    # # http://www.mathcs.emory.edu/~cheung/Courses/561/Syllabus/6-Fortran/array1.html
    # C <- matrix(C,nrow=5)
    # 
    # # CD = ZENITH ANGLE (IN DEGREES) FOR WHICH C VALUES ARE AVAILABLE.
    # CD <- c(0,37,48.19,60,70,75,80)
    # 
    # # Changes made Sep 98 - Heidi Maass  
    # for (K in 2:7) {
    #     if (ZEND < CD[K]) {
    #         KK <- K
    #         break
    #     }
    # }
    # 
    # if (abs(ZEND-80) < 1e-8) {KK <- 7}
    # 
    # FRCTN <- (ZEND - CD[KK-1]) / (CD[KK] - CD[KK-1])
    # 
    # C1 <- rep(NA,5)
    # for (L in 1:5) {
    #     C1[L] <- (C[L,KK] - C[L,KK-1]) * FRCTN + C[L,KK-1]
    # }
    # 
    # C2 <- rep(NA,24)
    # C2[1] <- C1[1]
    # L <- 1
    # for (L1 in 2:4) {
    #     CINC <- (C1[L1-1] - C1[L1])/5
    #     for (L2 in 1:5) {
    #         L <- L+1
    #         C2[L] <- C2[L-1] - CINC
    #     }
    # }
    # 
    # for (L1 in 17:24) {
    #     WLINC <- (LAMBDA[L1] - LAMBDA[L1-1]) / 160
    #     CINC <- (C1[5] - C1[4]) * WLINC
    #     L <- L+1
    #     C2[L] <- C2[L-1] + CINC
    # }
    
    
    #***************************************************************************
    # DIFFUSE HORIZONTAL IRRADIANCE FOR 24 WAVELENGTHS
    
    XX <- EXTER * cos(ZEN) * TO * TU * TW
    # Rayleigh scattering component
    R <- XX * TA*(1-TR)*0.5
    # Aerosol scattering component
    A <- XX * TR*(1-TA)*0.928*0.82
    # Component accounting for multiple reflection of irradiance between ground and air
    # G <- (DIR * cos(ZEN) + (R+A) * C2) * 0.05 * ROS / (1 - 0.05 * ROS)
    G <- (DIR * cos(ZEN) + (R+A)) * 0.05 * ROS / (1 - 0.05 * ROS)
    # Total, including correction factor for sum of R and A
    # DIF <- (R+A) * C2 + G
    DIF <- (R+A) + G
    
    #***************************************************************************
    # INTERPOLATION TO OTHER WAVELENGTHS wl
    
    DIRECT <- approx(LAMPTH,DIR,xout=wl)$y
    DIFFUSE <- approx(LAMPTH,DIF,xout=wl)$y
    
    return(list(DIRECT=DIRECT,DIFFUSE=DIFFUSE))
    
}
