Isubsurface_laval <- function(chlz,ZEN,Z,lambda,I0,Zmax,awat) {
    
    # Irradiance I0 input units: uEINSTEINS/M^2/S/NM
    
    # Constants for CDOM absorption and NAP (detritus) absorption,
    # from Matsuoka et al. equations (2007, 2011), based on chlorophyll
    Acdom <- 0.0858 # (2007)
    Ecdom <- 0.167 # (2007)
    Scdom <- 0.0180 # (2011)
    Anap <- 0.0131
    Enap <- 0.528
    Snap <- 0.0104 # (2011)
    
    # Water absorption
    awat2 <- awat[2,]
    
    # Get surface planar irradiance.
    Iz_planar <- I0
    
    PUR_z <- rep(0,length(Z))
    PAR_z <- rep(0,length(Z))
    
    for (i in 1:length(Z)) {
        
        #-----------------------------------------------------------------------
        # Absorption and backscattering for depth i
        
        # CDOM absorption (based on chla)
        acdom443_z <- Acdom*chlz[i]^Ecdom # Matsuoka et al. (2007) has tables with Acdom, Ecdom, Scdom, and NAP values
        acdom_z <- ays(acdom443_z,lambda=lambda,Scdom)
        
        # NAP absorption
        anap443_z <- Anap*chlz[i]^Enap
        anaps_z <- ays(anap443_z,lambda=lambda,Snap)
        
        # Phytoplankton absorption
        aph440_z <- 0.0293*chlz[i]^0.557 # Matsuoka et al. (2011)
        aph_z <- alp(aph440_z,alpbet) # vector of values for each wavelenth
        aph_z_norm <- aph_z/aph440_z # normalized
        
        # Total water absorption at this depth for each wavelength. Matsuoka et al. 2011
        atot_Matsuoka <- awat2 + acdom_z + anaps_z + aph_z
        
        # Total backscattering
        bp_555 <- 0.004*chlz[i]^0.357 # Wang et al. (2005)
        bb_total <- bb(bp_555,lambda=lambda,bw_550=0.00193)
        
        
        #-----------------------------------------------------------------------
        # Kd, Iz, PAR, and PUR for this depth
        
        # Compute Kd (coefficient of attenuation)
        # Lee et al 2005a, equation 11: Approximation of Kd (mean Kd) in euphotic zone
        Kd_Matsuoka <- (1+0.005*ZEN)*atot_Matsuoka+4.18*(1-0.52*exp(-10.8*atot_Matsuoka))*bb_total
        
        # Light attenuation according to Kd_Matsuoka_Sal
        # Need to use previous depth Iz_planar instead of I0 (surface planar Iz)
        # to account for changes in the attenuation at different depths, if
        # profile is inhomogeneous.
        Iz_planar <- Iz_planar*exp(-Kd_Matsuoka)
        
        # Morel 1991(?): Convert planar irradiance to scalar irradiance and remove NA values.
        # Iz_planar = |cos| * Iz_scalar
        # |cos| value from this Kd approximation: Kd = (a+b)/cos(solar_zenith_angle))
        Iz_scalar <- Iz_planar*(Kd_Matsuoka/(atot_Matsuoka+bb_total))
        Iz_scalar[is.na(Iz_scalar)] <- 0
        
        # Compute PAR and PUR at this depth by integrating across all wavelengths
        # PAR = photosynthetically active radiation
        # PUR = photosynthetically usable radiation, computed by multiplying
        #       irradiance by phytoplankton absorption normalized to phytoplankton
        #       absorption at 440nm.
        PUR_z[i] <- trapz(lambda,aph_z_norm*Iz_scalar)
        PAR_z[i] <- trapz(lambda,Iz_scalar)
        
        #-----------------------------------------------------------------------
        # End calculations if we've reached euphotic depth (when PAR <= 1% surface PAR)
        if ((i > 1) & (PAR_z[i] <= 0.01*PAR_z[1])) {
            # Return PUR, PAR, euphotic depth, and index of euphotic depth
            return(list(PUR_z,PAR_z,Z[i],i))
        }
        
    }
    
}
