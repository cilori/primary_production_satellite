# ORIGINAL BIO MODEL
# NEAREST-NEIGHBOUR METHOD, PLATT ET AL 2008
# WRITTEN IN FORTRAN BY GEORGE N WHITE III
# TRANSLATED TO R BY STEPHANIE CLAY, APRIL 2019
#
#
# SEE NEAR THE BOTTOM OF THE SCRIPT FOR EXAMPLE
#
# This script contains dwcpn(), the main primary production function for the BIO
# model, and dwcpn_C(), which is the same function with the exception that it uses
# bird() and sam_penguin() functions translated to C++ instead of R (this is faster).
#
#*******************************************************************************

# source("R/bird.R")
# source("R/sam_penguin.R") # exp_m600 and yellow_substances are used in sam_penguin
# source("R/yellow_substances.R")
# source("R/exp_m600.R")
# source("R/I0_corr_season.R")
# source("R/I0_corr_cloud_BIO.R")
# source("R/I0_corr_refl_loss.R")
# source("R/daylat.R")
# source("R/zenith_bio.R")
# library(caTools)

dwcpn <- function(input_df, uniform_profile=TRUE, bp=NULL, gauss=TRUE, cloud_corr=TRUE) {
    
    # INPUT:
    #   input_df = numeric vector or dataframe row containing the following values in this order:
    #           --doy (day of year)
    #           --longitude
    #           --latitude
    #           --bathymetry
    #           --chlorophyll-a
    #           --photosynthetically active radiation (PAR)
    #           --total cloud cover (between 0 and 1)
    #           --alphaB (PI curve parameter)
    #           --PmB (PI curve parameter)
    #           If uniform_profile=FALSE and gauss=TRUE, the following biomass profile parameters:
    #               --z_m (max depth)
    #               --B_0 (background chlorophyll-a)
    #               --h (related to the height of the curve)
    #               --sigma (related to the width of the curve)
    #   uniform_profile = TRUE if biomass profile is uniform, FALSE otherwise
    #   bp = optional numeric matrix with column 1 = depths and column 2 = chla,
    #        to be used only if you want a nonuniform profile but NOT a gaussian profile
    #        (i.e. if uniform_profile=FALSE and gauss=FALSE). This is used to interpolate
    #        values between known measurements to every depth in the profile.
    #        example format:     bp <- cbind(seq(1,100,by=10),   # depths (metres)
    #                                        rep(1,10))          # chla
    #   gauss = TRUE to model a gaussian biomass profile (if uniform_profile = FALSE)
    #   cloud_corr = TRUE to do cloud correction of incoming radiance
    
    
    # OUTPUT:
        # Primary production
        # Chlorophyll-a integrated over the water column
        # PAR integrated over water column
        # Surface PAR
        # Spectral irradiance XX integrated over water column
        # Euphotic depth, integrated over time
        # Yellow substances, integrated over depth and time (calculated in sam_penguin)
    
    
    #***************************************************************************
    # DEFINE VARIABLES
    
    LAMBDA = seq(400,700,by=5) # Vector of wavelengths
    dz = 0.5                   # Depth interval between computations (metres)
    Zmin = 0                   # Min depth
    Zmax = 250                 # Max depth
    Z = seq(Zmin,Zmax,by=dz)   # Vector of depths
    Rsize = length(LAMBDA)     # Number of wavelengths (400 to 700 nm, 5nm increments)
    nstep = length(Z)          # Length of depths vector, where dz=0.5m and Zmax=250m
    
    input_df <- as.numeric(input_df)
    doy = input_df[1]              # Day of year
    long = input_df[2]             # Longitude
    lat = input_df[3]              # Latitude
    BotZ = input_df[4]             # Bathymetry
    PARsat = input_df[6]           # Satellite PAR for scaling bird PAR
    if (cloud_corr) {
        Cloud = input_df[7]        # Fraction of cloud cover (between 0 and 1)
    }
    alphaB = input_df[8]           # Initial slope of PI curve
    P_mB = input_df[9]             # Assimilation number (maximum photosynthetic rate)
    
    
    #***************************************************************************
    # CHLA PROFILE ####
    if (uniform_profile) {
        chl = input_df[5]
        chlz = rep(chl,length(Z))
    } else {
        # Different methods for creating chla depth profile
        #   1. shifted gaussian formula
        #       (if parameters missing, return NA)
        #   2. profile given by user (interpolation to missing points)
        #       (if bp matrix missing, return NA)
        
        if (gauss) { # use shifted gaussian
            z_m = input_df[10]
            B_0 = input_df[11]
            h = input_df[12]
            sigma = input_df[13]
            chlz = shifted_gaussian(tv=Z, B0=B_0, h=h, sigma=sigma, tmax=z_m)
            exp_term = 0.5 * ((Z - z_m)/sigma)^2
            chlz[exp_term > 675] = B_0 #Avoid values too low?
        } else { # use interpolation
            chlz <- approx(bp[,1], bp[,2], xout=Z, method="linear",rule=2)$y
        }
        
        if (all(is.na(chlz) | (chlz < 0))) {
            return(c("PP_BIO"=NA,
                     "Chlz_BIO"=NA,
                     "PAR_BIO"=NA,
                     "PAR0_BIO"=NA,
                     "XX_BIO"=NA,
                     "EuphoD_BIO"=NA,
                     "yelsub_BIO"=NA))
        }
    }
    
    
    #***************************************************************************
    # TIME AND ZENITH ANGLE STUFF ####
    
    # Compute variables based on day and latitude, used in calculation of time of
    # dawn and solar zenith angle.
    Day = doy-1
    dl_vars = daylat(day=Day,lat=lat) # Jan1 = Day 0
    delta = dl_vars$delta
    phi = dl_vars$phi
    phidel = dl_vars$phidel
    
    # Check for 24-hour darkness or sunlight
    if (phidel < -1) { # 24hr sunlight
        phidel = -1
    } else if (phidel > 1) { # 24hr darkness
        return(c("PP_BIO"=0,
                 "Chlz_BIO"=trapz(Z,chlz),
                 "PAR_BIO"=0,
                 "PAR0_BIO"=0,
                 "XX_BIO"=0,
                 "EuphoD_BIO"=0,
                 "yelsub_BIO"=0))
    }
    
    # Compute dawn
    Dawn = 12 - acos(phidel)*(180/pi)/15
    T_start = 0                             #Calculation start index (i.e. first index where ZEND < 80)
    Delta_T = (12-Dawn)/12                  #Compute time step
    T_seq = seq(Dawn+Delta_T,12,by=Delta_T) #Sequence of time intervals
    
    
    #***************************************************************************
    # TIME LOOP PART 1 ####
    # Dawn to noon, doubled at the end
    
    EuphoD = rep(0,12)      #Euphotic depth
    PAR0 = rep(0,12)        #Surface solar irradiance
    # Variables below are integrated over depth within the time loop:
    PP = rep(0,12)          #Carbon primary production, nonuniform profile
    PAR = rep(0,12)         #Photosynthetically active radiation (irradiance integrated over wavelength at each depth)
    XX = rep(0,12)          #Irradiance used in PP formula
    I0 = rep(0,12)          #Surface irradiance
    yelsub = rep(0,12)      #yellow substances in the water column (for testing)
    
    ZEN_hr = rep(0,12)
    DIR = matrix(0,nrow=Rsize,ncol=12)
    DIF = matrix(0,nrow=Rsize,ncol=12)
    
    
    for (T_idx in 1:12) {

        # Get next time interval
        Time = T_seq[T_idx]

        #***********************************************************************
        # ZENITH ANGLE ####

        # Compute zenith angle of sun in radians and convert to degrees.
        ZEN = zenith_bio(Time, long, delta, phi)
        ZEND = ZEN * (180/pi)

        ZEN_hr[T_idx] = ZEN

        # Restrict to [ 0 < zenith < 80 ] because bird_revised can't accept ZEND > 80
        # NOTE: in Fortran script, looks like ZEND > 80 was converted to ZEND = 79 and used
        if (ZEND < 0 | ZEND > 80) {next} # Skip to next time step
        
        # T_start is T_idx value when ZenD becomes < 80 degrees
        if (T_start==0) {T_start = T_idx}
        
        #***********************************************************************
        # SURFACE IRRADIANCE ####

        # Compute clear sky direct and diffuse components of spectral irradiance
        # at sea level, for a given zenith angle (zenith < 80 degrees).
        # UNITS: watts * metre^-2 * micrometre^-1
        DIR_DIF = bird(ZEN=ZEN,wl=seq(0.4,0.7,by=0.005))
        I0_direct = DIR_DIF[[1]]
        I0_diffuse = DIR_DIF[[2]]
        
        # CORRECT FOR SEASONAL VARIATION IN SOLAR ENERGY
        I0_season <- I0_corr_season(I0_direct,I0_diffuse,Day)
        I0_direct <- I0_season$I0_direct
        I0_diffuse <- I0_season$I0_diffuse
        
        if (cloud_corr) {
            # CLOUD COVER CORRECTION
            I0_cloud <- I0_corr_cloud_BIO(I0_direct,I0_diffuse,ZEN,Cloud)
            I0_direct <- I0_cloud$I0_direct
            I0_diffuse <- I0_cloud$I0_diffuse
        }

        DIR[,T_idx] = I0_direct
        DIF[,T_idx] = I0_diffuse
        
        
        # If you don't want to correct by satellite PAR, comment from this star
        # line down to the star line at the start of time loop part 2 (below)
        #***********************************************************************

    }
    
    
    #***************************************************************************
    # ADJUST I0_DIRECT AND I0_DIFFUSE BASED ON SATELLITE PAR ####
    
    temp_par = (DIR + DIF)
    
    # Convert units of DIR and DIF: watts/m^2/micron --> Einsteins/m^2/hr/nm
    # (Integrating over wavelength removes nm, and then integrating over day
    # gives the full PAR per day, same as NASA)
    # MAKE SURE YOU'RE USING FILES DIRECTLY FROM NASA THAT HAVE PAR IN UNITS OF
    # EINSTEIN/M^2/DAY
    # Similar conversion steps are outlined in the time loop below.
    unit_scale = LAMBDA * 36 / (6.6261 * 2.9979 * 6.022 * 10^7)
    temp_par = temp_par * unit_scale
    temp_lam = LAMBDA
    
    # Integrate over wavelength
    temp_par = as.numeric(sapply(1:ncol(temp_par),function(i) {trapz(temp_lam,temp_par[,i])}))
    
    # Integrate over time using the method described near the bottom of this
    # script, which is dependent on the first time interval where ZEN > 80
    # degrees, given by the index T_start.
    time_interval = seq(Dawn + Delta_T, 12, by=Delta_T)
    t1 = c(Dawn,time_interval[T_start])
    t2 = time_interval[T_start:12]
    temp_par = (trapz(t1,c(0,temp_par[T_start])) + trapz(t2,temp_par[T_start:12])) * 2
    
    # Scale direct and diffuse values by a factor based on daily satellite PAR
    # and daily "bird" surface irradiance.
    scale_factor = PARsat/temp_par
    DIR = DIR * scale_factor
    DIF = DIF * scale_factor
    
    
    # TIME LOOP PART 2 ####
    # Dawn to noon, using updated direct and diffuse
    for (T_idx in 1:12) {

        #***********************************************************************
        
        ZEN = ZEN_hr[T_idx]
        ZEND = ZEN * (180/pi)
        if (ZEND < 0 | ZEND > 80) {next}
        
        # CALCULATE REFLECTION LOSSES AT AIR-SEA INTERFACE
        I0_refl <- I0_corr_refl_loss(DIR[,T_idx],DIF[,T_idx],ZEN)
        I0_direct <- I0_refl$I0_direct
        I0_diffuse <- I0_refl$I0_diffuse
        ZENW <- I0_refl$ZENW
        
        conversion_factor = LAMBDA * 36 / (19.87 * 6.022 * 10^7)

        # GET SCALAR IRRADIANCE ABOVE WATER, in Einsteins/m^2/hr/nm
        I0_direct = I0_direct * conversion_factor
        I0_diffuse = I0_diffuse * conversion_factor
        
        
        #***********************************************************************
        # SUBSURFACE IRRADIANCE, PHOTIC DEPTH ####
        
        sampenguin = sam_penguin(chl=chlz,
                                 I0_direct=I0_direct,
                                 I0_diffuse=I0_diffuse,
                                 alphaB=alphaB,
                                 Z=Z,
                                 ZENW=ZENW,
                                 LAMBDA=LAMBDA)
        PARz = sampenguin$PARz            # PAR at each depth
        XXz = sampenguin$XXz              # irradiance at each depth
        EuphoD[T_idx] = sampenguin$EuphoD # photic depth
        id = sampenguin$id                # index of photic depth
        yelsubz = sampenguin$yelsubz      # yellow substances
        
        # Surface PAR for this time index (T_idx)
        PAR0[T_idx] = PARz[1]

        # Primary production at each depth for this time index (T_idx)
        Pz = chlz * XXz/sqrt(1 + (XXz/P_mB)^2)

        if(EuphoD[T_idx] > BotZ) {
            # Set photic depth equal to bathymetry at this pixel.
            id = BotZ * 2 + 1 # since we're doing every half metre depth
            EuphoD[T_idx] = BotZ
        }
        
        
        #***********************************************************************
        # INTEGRATION OVER DEPTH ####
        # Integrating primary production, PAR, XX (spectral irradiance) down to
        # photic depth for this time index
        
        PP[T_idx] = trapz(Z[1:id], Pz[1:id])
        PAR[T_idx] = trapz(Z[1:id], PARz[1:id])
        XX[T_idx] = trapz(Z[1:id], XXz[1:id])
        yelsub[T_idx] = trapz(Z[1:id], yelsubz[1:id])
        
    }
    
    
    if (T_start==0) {
        # All ZEND > 80, irradiance can't be computed with Bird
        return(c("PP_BIO"=NA,
                 "Chlz_BIO"=trapz(Z,chlz),
                 "PAR_BIO"=NA,
                 "PAR0_BIO"=NA,
                 "XX_BIO"=NA,
                 "EuphoD_BIO"=NA,
                 "yelsub_BIO"=NA))
    }
    
    
    #***************************************************************************
    # INTEGRATION OVER TIME ####
    # Integrating PP, PAR, PAR0, XX (irradiance)
    
    # Integrate from dawn to "start time" given by the time at index T_start,
    # where T_start is the first time index where ZEND < 80 degrees. All values
    # are assumed to = 0 at dawn. Then integrate over the remaining time
    # intervals until noon, and add it to the first number.
    # Multiply answer by 2 since we have only integrated over half the day.
    time_interval = seq(Dawn + Delta_T, 12, by=Delta_T)
    t1 = c(Dawn,time_interval[T_start])
    t2 = time_interval[T_start:12]
    PP_BIO   = (trapz(t1,c(0,PP[T_start]))   + trapz(t2,PP[T_start:12]))   * 2
    PAR_BIO  = (trapz(t1,c(0,PAR[T_start]))  + trapz(t2,PAR[T_start:12]))  * 2
    PAR0_BIO = (trapz(t1,c(0,PAR0[T_start])) + trapz(t2,PAR0[T_start:12])) * 2
    XX_BIO   = (trapz(t1,c(0,XX[T_start]))   + trapz(t2,XX[T_start:12]))   * 2
    EuphoD_BIO   = median(EuphoD, na.rm=TRUE)
    yelsub_BIO   = (trapz(t1,c(0,yelsub[T_start]))   + trapz(t2,yelsub[T_start:12]))   * 2
    
    # DAILY VALUES
    return(c("PP_BIO"=PP_BIO,
             "Chlz_BIO"=trapz(Z,chlz),
             "PAR_BIO"=PAR_BIO,
             "PAR0_BIO"=PAR0_BIO,
             "XX_BIO"=XX_BIO,
             "EuphoD_BIO"=EuphoD_BIO,
             "yelsub_BIO"=yelsub_BIO))
    
}



#*******************************************************************************
# This function is identical to dwcpn() above, but it uses a C++ version of the
# bird() and sam_penguin() functions for speed.
# Note: bird() calculates direct and diffuse components of surface irradiance for the selected wavelengths.
#       sam_penguin() calculates underwater irradiance down to the euphotic depth.

dwcpn_C <- function(input_df, uniform_profile=TRUE, bp=NULL, gauss=TRUE, cloud_corr=TRUE) {
    
    #***************************************************************************
    # DEFINE VARIABLES
    
    LAMBDA = seq(400,700,by=5) # Vector of wavelengths
    dz = 0.5                   # Depth interval between computations (metres)
    Zmin = 0                   # Min depth
    Zmax = 250                 # Max depth
    Z = seq(Zmin,Zmax,by=dz)   # Vector of depths
    Rsize = length(LAMBDA)     # Number of wavelengths (400 to 700 nm, 5nm increments)
    nstep = length(Z)          # Length of depths vector, where dz=0.5m and Zmax=250m
    
    input_df <- as.numeric(input_df)
    doy = input_df[1]              # Day of year
    long = input_df[2]             # Longitude
    lat = input_df[3]              # Latitude
    BotZ = input_df[4]             # Bathymetry
    PARsat = input_df[6]           # Satellite PAR for scaling bird PAR
    if (cloud_corr) {
        Cloud = input_df[7]        # Fraction of cloud cover (between 0 and 1)
    }
    alphaB = input_df[8]           # Initial slope of PI curve
    P_mB = input_df[9]             # Assimilation number (maximum photosynthetic rate)
    
    
    #***************************************************************************
    # CHLA PROFILE ####
    if (uniform_profile) {
        chl = input_df[5]
        chlz = rep(chl,length(Z))
    } else {
        # Different methods for creating chla depth profile
        #   1. shifted gaussian formula
        #       (if parameters missing, return NA)
        #   2. profile given by user (interpolation to missing points)
        #       (if bp matrix missing, return NA)
        
        if (gauss) { # use shifted gaussian
            z_m = input_df[10]
            B_0 = input_df[11]
            h = input_df[12]
            sigma = input_df[13]
            chlz = shifted_gaussian(tv=Z, B0=B_0, h=h, sigma=sigma, tmax=z_m)
            exp_term = 0.5 * ((Z - z_m)/sigma)^2
            chlz[exp_term > 675] = B_0 #Avoid values too low?
        } else { # use interpolation
            chlz <- approx(bp[,1], bp[,2], xout=Z, method="linear",rule=2)$y
        }
        
        if (all(is.na(chlz) | (chlz < 0))) {
            return(c("PP_BIO"=NA,
                     "Chlz_BIO"=NA,
                     "PAR_BIO"=NA,
                     "PAR0_BIO"=NA,
                     "XX_BIO"=NA,
                     "EuphoD_BIO"=NA,
                     "yelsub_BIO"=NA))
        }
    }
    
    
    #***************************************************************************
    # TIME AND ZENITH ANGLE STUFF ####
    
    # Compute variables based on day and latitude, used in calculation of time of
    # dawn and solar zenith angle.
    Day = doy-1
    dl_vars = daylat(day=Day,lat=lat) # Jan1 = Day 0
    delta = dl_vars$delta
    phi = dl_vars$phi
    phidel = dl_vars$phidel
    
    # Check for 24-hour darkness or sunlight
    if (phidel < -1) { # 24hr sunlight
        phidel = -1
    } else if (phidel > 1) { # 24hr darkness
        return(c("PP_BIO"=0,
                 "Chlz_BIO"=trapz(Z,chlz),
                 "PAR_BIO"=0,
                 "PAR0_BIO"=0,
                 "XX_BIO"=0,
                 "EuphoD_BIO"=0,
                 "yelsub_BIO"=0))
    }
    
    # Compute dawn
    Dawn = 12 - acos(phidel)*(180/pi)/15
    T_start = 0                             #Calculation start index (i.e. first index where ZEND < 80)
    Delta_T = (12-Dawn)/12                  #Compute time step
    T_seq = seq(Dawn+Delta_T,12,by=Delta_T) #Sequence of time intervals
    
    
    #***************************************************************************
    # TIME LOOP PART 1 ####
    # Dawn to noon, doubled at the end
    
    EuphoD = rep(0,12)      #Euphotic depth
    PAR0 = rep(0,12)        #Surface solar irradiance
    # Variables below are integrated over depth within the time loop:
    PP = rep(0,12)          #Carbon primary production, nonuniform profile
    PAR = rep(0,12)         #Photosynthetically active radiation (irradiance integrated over wavelength at each depth)
    XX = rep(0,12)          #Irradiance used in PP formula
    I0 = rep(0,12)          #Surface irradiance
    yelsub = rep(0,12)      #yellow substances in the water column (for testing)
    
    ZEN_hr = rep(0,12)
    DIR = matrix(0,nrow=Rsize,ncol=12)
    DIF = matrix(0,nrow=Rsize,ncol=12)
    
    
    for (T_idx in 1:12) {
        
        # Get next time interval
        Time = T_seq[T_idx]
        
        #***********************************************************************
        # ZENITH ANGLE ####
        
        # Compute zenith angle of sun in radians and convert to degrees.
        ZEN = zenith_bio(Time, long, delta, phi)
        ZEND = ZEN * (180/pi)
        
        ZEN_hr[T_idx] = ZEN
        
        # Restrict to [ 0 < zenith < 80 ] because bird_revised can't accept ZEND > 80
        # NOTE: in Fortran script, looks like ZEND > 80 was converted to ZEND = 79 and used
        if (ZEND < 0 | ZEND > 80) {next} # Skip to next time step
        
        # T_start is T_idx value when ZenD becomes < 80 degrees
        if (T_start==0) {T_start = T_idx}
        
        #***********************************************************************
        # SURFACE IRRADIANCE ####
        
        # Compute clear sky direct and diffuse components of spectral irradiance
        # at sea level, for a given zenith angle (zenith < 80 degrees).
        # UNITS: watts * metre^-2 * micrometre^-1
        DIR_DIF = bird_C(ZEN=ZEN,wl=seq(0.4,0.7,by=0.005))
        I0_direct = DIR_DIF[[1]]
        I0_diffuse = DIR_DIF[[2]]
        
        # CORRECT FOR SEASONAL VARIATION IN SOLAR ENERGY
        I0_season <- I0_corr_season(I0_direct,I0_diffuse,Day)
        I0_direct <- I0_season$I0_direct
        I0_diffuse <- I0_season$I0_diffuse
        
        if (cloud_corr) {
            # CLOUD COVER CORRECTION
            I0_cloud <- I0_corr_cloud_BIO(I0_direct,I0_diffuse,ZEN,Cloud)
            I0_direct <- I0_cloud$I0_direct
            I0_diffuse <- I0_cloud$I0_diffuse
        }
        
        DIR[,T_idx] = I0_direct
        DIF[,T_idx] = I0_diffuse
        
        
        # If you don't want to correct by satellite PAR, comment from this star
        # line down to the star line at the start of time loop part 2 (below)
        #***********************************************************************
        
    }
    
    
    #***************************************************************************
    # ADJUST I0_DIRECT AND I0_DIFFUSE BASED ON DAILY SATELLITE PAR ####
    
    temp_par = (DIR + DIF)
    
    # Convert units of DIR and DIF: watts/m^2/micron --> Einsteins/m^2/hr/nm
    # (Integrating over wavelength removes nm, and then integrating over day
    # gives the full PAR per day, same as NASA)
    # MAKE SURE YOU'RE USING FILES DIRECTLY FROM NASA THAT HAVE PAR IN UNITS OF
    # EINSTEIN/M^2/DAY
    # Similar conversion steps are outlined in the time loop below.
    unit_scale = LAMBDA * 36 / (6.6261 * 2.9979 * 6.022 * 10^7)
    temp_par = temp_par * unit_scale
    temp_lam = LAMBDA
    
    # Integrate over wavelength
    temp_par = as.numeric(sapply(1:ncol(temp_par),function(i) {trapz(temp_lam,temp_par[,i])}))
    
    # Integrate over time using the method described near the bottom of this
    # script, which is dependent on the first time interval where ZEN > 80
    # degrees, given by the index T_start.
    time_interval = seq(Dawn + Delta_T, 12, by=Delta_T)
    t1 = c(Dawn,time_interval[T_start])
    t2 = time_interval[T_start:12]
    temp_par = (trapz(t1,c(0,temp_par[T_start])) + trapz(t2,temp_par[T_start:12])) * 2
    
    # Scale direct and diffuse values by a factor based on daily satellite PAR
    # and daily "bird" surface irradiance.
    scale_factor = PARsat/temp_par
    DIR = DIR * scale_factor
    DIF = DIF * scale_factor
    
    
    # TIME LOOP PART 2 ####
    # Dawn to noon, using updated direct and diffuse
    for (T_idx in 1:12) {
        
        #***********************************************************************
        
        ZEN = ZEN_hr[T_idx]
        ZEND = ZEN * (180/pi)
        if (ZEND < 0 | ZEND > 80) {next}
        
        # CALCULATE REFLECTION LOSSES AT AIR-SEA INTERFACE
        I0_refl <- I0_corr_refl_loss(DIR[,T_idx],DIF[,T_idx],ZEN)
        I0_direct <- I0_refl$I0_direct
        I0_diffuse <- I0_refl$I0_diffuse
        ZENW <- I0_refl$ZENW
        
        conversion_factor = LAMBDA * 36 / (19.87 * 6.022 * 10^7)
        
        # GET SCALAR IRRADIANCE ABOVE WATER, in Einsteins/m^2/hr/nm
        I0_direct = I0_direct * conversion_factor
        I0_diffuse = I0_diffuse * conversion_factor
        
        
        #***********************************************************************
        # SUBSURFACE IRRADIANCE, PHOTIC DEPTH ####
        
        sampenguin = sam_penguin_C(chl=chlz,
                                 I0_direct=I0_direct,
                                 I0_diffuse=I0_diffuse,
                                 alphaB=alphaB,
                                 Z=Z,
                                 ZENW=ZENW,
                                 LAMBDA=LAMBDA)
        PARz = sampenguin$PARz            # PAR at each depth
        XXz = sampenguin$XXz              # irradiance at each depth
        EuphoD[T_idx] = sampenguin$EuphoD # photic depth
        id = sampenguin$id                # index of photic depth
        yelsubz = sampenguin$yelsubz      # yellow substances
        
        # Surface PAR for this time index (T_idx)
        PAR0[T_idx] = PARz[1]
        
        # Primary production at each depth for this time index (T_idx)
        Pz = chlz * XXz/sqrt(1 + (XXz/P_mB)^2)
        
        if(EuphoD[T_idx] > BotZ) {
            # Set photic depth equal to bathymetry at this pixel.
            id = BotZ * 2 + 1 # since we're doing every half metre depth
            EuphoD[T_idx] = BotZ
        }
        
        
        #***********************************************************************
        # INTEGRATION OVER DEPTH ####
        # Integrating primary production, PAR, XX (spectral irradiance) down to
        # photic depth for this time index
        
        PP[T_idx] = trapz(Z[1:id],Pz[1:id])
        PAR[T_idx] = trapz(Z[1:id],PARz[1:id])
        XX[T_idx] = trapz(Z[1:id],XXz[1:id])
        yelsub[T_idx] = trapz(Z[1:id],yelsubz[1:id])
        
    }
    
    
    if (T_start==0) {
        # All ZEND > 80, irradiance can't be computed with Bird
        return(c("PP_BIO"=NA,
                 "Chlz_BIO"=trapz(Z,chlz),
                 "PAR_BIO"=NA,
                 "PAR0_BIO"=NA,
                 "XX_BIO"=NA,
                 "EuphoD_BIO"=NA,
                 "yelsub_BIO"=NA))
    }
    
    
    #***************************************************************************
    # INTEGRATION OVER TIME ####
    # Integrating PP, PAR, PAR0, XX (irradiance)
    
    # Integrate from dawn to "start time" given by the time at index T_start,
    # where T_start is the first time index where ZEND < 80 degrees. All values
    # are assumed to = 0 at dawn. Then integrate over the remaining time
    # intervals until noon, and add it to the first number.
    # Multiply answer by 2 since we have only integrated over half the day.
    time_interval = seq(Dawn + Delta_T, 12, by=Delta_T)
    t1 = c(Dawn,time_interval[T_start])
    t2 = time_interval[T_start:12]
    PP_BIO   = (trapz(t1,c(0,PP[T_start]))   + trapz(t2,PP[T_start:12]))   * 2
    PAR_BIO  = (trapz(t1,c(0,PAR[T_start]))  + trapz(t2,PAR[T_start:12]))  * 2
    PAR0_BIO = (trapz(t1,c(0,PAR0[T_start])) + trapz(t2,PAR0[T_start:12])) * 2
    XX_BIO   = (trapz(t1,c(0,XX[T_start]))   + trapz(t2,XX[T_start:12]))   * 2
    EuphoD_BIO   = median(EuphoD, na.rm=TRUE)
    yelsub_BIO   = (trapz(t1,c(0,yelsub[T_start]))   + trapz(t2,yelsub[T_start:12]))   * 2
    
    # DAILY VALUES
    return(c("PP_BIO"=PP_BIO,
             "Chlz_BIO"=trapz(Z,chlz),
             "PAR_BIO"=PAR_BIO,
             "PAR0_BIO"=PAR0_BIO,
             "XX_BIO"=XX_BIO,
             "EuphoD_BIO"=EuphoD_BIO,
             "yelsub_BIO"=yelsub_BIO))
    
}


#*******************************************************************************

# (See top of function "dwcpn" for details on the input variables)

# EXAMPLE ####

# # Using in situ data and constant PI and biomass profile parameters:
# uniform_profile <- TRUE # constant biomass profile
# cloud_corr <- FALSE # no cloud correction
# input_df <- data.frame(doy = c(182, 182, 182), # day of year
#                        longitude = c(-63.99128, -63.93519, -63.87911),
#                        latitude = c(42.02083, 42.02083, 42.02083),
#                        bathymetry = c(2532, 2518, 2496), # see note
#                        chla = c(0.31711, 0.31392, 0.31575),
#                        par = c(49.48521, 50.13367, 50.45139),
#                        cloud = rep(NA, 3), # dummy column for cloud fraction since it's not being used
#                        alphaB = rep(1.2, 3),
#                        PmB = rep(1.69, 3),
#                        stringsAsFactors=FALSE)
# # RUN THE CODE
# result <- apply(X=input_df,
#                 MARGIN=1,
#                 FUN=dwcpn,
#                 uniform_profile=uniform_profile,
#                 cloud_corr=cloud_corr)
# print(t(result))


#*******************************************************************************
# NOTES ####

# Steps to get conversion factor in the "TIME LOOP PART 2" section of the code:
#
# CONVERT I0 UNITS FROM WATTS/M^2/MICRON TO EINSTEINS/M^2/HR/NM
# (FOR USE IN "PAR", NOT IRRADIANCE XX - it's converted back to the right
# units in the sam_penguin function before using it in the formula)

# https://en.wikipedia.org/wiki/Einstein_(unit)
# "...Photosynthetically active radiation (PAR) was formerly often
# reported in microeinsteins per second per square meter."
#
# Einstein = 1 mole (6.022 * 10^23) of photons
# Micron = micrometre = 10^-6 m
# Watt = 1 Joule/second = 1 kg * m^2 / s^3
#
# I0_direct and I0_diffuse units (micron = micrometre):
#       watts / metre^2 / micron
#
#
# STEPS TO CONVERT UNITS:
#
# 1. Convert watts/m^2/micron to watts/m^2/nm
#       I0_direct = I0_direct / (10^3)
#
# 2. Convert watts/m^2/nm to photons/m^2/s/nm
#    Source: Ocean Optics Web Book: Light and Radiometry: Geometrical Radiometry
#    To convert spectral total scalar irradiance (I0_direct) to number of
#    photons per second:
#       I0_direct = I0_direct * LAMBDA * 10^(-9) / (h * c)
#       h = 6.6261 * 10^(-34) (UNITS: J * s)
#       c = 2.9979 * 10^8 (UNITS: m / s)
#    Note that lambda is currently in units of nanometres, so a factor of
#    10^(-9) is included in the equation above to convert it to metres,
#    as c is also in metres.
#
# 3. Convert photons/m^2/s/nm to einsteins/m^2/s/nm
#    einstein = 6.022 * 10^23 photons
#       I0_direct = I0_direct / (6.022 * 10^23)
#
# 4. Convert einsteins/m^2/s/nm to einsteins/m^2/hr/nm
#       I0_direct = I0_direct * (36 * 10^2)
#
# ALL STEPS COMBINED - multiply I0_direct and I0_diffuse each by:
#   LAMBDA * (36 / (6.6261 * 2.9979 * 6.022)) * (10^(-9) * 10^2 / ((10^(-34)) * (10^8) * (10^3) * (10^23)))
# = LAMBDA * 36 / (6.6261 * 2.9979 * 6.022 * 10^7)
# where LAMBDA is in nanometres
