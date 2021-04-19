
# REFERENCES:
#=============================

# DESCRIPTION OF MODEL
# ATBD Belanger and Babin 2011

# Ardyna, Mathieu & Babin, M & Gosselin, Michel & Devred, Emmanuel & Belanger, Simon & Matsuoka, Atsushi & Tremblay, J.-E. (2013). Parameterization of vertical chlorophyll a in the Arctic Ocean: Impact of the subsurface chlorophyll maximum on regional, seasonal, and annual primary production estimates. Biogeosciences. 10. 4383-4404. 10.5194/bg-10-4383-2013. 

# ARRIGO, KEVIN R. & Sullivan, Cornelius. (1994). A high resolution bio-optical model of microalgal growth: Tests using sea-ice algal community time-series data. Limnology and Oceanography - LIMNOL OCEANOGR. 39. 609-631. 10.4319/lo.1994.39.3.0609.

# Huot, Yannick & Babin, M & Bruyant, F. (2013). Photosynthetic parameters in the Beaufort Sea in relation to the phytoplankton community structure. Biogeosciences Discussions. 10. 10.5194/bgd-10-1551-2013. 

# Lee, Zhong-Ping & Du, Ke-Ping & Arnone, Robert. (2005). A Model for the Diffuse Attenuation Coefficient of Downwelling Irradiance. J. Geophys. Res. 110. 10.1029/2004JC002275. 

# Matsuoka, Atsushi & Hill, Victoria & Huot, Yannick & Babin, Marcel & Bricaud, Annick & Matsuoka, Citation. (2011). Seasonal variability in the light absorption properties of western Arctic waters: Parameterization of the individual components of absorption for ocean color applications. Journal of Geophysical Research. 116. 10.1029/2009JC005594. 

# Matsuoka, Atsushi & Huot, Yannick & Shimada, Koji & Saitoh, Sei-Ichi & Babin, Marcel. (2007). Bio-optical characteristics of the Western Arctic Ocean: Implications for ocean color algorithms. Canadian journal of remote sensing. 33. 10.5589/m07-059. 

# Morel, Andre. (1991). Light and marine photosynthesis: A spectral model with geochemical and climatological implications. Progress In Oceanography. 26. 263-306. 10.1016/0079-6611(91)90004-6. 

# Reynolds, A., Rick & Stramski, Dariusz & Greg Mitchell, B. (2001). A chlorophyll-dependent semianalytical model derived from field measurements of absorption and backscattering coefficients within the Southern Ocean. Journal of Geophysical Research. 106. 7125-7138. 10.1029/1999JC000311. 

# Srikanth, Ayyala Somayajula & Devred, Emmanuel & Belanger, Simon & Antoine, David & vellucci, vincenzo & Babin, Marcel. (2018). Evaluation of sea-surface photosynthetically available radiation algorithms under various sky conditions and solar elevations. Applied Optics. 27. 10.1364/AO.57.003088. 

# Wang, Jian & F Cota, Glenn & A Ruble, David. (2005). Absorption and backscattering in the Beaufort and Chukchi Seas. J. Geophys. Res. 110. 10.1029/2002JC001653. 


################################################################################
# EXTRA FUNCTIONS

# Compute absorption of CDOM, NAP
ays <- function(x, lambda=lambda, slope=slope) {
    x * exp(-slope*(lambda - 440))
}

# Compute absorption of phytoplankton
alp <- function (aphyto_z_440, alpbet) {
    alpbet[,1]*aphyto_z_440^alpbet[,2]   
}

# Wang et al 2005, equations 10, 12, 13 & Reynolds et al (2001)
bb <- function(bp_555, lambda=lambda, bw_550) {
    (bw_550+bp_555)*((555/lambda)^(-2.348*log10(bp_555)-4.353))
}


################################################################################

# VARIABLES

# UNITS: ARDYNA ET AL 2013, TABLE 1


#----------------------------------------------------------
# Get table of Ed (spectral irradiance) based on lambda (wavelength), theta (solz,
# solar zenith angle), ozone (O3), and taucl (total cloud fraction). 2 options:

# OPTION 1
# O3=300DU i.e. Dobson Units; 300 is the average value under standard conditions
# taucl=0, clear sky
lut <- read.table("data/LAVAL_Ed0moins_LUT.dat")
lambda <- seq(290,700,5)
thetas <- seq(0,90,5)
O3 <- seq(100,550,50)
taucl <- c(0,1,2,4,8,16,32,64)
# Reduce lut by selected wavelengths (400 to 700nm)
wl_idx <- 23:83
lut <- lut[,wl_idx]
lambda <- lambda[wl_idx]
# For each combination of solar zenith angle, O3, and cloud fraction (taucl), get the LUT value of Ed0moins.
pyramide <- array(NA,c(length(lambda),19,10,8)) # lambda, thetas, O3, taucl
ct <- 1
for (i in 1:19) { # theta loop (sun angle)
    for (j in 1:10) { # O3 loop
        for (k in 1:8) { # taucl loop
            pyramide[,i,j,k] <- as.numeric(lut[ct,])
            ct <- ct+1
        }
    }
}
# O3 index 5 for 300DU, taucl index 1 for 0 cloud fraction
Ed0moins <- t(pyramide[,,5,1])

# # OPTION 2
# # Use Ed0moins at O3=300 and taucl=5 (credit: Simon Lambert?)
# # NOTE: THIS CURRENTLY WILL NOT WORK because the function PP_Zeu_Zpd adjusts the
# # spectral irradiance assuming it was retrieved for a clear sky (Srikanth et al 2018).
# load("Ed0moins_thetas_wv.Rdata") # contains variable Ed0moins, matrix of 19 thetas x 61 wavelengths
#----------------------------------------------------------


# Water absorption
awtable <- read.table(file="data/LAVAL_Absorption_Water.csv", sep=";")
Aw <- spline(awtable$V1, awtable$V2, xout=lambda, method="natural") # interpolate to get values for all necessary wavelengths lambda
awat <- t(cbind(Aw$x,Aw$y)) # create a new table containing all selected wavelengths and their corresponding Aw

# Morel coefficients (2001)
awtable_Morel <- read.table(file="data/LAVAL_Tableau_Water_Morel & Moriterano.csv", sep=";", header=T)
Kw <- awtable_Morel$'Kw' # attenuation coef of water using Morel's method (2001)
E <- awtable_Morel$'e'
X <- awtable_Morel$'X'

# Model parameters for Kd with different solar zenith angles (ATBD Belanger and Babin 2011, equation 4)
Kd_m <- matrix(c(1.044,1.108,1.32,4.173,4.245,4.12,0.53,0.526,0.504,11.157,10.942,10.304),nrow=4,ncol=3,byrow=T)

# Phytoplankton absorption
aphytable <- read.table(file="data/LAVAL_Absorption_phyto.csv", sep=";") # Matsuoka et al. (2011)
alpbet <- cbind(aphytable$V2,aphytable$V3)


################################################################################

# source("R/zenith_laval.R")
# source("R/Photoperiod_time.R") # this uses days_to_J2000, hour_angle, and sun_lon
# source("R/days_to_J2000.R")
# source("R/hour_angle.R")
# source("R/sun_lon.R")
# source("R/Isurface_laval.R")
# source("R/Isubsurface_laval.R")
# source("R/I0_corr_cloud_LAVAL.R")

LAVAL_PP <- function(input_df,doy,year,lambda,thetas,awat,eklaval=FALSE,uniform_profile=TRUE,bp=NULL,gauss=TRUE,cloud_corr=TRUE) {
    
    #---------------------------------------------------------------------------
    # DESCRIPTION
    
    # Get primary production (PP) in units of mgC m^-2 d^-1
    # (d = day)
    # If profile is not uniform, you need either profile parameters, or bp vector.
    
    # MODELE ECRIT PAR MATTHIEU ARDYNA (ARDYNA ET AL 2013),
    # MODIFIE PAR LAURENT OZIEL ET EMMANUEL DEVRED (2016-2017) pour le PROJET IGS
    # Modified by Stephanie Clay, early 2019
    
    
    #---------------------------------------------------------------------------
    # DEFINE VARIABLES
    
    # input_df         = a vector containing the data and parameters:
    #             lon, lat, chl, tcc, alphab, pmb, and (optionally) zm, B0, h, sigma
    # doy       = day of year, single numeric value
    # bp        = a matrix where column1 = depth and column2 = chla at that depth
    # lambda    = vector of wavelengths (400 to 700 at 5nm intervals)
    
    dz <- 1                     # Depth interval between computations (metres)
    Zmin <- 1                   # Min depth
    Zmax <- 100                 # Max depth
    Z <- seq(Zmin,Zmax,by=dz)   # Vector of depths
    
    input_df <- as.numeric(input_df)
    xlon <- input_df[1]
    ylat <- input_df[2]
    if (cloud_corr) {
        TCC <- input_df[4]  # total cloud cover
    }
    alphaB <- input_df[5]   # Initial slope of curve (not used in Laval model, but BIO model for comparison)
    P_mB <- input_df[6]     # Maximum photosynthetic rate
    
    
    #---------------------------------------------------------------------------
    # CHLA PROFILE  (shifted gaussian)
    
    if (uniform_profile) {
        chl <- input_df[3]
        chlz <- rep(chl,length(Z))
    } else {
        # Different methods for creating chla depth profile
        #   1. gaussian formula
        #       (if parameters missing, return NA)
        #   2. profile given by user (interpolation to missing points)
        #       (if bp matrix missing, return NA)
        
        if (gauss) { # use shifted gaussian
            zm <- input_df[7]
            B0 <- input_df[8]
            h <- input_df[9]
            sigma <- input_df[10]
            chlz = shifted_gaussian(tv=Z, B0=B0, h=h, sigma=sigma, tmax=zm) # chla profile vector, length = Z
            exp_term = 0.5 * ((Z - zm)/sigma)^2
            chlz[exp_term > 675] = B0 #Avoid values too low?
        } else { # use interpolation
            chlz <- approx(bp[,1], bp[,2], xout=Z, method="linear",rule=2)$y
        }
        
        # if bad chla values for the entire profile, return NA
        if (all(is.na(chlz) | (chlz < 0))) {
            return(c("PP_LAVAL"=NA,
                     "Chlz_LAVAL"=NA,
                     "PAR_LAVAL"=NA,
                     "PAR0_LAVAL"=NA,
                     "XX_LAVAL"=NA,
                     "PUR_LAVAL"=NA,
                     "Ek_LAVAL"=NA,
                     "Ek_BIO"=NA,
                     "I0_LAVAL"=NA))
        }
        
    }
    
    
    #---------------------------------------------------------------------------
    # TIME AND ZENITH ANGLE STUFF
    
    # Compute variables based on day and latitude, used in calculation of time of
    # dawn and solar zenith angle.
    Day <- doy-1
    dl_vars = daylat(day=Day,lat=ylat) # Jan1 = Day 0
    theta = dl_vars[[1]]
    delta = dl_vars[[2]]
    phi = dl_vars[[3]]
    phidel = dl_vars[[4]]
    
    # Check for 24-hour darkness or sunlight
    if (phidel < -1) { # 24hr sunlight
        phidel = -1
    } else if (phidel > 1) { # 24hr darkness
        return(c("PP_LAVAL"=0,
                 "Chlz_LAVAL"=trapz(Z,chlz),
                 "PAR_LAVAL"=0,
                 "PAR0_LAVAL"=0,
                 "XX_LAVAL"=0,
                 "PUR_LAVAL"=0,
                 "Ek_LAVAL"=0,
                 "Ek_BIO"=0,
                 "I0_LAVAL"=0))
    }
    
    
    # VECTOR OF TIME STEPS
    
    # LAVAL METHOD
    Delta_T <- 3
    T_seq <- seq(0,24,by=Delta_T)
    
    # # BIO METHOD
    # Dawn <- 12 - acos(phidel)*(180/pi)/15   #Compute dawn
    # Delta_T <- (12-Dawn)/12                 #Compute time step from Dawn to Noon
    # T_seq <- seq(Dawn+Delta_T,12,by=Delta_T)#Sequence of time intervals
    
    
    #---------------------------------------------------------------------------
    
    # Collect PUR and PAR at each depth (rows) and hour interval (columns), and
    # euphotic depth at each hour interval
    PUR_z_time <- matrix(0,nrow=length(Z),ncol=length(T_seq))
    PAR_z_time <- matrix(0,nrow=length(Z),ncol=length(T_seq))
    XX_z_time <- matrix(0,nrow=length(Z),ncol=length(T_seq))
    EuphoD <- rep(0,length(T_seq))
    I0_time <- rep(0,length(T_seq))
    
    # TIME LOOP
    for (T_idx in 1:length(T_seq)) {
        
        Time <- T_seq[T_idx]
        
        
        #-----------------------------------------------------------------------
        # ZENITH ANGLE
        
        # # Solar zenith angle calculation.
        # ZEN <- zenith_bio(TIME=Time, LONG=xlon, delta=delta, phi=phi)
        # ZEND <- ZEN * 180/pi
        ZEND <- zenith_laval(doy=doy, hour=Time, xcoor=xlon, ycoor=ylat)$szendeg
        
        
        #-----------------------------------------------------------------------
        # SURFACE IRRADIANCE

        # GET SURFACE IRRADIANCE FOR THIS TIME INTERVAL AND SOLAR ZENITH ANGLE
        # NOTE: The look-up table used for this, Ed0moins, contains input for
        #       solz, lambda, O3, and taucl (cloud fraction). In the main script,
        #       it's subsetted to use only the values for O3=300DU and taucl=0:
        #       Ed0moins <- t(pyramide[,,5,1])

        # ATBD Belanger and Babin 2011:
        # I0 = I0_diffuse * 0.934 + Ed_direct * (1-ref)
        # ref = specular reflection coefficient (function of solar zenith angle
        #       and refraction indices of air and water)
        I0 <- Isurface_laval(LUT=Ed0moins,thetas=thetas,solz=ZEND)

        #---------------------------------------------------------------------------
        # CONVERT I0 UNITS FROM uE/m^2/s ==> E m-2 d-1

        #**************************************************************************************
        # DOUBLE-CHECK THIS
        # (originally I0*3600*24/1000, but to get from uE to E, shouldn't it be 
        #  divided by a factor of 10^6, not 10^3?)
        #**************************************************************************************
        I0 <- I0*3600*24/10^6  # convert units from uE m-2 s-1 ==> E m-2 d-1
        #I0 <- I0*3600/10^6  # convert units from uE m-2 s-1 ==> E m-2 hr-1
        #I0 <- I0/(10^6) # convert units from uE m-2 s-1 ==> E m-2 s-1
        
        
        
        # #***********************************************************************
        # # Using BIO v1 model functions temporarily to compare
        # # Note that surface irradiance must also be scaled according to PAR
        #
        # # bird requires ZEND < 80 - check here
        # if (ZEND > 80) {
        #     PUR_z_time[,T_idx] <- rep(0,length(Z))
        #     PAR_z_time[,T_idx] <- rep(0,length(Z))
        #     XX_z_time[,T_idx] <- rep(0,length(Z))
        #     EuphoD[T_idx] <- 1
        #     I0_time[T_idx] <- 0
        #     next
        # }
        # 
        # # bird surface irradiance test:
        # # NOTE: THIS REQUIRES ZEN IN RADIANS
        # # output irradiance in WATTS/M^2/MICRON
        # DIR_DIF = bird(ZEN=ZEN,wl=seq(0.4,0.7,by=0.005))
        # I0_direct = DIR_DIF[[1]]
        # I0_diffuse = DIR_DIF[[2]]
        # 
        # # CORRECT FOR SEASONAL VARIATION IN SOLAR ENERGY
        # I0_season <- I0_corr_season(I0_direct,I0_diffuse,Day)
        # I0_direct <- I0_season$I0_direct
        # I0_diffuse <- I0_season$I0_diffuse
        # 
        # # CLOUD COVER CORRECTION
        # I0_cloud <- I0_corr_cloud_BIO(I0_direct,I0_diffuse,ZEN,TCC)
        # I0_direct <- I0_cloud$I0_direct
        # I0_diffuse <- I0_cloud$I0_diffuse
        # 
        # # CALCULATE REFLECTION LOSSES AT AIR-SEA INTERFACE
        # I0_refl <- I0_corr_refl_loss(I0_direct,I0_diffuse,ZEN)
        # I0_direct <- I0_refl$I0_direct
        # I0_diffuse <- I0_refl$I0_diffuse
        # 
        # # TO OUTPUT TO NETCDF FOR COMPARISON
        # 
        # # convert to EINSTEINS/M^2/HR/NM (see model_BIO.R for detailed steps)
        # conversion_factor = lambda * 36 / (19.87 * 6.022 * 10^7)
        # # GET SCALAR IRRADIANCE ABOVE WATER, in Einsteins/m^2/hr/nm
        # I0_direct = I0_direct * conversion_factor
        # I0_diffuse = I0_diffuse * conversion_factor
        # I0 <- I0_direct * cos(ZEN) + I0_diffuse
        # I0_time[T_idx] <- I0        
        # 
        # # FOR THE SUBSURFACE IRRADIANCE CALCULATION WITH THE LAVAL MODEL
        # 
        # # bio conversion gives this
        # # EINSTEINS/M^2/HR/NM
        # # for laval, we need uEINSTEINS/M^2/S/NM
        # conversion_factor = 10^6 / 3600
        # I0_direct = I0_direct * conversion_factor
        # I0_diffuse = I0_diffuse * conversion_factor
        # I0 <- I0_direct * cos(ZEN) + I0_diffuse        
        # 
        # #***********************************************************************
        
        
        
        if (cloud_corr) {
            # CLOUD COVER CORRECTION
            I0 <- I0_corr_cloud_LAVAL(I0,TCC)
            # I0 <- I0_corr_cloud_BIO_total(I0,ZEN,TCC)
        }
        
        # Check - if surface irradiance = 0, set PAR and PUR to 0 for this time
        # step, and Euphotic depth to 1 (i.e. surface), then move to next time step.
        if (all(I0==0)) {
            PUR_z_time[,T_idx] <- rep(0,length(Z))
            PAR_z_time[,T_idx] <- rep(0,length(Z))
            XX_z_time[,T_idx] <- rep(0,length(Z))
            EuphoD[T_idx] <- 1
            I0_time[T_idx] <- 0
            next
        }
        
        
        #-----------------------------------------------------------------------
        # SUBSURFACE IRRADIANCE, PHOTIC DEPTH
        
        Isur_sub <- Isubsurface_laval(chlz,ZEND,Z,lambda,I0,Zmax,awat)
        PUR_z_time[,T_idx] <- Isur_sub[[1]]
        PAR_z_time[,T_idx] <- Isur_sub[[2]]
        XX_z_time[,T_idx] <- Isur_sub[[3]]
        EuphoD[T_idx] <- Isur_sub[[4]]
        I0_time[T_idx] <- trapz(lambda,I0)
        
    }
    
    ############################################################################
    # Compute daily mean PUR and Ek for each depth
    
    # Get daylength and convert units to same as used for I0????
    date_obj <- as.Date(doy,origin=paste0(year,"-01-01"))
    m <- as.numeric(format(date_obj,"%m"))
    d <- as.numeric(format(date_obj,"%d"))
    daylen <- Photoperiod_time(latit=ylat,longit=xlon,y=year,m=m,day=d,tzone=-4)$daylen # Photoperiod_time daylen is in hours
    
    # Get mean PUR by integrating PUR across time at each depth, then dividing by
    # the day length (available daylight). ATBD Belanger and Babin 2011, equation 10
    PUR_mean_z <- rep(0,length(Z))
    for (j in 1:length(Z)) {
        PUR_mean_z[j] <- trapz(T_seq,PUR_z_time[j,])
    }
    PUR_mean_z <- PUR_mean_z/daylen
    
    # Original Ek model from Arrigo and Sullivan (1994), equation 13
    # Updated coefficients from Huot et al. (2013), combined and reported in
    # Ardyna et al. (2013), equation 7
    Ek_z_LAVAL <- 25.7/(1+2.2*exp(-0.336*PUR_mean_z))
    
    # FOR TESTING, not necessarily accurate
    # ***NOTE: This is a single value, not a vector of values for each depth like Ek_z_Laval
    Ek_BIO <- P_mB/alphaB
    
    
    ############################################################################
    # Calculate PP for each depth and time interval, then integrate over depth and time.
    # Also do PAR, PUR, and XX (not necessary, just good for comparisons to other models).
    
    PP_time <- rep(0,length(T_seq))
    PAR_time <- rep(0,length(T_seq))
    PUR_time <- rep(0,length(T_seq))
    XX_time <- rep(0,length(T_seq))
    
    for (T_idx in 1:length(T_seq)) {
        
        # Get PP at each depth for this hour (ATBD Belanger and Babin 2011, equation 5)
        # UNITS:
        #       chlz    --> mg Chla / m^3
        #       P_mB    --> mg C / mg Chla / h
        #       PUR     --> umol photons / m^2 / s
        #       Ek      --> umol photons / m^2 / s
        #       PP      --> mg C / m^2 / day
        # (integrate over depth: m^3 ----> m^2
        #  integrate over time: h ----> day)
        
        PAR_z <- PAR_z_time[,T_idx]
        PUR_z <- PUR_z_time[,T_idx]
        XX_z <- XX_z_time[,T_idx]
        
        # units of PUR and Ek cancel, so they aren't affected when you integrate
        # over the full day by hour?
        if (eklaval) {
            PP_z <- (chlz*P_mB)*(1-exp(-PUR_z/Ek_z_LAVAL))
        } else {
            PP_z <- (chlz*P_mB)*(1-exp(-PUR_z/Ek_BIO))
        }
        
        
        # # BIO METHOD FOR TESTING
        # PP_z <- chlz * XX_z/sqrt(1 + (XX_z/P_mB)^2)
        
        
        # Integrate with respect to depth to get total PP in the water column
        # down to the euphotic depth for this hour interval.
        if (all(PP_z==0)) {
            # Check if all PP=0 (if they do, trapz will return NA, so this must be considered separately)
            PP_time[T_idx] <- 0
            PAR_time[T_idx] <- 0
            PUR_time[T_idx] <- 0
            XX_time[T_idx] <- 0
        } else {
            euph_layer <- 1:EuphoD[T_idx]
            PP_time[T_idx] <- trapz(euph_layer,PP_z[euph_layer])
            PAR_time[T_idx] <- trapz(euph_layer,PAR_z[euph_layer])
            PUR_time[T_idx] <- trapz(euph_layer,PUR_z[euph_layer])
            XX_time[T_idx] <- trapz(euph_layer,XX_z[euph_layer])
        }
        
    }
    
    # Integrate with respect to time to get total PP for this day.
    PP_LAVAL <- trapz(T_seq,PP_time)
    PAR_LAVAL <- trapz(T_seq,PAR_time)
    PAR0_LAVAL <- trapz(T_seq,PAR_z_time[1,])
    PUR_LAVAL <- trapz(T_seq,PUR_time)
    XX_LAVAL <- trapz(T_seq,XX_time)
    I0_LAVAL <- trapz(T_seq,I0_time)
    
    
    ############################################################################
    # DAILY VALUES
    
    # MULTIPLY BY 2 IF INTEGRATING ONLY UNTIL NOON (when comparing to BIO v1 model)
    return(c("PP_LAVAL"=PP_LAVAL,#*2,
             "Chlz_LAVAL"=trapz(Z,chlz),#*2,
             "PAR_LAVAL"=PAR_LAVAL,#*2,
             "PAR0_LAVAL"=PAR0_LAVAL,#*2,
             "XX_LAVAL"=XX_LAVAL,#*2,
             "PUR_LAVAL"=PUR_LAVAL,#*2,
             "Ek_LAVAL"=trapz(Z,Ek_z_LAVAL),
             "Ek_BIO"=Ek_BIO,
             "I0_LAVAL"=I0_LAVAL))#*2)))
    
}
