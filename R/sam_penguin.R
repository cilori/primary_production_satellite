# SUBSURFACE IRRADIANCE (sam_penguin)

sam_penguin <- function(chl,I0_direct,I0_diffuse,alphaB,Z,ZENW,LAMBDA) {
    
    # Modified from raman??????????????
    # Output: PAR, irradiance XX, euphotic depth, Z index of euphotic depth
    
    
    # INPUT         DESCRIPTION                         UNITS
    #***************************************************************************
    # chl           chla profile                        mg/m^3
    # I0_direct     direct irradiance at surface        watts/micrometer
    # I0_diffuse    diffuse irradiance at surface       watts/micrometer
    # alphaB        initial slope of PI curve           mgC/mgChl/(W/m^2)/h?????????Platt1991
    # Z             vector of depths                    meters
    # ZENW          solar zenith angle, below water     radians
    # LAMBDA        vector of wavelengths               nanometers
    #
    #***************************************************************************
    #
    # A        = TOTAL ABSORPTION COEFFICIENT
    # AC440    = ABSOLUTE VALUE OF ABSORPTION ASSOCIATED WITH
    #            CHLOROPHYLL, AT 440 NM.
    # AC550    = ABSOLUTE VALUE OF ABSORPTION ASSOCIATED WITH
    #            CHLOROPHYLL, AT 550 NM.
    # AW       = ABSORPTION COEFFICIENT OF PURE SEAWATER
    # AY       = ABSORPTION COEFFICIENT OF YELLOW SUBSTANCES
    # AY440    = ABSORPTION COEFFICIENT OF YELLOW SUBSTANCES,
    #            AT 440 NM
    # BB       = TOTAL BACKSCATTERING COEFFICIENT
    # BC550    = SCATTERING COEFFICIENT ASSOCIATED WITH CHLOROPHYLL,
    #            AT 550 NM.
    # BW500    = SCATTERING COEFFICIENT OF PURE SEA WATER AT 500 NM
    # K        = DIFFUSE ATTENUATION COEFFICIENT FOR DOWNWELLING LIGHT
    #
    #***************************************************************************
    # CONSTANTS
    
    nstep = length(Z)
    dz = Z[2]-Z[1]
    
    # ATLANTIC ZONE ("! Atl Zone - Nov 24, 1999 ")
    pc1 = c(0.0744, 0.0773, 0.0778, 0.0776, 0.0773, 0.0745, 0.0746, 0.0748, 0.0723, 0.0668, 0.0612, 0.0571, 0.0546, 0.0516, 0.0488, 0.0446, 0.0401, 0.036, 0.0322, 0.0283, 0.0247, 0.0217, 0.019, 0.017, 0.0151, 0.0134, 0.012, 0.0106, 0.0093, 0.0081, 0.007, 0.0061, 0.0056, 0.0056, 0.0064, 0.007, 0.0089, 0.0104, 0.0103, 0.0118, 0.0105, 0.0102, 0.0113, 0.0112, 0.0116, 0.0119, 0.0128, 0.0135, 0.0116, 0.0103, 0.0104, 0.012, 0.0168, 0.0224, 0.0277, 0.0296, 0.0269, 0.0207, 0.0135, 0.0082, 0.0072)
    pc2 = c(0.0104, 0.0105, 0.0112, 0.0115, 0.0114, 0.0117, 0.0119, 0.0124, 0.0125, 0.0121, 0.0114, 0.0109, 0.0106, 0.0103, 0.0096, 0.009, 0.0084, 0.0081, 0.0079, 0.0078, 0.0077, 0.0073, 0.007, 0.0066, 0.0063, 0.0061, 0.0058, 0.0056, 0.0053, 0.005, 0.0046, 0.0042, 0.0037, 0.0034, 0.0031, 0.003, 0.0029, 0.0028, 0.0028, 0.0025, 0.0024, 0.0025, 0.0027, 0.003, 0.0033, 0.0034, 0.0035, 0.0035, 0.0036, 0.0036, 0.0037, 0.0043, 0.0056, 0.0074, 0.0089, 0.0092, 0.0081, 0.0061, 0.0039, 0.0023, 0.0011)
    rate = c(0.5582, 0.5706, 0.611, 0.6465, 0.6668, 0.7038, 0.7335, 0.7657, 0.7975, 0.8253, 0.8571, 0.8993, 0.9557, 1.0389, 1.0929, 1.1765, 1.259, 1.3455, 1.3954, 1.4338, 1.4485, 1.4396, 1.434, 1.415, 1.4161, 1.399, 1.365, 1.3307, 1.3254, 1.293, 1.2783, 1.2116, 1.0704, 0.8557, 0.6616, 0.5854, 0.47, 0.429, 0.4673, 0.4162, 0.4494, 0.4422, 0.3926, 0.4014, 0.4103, 0.4147, 0.4125, 0.4179, 0.5266, 0.6462, 0.6428, 0.5651, 0.4828, 0.5287, 0.5984, 0.6433, 0.6482, 0.5811, 0.4876, 0.4428, 0.3376)
    
    # # GLOBAL ("!all cruises-764 samples")
    # pc1 = c(0.047944, 0.051486, 0.054383, 0.056347, 0.056745, 0.057034, 0.058144, 0.05964, 0.058356, 0.054606, 0.050313, 0.048452, 0.048273, 0.048178, 0.047131, 0.044843, 0.04204, 0.03852, 0.035078, 0.031022, 0.026964, 0.023036, 0.019536, 0.016517, 0.014005, 0.01183, 0.009916, 0.008453, 0.007296, 0.00634, 0.00538, 0.00445, 0.003657, 0.003186, 0.002969, 0.003008, 0.003264, 0.003673, 0.004131, 0.004329, 0.00422, 0.003967, 0.003673, 0.003584, 0.003707, 0.003754, 0.004042, 0.004425, 0.005037, 0.005738, 0.005899, 0.005857, 0.006762, 0.00968, 0.013426, 0.015434, 0.01451, 0.009734, 0.00483, 0.002587, 0.00204)
    # pc2 = c(0.020423, 0.020797, 0.02153, 0.021775, 0.021679, 0.021509, 0.021798, 0.022315, 0.022192, 0.021124, 0.019786, 0.018766, 0.017943, 0.017007, 0.015666, 0.014284, 0.01308, 0.012464, 0.012088, 0.011885, 0.011534, 0.011123, 0.010719, 0.0103, 0.009882, 0.009532, 0.009208, 0.008798, 0.00835, 0.007807, 0.007167, 0.006493, 0.005868, 0.005399, 0.005148, 0.005034, 0.005074, 0.00512, 0.005014, 0.004772, 0.004541, 0.004633, 0.005071, 0.005604, 0.005975, 0.006203, 0.00639, 0.006543, 0.006368, 0.005994, 0.00613, 0.007297, 0.00978, 0.013016, 0.015776, 0.016276, 0.014089, 0.010754, 0.007146, 0.004245, 0.00251)
    # rate = c(1.205002, 1.220828, 1.277448, 1.342184, 1.406491, 1.465066, 1.538236, 1.606021, 1.681225, 1.731006, 1.776641, 1.777083, 1.774193, 1.778021, 1.798992, 1.825786, 1.834832, 1.870123, 1.856917, 1.842261, 1.808863, 1.802813, 1.824051, 1.844891, 1.853456, 1.880783, 1.939347, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1.992925, 1.995062, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1.921015, 1.784356, 1.695164, 1.534441, 1.242857)
    
    AW = c(0.00663, 0.00530, 0.00473, 0.00444, 0.00454, 0.00478, 0.00495, 0.00530,
           0.00635, 0.00751, 0.00922, 0.00962, 0.00979, 0.01011, 0.0106, 0.0114,
           0.0127, 0.0136, 0.0150, 0.0173, 0.0204, 0.0256, 0.0325, 0.0396, 0.0409,
           0.0417, 0.0434, 0.0452, 0.0474, 0.0511, 0.0565, 0.0596, 0.0619, 0.0642,
           0.0695, 0.0772, 0.0896, 0.1100, 0.1351, 0.1672, 0.2224, 0.2577, 0.2644,
           0.2678, 0.2755, 0.2834, 0.2916, 0.3012, 0.3108, 0.325, 0.340, 0.371,
           0.410, 0.429, 0.439, 0.448, 0.465, 0.486, 0.516, 0.559, 0.624)
    
    BW500 = 0.00288
    BR488 = 0.00027
    
    BW = BW500*(LAMBDA/500)^(-4.3)
    BBR = 0.5*BR488*(LAMBDA/488)^(-5.3)
    AY <- sapply(-0.014*(LAMBDA-440),exp_m600)
    
    
    #***************************************************************************
    # COMPUTE IRRADIANCE AND MU_D JUST BELOW SURFACE
    
    # Helpful material:
    #
    # Diagram of above water and below water radiance: Kirk 1983, page 40
    #
    # http://www.oceanopticsbook.info/view/light_and_radiometry/geometrical_radiometry
    # Ed = spectral downward plane irradiance
    # Eod = spectral downward scalar irradiance
    #
    # http://www.oceanopticsbook.info/view/overview_of_optical_oceanography/apparent_optical_properties
    # mu_d = mean cosine of downwelling irradiance = Ed/Eod
    #
    # Sathyendranath and Platt 1998
    #
    # Kirk 1983, section 1.3
    
    
    # Get total scalar and planar irradiance Iz
    Iz_scalar = I0_direct + I0_diffuse
    Iz_planar = I0_direct*cos(ZENW) + I0_diffuse # should I0_diffuse be multiplied by 0.831 here?
    mu_d = (I0_direct*cos(ZENW) + I0_diffuse*0.831)/Iz_scalar
    
    PARz = rep(0,nstep) # PAR at each depth
    XXz = rep(0,nstep)  # irradiance at each depth
    yelsubz = rep(0, nstep) # yellow substances at each depth (for testing)
    
    for (i in 1:nstep) {
        
        # Chlorophyll at depth i
        chlz <- chl[i]
        
        #***********************************************************************
        # Absorption and backscattering for depth i
        
        # Devred 2006
        AC = pc1*(1-sapply(-rate*chlz,exp_m600)) + pc2*chlz
        AC440 = AC[9]
        yelsubz[i] = yellow_substances(chlz)
        if (is.na(yelsubz[i])) {
            
            # THIS ONLY HAPPENS IF CHL < 0
            # HOW TO HANDLE THIS ERROR?
            
            next
        }
        AY440 = yelsubz[i] * AC440
        ACmean = mean(AC)
        A = AW + AC + AY440*AY + 2*BBR
        
        # Loisel and Morel 1998
        bbtilda = (0.78+0.42*(-log10(chlz)))*0.01
        if (bbtilda < 0.0005) {bbtilda = 0.0005
        } else if (bbtilda > 0.01) {bbtilda = 0.01}
        BC660 = 0.407*chlz^0.795
        BC = BC660*(660/LAMBDA)^(-log10(chlz))
        BC[BC < 0] = 0
        BB = BC*bbtilda + BW*0.5
        
        
        #***********************************************************************
        # PAR and XX for this depth, new K and Iz for next depth
        
        # Compute PAR for this depth (integration over wavelength)
        PARz[i] = sum(Iz_scalar * 5) # Riemann sum with left sides (used originally)
        #PARz[i] = trapz(LAMBDA, Iz_scalar) # trapezoidal rule
        
        # PAR units:
        #   BIO:    E/m^2/hr
        #   Laval:  uE/m^2/s
        
        # Convert BIO units to Laval units for easier comparison:
        PAR_conversion_factor <- (10^6)/3600
        PARz[i] <- PARz[i] * PAR_conversion_factor
        
        
        # CONVERT EINSTEINS/M^2/HR/NM --> W/M^2 FOR IRRADIANCE XX USED IN PP MODEL
        #
        # http://strang.smhi.se/extraction/units-conversion.html
        # Morel and Smith (1974):
        # 1 W of PAR = 2.77 A 10^18 quanta/second for marine atmospheres above
        # the water surface with Sun altitudes above 22 degrees (solar zenith
        # angle smaller than 68 degrees).
        #
        # 1. EINSTEINS/M^2/HR/NM --> EINSTEINS/M^2/S/NM
        #       Iz = Iz / (36 * 10^2)
        # 2. EINSTEINS/M^2/S/NM --> PHOTONS/M^2/S/NM
        #       Iz = Iz * (6.022 * 10^23)
        # 3. PHOTONS/M^2/S/NM --> W/M^2/NM
        #       Iz = Iz / (2.77 * 10^18)
        # ALL COMBINED: = 6.022 * 10^23 / ((36 * 10^2) * (2.77 * 10^18))
        #               = 6022 / (2.77 * 36)
        fXX_conversion_factor <- 6022 / (2.77 * 36)
        # NOTE: AC/ACmean and 1/mu_d are used to normalize the function?????????????????
        fXX = (alphaB * Iz_planar) * (AC/ACmean) * (1/mu_d) * fXX_conversion_factor # Vectors: AC, Iz, mu_d
        XXz[i] = sum(fXX * 5) # Riemann sum with left sides (used originally)
        # XXz[i] = trapz(LAMBDA,fXX) # trapezoidal rule
        
        # Find K, diffuse attenuation coefficient, for this layer at each
        # wavelength, and use it to update Iz at this depth for each wavelength.
        K = (A+BB)/mu_d
        Iz_planar = Iz_planar*sapply(-K*dz,exp_m600)
        Iz_scalar = Iz_scalar*sapply(-K*dz,exp_m600)
        
        # End calculations if we've reached euphotic depth (when PAR <= 1% surface PAR)
        if ((i > 1) & (PARz[i] < 0.01*PARz[1])) {
            # Return PAR, XX, euphotic depth, index of euphotic depth
            return(list(PARz=PARz, XXz=XXz, EuphoD=Z[i], id=i, yelsubz=yelsubz))
        }
        
    }
    
    return(list(PARz=PARz, XXz=XXz, EuphoD=Z[nstep], id=nstep, yelsubz=yelsubz))
    
}
