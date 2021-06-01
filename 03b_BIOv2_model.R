
#*******************************************************************************
# DATA ####

# We define lambda for the visible spectrum (i.e. PAR)
lambda = 400:700

# Absorption:
# Pure seawater is taken from Pope and Fry 1997
awp <- read.table("data/pope97.dat",skip=6) # absorption of pure sea water according to Pope and Fry (1997)
res <- approx(awp[,1],awp[,2],lambda) # asw is absorption by pure seawater at wavelengths of interest
asw <- res$y


#*******************************************************************************
# R FUNCTION ####

# # test data (get remaining variables from PAR_resolved():
# chlpix = 8
# alphaB = 0.033 * 1000 # those have to be defined before in mg C (mg chl m-3)-1 s-1
# PBm = 2.5 # Those have to be defined before

# chlpix = numeric vector of chl in water column
pp_BIO_v2 <- function(chlpix, alphaB, PBm, Eqdifw, Eqdirw, zendR, zendw, asw) {

    dz = 0.5                   # Depth interval between computations (metres)
    Zmin = 0.0                 # Min depth, we have computed the light just under the water, so the next
    # computation occurs at 0.5m 
    Zmax = 120                 # Max depth
    Z = seq(Zmin,Zmax,by=dz)   # Vector of depths
    Rsize = length(lambda)     # Number of wavelengths (400 to 700 nm, 5nm increments)
    nstep = length(Z)          # Length of depths vector, where dz=0.5m and Zmax=250m
    xhr = 1:23                 # hours
    
    
    #***********************
    # fix chl profile
    if (length(chlpix)==1) {
        chlpix <- rep(chlpix, nstep)
    }
    # # tmax = depth at max chla
    # chlpix <- shifted_gaussian(tv=depths, B0=B0, beta=0, h=h, sigma=sigma, tmax=tmax)
    
    # Here is the table where we store the light at each depth (i.e., 50cm), hour and lambda
    E0dir = array(NaN,c(23,nstep,Rsize))
    E0dif = array(NaN,c(23,nstep,Rsize))
    # PP = matrix(0,23,nstep)

    #Here we just initialise the light at the top of the water column from the atmospheric model
    E0dif[,1,] = Eqdifw
    E0dir[,1,] = Eqdirw
    
    
    ###############################
    # BIO-OPTICS INFORMATION ######
    ###############################
    
    # Some information on "constant" IOPs such as absorption and backscattering from pure seawater
    # Seawater backscattering 
    BW500 = 0.00288
    BW = BW500*(lambda/500)^(-4.3)
    
    # for phytoplankton absorption we used the model from Devred et al. 2006 !!! THIS NEED TO BE UPDATED
    # AC = pc1*(1-exp(-rate*chlz) + pc2*chlz with:
    # ATLANTIC ZONE
    lamc = seq(400,700,5)
    pc1 = c(0.0744, 0.0773, 0.0778, 0.0776, 0.0773, 0.0745, 0.0746, 0.0748, 0.0723, 0.0668, 0.0612, 0.0571, 0.0546, 0.0516, 0.0488, 0.0446, 0.0401, 0.036, 0.0322, 0.0283, 0.0247, 0.0217, 0.019, 0.017, 0.0151, 0.0134, 0.012, 0.0106, 0.0093, 0.0081, 0.007, 0.0061, 0.0056, 0.0056, 0.0064, 0.007, 0.0089, 0.0104, 0.0103, 0.0118, 0.0105, 0.0102, 0.0113, 0.0112, 0.0116, 0.0119, 0.0128, 0.0135, 0.0116, 0.0103, 0.0104, 0.012, 0.0168, 0.0224, 0.0277, 0.0296, 0.0269, 0.0207, 0.0135, 0.0082, 0.0072)
    pc2 = c(0.0104, 0.0105, 0.0112, 0.0115, 0.0114, 0.0117, 0.0119, 0.0124, 0.0125, 0.0121, 0.0114, 0.0109, 0.0106, 0.0103, 0.0096, 0.009, 0.0084, 0.0081, 0.0079, 0.0078, 0.0077, 0.0073, 0.007, 0.0066, 0.0063, 0.0061, 0.0058, 0.0056, 0.0053, 0.005, 0.0046, 0.0042, 0.0037, 0.0034, 0.0031, 0.003, 0.0029, 0.0028, 0.0028, 0.0025, 0.0024, 0.0025, 0.0027, 0.003, 0.0033, 0.0034, 0.0035, 0.0035, 0.0036, 0.0036, 0.0037, 0.0043, 0.0056, 0.0074, 0.0089, 0.0092, 0.0081, 0.0061, 0.0039, 0.0023, 0.0011)
    rate = c(0.5582, 0.5706, 0.611, 0.6465, 0.6668, 0.7038, 0.7335, 0.7657, 0.7975, 0.8253, 0.8571, 0.8993, 0.9557, 1.0389, 1.0929, 1.1765, 1.259, 1.3455, 1.3954, 1.4338, 1.4485, 1.4396, 1.434, 1.415, 1.4161, 1.399, 1.365, 1.3307, 1.3254, 1.293, 1.2783, 1.2116, 1.0704, 0.8557, 0.6616, 0.5854, 0.47, 0.429, 0.4673, 0.4162, 0.4494, 0.4422, 0.3926, 0.4014, 0.4103, 0.4147, 0.4125, 0.4179, 0.5266, 0.6462, 0.6428, 0.5651, 0.4828, 0.5287, 0.5984, 0.6433, 0.6482, 0.5811, 0.4876, 0.4428, 0.3376)
    pc1f = approx(lamc,pc1,lambda)$y
    pc2f = approx(lamc,pc2,lambda)$y
    ratef = approx(lamc,rate,lambda)$y
    
    # the slope of the exponential decrease is from Bricaud et al. 1981 Limnology and Oceanography,
    # !!!! this also needs to be udpated !!!!
    Sy = 0.014
    # fraction of ays443 to aph443 = 0.35 !!!!!! This has to be checked as well !!!!!!!!!!!!!
    fracys = 0.35
    
    day_hrs <- which(zendR <= 90)
    
    PP <- underwater_irradiance(E0dir, E0dif, day_hrs, nstep, chlpix, lambda, pc1f, pc2f, ratef, fracys, Sy, asw, zendw, alphaB, PBm, BW)
    
    # /2 because doing measurements every half metre
    finalPP <- sum(PP)/2
    
    return(list(xhr=xhr, Z=Z, PPgrid=PP, PP=finalPP))
    
}


underwater_irradiance <- function(E0dir, E0dif, day_hrs, nstep, chlpix, lambda, pc1f, pc2f, ratef, fracys, Sy, asw, zendw, alphaB, PBm, BW) {
    
    PP <- matrix(0,23,nstep)
    
    for (it in day_hrs) # loop on the time of day (for daylight zendR)
    {
        
        for (iz in 2:nstep) # loop on depth, we start at the second step as the first one is Eqdw
        {
            
            chlz = chlpix[iz]
            
            # To compute the underwater light field, we need to know total absorption and total backscattering
            # a = aw + aphy + acdom, bb = bbw + bbp, this as to be define for each wavelength
            
            # backscattering is modelled according to Loisel and Morel 1998
            bbtilda = (0.78+0.42*(-log10(chlz)))*0.01 # bbtilda et le rapport entre bb et b.
            if (bbtilda < 0.0005) {bbtilda = 0.0005
            } else if (bbtilda > 0.01) {bbtilda = 0.01}
            BC660 = 0.407*chlz^0.795 # ici, bc corresponds to b_chl a 660 nm (scattering by chlorophyll-a)
            BC = BC660*(660/lambda)^(-log10(chlz)) # BC spectral dependence depends on chlz,
            # this means that at chl = 1 tere is no spectral dependence.This is weird and has to be checked
            BC[BC < 0] = 0
            bbt = BC*bbtilda + BW*0.5 # Here, we get total backscattering.
            
            #phytoplankton absorption ### !!!!! THIS HAS TO BE UPDATED.
            aphy = pc1f*(1-exp(-ratef*chlz)) + pc2f*chlz
            
            # Yellow substances and detritus absorption is derived from chl
            ays = fracys * aphy[44] * exp(-Sy * (lambda -443))
            
            # Now we can compute total absorption
            at = asw + ays + aphy
            
            # Here we can start to compute attenuation at each depth for each wavelength
            Kddir = (at + bbt) / cos(zendw[it])
            Kddif = (at + bbt) / 0.83
            
            E0dir[it,iz,] = E0dir[it,(iz-1),] * exp(-Kddir*0.5)
            E0dif[it,iz,] = E0dif[it,(iz-1),] * exp(-Kddif*0.5)
            
            # Here we compute the photosynthetic action spectrum according to Kyewangala et al. 1997
            apbar = mean(aphy)
            alphaBc = alphaB * aphy / apbar
            
            PIz = 1/cos(zendw[it]) * sum(alphaBc *  E0dir[it,iz,]) + 1.20 * sum(alphaBc * E0dif[it,iz,])
            PP[it,iz] = chlz * PIz/sqrt(1 + (PIz/PBm)^2)
            
            # Here we need to compute PP at 0 as the loop starts at 2:
            if (iz == 2)
            {
                PIz = 1/cos(zendw[it]) * sum(alphaBc *  E0dir[it,1,]) + 1.20 * sum(alphaBc * E0dif[it,1,])
                PP[it,1] = chlz * PIz/sqrt(1 + (PIz/PBm)^2)
            }
            
        } # end loop on depth
        
    } # end loop on time of day (xhr)
    
    return(PP)
    
}



#*******************************************************************************
# C++ FUNCTION ####

pp_BIO_v2_c <- function(chlpix, alphaB, PBm, Eqdifw, Eqdirw, zendR, zendw, asw) {
    
    dz = 0.5                   # Depth interval between computations (metres)
    Zmin = 0.0                 # Min depth, we have computed the light just under the water, so the next
    # computation occurs at 0.5m 
    Zmax = 120                 # Max depth
    Z = seq(Zmin,Zmax,by=dz)   # Vector of depths
    Rsize = length(lambda)     # Number of wavelengths (400 to 700 nm, 5nm increments)
    nstep = length(Z)          # Length of depths vector, where dz=0.5m and Zmax=250m
    xhr = 1:23                 # hours
    
    
    #***********************
    # fix chl profile
    if (length(chlpix)==1) {
        chlpix <- rep(chlpix, nstep)
    }
    # # tmax = depth at max chla
    # chlpix <- shifted_gaussian(tv=depths, B0=B0, beta=0, h=h, sigma=sigma, tmax=tmax)
    
    uniform_profile = length(unique(chlpix))==1
    
    
    ###############################
    # BIO-OPTICS INFORMATION ######
    ###############################
    
    # Some information on "constant" IOPs such as absorption and backscattering from pure seawater
    # Seawater backscattering 
    BW500 = 0.00288
    BW = BW500*(lambda/500)^(-4.3)
    
    # for phytoplankton absorption we used the model from Devred et al. 2006 !!! THIS NEED TO BE UPDATED
    # AC = pc1*(1-exp(-rate*chlz) + pc2*chlz with:
    # ATLANTIC ZONE
    lamc = seq(400,700,5)
    pc1 = c(0.0744, 0.0773, 0.0778, 0.0776, 0.0773, 0.0745, 0.0746, 0.0748, 0.0723, 0.0668, 0.0612, 0.0571, 0.0546, 0.0516, 0.0488, 0.0446, 0.0401, 0.036, 0.0322, 0.0283, 0.0247, 0.0217, 0.019, 0.017, 0.0151, 0.0134, 0.012, 0.0106, 0.0093, 0.0081, 0.007, 0.0061, 0.0056, 0.0056, 0.0064, 0.007, 0.0089, 0.0104, 0.0103, 0.0118, 0.0105, 0.0102, 0.0113, 0.0112, 0.0116, 0.0119, 0.0128, 0.0135, 0.0116, 0.0103, 0.0104, 0.012, 0.0168, 0.0224, 0.0277, 0.0296, 0.0269, 0.0207, 0.0135, 0.0082, 0.0072)
    pc2 = c(0.0104, 0.0105, 0.0112, 0.0115, 0.0114, 0.0117, 0.0119, 0.0124, 0.0125, 0.0121, 0.0114, 0.0109, 0.0106, 0.0103, 0.0096, 0.009, 0.0084, 0.0081, 0.0079, 0.0078, 0.0077, 0.0073, 0.007, 0.0066, 0.0063, 0.0061, 0.0058, 0.0056, 0.0053, 0.005, 0.0046, 0.0042, 0.0037, 0.0034, 0.0031, 0.003, 0.0029, 0.0028, 0.0028, 0.0025, 0.0024, 0.0025, 0.0027, 0.003, 0.0033, 0.0034, 0.0035, 0.0035, 0.0036, 0.0036, 0.0037, 0.0043, 0.0056, 0.0074, 0.0089, 0.0092, 0.0081, 0.0061, 0.0039, 0.0023, 0.0011)
    rate = c(0.5582, 0.5706, 0.611, 0.6465, 0.6668, 0.7038, 0.7335, 0.7657, 0.7975, 0.8253, 0.8571, 0.8993, 0.9557, 1.0389, 1.0929, 1.1765, 1.259, 1.3455, 1.3954, 1.4338, 1.4485, 1.4396, 1.434, 1.415, 1.4161, 1.399, 1.365, 1.3307, 1.3254, 1.293, 1.2783, 1.2116, 1.0704, 0.8557, 0.6616, 0.5854, 0.47, 0.429, 0.4673, 0.4162, 0.4494, 0.4422, 0.3926, 0.4014, 0.4103, 0.4147, 0.4125, 0.4179, 0.5266, 0.6462, 0.6428, 0.5651, 0.4828, 0.5287, 0.5984, 0.6433, 0.6482, 0.5811, 0.4876, 0.4428, 0.3376)
    
    # this uses a C++ function, which is generally veeeeeeeeeery slightly faster
    # (but the speed increase adds up when processing millions of pixels)
    pc1f = higher_res(lamc,lambda,pc1)
    pc2f = higher_res(lamc,lambda,pc2)
    ratef = higher_res(lamc,lambda,rate)
    
    # the slope of the exponential decrease is from Bricaud et al. 1981 Limnology and Oceanography,
    # !!!! this also needs to be udpated !!!!
    Sy = 0.014
    # fraction of ays443 to aph443 = 0.35 !!!!!! This has to be checked as well !!!!!!!!!!!!!
    fracys = 0.35
    
    day_hrs <- which(zendR <= 90)
    
    if (uniform_profile) {
        PP <- underwater_irradiance_c(Eqdirw, Eqdifw, chlpix, lambda, pc1f, pc2f, ratef, fracys, Sy, asw, zendR, zendw, alphaB, PBm, BW, min(day_hrs)-1, max(day_hrs)-1)
    } else {
        PP <- underwater_irradiance_NU_c(Eqdirw, Eqdifw, chlpix, lambda, pc1f, pc2f, ratef, fracys, Sy, asw, zendR, zendw, alphaB, PBm, BW, min(day_hrs)-1, max(day_hrs)-1)
    }
    
    # /2 because doing measurements every half metre
    finalPP <- sum(PP)/2
    
    return(list(xhr=xhr, Z=Z, PPgrid=PP, PP=finalPP))
    
}

