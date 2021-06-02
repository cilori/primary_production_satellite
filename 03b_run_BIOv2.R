library(Rcpp)
library(RcppArmadillo)
source("03b_BIOv2_model.R")
source("R/PAR_resolved.R")
sourceCpp("underwater_irradiance.cpp")
sourceCpp("surface_irradiance.cpp")
# library(oceancolouR) # for shifted_gaussian (chl profile)


#*******************************************************************************

# PRIMARY PRODUCTION BIO MODEL VERSION 2

# Step 1: calculate surface irradiance (use Gregg-Carder model)
# Step 2: get chl profile
# Step 3: calculate irradiance in the water column
# Step 4: calculate PP
# Step 5: integrate



# COMMENTS

# In 03b_run_BIOv2.R
    # Underwater field as a function on Chl 
    # For now, homogeneous profile, we can make it more complicated later

# In PAR_resolved.R
    # PERHAPS WE SHOULD KEEP DIRECT AND DIFFUSE ALL ALONG

# In 03b_BIOv2_model.R
    # CHECK/UPDATE:
    # - for phytoplankton absorption we used the model from Devred et al. 2006
    #       AC = pc1*(1-exp(-rate*chlz) + pc2*chlz with:
    # - the slope of the exponential decrease is from Bricaud et al. 1981 Limnology and Oceanography,
    #       Sy = 0.014
    # - fraction of ays443 to aph443 = 0.35
    #       fracys = 0.35
    # - phytoplankton absorption
    #       aphy = pc1f*(1-exp(-ratef*chlz)) + pc2f*chlz
    # - BC = BC660*(660/lambda)^(-log10(chlz)) # BC spectral dependence depends on chlz,
    #   # this means that at chl = 1 tere is no spectral dependence.This is weird and has to be checked


#*******************************************************************************

# VARIABLES FOR SURFACE PAR
latipxl = 50.
daypxl = 173
yearpxl = 2018
SatPAR = 24

# VARIABLES FOR UNDERWATER IRRADIANCE AND PP
chlpix = 8
alphaB = 0.033 * 1000 # those have to be defined before in mg C (mg chl m-3)-1 s-1
PBm = 2.5 # Those have to be defined before

uniform_profile = TRUE



#*******************************************************************************
# GET SURFACE PAR


ptm <- Sys.time()
presolved <- PAR_resolved(latipxl = latipxl,
                          daypxl = daypxl,
                          yearpxl = yearpxl,
                          SatPAR = SatPAR)
Eqdifw <- presolved$Eqdifw
Eqdirw <- presolved$Eqdirw
# Eqdw <- presolved$Eqdw
zendR <- presolved$zendR
zendw <- presolved$zendw
print(Sys.time() - ptm)


#*******************************************************************************
# GET UNDERWATER IRRADIANCE AND PP


ptm <- Sys.time()
if (uniform_profile) {
    result <- pp_BIO_v2_c(chlpix = chlpix,
                          alphaB = alphaB,
                          PBm = PBm,
                          Eqdifw = Eqdifw,
                          Eqdirw = Eqdirw,
                          zendR = zendR,
                          zendw = zendw,
                          asw = asw)
} else {
    result <- pp_BIO_v2_NU_c(chlpix = chlpix,
                          alphaB = alphaB,
                          PBm = PBm,
                          Eqdifw = Eqdifw,
                          Eqdirw = Eqdirw,
                          zendR = zendR,
                          zendw = zendw,
                          asw = asw)
}

print(Sys.time() - ptm)
print(result$PP)




stop()



#*******************************************************************************
# PRINT RESULTS

PP <- result$PP
print(PP)

image.plot(result$xhr, result$Z[1:101], result$PPgrid[,101:1],
           xlab = "Time (h)", ylab = "depth (m)",axes = F)
axis(1)
axis(2,at = seq(0,50,5),labels = seq(50,0,-5))
box()


#*******************************************************************************






stop()


# # do some speedy stuff
# library(memoise)
# library(compiler)
# # pp_BIO_v2 <- memoise(cmpfun(pp_BIO_v2))
# # PAR_resolved <- memoise(cmpfun(PAR_resolved))
# # # IT MAKES NO DIFFERENCE


# speed test with model v1 - original model is faster (mostly because of the lower
# spectral resolution - 5nm instead of 1nm), need to speed up the new version
source("03a_BIOv1_model.R")
source("R/bird.R")
source("R/sam_penguin.R") # exp_m600 and yellow_substances are used in sam_penguin
source("R/yellow_substances.R")
source("R/exp_m600.R")
source("R/I0_corr_season.R")
source("R/I0_corr_cloud_BIO.R")
source("R/I0_corr_refl_loss.R")
source("R/daylat.R")
source("R/zenith_bio.R")
library(caTools)
profvis::profvis({
    resultv1 <- dwcpn(input_df=data.frame(doy=daypxl, lon=-60, lat=latipxl, bathy=500, chl=chlpix,
                                          par=40, tcc=0.3, alphab=1, pmb=1.69),
                      uniform_profile=TRUE,
                      gauss=FALSE)
})


