
source("03b_BIOv2_model.R")
source("R/PAR_resolved.R")

#*******************************************************************************
# COMMENTS

# Underwater field as a function on Chl 
# For now, homogeneous profile, we can make it more complicated later


# In PAR_resolved.R
    ##### PERHAPS WE SHOULD KEEP DIRECT AND DIFFUSE ALL ALONG

# In model_BIO_v2.R
    # CHECK/UPDATE:
    # - for phytoplankton absorption we used the model from Devred et al. 2006 !!! THIS NEED TO BE UPDATED
    #       AC = pc1*(1-exp(-rate*chlz) + pc2*chlz with:
    # - the slope of the exponential decrease is from Bricaud et al. 1981 Limnology and Oceanography,
    #       Sy = 0.014
    # - fraction of ays443 to aph443 = 0.35
    #       fracys = 0.35
    # - phytoplankton absorption
    #       aphy = pc1f*(1-exp(-ratef*chlz)) + pc2f*chlz


#*******************************************************************************

# profvis::profvis({
ptm <- Sys.time()

# GET PAR
latipxl = 50.
daypxl = 173
yearpxl = 2018
presolved <- PAR_resolved(latipxl = latipxl,
                          daypxl = daypxl,
                          yearpxl = yearpxl)
Eqdifw <- presolved$Eqdifw
Eqdirw <- presolved$Eqdirw
Eqdw <- presolved$Eqdw

# GET PP
chlpix = 8
alphaB = 0.033 * 1000 # those have to be defined before in mg C (mg chl m-3)-1 s-1
PBm = 2.5 # Those have to be defined before
result <- pp_BIO_v2(chlpix = chlpix,
                    alphaB = alphaB,
                    PBm = PBm,
                    Eqdifw = Eqdifw,
                    Eqdirw = Eqdirw,
                    Eqdw = Eqdw)

# })
print(Sys.time() - ptm)




PP <- result$PP

print(PP)



# AFTER FORMATTING:
# speed check - is it faster to read data and format it from iwthin a function (pope),
# or define it ahead of time and pass it into the function?








# # speed test with model v1 - original model is faster, need to speed up the new version
# source("03a_BIOv1_model.R")
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
# profvis::profvis({
#     resultv1 <- dwcpn(input_df=data.frame(doy=daypxl, lon=-60, lat=latipxl, bathy=500, chl=chlpix,
#                                           par=40, tcc=0.3, alphab=1, pmb=1.69),
#                       uniform_profile=TRUE,
#                       gauss=FALSE)
# })


