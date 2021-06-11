library(Rcpp)
library(RcppArmadillo)
library(oceancolouR)
library(fst)
library(dplyr)
library(lubridate)
library(pbapply)
source("03b_BIOv2_model.R")
source("R/PAR_resolved.R")
sourceCpp("underwater_irradiance.cpp")
sourceCpp("surface_irradiance.cpp")


#*******************************************************************************

# PRIMARY PRODUCTION BIO MODEL VERSION 2

# Step 1: calculate surface irradiance (use Gregg-Carder model)
# Step 2: get chl profile
# Step 3: calculate irradiance in the water column
# Step 4: calculate PP
# Step 5: integrate



#*******************************************************************************
# SET VARIABLES

# FOR TESTING - VARIABLES FOR SURFACE PAR
latipxl = 50.
daypxl = 173
yearpxl = 2018
SatPAR = 24
#  FOR TESTING - VARIABLES FOR UNDERWATER IRRADIANCE AND PP
chlpix = 8
alphaB = 0.033 * 1000 # those have to be defined before in mg C (mg chl m-3)-1 s-1
PBm = 2.5 # Those have to be defined before
uniform_profile = TRUE
# merge into dataframe format
all_data <- data.frame(lat=latipxl,
                       doy=daypxl,
                       par=SatPAR,
                       chl=chlpix,
                       bin=1,
                       stringsAsFactors = FALSE)


# yearpxl = 2015
# uniform_profile = TRUE
# alphaB = 0.033 * 1000 # those have to be defined before in mg C (mg chl m-3)-1 s-1
# PBm = 2.5 # Those have to be defined before
# 
# 
# 
# #*******************************************************************************
# # COLLECT DATA
# 
# cat("Loading data...\n\n")
# 
# # CHL_POLY4 (MODIS, 8day, filled gaps)
# chl_path <- "../data_gaps/filled_files/"
# chl_files <- list.files(chl_path, pattern=".fst")
# chl_file <- chl_files[grepl(yearpxl, chl_files) & grepl("_CHL_OCX_", chl_files)]
# chl_data <- read_fst(file.path(chl_path, chl_file)) %>%
#     dplyr::filter(is.finite(var_imputed))
# 
# # Daily MODIS PAR (to be averaged to 8day)
# data("pancan_bins_4km")
# data("pancan_lons_4km")
# data("pancan_lats_4km")
# par_path <- "/mnt/data3/claysa/MODIS/PAR/PANCAN/annual_fst/"
# # par_path <- chl_path # filled files
# par_files <- list.files(par_path, pattern=".fst")
# par_file <- par_files[grepl(yearpxl, par_files) & grepl("_PAR_", par_files)]
# par_data <- read_fst(file.path(par_path, par_file))
# par_data <- matrix(par_data$var, nrow=num_pix[["PANCAN"]][["4km"]])
# par_data_8day <- avg_columns(par_data, year=yearpxl)
# par_data_8day <- data.frame(bin=rep(pancan_bins_4km,46),
#                             lon=rep(pancan_lons_4km,46),
#                             lat=rep(pancan_lats_4km,46),
#                             week=rep(1:46,each=num_pix[["PANCAN"]][["4km"]]),
#                             par=as.numeric(par_data_8day),
#                             stringsAsFactors = FALSE) %>%
#     dplyr::mutate(time=yearpxl+(week-1)/46,
#                   doy=yday(week8_date(yearpxl,week))) %>%
#     dplyr::filter(is.finite(par))
# 
# all_data <- dplyr::left_join(chl_data, par_data_8day, by=c("bin", "time")) %>%
#     dplyr::rename(chl=var_imputed) %>%
#     dplyr::filter(is.finite(par) & is.finite(chl))
# 
# 
# 
# 
# # with the current LUT for PARday, you can only go to 80deg lat
# # just do a subset for testing
# all_data <- all_data %>% dplyr::filter(doy==97 & lat <= 80)
# 
# 
# 
# 
# metadata <- all_data %>% dplyr::select(bin, time, doy)
# all_data <- all_data %>% dplyr::select(lat, doy, par, chl, bin)
# 
# 
# 
# 
# 
# stop()








#*******************************************************************************
# CALCULATE PP

cat("Calculating primary production...\n\n")


full_PP_calc <- function(df) {
    
    latipxl <- df[1]
    daypxl <- df[2]
    SatPAR <- df[3]
    chlpix <- df[4]
    
    # GET SURFACE PAR
    presolved <- PAR_resolved(latipxl = latipxl,
                              daypxl = daypxl,
                              yearpxl = yearpxl,
                              SatPAR = SatPAR)
    Eqdifw <- presolved$Eqdifw
    Eqdirw <- presolved$Eqdirw
    zendR <- presolved$zendR
    zendw <- presolved$zendw
    
    # GET UNDERWATER IRRADIANCE AND PP
    result <- pp_BIO_v2_c(chlpix = chlpix,
                          alphaB = alphaB,
                          PBm = PBm,
                          Eqdifw = Eqdifw,
                          Eqdirw = Eqdirw,
                          zendR = zendR,
                          zendw = zendw,
                          asw = asw)
    
    return(result$PP)
    
}

full_PP_calc_NU <- function(df) {
    
    latipxl <- df[1]
    daypxl <- df[2]
    SatPAR <- df[3]
    chlpix <- df[4]
    
    # GET SURFACE PAR
    presolved <- PAR_resolved(latipxl = latipxl,
                              daypxl = daypxl,
                              yearpxl = yearpxl,
                              SatPAR = SatPAR)
    Eqdifw <- presolved$Eqdifw
    Eqdirw <- presolved$Eqdirw
    zendR <- presolved$zendR
    zendw <- presolved$zendw
    
    # GET UNDERWATER IRRADIANCE AND PP
    result <- pp_BIO_v2_NU_c(chlpix = chlpix,
                             alphaB = alphaB,
                             PBm = PBm,
                             Eqdifw = Eqdifw,
                             Eqdirw = Eqdirw,
                             zendR = zendR,
                             zendw = zendw,
                             asw = asw)
    
    return(result$PP)
    
}



if (uniform_profile) {
    result <- pbapply(all_data,
                      MARGIN=1,
                      FUN=full_PP_calc)
} else {
    result <- pbapply(all_data,
                      MARGIN=1,
                      FUN=full_PP_calc_NU)
}






print(result)
stop()




all_data <- dplyr::left_join(metadata,
                             all_data %>% dplyr::mutate(PP=as.numeric(result)),
                             by=c("bin", "doy")) %>%
    dplyr::select(bin, time, doy, par, chl, PP)

write_fst(all_data, path=paste0("03_PP_output/8day/BIOv2_NWA_8day_CHL-OCX_weekdoy097_", yearpxl, ".fst"), compress=100)




# PP_final <- rep(NaN, nrow(all_data))
# 
# for (row in 1:1000) {#nrow(all_data)) {
#     
#     latipxl <- all_data[row,"lat"]
#     daypxl <- all_data[row,"doy"]
#     SatPAR <- all_data[row,"par"]
#     chlpix <- all_data[row,"chl"]
#     
#     
#     # GET SURFACE PAR
#     
#     # ptm <- Sys.time()
#     presolved <- PAR_resolved(latipxl = latipxl,
#                               daypxl = daypxl,
#                               yearpxl = yearpxl,
#                               SatPAR = SatPAR)
#     Eqdifw <- presolved$Eqdifw
#     Eqdirw <- presolved$Eqdirw
#     zendR <- presolved$zendR
#     zendw <- presolved$zendw
#     # print(Sys.time() - ptm)
#     
#     
#     # GET UNDERWATER IRRADIANCE AND PP
#     
#     # ptm <- Sys.time()
#     if (uniform_profile) {
#         result <- pp_BIO_v2_c(chlpix = chlpix,
#                               alphaB = alphaB,
#                               PBm = PBm,
#                               Eqdifw = Eqdifw,
#                               Eqdirw = Eqdirw,
#                               zendR = zendR,
#                               zendw = zendw,
#                               asw = asw)
#     } else {
#         result <- pp_BIO_v2_NU_c(chlpix = chlpix,
#                                  alphaB = alphaB,
#                                  PBm = PBm,
#                                  Eqdifw = Eqdifw,
#                                  Eqdirw = Eqdirw,
#                                  zendR = zendR,
#                                  zendw = zendw,
#                                  asw = asw)
#     }
#     
#     # print(Sys.time() - ptm)
#     # print(result$PP)
#     
#     PP_final[row] <- result$PP
#     
# }







# #*******************************************************************************
# # PRINT RESULTS FOR A SINGLE PIXEL
# 
# PP <- result$PP
# print(PP)
# 
# image.plot(result$xhr, result$Z[1:101], result$PPgrid[,101:1],
#            xlab = "Time (h)", ylab = "depth (m)",axes = F)
# axis(1)
# axis(2,at = seq(0,50,5),labels = seq(50,0,-5))
# box()



# # do some speedy stuff
# library(memoise)
# library(compiler)
# # pp_BIO_v2 <- memoise(cmpfun(pp_BIO_v2))
# # PAR_resolved <- memoise(cmpfun(PAR_resolved))
# # # IT MAKES NO DIFFERENCE


# # speed test with model v1 - original model is faster (mostly because of the lower
# # spectral resolution - 5nm instead of 1nm), need to speed up the new version
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


