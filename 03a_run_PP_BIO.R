cat("\014")     # Clear console
rm(list=ls())   # Clear global environment

# INPUT:
#   csv files containing satellite data, PI curve parameters, and (optionally)
#   biomass profile parameters, formatted using 02a_format_PP_input.R (one for
#   each composite, either 8day or monthly).
#   NOTES:
#       CLOUD COVER MUST BE BETWEEN 0 AND 1
#       NO NA VALUES

# OUTPUT:
#   netcdf containing binned grid(s) of data values, on an NWA (Northwest Atlantic)
#   sized grid (i.e. flat, length 295425, with NA for missing data values -- this
#   makes the file easier to read later by using the pre-defined NWA bin numbers)
# data layers:
#       Primary production
#       Chla, integrated over depth and time
#       PAR (photosynthetically active radiation in the water column, integrated over depth and time)
#       PAR0 (photosynthetically active radiation at the surface, integrated over time)
#       "XX_BIO" spectral irradiance integrated over depth and time
#               fXX <- (alphaB * Iz_scalar) * (AC/ACmean) * (1/mu_d) * fXX_conversion_factor
#               XX_z[i] <- sum(fXX * 5)
#               XX = XX_z integrated over depth down to euphotic depth
#       Euphotic depth, integrated over time
#       Yellow substances, integrated over depth and time


#*******************************************************************************
# TROUBLESHOOTING AND TIPS

# If you get an error like this, check the input files for the selected year, month,
# and composite (8day or monthly) actually exist, in the "input" subdirectory.
#   Error in fread(paste0(input_path, file_main)) : 
#       File [directory_name] is a directory. Not yet implemented.


#*******************************************************************************
# LIBRARIES, CUSTOM FUNCTIONS

library(ncdf4)          # to write outputto netcdf files
library(caTools)        # for integration across wavelength/depth/time using "trapz" function
library(data.table)     # for reading large csv input files with fast "fread" function
library(dplyr)          # for organizing input and output
library(oceancolouR)    # for various custom functions
library(pbapply)        # for progress bars
library(stringr)

source("PP_BIO/bird.R")
source("PP_BIO/sam_penguin.R")
source("PP_BIO/irradiance_corrections.R")
source("PP_BIO/zenith.R")
source("PP_BIO/model_BIO.R")

# For speeding up code:
# library(profvis)    # speed tests
library(compiler)   # to compile functions
library(Rcpp)       # to incorporate C++ in the code
# library(parallel)   # for parallel processing


#*******************************************************************************
# VARIABLES THAT CAN BE CHANGED

# Years to process (numeric vector)
y_list <- 2018

# Months to process (numeric vector)
m_list <- 7

# Length of composites used for PP calculation (8day or monthly)
interval <- "8day"

# Chla profile uniform or not?
uniform_profile <- TRUE

# If nonuniform, use Gaussian curve?
# (Must have properly named and formatted csv files containing parameters, created using format_input.R)
gauss <- TRUE

# If nonuniform and not Gaussian, use this input matrix (values can be changed).
# Values in between these depths will be interpolated.
bp <- cbind(seq(1,100,by=10),   # depths (metres)
            rep(1,10))          # chla

# Use a single value for each PI parameter? (i.e. the same value for every pixel in the image)
use_constant_alphab <- FALSE
use_constant_pmb <- FALSE

# If using a constant alphab or pmb for the whole region, how should it be calculated?
#       OPTIONS: predefined (alphab=1, pmb=1.69), mean, or median
#       Note that if the image is subsetted below, the mean or median will be calculated using only the values inside those boundaries
alphab_method <- "predefined"
pmb_method <- "predefined"

# Output filename format = BIO_NWA_[interval]_[date]_[nc_out].nc
# [date] = YYYYMM or YYYYDDD for monthly or 8day respectively.
nc_out <- "alphab-pmb-variable_uniform-profile"

# Do a subsetted region (for testing)
subset <- TRUE
lat_bounds <- c(47,55)
lon_bounds <- c(-53,-43)

# # Number of clusters to use in parallel processing
# num_cl <- 10 # detectCores()-1, or manually: 10 for hecla, 3 for Windows, 12 for Loki


#*******************************************************************************
# MAIN CODE

# for writing to output
data("nwa_bins_4km")

# For 8day interval
mvec <- as.numeric(sapply(1:46,function(i) format(as.Date((8*0:45)[i],origin=paste0("2001-01-01")),"%m")))
jvec <- str_pad(((8*0:45)+1),width=3,side="left",pad="0")

# Path of input formatted csv files
input_path <- paste0("02_PP_input/", interval, "/")

# Path for output files
output_path <- paste0("03_PP_output/",interval,"/")

files <- list.files(input_path, pattern=".csv")
file_dates <- sapply(1:length(files), function(i) {strsplit(files[i], "_")[[1]][1]})

#-----------------------------
# Loop through years
for (y in y_list) {
    
    #-----------------------------
    # Loop through months
    for (m in m_list) {
        
        # Get a list of formatted dates to match to files to process for this
        # month, based on choice of 8day or monthly composites.
        if (interval=="monthly") {
            # First day of the month
            jvec_sub <- as.numeric(format(as.Date(paste0(y,str_pad(m,width=2,side="left",pad="0"),"01"),format="%Y%m%d"),"%j"))
            input_datestrs <- paste0(y,str_pad(m,width=2,side="left",pad="0"))
        } else if (interval=="8day") {
            # First day of each 8day interval
            jvec_sub <- jvec[mvec %in% m]
            input_datestrs <- paste0(y,jvec_sub)
        }
        
        #-----------------------------
        # Loop through files (1 if monthly, 3-4 if 8day)
        for (di in 1:length(input_datestrs)) {
            
            # Get day of year (first day of month, or first day of 8day interval)
            doy <- as.numeric(jvec_sub[di])
            
            #*******************************************************************
            # GET DATA
            
            file <- files[which(file_dates==input_datestrs[di])]
            
            if (length(file)==0) {next}
            
            vars_to_select <- c("doy", "lon", "lat", "bathy", "chl", "par", "tcc", "alphab", "pmb")
            
            # bin, lon, lat, bathy, chl, sst, par, tcc
            file_main <- file[which(grepl("satdata", file))]
            input <- fread(paste0(input_path, file_main)) %>% as.data.frame()
            
            # add alphab and pmb (PI curve parameter) columns to input
            file_pi <- file[which(grepl("PIparams", file))]
            input <- bind_cols(input, fread(paste0(input_path, file_pi)) %>% as.data.frame())
            
            # if non-uniform profile and gaussian curve, add BP parameter columns to input
            if (!uniform_profile & gauss) {
                file_bp <- file[which(grepl("BPparams", file))]
                input <- bind_cols(input, fread(paste0(input_path, file_bp)) %>% as.data.frame())
                vars_to_select <- c(vars_to_select, "zm", "B0", "h", "sigma")
            }
            
            # subset (better for testing)
            if (subset) {
                input <- input %>%
                    dplyr::filter(between(lat, lat_bounds[1], lat_bounds[2]),
                                  between(lon, lon_bounds[1], lon_bounds[2]))
            }
            
            # if you're using one pair of PI parameters for every pixel, set them
            # all to either the mean value of the region to process, the median value, 
            # or another predefined value (alphab=1 and PmB=1.69)
            if (use_constant_alphab) {
                if (alphab_method=="median") {
                    input$alphab <- median(input$alphab)
                } else if (alphab_method=="mean") {
                    input$alphab <- mean(input$alphab)
                } else if (alphab_method=="predefined") {
                    input$alphab <- 1 # Laval alphab = 1
                }
            }
            if (use_constant_pmb) {
                if (pmb_method=="median") {
                    input$pmb <- median(input$pmb)
                } else if (pmb_method=="mean") {
                    input$pmb <- mean(input$pmb)
                } else if (pmb_method=="predefined") {
                    input$pmb <- 1.69 # Laval pmb = 1.69
                }
            }
            
            # grab only the columns needed for the BIO model
            input_subset <- input %>%
                dplyr::mutate(doy=doy) %>%
                dplyr::select(all_of(vars_to_select)) %>%
                as.data.frame()
            
            
            #*******************************************************************
            # CALCULATE PP FOR EACH VALID INDEX
            
            
            # # OPTION 1: USE FUNCTION COMPILING AND PARALELL PROCESSING
            # # Compile the function
            # cmp_dwcpn <- cmpfun(dwcpn)
            # # Create clusters for parallel processing.
            # cl <- makeCluster(num_cl)
            # # Load necessary variables and libraries into cluster.
            # clusterExport(cl, c(# variables
            #                     "input_subset","uniform_profile","gauss","bp","doy",
            #                     # custom functions
            #                     "cmp_dwcpn","bird","exp_m600","zenith_bio","sam_penguin",
            #                     "daylat","I0_corr_season","I0_corr_cloud_BIO",
            #                     "I0_corr_refl_loss","yellow_substances"))
            # clusterEvalQ(cl, c(library(caTools), library(dplyr)))
            # # Get results using pblapply to include progress bar.
            # result <- pbapply(input_subset[1:500,],
            #                   MARGIN=1,
            #                   FUN=cmp_dwcpn,
            #                   uniform_profile=uniform_profile,
            #                   bp=bp,
            #                   gauss=gauss,
            #                   cl=cl)
            # # Stop parallel processing and return processing power to other operations.
            # stopCluster(cl)
            # result <- t(result) %>% as.data.frame
            
            
            # OPTION 2: USE FUNCTION COMPILING AND RCPP FOR C++ INSTEAD (no parallel processing, faster)
            sourceCpp("PP_BIO/bio_pp_c.cpp")    # source the C++ script
            cmp_dwcpn_C <- cmpfun(dwcpn_C)      # compile the C++ function
            result <- pbapply(input_subset,
                              MARGIN=1,
                              FUN=cmp_dwcpn_C,
                              uniform_profile=uniform_profile,
                              bp=bp,
                              gauss=gauss)
            result <- result %>% t() %>% as.data.frame
            
            
            #*******************************************************************
            # ERROR CHECKING AND SPEED TESTS
            
            # # find the index creating an error
            # for (i in 1:nrow(input_subset)) {
            #     result <- dwcpn(input_df=input_subset[i,],
            #                     uniform_profile=uniform_profile,
            #                     bp=bp,
            #                     gauss=gauss)
            # }

            # # check if compiling the function makes it faster
            # cmp_dwcpn <- cmpfun(dwcpn)
            # speed_test <- profvis({
            #     result1 <- dwcpn(input_df=input_subset[1,],
            #                      uniform_profile=uniform_profile,
            #                      bp=bp,
            #                      gauss=gauss)
            #     result2 <- cmp_dwcpn(input_df=input_subset[1,],
            #                          uniform_profile=uniform_profile,
            #                          bp=bp,
            #                          gauss=gauss)
            # })
            
            # # try using C++ with the Rcpp package for speed
            # sourceCpp("PP_BIO/bio_pp_c.cpp")
            # cmp_dwcpn <- cmpfun(dwcpn)
            # cmp_dwcpn_C <- cmpfun(dwcpn_C)
            # speed_test <- profvis({
            #     result1 <- apply(input_subset[1:100,],
            #                      MARGIN=1,
            #                      FUN=cmp_dwcpn,
            #                      uniform_profile=uniform_profile,
            #                      bp=bp,
            #                      gauss=gauss)
            #     result2 <- apply(input_subset[1:100,],
            #                      MARGIN=1,
            #                      FUN=cmp_dwcpn_C,
            #                      uniform_profile=uniform_profile,
            #                      bp=bp,
            #                      gauss=gauss)
            # })
            # # check to make sure the results are the same -- if so, the expression below should be TRUE
            # # (note that they might not be 100% identical due to floating-point differences in R and C++)
            # all(abs(t(result1)-t(result2)) < 1e-10)
            
            
            #*******************************************************************
            # WRITE OUTPUT TO FILES
            
            output_name <- paste0(output_path, "BIO_NWA_", interval, "_", input_datestrs[di], "_", nc_out)
            if (subset) {
                output_name <- paste0(output_name, "_", paste0(lat_bounds,collapse="-"),"N_", paste0(lon_bounds,collapse="-"),"W")
            }
            cat(paste0("Writing results to ",output_name,"\n"))
            
            # format output, remove invalid data, place on NWA grid
            output <- dplyr::left_join(data.frame(bin=nwa_bins_4km, stringsAsFactors = FALSE),
                                       dplyr::bind_cols(input, result) %>%
                                           dplyr::select(bin, PP_BIO) %>%
                                           dplyr::filter(PP_BIO >= 0),
                                       by="bin") %>%
                      dplyr::select(PP_BIO)
            
            # create dimensions
            dim_bindata <- ncdim_def(name="binDataDim",
                                     units = "",
                                     vals = 1:nrow(output),
                                     unlim = FALSE,
                                     create_dimvar = FALSE)
            dim_time <- ncdim_def(name = "time",
                                  units = "seconds since 1970-01-01 00:00:00.000 +0000",
                                  vals = as.numeric(as.POSIXct(as.Date(doy, origin = paste0(y, "-01-01")))),
                                  unlim = TRUE,
                                  calendar = "standard")
            # create output variable
            output_var <- ncvar_def(name="PP",
                                    units="mg C / m^2 / day",
                                    dim=list(dim_time, dim_bindata),
                                    missval=NA,
                                    longname="Primary Production using BIO model")
            # create new output netcdf
            ncout <- nc_create(paste0(output_name, ".nc"), output_var, force_v4=TRUE)
            # put variable in file
            ncvar_put(ncout, output_var, vals=output$PP_BIO)
            # close file
            nc_close(ncout)
            
            # write the rows with valid data (in the subsetted region, if applicable),
            # and all the columns (all input and output) to Rdata file
            rda_output <- dplyr::bind_cols(input, result)
            save(rda_output, file=paste0(output_name, ".rda"), compress=TRUE)
            
            
        } # list of files for this month
        
    } # list of months
    
} # list of years

