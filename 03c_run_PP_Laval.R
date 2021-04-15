cat("\014")     # Clear console
rm(list=ls())   # Clear global environment



####################
# 2021-04-15 WARNING:
# THIS HAS NOT BEEN THOROUGHLY TESTED IN SEVERAL MONTHS (beyond minimal testing
# today to make sure it runs through without R errors)
# Previously is was being tested with the BIO v1 model, and might still contain
# temporary blocks of code in this script or model_Laval.R that were used to try
# to calculate and compare the same variables between models. Use with caution
####################



# Stephanie.Clay@dfo-mpo.gc.ca
# 10 JAN 2019

# This script uses input created by format_input.R to compute primary production
# using the Laval model. For references and description, see model_Laval.R

# INPUT:
#   csv files containing satellite data (and, optionally, parameters) formatted
#   using format_input.R (one for each composite, either 8day or monthly).
#   NOTES:
#       CLOUD COVER MUST BE BETWEEN 0 AND 1
#       NO NA VALUES

# OUTPUT:
# csv files containing output where:
#           rows = pixels in the selected region,
#           columns = variables to output (PP, and any others)
# netcdf files containing the same variables as csv, but gridded to view in map form
#   Variables in NetCDF are listed in this order:
#       Primary production
#       Total integrated chla in water column
#       PAR (photosynthetically active radiation in the water column)
#       PAR0 (photosynthetically active radiation at the surface)
#       PUR (photosynthetically usable radiation in the water column)
#       Ek (irradiance at the point of maximum photosynthetic rate)
#       I0 (irradiance at the)


#*******************************************************************************
# TROUBLESHOOTING AND TIPS

# Parallel processing progress bar coding tips:
# https://stackoverflow.com/questions/40684721/how-to-show-the-progress-of-code-in-parallel-computation-in-r

# If you get this error, check the input files for the selected year, month,
# and composite (8day or monthly) actually exist, in the "input" subdirectory.
#   Error in fread(paste0(input_path, file_main)) : 
#       File 'input/' is a directory. Not yet implemented.


#*******************************************************************************
# LIBRARIES, CUSTOM FUNCTIONS

library(ncdf4)
library(raster)
library(parallel)
library(pbapply)
library(data.table)
library(dplyr)
library(oceancolouR)
library(stringr)

source("PP_LAVAL/model_Laval.R")
source("PP_LAVAL/Laval_zenith.R")
source("PP_LAVAL/Laval_daylight.R")
source("PP_LAVAL/Laval_irradiance_functions.R")


#*******************************************************************************
# VARIABLES THAT CAN BE CHANGED

# Years to process (numeric vector)
y_list <- c(2018)

# Months to process (numeric vector)
m_list <- 7

# Length of composites used for PP calculation (8day or monthly)
interval <- "monthly"

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
use_constant_alphab <- TRUE
use_constant_pmb <- TRUE

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
lat_bounds <- c(52,56)
lon_bounds <- c(-60,-55)

# Number of clusters to use in parallel processing (pretty much necessary when
# processing huge rasters of data).
num_cl <- 12


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
            
            vars_to_select <- c("lon", "lat", "chl", "tcc", "alphab", "pmb")
            
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
            
            # Create clusters for parallel processing.
            cl <- makeCluster(min(detectCores()-1, num_cl))
            # Load necessary variables and libraries into cluster.
            clusterExport(cl, c("doy","y","thetas","awat","ays","alp","bb","daylat",
                                "alpbet","Photoperiod_time","zenith_laval",
                                "Ed0moins","input_subset","FNday","FNrange",
                                "f0","f1","FNsun","Isurface_laval","Isubsurface_laval",
                                "uniform_profile","gauss","bp","I0_corr_cloud_LAVAL"))
            clusterEvalQ(cl, c(library(caTools), library(stringr)))
            # Run parallel processing on input_subset to compute results.
            result <- pbapply(input_subset,
                              MARGIN=1,
                              FUN=LAVAL_PP,
                              doy=doy,
                              year=y,
                              lambda=lambda,
                              thetas=thetas,
                              awat=awat,
                              uniform_profile=uniform_profile,
                              gauss=gauss,
                              cl=cl)
            # Stop parallel processing and return processing power to other operations.
            stopCluster(cl)
            
            
            # # FOR TESTING INDIVIDUAL PIXELS INSTEAD (troubleshooting)
            # result <- list()
            # for (pp_idx in 1:nrow(input_subset)) {
            #     result[[pp_idx]] <- LAVAL_PP(input_df=input_subset[pp_idx,],
            #                                  doy=doy,
            #                                  year=y,
            #                                  lambda=lambda,
            #                                  thetas=thetas,
            #                                  awat=awat,
            #                                  uniform_profile=uniform_profile,
            #                                  bp=bp,
            #                                  gauss=gauss)
            # }
            
            result <- result %>% t() %>% as.data.frame
            
            
            #*******************************************************************
            # WRITE OUTPUT TO FILES
            
            output_name <- paste0(output_path, "LAVAL_NWA_", interval, "_", input_datestrs[di], "_", nc_out)
            if (subset) {
                output_name <- paste0(output_name, "_", paste0(lat_bounds,collapse="-"),"N_", paste0(lon_bounds,collapse="-"),"W")
            }
            cat(paste0("Writing results to ",output_name,"\n"))
            
            
            # format output, remove invalid data, place on NWA grid
            output <- dplyr::left_join(data.frame(bin=nwa_bins_4km, stringsAsFactors = FALSE),
                                       dplyr::bind_cols(input, result) %>%
                                           dplyr::select(bin, PP_LAVAL) %>%
                                           dplyr::filter(PP_LAVAL >= 0),
                                       by="bin") %>%
                dplyr::select(PP_LAVAL)
            
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
                                    longname="Primary Production using Laval model")
            # create new output netcdf
            ncout <- nc_create(paste0(output_name, ".nc"), output_var, force_v4=TRUE)
            # put variable in file
            ncvar_put(ncout, output_var, vals=output$PP_LAVAL)
            # close file
            nc_close(ncout)
            
            # write the rows with valid data (in the subsetted region, if applicable),
            # and all the columns (all input and output) to Rdata file
            rda_output <- dplyr::bind_cols(input, result)
            save(rda_output, file=paste0(output_name, ".rda"), compress=TRUE)
            
            
        } # list of files for this month
        
    } # list of months
    
} # list of years
