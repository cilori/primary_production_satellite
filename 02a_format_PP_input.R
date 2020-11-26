library(dplyr)
library(data.table)
library(oceancolouR)
library(ncdf4)
library(stringr)
source("02b_get_param.R") # to get profile and PI curve parameters


#===============================================================================
# DESCRIPTION
#===============================================================================
#
# Create standard input for later use in BIO or Laval primary production scripts.
#
# Output files (stored in ../data/formatted/):
#       1. sat_data_YYYYMM_formatted_interval_created[created date].csv
#               Columns: bin, lon, lat, bathy, chl, sst, par, tcc, yelsub
#       2. parameters_BP_YYYYMM_formatted_inteval_created[created date].csv
#               Columns: alphab, pmb
#       3. parameters_PI_YYYYMM_formatted_inteval_created[created date].csv
#               Columns: zm, B0, h, sigma
# 
# bin, lon, lat, bathy = bin number, latitude, longitude, bathymetry
# chl = (satellite) chlorophyll-a
# sst = (satellite) sea surface temperature
# par = (satellite) photosynthetically active radiation
# tcc = (satellite) total cloud cover
# alphab = (PI curve parameter) initial slope normalized to biomass concentration
# pmb = (PI curve parameter) assimiliation number, i.e. maximum photosynthetic rate, normalized to biomass concentration
# zm, B0, h, sigma = chlorophyll-a profile parameters


#*******************************************************************************
# VARIABLES TO CHANGE

# NOTES:
#   Input file boundaries:  PANCAN (for CHL/PAR/SST), and NWA (for TCC)
#   Output file boundaries: NWA

# Years to process (numeric vector)
y_list <- 2018

# Months to process (numeric vector)
m_list <- 7

# Options: 8day or monthly (i.e. get data for weekly PP processing, or monthly?)
interval <- "monthly"

# VIIRS-SNPP or MODIS
sensor <- "VIIRS-SNPP"

# base input path to chl, par, and sst files (not including sensor/variable/region/year subfolders)
# input_path <- "/mnt/data3/claysa/"      # using pancan dataset
input_path <- "01c_other_variables/"    # using your own dataset

# number of clusters for parallel processing
num_cl <- 10


#*******************************************************************************

input_path <- paste0(input_path, sensor, "/")

# function to read and format data for daily CHL, PAR, SST netCDF files
get_data <- function(file, variable) {
    nc <- nc_open(file)
    dat <- ncvar_get(nc, variable)
    nc_close(nc)
    return(matrix(dat, ncol=1))
}

data("nwa_lats_4km")
data("nwa_lons_4km")
data("nwa_bins_4km")
data("nwa_bath_4km")
data("pancan_bins_4km")

# (nn_pre.f90: issue warnings for extreme values, but still use them.  A large
# number of warnings might indicate a data format problem.)

# Reference parameters for nonuniform biomass profile
bp_par <- read.csv("02a_LUT_cp-20141117.csv") %>%
    mutate(dates=as.Date(paste0(year,"-",mon,"-",day))) %>%
    filter(chlsur>=0,!is.na(chlsur),!is.na(tempsur))#,
           #between(sigma,0,250),between(zm,-100,200),between(rho,0,1))

# Reference parameters for PI curve
pi_par <- read.csv("02a_LUT_pi-20141113.csv") %>%
    mutate(dates = as.Date(paste0(year,"-",mon,"-",day)),
           dn = format(as.Date(dates,origin="1970-01-01"),"%j")) %>%
    filter(alpha >= 0, !is.na(pm), !is.na(alpha))#,
           #between(alpha,0,0.6),between(PmB,0,15))

# Path for output files
output_path <- get_dir(paste0("02_PP_input/", interval, "/"))

# End of the output name (after the date which is automatically used, before the extension .csv)
output_name_end <- paste0(interval, "_NWA_", sensor, "_created",Sys.Date())

# Make a vector of months and days at 8-day intervals for weekly images to extract month-day combinations.
# NOTE ABOUT FORMATS:
#   VIIRS CHL/SST/PAR: YYYYDDD, 8 day intervals, not accounting for leap years
#   MODIS TCC: YYYY-MM-DD, not accounting for leap years
mvec <- as.numeric(sapply(1:46,function(i) format(as.Date((8*0:45)[i],origin=paste0("2001-01-01")),"%m")))
dvec <- sapply(1:46,function(i) format(as.Date((8*0:45)[i],origin=paste0("2001-01-01")),"%d"))
jvec <- str_pad(8*(0:45)+1, width=3, side="left", pad="0")


#-----------------------------
# Loop through years
for (y in y_list) {
    
    all_chl_files <- list.files(paste0(input_path, "CHL_OCX/PANCAN/", y))
    all_par_files <- list.files(paste0(input_path, "PAR/PANCAN/", y))
    all_sst_files <- list.files(paste0(input_path, "SST/PANCAN/", y))
    
    if (length(all_chl_files)==0 | length(all_par_files)==0 | length(all_sst_files)==0) {next}
    
    chl_doys <- as.numeric(sapply(all_chl_files, substr, start=6, stop=8))
    par_doys <- as.numeric(sapply(all_par_files, substr, start=6, stop=8))
    sst_doys <- as.numeric(sapply(all_sst_files, substr, start=6, stop=8))
    
    #-----------------------------
    # Loop through months
    for (m in m_list) {
        
        # Get a list of days of the year for this month
        if (interval=="monthly") {
            dvec_sub <- "01"
            jvec_sub <- jvec[mvec %in% m][1]
        } else if (interval=="8day") {
            dvec_sub <- dvec[mvec %in% m]
            jvec_sub <- jvec[mvec %in% m]
        }
        
        #-----------------------------
        # Loop through indices (1 index if monthly, 3-4 if 8day)
        for (id in 1:length(dvec_sub)) {
            
            day <- dvec_sub[id]
            doy <- jvec_sub[id]
            
            # get tcc filename and check if it exists
            # NOTE: tcc filenames use YYYYMMDD, but are not affected by leap years
            tcc_file <- ifelse(interval=="8day",
                               paste0("01b_formatted_TCC/", interval, "/MYDAL2_E_CLD_FR_",y,"-",str_pad(m,width=2,side="left",pad="0"),"-",day,".csv"),
                               paste0("01b_formatted_TCC/", interval, "/MYDAL2_M_CLD_FR_",y,"-",str_pad(m,width=2,side="left",pad="0"),".csv"))
            
            if (!file.exists(tcc_file)) {next}
            
            # get chl, par, and sst filenames and check if they exist
            if (interval=="8day") {
                doy_vec <- days_vector(year=y, week=which(jvec==doy))
            } else if (interval=="monthly") {
                doy_vec <- days_vector(year=y, month=m)
            }
            
            # pick the indices of the files for these days
            chl_file_inds <- which(chl_doys %in% doy_vec)
            par_file_inds <- which(par_doys %in% doy_vec)
            sst_file_inds <- which(sst_doys %in% doy_vec)
            
            # any daily input files missing?
            if (length(chl_file_inds) != length(doy_vec) | length(par_file_inds) != length(doy_vec) | length(sst_file_inds) != length(doy_vec)) {next}
            
            # grab the filenames at these indices
            chl_files <- all_chl_files[chl_file_inds]
            par_files <- all_par_files[par_file_inds]
            sst_files <- all_sst_files[sst_file_inds]
            
            # GET CHL, PAR, SST DATA (by writing daily data to matrix and condensing daily data into 8day or monthly)
            chl_mat <- get_data(paste0(input_path, "CHL_OCX/PANCAN/", y, "/", chl_files[1]), "chlor_a")
            par_mat <- get_data(paste0(input_path, "PAR/PANCAN/", y, "/", par_files[1]), "par")
            sst_mat <- get_data(paste0(input_path, "SST/PANCAN/", y, "/", sst_files[1]), "sst")
            for (d in 2:length(doy_vec)) {
                chl_mat <- cbind(chl_mat, get_data(paste0(input_path, "CHL_OCX/PANCAN/", y, "/", chl_files[d]), "chlor_a"))
                par_mat <- cbind(par_mat, get_data(paste0(input_path, "PAR/PANCAN/", y, "/", par_files[d]), "par"))
                sst_mat <- cbind(sst_mat, get_data(paste0(input_path, "SST/PANCAN/", y, "/", sst_files[d]), "sst"))
            }
            chl <- rowMeans(chl_mat, na.rm=TRUE)
            par <- rowMeans(par_mat, na.rm=TRUE)
            sst <- rowMeans(sst_mat, na.rm=TRUE)
            
            # Merge data by bin number and add TCC to it
            sat_vars <- left_join(x=data.frame(bin=nwa_bins_4km,
                                               lon=nwa_lons_4km,
                                               lat=nwa_lats_4km,
                                               bathy=nwa_bath_4km,
                                               stringsAsFactors = FALSE),
                                  y=data.frame(bin=pancan_bins_4km, chl=chl, sst=sst, par=par, stringsAsFactors = FALSE),
                                  by="bin")
            sat_vars <- left_join(x=sat_vars, y=fread(tcc_file), by="bin") %>%
                dplyr::mutate(bathy=abs(bathy)) %>%
                dplyr::filter(is.finite(chl) & is.finite(sst) & is.finite(par) & is.finite(tcc))# & between(tcc,0,1) & between(chl,0.0095,64.0005))
            
            
            #-------------------------------------------------------------------
            # GET PARAMETERS - "valid_chl" SHOULD BE LOGGED CHLA
            
            # This is using the current method coded in R, which is slightly
            # different from the one in the original complex Fortran BIO code.
            # How? We may never know
            
            cat("Computing biomass profile and PI curve parameters...\n")
            
            valid_sst <- sat_vars$sst
            valid_chl <- log(sat_vars$chl)
            
            # Biomass profile
            bp_param <- get_param(logchla=valid_chl,
                                  sst=valid_sst,
                                  month=m,
                                  day=as.numeric(doy),
                                  ref_tab=bp_par,
                                  param_type="bp",
                                  num_cl=num_cl)
            
            # PI curve
            pi_param <- get_param(logchla=valid_chl,
                                  sst=valid_sst,
                                  month=m,
                                  day=as.numeric(doy),
                                  ref_tab=pi_par,
                                  param_type="pi",
                                  num_cl=num_cl)
            
            
            #-------------------------------------------------------------------
            # FORMAT AND WRITE OUTPUT TO .CSV FILE
            
            datestr <- ifelse(interval=="8day", paste0(y, doy), paste0(y, str_pad(m,width=2,side="left",pad="0")))
            
            cat(paste0("Writing output for ", datestr, "...\n"))
            
            # Write satellite variables, biomass profile parameters, and PI curve parameters
            fwrite(sat_vars, file = paste0(output_path, datestr, "_satdata_", output_name_end, ".csv"))
            fwrite(bp_param[,c(3,4,1,2)], file = paste0(output_path, datestr, "_BPparams_", output_name_end, ".csv"))
            fwrite(pi_param[,c(2,1)], file = paste0(output_path, datestr, "_PIparams_", output_name_end, ".csv"))
            
        }
        
    }

}

