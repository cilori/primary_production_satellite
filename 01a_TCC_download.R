# Select years and months and the composite length ("interval") below,
# then run this script to get a .txt file containing a list of full TCC
# (total cloud cover) filenames with their addresses for downloading.
# File list is in the format YYYY_get_tcc.txt and is written to the folder
# given by the "download_folder" variable below.

# Then run this line from the Linux command line within the folder
# containing the .txt list file:
# wget --random-wait --content-disposition -i [filename here, including .txt extension]

# TCC files: MODIS-Aqua, global PNG, 1800x3600

# 2019-10-09 tried running the wget line with system() and it was lagging terribly,
# don't know why. Possibly the loop starts the wget command and immediately goes
# to the next loop without waiting for it to finish, and starts a new one, so
# they're running simultaneously.

library(stringr)

years <- 2018
months <- 1:12
interval <- "8day" # 8day or monthly


#===============================================================================
# MAIN CODE
#===============================================================================

download_folder <- paste0("01a_downloaded_TCC/", interval, "/")

# Base url (without the final forward slash)
baseurl="https://neo.sci.gsfc.nasa.gov/archive/gs/"

# Make a vector of months and days at 8-day intervals for weekly images to extract month-day
# combinations. NOTE: MODIS TCC files do not account for leap years in their filenames.
mvec <- as.numeric(sapply(1:46,function(i) format(as.Date((8*0:45)[i],origin=paste0("2001-01-01")),"%m")))
dvec <- as.numeric(sapply(1:46,function(i) format(as.Date((8*0:45)[i],origin=paste0("2001-01-01")),"%d")))

# Make a list of file names and addresses based on selected years, time intervals, and months.
for (y in years) {
    filename_list <- c()
    for (m in months) {
        mname <- str_pad(m,width=2,side="left",pad="0")
        if (interval=="8day") {
            days <- dvec[mvec %in% m]
            days <- sapply(1:length(days),function(i) str_pad(days[i],width=2,side="left",pad="0"))
            datelist <- paste0(y,"-",mname,"-",days)
            filename_list <- rbind(filename_list,cbind(paste0(baseurl,"MYDAL2_E_CLD_FR/MYDAL2_E_CLD_FR_",datelist,".PNG")))
        } else if (interval=="monthly") {
            filename_list <- rbind(filename_list,cbind(paste0(baseurl,"MYDAL2_M_CLD_FR/MYDAL2_M_CLD_FR_",y,"-",mname,".PNG")))
        }
    }
    gettcc_fname <- paste0(download_folder,y,"_get_tcc.txt")
    lapply(filename_list, write, gettcc_fname, append=T, ncolumns=1)
}
