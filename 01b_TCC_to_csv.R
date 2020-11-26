# Read TCC (total cloud cover) PNG file, generate a raster of bins the same size,
# and assign each png pixel a bin number. Save as a csv file.

library(data.table)
library(dplyr)
library(png)
library(raster)
library(oceancolouR)
library(stringr)

interval <- "8day"      # 8day or monthly
year <- c(2018)         # numeric vector of years
months <- 7             # numeric vector of months


#*******************************************************************************

input_path <- paste0("01a_downloaded_TCC/",interval,"/")
output_path <- paste0("01b_formatted_TCC/",interval,"/")

if (interval=="8day") {
    name1 <- "MYDAL2_E_CLD_FR_"
    # Make a vector of months and days at 8-day intervals for weekly images to extract month-day combinations.
    # NOTE: MODIS TCC files do not account for leap years in their filenames.
    mvec <- as.numeric(sapply(1:46,function(i) format(as.Date((8*0:45)[i],origin=paste0("2001-01-01")),"%m")))
    dates <- sapply(1:46,function(i) paste0(year,format(as.Date((8*0:45)[i],origin=paste0("2001-01-01")),"-%m-%d")))
    dates <- dates[mvec %in% months]
} else if (interval=="monthly") {
    name1 <- "MYDAL2_M_CLD_FR_"
    dates <- paste0(year,"-",sapply(months,function(i) str_pad(i,width=2,side="left",pad="0")))
}

bins.rl <- gen_bin_grid(gen_start_bin())

data("nwa_bins_4km", package="oceancolouR")
data("nwa_lats_4km", package="oceancolouR")
data("nwa_lons_4km", package="oceancolouR")

for (i in 1:length(dates)) {
    
    basename <- paste0(name1,dates[i])
    input_file <- paste0(basename,".PNG")
    output_file <- paste0(basename,".csv")
    
    if (!file.exists(paste0(input_path,input_file))) {next}
    
    cat(paste0("Creating ", output_path, output_file, "...\n"))
    
    # Load TCC file
    tcc <- raster(readPNG(paste0(input_path,input_file)))
    
    # Match projection and extent from bins raster to TCC raster
    crs(tcc) <- crs(bins.rl)
    extent(tcc) <- extent(bins.rl)
    
    # Resample to match the size of the sst, chl, and par grids, using the nearest neighbour method (ngb)
    new_tcc <- resample(x=tcc, y=bins.rl, method="ngb")
    
    # Extract values matching the lats/lons of the NWA region
    new_tcc <- extract(x=new_tcc, y=data.frame(x=nwa_lons_4km, y=nwa_lats_4km, stringsAsFactors = FALSE))
    
    # Create output dataframe with NWA bins and extracted tcc, and arrange by bin number
    output <- data.frame(bin=nwa_bins_4km, tcc=new_tcc, stringsAsFactors = FALSE) %>% arrange(bin)
    
    # Write to output file
    fwrite(output, paste0(output_path, output_file), row.names=FALSE, quote=FALSE)
    
}
