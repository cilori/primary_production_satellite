# Primary production

*00_background* subfolder contains notes / images / research papers in pdf format for reference. Sample input and output data can be found in folders *01a*, *01b*, *01c*, *02*, *03*, and *04*, using a region in the NWA for July 2018.  

--------------------------------------------------------------------------------

## INSTRUCTIONS

### STEP 1A: Download MODIS TCC (total cloud cover)

- Edit variables in __01a_TCC_download.R__ and run it to create a list of TCC files to download, written to txt files in the `01a_downloaded_TCC/interval` subfolder, where *interval* is either *8day* or *monthly*, and the txt files are in the naming format __*year*_get_tcc.txt__.  
- Navigate to the folder containing the txt files, and type on the command line: `wget --random-wait --content-disposition -i [filename here, including .txt extension]`  
- Convert TCC .png files to .csv by editing the variables in the __01b_TCC_to_csv.R__ script and running it.  


### STEP 1B: Download VIIRS-SNPP CHL, PAR, SST

**If you have access to the PanCanadian dataset (located on _hecla_ as of 2020-11-20), skip this step**  

- Download daily VIIRS-SNPP files from the NASA ocean colour website (*note: you must download daily instead of 8day files for consistency with the formatting script*):  
This can be done using the *01a_download_daily_L3b.sh* from the panCan_processing repo (https://github.com/BIO-RSG/panCan_processing)  

- Subset the files to the PanCanadian grid (panCan_processing: *02a_subset_L3b_to_panCan_grid.py*)  

- Transfer the files to the `01c_other_variables` subfolder in the primary production repo, sorted into subfolders following the convention `sensor/variable/PANCAN/year`, where sensor is MODIS or VIIRS-SNPP, variable is SST, PAR, or CHL_OCX, and year is a 4-digit number.  

- When running `02a_format_PP_input.R`, make sure to adjust the input_path variable to the `01c_other_variables` folder.  

__NOTE 1:__ You need an Earthdata account to download these (see top of *01a_download_daily_L3b.sh* script for instructions)

__NOTE 2:__ In each script, don't forget to adjust input/output paths to the appropriate directories.  


### STEP 2: Format input

Format input data (CHL, PAR, SST, TCC) by editing the variables in __02a_format_PP_input.R__ and running it. Note that CHL, PAR, and SST are from the PanCanadian dataset.  
This merges the input data for each variable, with their bin numbers, latitudes, longitudes, and bathymetry. CHL, SST, and time of year are used as input to the function in __02b_get_param.R__, which retrieves corresponding BP (biomass profile) and PI (photosynthesis-irradiance) curve parameters, based on the look-up tables in *02a_LUT_cp-20141117.csv* and *02a_LUT_pi-20141113.csv*.  
All the data is then combined into a format that can be used as input to the primary production scripts, pixels with missing input data are removed, and the table is written to csv.  


### STEP 3: Calculate PP

Edit the variables in __03a_run_PP_BIO.R__ or __03b_run_PP_Laval.R__ and run them to create primary production files using either the BIO or Laval models, which will be written to the `03_PP_output/interval` subfolder, where *interval* is either *8day* or *monthly*. These "run" scripts take the user input and read and organize the formatted csv input, then call the main function for the selected model (*dwcpn* for the BIO model, located in `PP_BIO/model_BIO.R`, or *LAVAL_PP*, located in `PP_Laval/model_Laval.R`).  
Two files are created:  

- **NetCDF (.nc)**: This contains only one variable, the primary production in 4km-resolution, in the same size as the NWA (Northwest Atlantic) binned grid so that the values can be easily matched to NWA bin numbers. To view an image, you can run the following code, replacing *PP_filename* with your filename, including path and extension (.nc):  
```{r}
library(ncdf4)
library(raster)
library(oceancolouR)
nc <- nc_open(PP_filename)
PP <- ncvar_get(nc, "PP")
nc_close(nc)
data("nwa_bins_4km")
PP_rast <- var_to_rast(data.frame(bin=nwa_bins_4km, PP=PP, stringsAsFactors=FALSE), ext=c(-95,-42,39,82))
spplot(PP_rast)
```
- **Rdata (.rda)**: This contains the input data, primary production, and other variables outputted by the model function (for example, surface PAR). To view the first few rows of the table, run the following code:  
```{r}
load(PP_filename)
head(rda_output)
```

To view another variable produced with primary production and stored in the rda file (for example, surface PAR), you can use this code:  
```{r}
library(tidyr)
par0_df <- rda_output %>% dplyr::select(bin, PAR0_BIO)
par0_rast <- var_to_rast(df=par0_df, ext=c(-95,-42,39,82))
spplot(par0_rast)
```

Note that the *bird* and *sam_penguin* functions in the BIO model have been converted to C++ for the sake of speed, and are run by using the Rcpp package. (No changes needed in the code, just adjust variables at the top of __03a_run_PP_BIO.R__ as usual). *Bird* is the function that calculates the surface irradiance at all wavelengths for direct and diffuse irradiance separately, and *sam_penguin* calculates the underwater irradiance in the water column down to the euphotic depth.  


### STEP 4: Examine output

Examine the images in Rmarkdown summary files. You can also compare them to other models, such as the VGPM model, found here:  
`http://orca.science.oregonstate.edu/2160.by.4320.8day.hdf.vgpm.v.chl.v.sst.php `  


--------------------------------------------------------------------------------

## Data sources

**CHL, PAR, SST**:  

    PanCanadian dataset (see panCan_processing repository for details)  
    Level-3 binned (L3b) files from NASA  
    Link: https://oceandata.sci.gsfc.nasa.gov/  
    
**TCC**:  

    MODIS-Aqua 8day cloud fraction  
    Link: https://neo.sci.gsfc.nasa.gov/view.php?datasetId=MYDAL2_E_CLD_FR&date=2020-10-01  

**PI curve parameters (alphaB, PmB)**:  

    02a_LUT_pi-20141113.csv  

**Biomass profile parameters (h, sigma, B0, zm)**:  

    02a_LUT_cp-20141117.csv  
    (Note: this is only used if you have set uniform_profile = FALSE in 03a_run_PP_BIO.R or 03b_run_PP_Laval.R)


--------------------------------------------------------------------------------

## BIO primary production model

Brief summary of the steps in the PP calculation in *model_BIO.R*:  

1. Get *chla* profile (uniform or non-uniform)  
2. Calculate angles and time variables  
3. Loop over time intervals/zenith angles. At each iteration:  
    a. "bird" function: Calculate direct and diffuse components of surface irradiance for all the wavelengths  
    b. Correct for seasonal variation and cloud cover  
4. Scale/correct irradiance based on satellite PAR  
5. Loop over time intervals/zenith angles. At each iteration:  
    a. Calculate reflection losses at air-sea interface  
    b. Convert values to get the right units for sam_penguin  
    c. "sam_penguin" function: Calculate subsurface irradiance down to euphotic depth (wavelength integration happens in here)  
    d. Integrate over depth  
6. Integrate over time  

Note there are two functions in *model_BIO.R*: *dwcpn()* and *dwcpn_C()*. They are identical except that dwcpn uses *bird()* and *sam_penguin()* for the surface and subsurface irradiance calculations respectively, while dwcpn_C uses *bird_C()* and *sam_penguin_C()*, which are C++ versions of the functions that are much faster.  


### References

Platt, Trevor & Sathyendranath, Shubha & Forget, Marie-Helene & White, George & Caverhill, Carla & Bouman, Heather & Devred, Emmanuel & Son, Seunghyun. (2008). Operational estimation of primary production at large geographical scales. Remote Sensing of Environment - REMOTE SENS ENVIRON. 112. 10.1016/j.rse.2007.11.018.

Prieur, Louis & Sathyendranath, Shubha. (1981). An optical classification of coastal and oceanic waters based on the specific spectral absorption curves of phytoplankton pigments, dissolved organic matter, and other particulate matters. Limnology and oceanography. 26. 671. 10.4319/lo.1981.26.4.0671.

1991 model:
Platt, Trevor & Caverhill, Carla & Sathyendranath, Shubha. (1991). Basin-scale estimates of oceanic primary production by remote sensing - The North Atlantic. J. Geophys. Res.. 96. 10.1029/91JC01118.
Sathyendranath, Shubha & Longhurst, Alan & Caverhill, Carla & Platt, Trevor. (1995). Regionally and seasonally differentiated primary production in the North Atlantic. Deep-Sea Research (Part I, Oceanographic Research Papers). 42. 1773-1802. 10.1016/0967-0637(95)00059-F.


--------------------------------------------------------------------------------

## Laval primary production model





