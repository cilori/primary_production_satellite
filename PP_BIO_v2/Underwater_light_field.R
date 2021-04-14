
# Underwater field as a function on Chl 

# For now, homogeneous profile, we can make it more complicated later

source("model_BIO_v2.R")
source("PAR_resolved.R")


chlpix = 8


PP <- pp_BIO_v2(chlpix = chlpix)$PP

print(PP)

