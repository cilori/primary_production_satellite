raman = 0
if (raman == 1) {lambda <- raman_wv(lambda_e)} else {lambda <- lambda_e}

# computation of average cosine

twopi = 2*pi
rt2pi = sqrt(twopi)
               
##################################################################
# Components
##################################################################
# PURE SEAWATER
#absorption by water

awp <- read.table("pope97.dat",skip=6) # absorption of pure sea water according to Pope and Fry (1997)

res <- approx(awp[,1],awp[,2],lambda) # asw is absorption by pure seawater at wavelengths of interest
asw <- res$y

# Scattering by water  
bw500 = 0.00288 #at 500 nm
bw = bw500*(lambda/500)^(-4.3)

##################################################################
# PHYTOPLANKTON
# backscattering by particles (phytoplankton en others, following
# Sathyendranath et al 2001)

physca <- function(bbp555,bbps) (bbp555*(555./lambda)^(bbps))

######################################################################
# YELLOW SUBSTANCES
#absorption by yellow substances

fays <- function(ay443,S)
{ay443*exp(-S*(lambda-443))}

#raman scattering at 488 nm
br488 = 0.00027
bbr=0.5*br488*(lambda/488)^(-5.3)

# s: ratio between backscattering and upward-scattering
s=1.

