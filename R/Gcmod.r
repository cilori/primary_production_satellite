
#*******************************************************************************
#  Atmospheric data

#  Opens and reads light data for the radiative transfer model
#  of Gregg and Carder 1990.  This is the data of Table 1 in the L&O
#  paper, with values added to extend the model to 800 nm.

data_atm <- read.table("data/gcirrad.txt",skip=5)
lam = rep(NaN,452)
Fobar = lam
oza = lam
ag = lam
aw = lam
oidx = seq(1,452,2)
eidx = seq(2,452,2)

lam[oidx] <- data_atm[,1]
lam[eidx] <- data_atm[,6]
Fobar[oidx] <- data_atm[,2]*10
Fobar[eidx] <- data_atm[,7]*10 #convert W/cm2/um to W/m2/nm
oza[oidx] <- data_atm[,3]
oza[eidx] <- data_atm[,8]
ag[oidx] <- data_atm[,4]
ag[eidx] <- data_atm[,9]
aw[oidx] <- data_atm[,5]
aw[eidx] <- data_atm[,10]


#*******************************************************************************
# Parameters to run atmospheric irradiance model

mxmu=10
mxphi=24
mxz=100
mxwave=90
mxnsky=10
mxcomp=10
mxnzvals=200
nlt=452

#rh = 0.8  # humidity default value
#am = 1 # aerosol model (1-10), 1 is default value (marine model)
#wv = 1.5 # precipitable water
#wsm = 4.0 # Average wind speed over 24 hours, default value (m/s)
#ws = 5.0 # instant wind speed (m/s)
#vis = 25  # visibility in km; enter 15 for default
#ro3 = -99 # ozone in Dobson Units, or enter -99 for climatological ozone
#     # values (computed from the day and Earth position)
#iblw = 0 # = 0 to compute above-surface values
#         # = 1 to compute below-surface values
#p0 = 29.92 # reference atmospheric pressure at sea level
#pres = 29.96 # Real time atmospheric pressure

rad = 180/pi
pi2 <- 2*pi


#*******************************************************************************
# FUNCTIONS

comp_sco3 <- function(ylat,xlon,jday,ro3) {
  rlat <- ylat/rad
  rlon <- xlon/rad
  to3 <- 235.0 + (150.0+40.0*sin(0.9865*(jday-30.0)) + 20.0*sin(3.0*(rlon)))*sin(1.28*rlat)*sin(1.28*rlat)
  if (ro3 > 0) {
    to3 <- ro3
  } else {
    to3 = 235.0 + (150.0+40.0*sin(0.9865*(jday-30.0)) + 20.0*sin(3.0*(rlon)))*sin(1.28*rlat)*sin(1.28*rlat)
  }
  sco3 = to3*1e-3
  return(sco3)
}


#  Correct for Earth-Sun distance
#  modify date if it is last day of a leap year
solar_irr <- function(Fobar,jday) {
  if(jday >365) {jday=365.0}
  Fo = Fobar * (1.0+1.67E-2*cos(pi2*(jday-3)/365.0))^2
  return(Fo)
}


# Computation of irradiance accross air-interface if needed
sfcrfl <- function(rad,theta,ws) {

  #  Computes surface reflectance for direct (rod) and diffuse (ros)
  #  components separately, as a function of theta, wind speed or stress.

  rn = 1.341      #index of refraction of pure seawater
  roair = 1.2E3   #density of air g/m3

  #  Foam and diffuse reflectance
  if (ws > 4.0) {
    if (ws < 7.0) {
      cn = 6.2E-4 + 1.56E-3/ws
      rof = roair*cn*2.2E-5*ws*ws - 4.0E-4
    } else {
      cn = 0.49E-3 + 0.065E-3*ws
      rof = (roair*cn*4.5E-5 - 4.0E-5)*ws*ws
    }
    rosps = 0.057
  } else {
    rof = 0.0
    rosps = 0.066
  }
  
  #  Direct
  #  Fresnel reflectance for theta < 40, ws < 2 m/s
  if (theta < 40.0 || ws < 2.0) {
    if (theta == 0.0) {
      rospd = 0.0211
    } else {
      rtheta = theta/rad
      sintr = sin(rtheta)/rn
      rthetar = asin(sintr)
      rmin = rtheta - rthetar
      rpls = rtheta + rthetar
      sinp = (sin(rmin)*sin(rmin))/(sin(rpls)*sin(rpls))
      tanp = (tan(rmin)*tan(rmin))/(tan(rpls)*tan(rpls))
      rospd = 0.5*(sinp + tanp)
    }
  } else { # Empirical fit otherwise
    a = 0.0253
    b = -7.14E-4*ws + 0.0618
    rospd = a*exp(b*(theta-40.0))
  }
  
  #  Reflectance totals
  rodaux = rospd + rof
  rosaux = rosps + rof
  list( rod = rodaux, ros = rosaux)
  
}


navaer <- function(rh,am,wsm,ws,vi) {
  
  #  Computes aerosol parameters according to a simplified version
  #  of the Navy marine aerosol model.
  
  ro <- c(0.03,0.24,2.0)
  r <- c(0.1,1.0,10.0)
  
  rlam = 0.55

  #  Relative humidity factor
  if (rh >= 100.0) rh = 99.9
  rnum = 2.0 - rh/100.0
  rden = 6.0*(1.0-rh/100.0)
  frh = (rnum/rden)**0.333

  #  Size distribution amplitude components
  a = c(2000.0*am*am,5.866*(wsm-2.2),0.01527*(ws-2.2)*0.05)
  if (a[2] < 0.5) a[2] = 0.5
  if (a[3] < 1.4E-5)a[3] = 1.4E-5

  #  Compute size distribution at three selected radii according to Navy method
  dndr <- rep(0,3)
  for (n in 1:3) {
    dndr[n] = 0.0
    for (i in 1:3) {
      rden = frh*ro[i]
      arg = log(r[n]/rden)*log(r[n]/rden)
      rval = a[i]*exp(-arg)/frh
      dndr[n] = dndr[n] + rval
    }
  }

  #  Least squares approximation

  sumx = sum(log10(r))
  sumy = sum(log10(dndr))
  sumxy = sum(log10(r)*log10(dndr))
  sumx2 = sum(log10(r)^2)
  
  gama = sumxy/sumx2
  rlogc = sumy/3.0 - gama*sumx/3.0
  alphaux = -(gama+3.0)

  #  Compute beta
  cext = 3.91/vis
  betaux = cext*rlam^alphaux

  #  Compute asymmetry parameter -- a function of alpha
  if (alphaux > 1.2) {
    asympaux = 0.65
  } else if (alphaux > 0.0) {
    asympaux = 0.82
  } else {
    asympaux = -0.14167*alphaux + 0.82
  }

  #  Single scattering albedo at 550; function of RH
  waux = (-0.0032*am + 0.972)*exp(3.06E-4*rh)

  list( alpha = alphaux, beta = betaux, asymp = asympaux, wa = waux)
  
}


atmodd <- function(iblw=iblw,rad=rad,lam=lam,theta=theta,oza=oza,
                   ag=ag,aw=aw,sco3=sco3,p=p,wv=wv,rh=rh,am=am,
                   wsm=wsm,ws=ws,vis=vis,Fo=Fo) {
  
  #  Model for atmospheric transmittance of solar irradiance through
  #  a maritime atmosphere.  Computes direct and diffuse separately.
  #  Includes water vapor and oxygen absorption.

  #  Compute atmospheric path lengths (air mass); pressure-corrected
  cosunz = cos(theta/rad)
  #  Modified March, 1994 according to Kasten and Young 1989.

  rex = -1.6364
  rtmp = (96.07995-theta)^rex
  rm = 1.0/(cosunz+0.50572*rtmp)
  rmp = p/p0*rm
  otmp = (cosunz*cosunz+44.0/6370.0)^0.5
  rmo = (1.0+22.0/6370.0)/otmp

  #  Obtain aerosol parameters; simplified Navy aerosol model
  res.aer <- navaer(rh,am,wsm,ws,vis)
  alpha <- res.aer$alpha
  beta <- res.aer$beta
  asymp <- res.aer$asymp
  wa <- res.aer$wa

  eta = -alpha
  #   Forward scattering probability
  alg = log(1.0-asymp)
  afs = alg*(1.459+alg*(.1595+alg*.4129))
  bfs = alg*(.0783+alg*(-.3824-alg*.5874))
  Fa = 1.0 - 0.5*exp((afs+bfs*cosunz)*cosunz)

  # Surface reflectance
  if (iblw == 1) {
    res.sfc <- sfcrfl(rad,theta,ws)
    rod <- res.sfc$rod
    ros <- res.sfc$ros
  } else {
    rod = 0.0
    ros = 0.0
  }

  # Compute spectral irradiance
  #   Rayleigh, by Bird's method
  rlam = lam*1.0E-3
  tr = 1.0/(115.6406*rlam^4 - 1.335*rlam^2)
  rtra = exp(-tr*rmp)
  #    Ozone
  to = oza*sco3   #optical thickness
  otra = exp(-to*rmo)   #transmittance
  #   Aerosols
  ta = beta*rlam^eta
  atra = exp(-ta*rm)
  taa = exp(-(1.0-wa)*ta*rm)
  tas = exp(-wa*ta*rm)
  #   Oxygen/gases
  gtmp = (1.0 + 118.3*ag*rmp)^0.45
  gtmp2 = -1.41*ag*rmp
  gtra = exp(gtmp2/gtmp)
  #   Water Vapor
  wtmp = (1.0+20.07*aw*wv*rm)^0.45
  wtmp2 = -0.2385*aw*wv*rm
  wtra = exp(wtmp2/wtmp)

  #  Direct irradiance
  Ediraux = Fo*cosunz*rtra*otra*atra*gtra*wtra*(1.0-rod)

  #   Diffuse irradiance
  dray = Fo*cosunz*gtra*wtra*otra*taa*0.5*(1.0-rtra^.95)
  daer = Fo*cosunz*gtra*wtra*otra*rtra^1.5*taa*Fa*(1.0-tas)

  #  Total diffuse
  Edifaux = (dray + daer)*(1.0-ros)

  Edaux = Ediraux + Edifaux

  list(Edir = Ediraux, Edif= Edifaux,Ed = Edaux)
  
}
