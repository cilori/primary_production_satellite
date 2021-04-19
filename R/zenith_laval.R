# Model used at Laval
# Sun position from doy, hour, ycoor (longi), xcoor (lati)
zenith_laval = function(doy=doy, hour=hour, xcoor=xcoor, ycoor=ycoor) {
    
    # long and lat are in degrees
    
    ltm = 0
    hr = floor(hour)
    minute = ( (hour - hr) * 60.)
    
    d2r     = pi/180.0
    r2d     = 1/d2r
    lsn     = 12.0+((ltm-xcoor)/15.0)
    latrad  = ycoor*d2r
    decrad  = 23.45*d2r*sin(d2r*360.*(284.+doy)/365.)
    decdeg  = decrad*r2d
    
    
    ha      = hr+(minute/60.0) #convert to floating
    hangle  = (lsn-ha)*60.0               #solrad is given for the hour preceeding the time given
    harad   = hangle*0.0043633            #convert hangle (in minutes) into radians
    #This is the same as multiplying the #hrs by 15 (deg./hr),
    #and convert d2r (*0.017453292)
    
    saltrad = asin((sin(latrad)*sin(decrad))+(cos(latrad)*cos(decrad)*cos(harad)))
    saltdeg = saltrad * r2d
    sazirad = asin(cos(decrad)*sin(harad)/cos(saltrad))
    sazideg = sazirad * r2d
    
    if (saltdeg < 0.0 | saltrad  > 180.0) {   # sun is below horizon
        saltdeg = 0.0
        saltrad = 0.0
        szendeg = 90.0
        szenrad = 90.0*d2r
        mass    = 1229**.5             # if solaralt=0 -> sin(0)=0
    }
    else {
        szendeg = 90-saltdeg
        szenrad = szendeg*d2r
        mass = (1229.0+(614.0*sin(saltrad))**2)**.5-(614*sin(saltrad))
    }
    
    return(list(szendeg=szendeg,szenrad=szenrad,airmass=mass))
}

