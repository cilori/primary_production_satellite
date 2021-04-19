sunang <- function(iday,hr,rad,xlon,ylat) {
    
    #     Computes sun azimuth and zenith angles for a given
    #     Julian day, GMT, latitude, and longitude.  This program
    #     is from the NMFS ELAS computer code.  Modified for
    #     standard coordinates (W long. negative), to correct
    #     for dateline problem, and to correct coefficients (taken
    #     from Iqbal, 1983, An Introduction to Solar Radiation).
    #     Watson Gregg, Research and Data Systems, Corp.
    #     Further modified to remove unused input (CDM, March 98)
    
    #     INPUT:
    #     iday = julian day
    #     hr = GMT in hours (e.g. 21.5 for 21 hours, 30 minutes GMT)
    #     rad = radian to degrees conversion factor
    #     ylat = latitude of pixel
    #     xlon = longitude of pixel
    #     OUTPUT:
    #     sunz = solar zenith angle in degrees
    #     suna = solar azimuth angle in degrees
    
    #     internal variables:
    #     sdec = solar declination angle in degrees
    #     thez = theta zero orbital position in degrees
    #     tc = time correction
    #     xha = solar hour angle in degrees
    
    #  Compute solar declination angle
    thez = 360.0*(iday-1)/365.0
    rthez = thez/rad
    sdec = 0.396372-22.91327*cos(rthez) + 4.02543*sin(rthez)- 0.387205*cos(2.0*rthez) + 0.051967*sin(2.0*rthez) - 0.154527*cos(3.0*rthez) + 0.084798*sin(3.0*rthez)
    rsdec = sdec/rad
    
    #  Time correction for solar hour angle, and solar hour angle
    tc = 0.004297 + 0.107029*cos(rthez) - 1.837877*sin(rthez) - 0.837378*cos(2.0*rthez) - 2.342824*sin(2.0*rthez)
    xha = (hr-12.0)*15.0 + xlon + tc
    if (xha >  180.0) {xha = xha-360.0}
    if (xha < -180.0) {xha = xha+360.0}
    rlat = ylat/rad
    #      rlon = xlon/rad
    rha = xha/rad
    
    #  Sun zenith
    costmp = sin(rlat)*sin(rsdec) + cos(rlat)*cos(rsdec)*cos(rha)
    eps = abs(costmp)
    if (eps > 1.1){
        'Error in acos argument in sun zenith'
        costmp} else { if (eps > 1.0) {
            if (costmp > 0.0)costmp=1.0
            if (costmp > 0.0)costmp=-1.0}}
    
    rsunz = acos(costmp)
    #  Sun azimuth
    sintmp = sin(abs(rha))*cos(rsdec)/sin(rsunz)
    rsuna = asin(sintmp)
    if(ylat > sdec) {rsuna = 180.0 / rad - rsuna}
    if(xha > 0.) {rsuna = 360.0/rad - rsuna}
    
    #  Convert to degrees
    suna = rsuna * rad
    
    sunz = rsunz*rad
    
    list( sunz = sunz, suna = suna)
}
