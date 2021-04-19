# CORRECT FOR SEASONAL VARIATION IN SOLAR ENERGY
# Input: surface direct and diffuse irradiance (each a vector of values across
# wavelengths), and the Day, a numeric ZERO-INDEXED day of year (0-364).
I0_corr_season <- function(I0_direct,I0_diffuse,Day) {
    
    # DayNo   = DAY NUMBER
    DayNo = c(0,3,31,42,59,78,90,93,120,133,151,170,181,183,206,212,243,265,273,
              277,304,306,334,355,365)
    
    # SolIr   = EXTRATERRESTRIAL SOLAR IRRADIANCE
    # Values are from Thekaekara 1977 (ref in JGR 1988 - s & p)
    SolIr = c(1399.,1399.,1393.,1389.,1378.,1364.,1355.,1353.,1332.,1324.,1316.,1310.,
              1309.,1309.,1312.,1313.,1329.,1344.,1350.,1353.,1347.,1375.,1392.,1398.,1399.)
    
    # Find nearest day index that is greater or equal to satellite Day.
    D = which(DayNo >= Day)[1]
    
    # Set value of solar irradiance according to SolIr data array
    if (D==1) {
        sol = SolIr[1]
    } else {
        sol = SolIr[D-1]-(SolIr[D-1]-SolIr[D]) * ((Day-DayNo[D-1])/(DayNo[D]-DayNo[D-1]))
    }
    
    # Adjust direct and diffuse components by fractional correction
    # Note: 1353. is max value (100%)
    I0_direct_corrected = I0_direct * sol/1353
    I0_diffuse_corrected = I0_diffuse * sol/1353
    
    return(list("I0_direct"=I0_direct_corrected,"I0_diffuse"=I0_diffuse_corrected))
    
}
