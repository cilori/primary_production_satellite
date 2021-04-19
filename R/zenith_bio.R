# Calculate zenith angle of the sun for a given location, date, and time
zenith_bio <- function(TIME, LONG, delta, phi) {
    
   # THE FORMULAE USED IN THE CALCULATIONS ARE TAKEN FROM:
   # G.W. PALTRIDGE AND C.M.R. PLATT, 1976, "RADIATIVE PROCESSES
   # IN METEOROLOGY AND CLIMATOLOGY", ELSEVIER, AMSTERDAM, P 60-63.
   # 
   # TIME (LOCAL) AS HOURS PLUS FRACION OF HOUR IS INPUT FROM THE MAIN SCRIPT
   # 
   # DELTA  = DECLINATION (error in delta equation is < 3', maximum error 35")
   # LONG   = LONGITUDE
   # PHI    = LATITUDE IN RADIANS
   # TH     = SUN'S HOUR ANGLE IN RADIANS, COUNTING FROM MID DAY,
   #          AND MEASURED WESTWARD FROM THE OBSERVERS'S MERIDIAN
   # TIME   = LOCAL APPARENT TIME, IN HOURS
   # ZEN    = ZENITH ANGLE
    
    # IF TIME IS IN GMT, THEN ACTIVATE THE NEXT LINE PLUS THE FOURTH LINE. 
    #TH = TIME + 4*LONG/60
    TH = TIME - 12
    TH = TH * (2*pi/24) 
    #TH = TH + EQ
    
    ZEN = sin(delta) * sin(phi) + cos(delta) * cos(phi) * cos(TH)
    if (ZEN < -1) ZEN = -1
    if (ZEN > 1) ZEN = 1
    ZEN = pi/2 - asin(ZEN)
    
    return(ZEN)

}
