# Get the days to J2000
# h is UT in decimal hours
# days_to_J2000 only works between 1901 to 2099 - see Meeus chapter 7
days_to_J2000 <- function(y, m, day, h) {
    y <- as.integer(y)
    m <- as.integer(m)
    day <- as.integer(day)
    luku <- as.integer(- 7 * (y + (m + 9)/12)/4 + 275*m/9 + day)
    # Typecasting needed for TClite on PC DOS at least, to avoid product overflow
    luku <- luku + as.integer(y*367)
    luku <- luku - 730531.5 + h/24.0
    return (luku)
}
