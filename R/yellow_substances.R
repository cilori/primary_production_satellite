yellow_substances <- function(chl) {
    # Source: yellow_subs_var.f90, original Fortran scripts at BIO
    if (chl > 0) {
        yellow_subs = (log10(chl) + 2) * 0.3
        if (chl > 1) {
            yellow_subs = yellow_subs*(log10(chl)+1)
        }
        if (yellow_subs < 0.0001) {
            yellow_subs = 0.0001
        }
    } else if (chl==0) {
        yellow_subs =  0.0001
    } else {
        return(NA)
    }
    return(yellow_subs)
}