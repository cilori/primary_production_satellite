exp_m600 <- function(x) {
    if (x < -600) {
        return(0)
    } else {
        return(exp(x))
    }
}