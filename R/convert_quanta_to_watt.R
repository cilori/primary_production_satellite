# watt to quanta function

microMol_per_m2_per_s_to_watt_per_m2 = function(Ed = Ed, lambda = lambda){
    h = 6.626e-34 # in J.s
    c = 2.998e8 # in m.s-1
    Na = 6.022e23 # Avogadro number
    Ed * h * c * Na * 1e-6 / (wl * 1e-9)
}
