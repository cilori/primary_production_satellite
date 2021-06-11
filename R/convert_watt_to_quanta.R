# watt to quanta function

Watt_perm_2_to_microMol_per_m2_per_s = function(Ed = Ed, lambda = lambda){
	h = 6.626e-34 # in J.s
	c = 2.998e8 # in m.s-1
	Na = 6.022e23 # Avogadro number
    Np = Ed * lambda * 1e-9/(h*c)
	Np / (Na * 1e-6)
}
