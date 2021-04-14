VGPM = function( chl = chl, irr =  irr , sst = sst, dayL = dayL) 
  {
  #C--------------------------------------------------------------------------*\
  
  #!Description:     opp_befa - computes daily primary productivity using
  #the Behrenfeld-Falkowski (BeFa) algorithm.  The BeFa
  #algorithm estimates productivity using surface chl
  #(mg m-3), surface irradiance (Einsteins m-2 d-1),
  #sea surface temperature (C), and day length (hours).
  #Pb_opt is modelled as a polynomial function of SST.
  
  #!Input Parameters:  
  # chl            Chlorophyll_a surface concentration in milligrams chlorophyl per cubic meter
  # irr            Photosynthetically available radiation in Einsteins per day per square meter
  #sst            Sea surface temperature in degrees Centigrade
  #dayL           Length day in decimal hours. 
  
  #Output Parameters: 
  # PP = Primary productivity in milligrams Carbon per square meter per day
  
  # References and Credits
  # 
  #    Behrenfeld,M.J; Falkowski,P.G.; 1997. Photosynthetic Rates Derived
  #    from Satellite-Based Chlorophyll Concentration.  Limnology and 
  #    Oceanography, Volume 42, Number 1
      


#   /* Calculate euphotic depth (z_eu) with Morel's Case I model.            */
#    /* Calculate chl_tot from Satellite Surface Chlorophyll Data.            */
    
    if (chl <  1.0) 
      {chl_tot = 38.0 * chl^0.425 }
  if (chl >=  1.0)
      {chl_tot = 40.2 * chl^0.507}
  
  
  z_eu = 200.0 * chl_tot^(-.293)
  
  if (z_eu <= 102.0) {z_eu = 568.2 * chl_tot^(-.746)}
  
  
#  /* Calculate the Pb_opt from satellite sea surface temperature (sst).    */
    
     pb_opt = 1.2956 + 2.749e-1*sst + 6.17e-2*sst^2 - 2.05e-2*sst^3 +
     2.462e-3*sst^4 - 1.348e-4*sst^5 + 3.4132e-6*sst^6 - 3.27e-8*sst^7
    if (sst < -10.0) {pb_opt = 0.00} 
    if (sst <  -1.0 & sst >= -10) {pb_opt = 1.13} 
    if (sst >  28.5) {pb_opt = 4.00}

  
  
#  /* calculate the irradiance function */
    
    irrFunc = 0.66125 * irr / ( irr + 4.1 )
  
  
#  /* Return the primary production calculation.                            */
    
    PP = pb_opt * chl * dayL * irrFunc * z_eu;
  
  return(list(PP = PP, z_eu = z_eu, pb_opt = pb_opt))
}