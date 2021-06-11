
library(ggplot2)
library(gganimate)
library(png)
library(gifski)
library(transformr)
library(oce)
library(insol)
library(stringr)
library(lubridate)
library(Rcpp)
library(patchwork)
library(dplyr)

source("R/bird.R")
source("R/I0_corr_season.R")
source("R/I0_corr_cloud_BIO.R")
source("R/I0_corr_refl_loss.R")
source("R/PAR_resolved.R")
source("R/Gcmod.r")
sourceCpp("underwater_irradiance.cpp")
sourceCpp("surface_irradiance.cpp")



# yearpxl = 2018
# hours <- 6:20
# 
# # days <- c(80,172,266,355) # test on equinoxes and solstices
# # lats <- c(40, 50, 60, 70, 80)
# days <- c(80,172) # test on equinoxes and solstices
# lats <- c(50)
# 
# 
# # use dataframe with expand.grid for all combinations of day/hour/lat
# totest <- expand.grid(days,hours,lats)
# 
# 
# irrad_df <- data.frame(matrix(nrow=dim(totest)[1]*301*2, ncol=9), stringsAsFactors = FALSE)
# colnames(irrad_df) <- c("Day", "Hour", "Latitude", "Solar_zenith_angle", "Model", "Wavelength", "Direct", "Diffuse", "Total")
# 
# bird_wl = seq(0.4,0.7,by=0.005)
# gc_wl = seq(0.4,0.7,by=0.001)
# 
# for (i in 1:nrow(totest)) {
#     
#     daypxl = totest[i,1]
#     hrpxl = totest[i,2]
#     latipxl = totest[i,3]
#     
#     # zenith angle, in radians (bird requires radians, gc requires degrees)
#     ZEN = (90 - sunAngle(t=str_replace(as.character(doyday(yearpxl,daypxl) + hours(hrpxl)),"AST","UTC"), longitude = 0, latitude = latipxl)$altitude) * pi/180
#     
#     coart_light <- NA
#     
#     # BIRD above-surface irradiance in w.m-2.um-1
#     irr.bird = bird(ZEN, bird_wl)
#     if (length(irr.bird$DIRECT)==1) {
#         bird_dir <- bird_dif <- bird_total <- rep(NA,301)
#     } else {
#         irr.seacor = I0_corr_season(I0_direct = irr.bird$DIRECT,I0_diffuse = irr.bird$DIFFUSE,Day = daypxl)
#         irr.cloudcorr = I0_corr_cloud_BIO(I0_direct = irr.seacor$I0_direct, I0_diffuse = irr.seacor$I0_diffuse, ZEN = ZEN,Cloud = 0)
#         bird_light = I0_corr_refl_loss(I0_direct = irr.cloudcorr$I0_direct, I0_diffuse = irr.cloudcorr$I0_diffuse, ZEN=ZEN)
#         bird_dir <- bird_light$I0_direct
#         bird_dir <- approx(x=bird_wl, y=bird_dir, xout=gc_wl)$y
#         bird_dif <- bird_light$I0_diffuse
#         bird_dif <- approx(x=bird_wl, y=bird_dif, xout=gc_wl)$y
#         bird_total <- bird_dir + bird_dif
#     }
#     
#     # GREGG-CARDER above-surface irradiance in w.m-2.nm-1
#     sco3 = comp_sco3(ylat = latipxl, xlon = 0,jday = daypxl,ro3 = 100)
#     Fo = solar_irr(Fobar,daypxl)
#     gc_light = atmodd_v2_c(rod=0, ros=0, rad=180/pi, lam=gc_wl*1000, theta=ZEN*180/pi, oza=oza, ag=ag, aw=aw, sco3=sco3, p=29.92, p0=29.92, wv=1.5, rh=0.8, am=1, wsm=0, ws=0, vis=40, Fo=Fo)
#     gc_dir <- gc_light$Edir * 1000
#     gc_dif <- gc_light$Edif * 1000
#     gc_total <- gc_light$Ed * 1000
#     
#     
#     inds <- ((301*2)*(i-1)+1):((301*2)*i)
#     irrad_df[inds,] <- data.frame(Day=daypxl, Hour=hrpxl, Latitude=latipxl, Solar_zenith_angle=ZEN*180/pi,
#                                   Model=c(rep("Bird",301), rep("Gregg-Carder",301)),
#                                   Wavelength=c(gc_wl, gc_wl),
#                                   Direct=c(bird_dir, gc_dir),
#                                   Diffuse=c(bird_dif, gc_dif),
#                                   Total=c(bird_total, gc_total),
#                                   stringsAsFactors = FALSE)
#     
# }
# 
# irrad_df$Day <- as.integer(irrad_df$Day)
# irrad_df$Hour <- as.integer(irrad_df$Hour)
# 
# for (day_to_test in days) {
#     for (lat_to_test in lats) {
#         test = irrad_df %>%
#             dplyr::filter(Day==day_to_test & Latitude==lat_to_test) %>%
#             tidyr::pivot_longer(Direct:Total, names_to="Type", values_to="Irradiance")
#         
#         p <- ggplot(test) +
#             geom_line(aes(x=Wavelength, y=Irradiance, color=Model)) +
#             theme_bw() +
#             labs(x="Wavelength (μm)", y="Above water surface irradiance (W.m-2.μm-1)") +
#             facet_grid(Type~., scales="free_y") +
#             theme(legend.position=c(0.8,1.05),
#                   legend.direction="horizontal") +
#             transition_time(Hour) +
#             ggtitle(paste0("Day: ", day_to_test, "\n",
#                            "Latitude: ", lat_to_test, "°N"),
#                     subtitle = paste0("Hour: {frame_time},   ",
#                                       "Solar zenith angle: {test %>% dplyr::filter(Hour == frame_time & Wavelength==0.4 & Model==\"Bird\" & Type==\"Direct\") %>% dplyr::select(Solar_zenith_angle) %>% round()}")) +
#             ease_aes('linear')
#         
#         animate(p, height = 800, width = 500)
#         anim_save(paste0("PARmodel_bird-vs-GC_day",day_to_test,"_lat",lat_to_test,".gif"))
#     }
# }





stop()
#*******************************************************************************

# USING ZENITH ANGLE INSTEAD OF DAY/LAT/HOUR - this is all you really need:


yearpxl = 2018
latipxl = 50
daypxl = 172

zeniths <- seq(10,90,by=10)

irrad_df <- data.frame(matrix(nrow=length(zeniths)*301*2, ncol=6), stringsAsFactors = FALSE)
colnames(irrad_df) <- c("Solar_zenith_angle", "Model", "Wavelength", "Direct", "Diffuse", "Total")

bird_wl = seq(0.4,0.7,by=0.005)
gc_wl = seq(0.4,0.7,by=0.001)

for (i in 1:length(zeniths)) {

    # zenith angle, in radians (bird requires radians, gc requires degrees)
    ZEN = zeniths[i] * pi/180

    coart_light <- NA

    # BIRD above-surface irradiance in w.m-2.um-1
    irr.bird = bird(ZEN, bird_wl)
    if (length(irr.bird$DIRECT)==1) {
        bird_dir <- bird_dif <- bird_total <- rep(NA,301)
    } else {
        irr.seacor = I0_corr_season(I0_direct = irr.bird$DIRECT,I0_diffuse = irr.bird$DIFFUSE,Day = daypxl)
        irr.cloudcorr = I0_corr_cloud_BIO(I0_direct = irr.seacor$I0_direct, I0_diffuse = irr.seacor$I0_diffuse, ZEN = ZEN,Cloud = 0)
        bird_light = I0_corr_refl_loss(I0_direct = irr.cloudcorr$I0_direct, I0_diffuse = irr.cloudcorr$I0_diffuse, ZEN=ZEN)
        bird_dir <- bird_light$I0_direct
        bird_dir <- approx(x=bird_wl, y=bird_dir, xout=gc_wl)$y
        bird_dif <- bird_light$I0_diffuse
        bird_dif <- approx(x=bird_wl, y=bird_dif, xout=gc_wl)$y
        bird_total <- bird_dir + bird_dif
    }

    # GREGG-CARDER above-surface irradiance in w.m-2.nm-1
    sco3 = comp_sco3(ylat = latipxl, xlon = 0,jday = daypxl,ro3 = 100)
    Fo = solar_irr(Fobar,daypxl)
    gc_light = atmodd_v2_c(rod=0, ros=0, rad=180/pi, lam=gc_wl*1000, theta=zeniths[i], oza=oza, ag=ag, aw=aw, sco3=sco3, p=29.92, p0=29.92, wv=1.5, rh=0.8, am=1, wsm=0, ws=0, vis=40, Fo=Fo)
    gc_dir <- gc_light$Edir * 1000
    gc_dif <- gc_light$Edif * 1000
    gc_total <- gc_light$Ed * 1000


    inds <- ((301*2)*(i-1)+1):((301*2)*i)
    irrad_df[inds,] <- data.frame(Solar_zenith_angle=zeniths[i],
                                  Model=c(rep("Bird",301), rep("Gregg-Carder",301)),
                                  Wavelength=c(gc_wl, gc_wl),
                                  Direct=c(bird_dir, gc_dir),
                                  Diffuse=c(bird_dif, gc_dif),
                                  Total=c(bird_total, gc_total),
                                  stringsAsFactors = FALSE)

}


test = irrad_df %>%
    dplyr::mutate(Sun_altitude = -Solar_zenith_angle + 90) %>%
    tidyr::pivot_longer(Direct:Total, names_to="Type", values_to="Irradiance")

p <- ggplot(test) +
    geom_line(aes(x=Wavelength, y=Irradiance, color=Model)) +
    theme_bw() +
    labs(x="Wavelength (μm)", y="Above water surface irradiance (W.m-2.μm-1)") +
    facet_grid(Type~., scales="free_y") +
    theme(legend.position=c(0.8,1.02),
          legend.direction="horizontal") +
    transition_states(Sun_altitude, state_length=c(0.8,rep(0.2,6),0.8,0.8)) +
    ggtitle("Sun altitude: {closest_state}") +
    ease_aes('cubic-in-out')

animate(p, height = 800, width = 500)
anim_save("PARmodel_bird-vs-GC_by-sun-altitude.gif")




stop()
#*******************************************************************************




yearpxl = 2018
days <- 1:365
hours <- seq(0,23,by=2)
lats <- c(40,50,60,70,80)

totest <- expand.grid(days,hours,lats)

sunalt_df <- data.frame(matrix(nrow=nrow(totest),ncol=4), stringsAsFactors = FALSE)
colnames(sunalt_df) <- c("Day", "Hour", "Latitude", "Sun_altitude")

for (i in 1:nrow(totest)) {
    daypxl = totest[i,1]
    hrpxl = totest[i,2]
    latipxl = totest[i,3]
    sunalt = sunAngle(t=str_replace(as.character(doyday(yearpxl,daypxl) + hours(hrpxl)),"AST","UTC"), longitude = 0, latitude = latipxl)$altitude
    sunalt_df[i,] <- c(daypxl, hrpxl, latipxl, sunalt)
}

test <- sunalt_df %>%
    dplyr::group_by(Day, Latitude) %>%
    dplyr::summarize(min_zen=min(Sun_altitude),
                     max_zen=max(Sun_altitude)) %>%
    dplyr::ungroup()

sunaltp <- ggplot(test) +
    geom_rect(aes(xmin=0, xmax=365, ymin=-90, ymax=0), fill="#222222", alpha=0.05) +
    geom_rect(aes(xmin=0, xmax=365, ymin=0, ymax=90), fill="white", alpha=0.05) +
    geom_text(aes(x=180,y=-85,label="DARK"),color="gold") +
    geom_text(aes(x=180,y=85,label="LIGHT"),color="gold") +
    geom_hline(yintercept=0, color="red") +
    geom_ribbon(aes(x=Day, ymin=min_zen, ymax=max_zen), fill="grey", alpha=0.8) +
    geom_text(aes(30, 0, label = "Horizon", vjust = -1), color="red") +
    theme_bw() +
    scale_y_continuous(expand = c(0,0), breaks=seq(-90,90,by=10), labels=seq(-90,90,by=10)) +
    scale_x_continuous(expand = c(0,0)) +
    labs(x="Day of year", y="Sun altitude") +
    theme(plot.title=element_text(size=14)) +
    transition_states(Latitude) +
    ggtitle("Latitude: {closest_state}°N") +
    ease_aes('cubic-in-out')


animate(sunaltp, height = 600, width = 350)
anim_save("sun_altitudes.gif")


