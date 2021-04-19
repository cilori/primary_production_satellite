/* =========================================================================
 prodprim.c
 
 author: Maxime Benoit-Gagne
 - Takuvik - Canada.
 David Dessailly
 - Laboratoire d'Oceanologie et Geoscience Wimereux-France.
 creation: November 2011.
 
 compiler:
 $ gcc -v
 Configured with: --prefix=/Applications/Xcode.app/Contents/Developer/usr --with-gxx-include-dir=/usr/include/c++/4.2.1
 Apple LLVM version 5.1 (clang-503.0.40) (based on LLVM 3.4svn)
 Target: x86_64-apple-darwin12.6.0
 Thread model: posix
 
 usage:
 ./prodprim
 	--rrs_type rrs_type
  --rrs_file rrs_file
  --atmosphere_type atmosphere_type
  [--atmosphere_file atmosphere_file]
  [--atmosphere_file00 atmosphere_file00
  --atmosphere_file03 atmosphere_file03
  --atmosphere_file06 atmosphere_file06
  --atmosphere_file09 atmosphere_file09
  --atmosphere_file12 atmosphere_file12
  --atmosphere_file15 atmosphere_file15
  --atmosphere_file18 atmosphere_file18
  --atmosphere_file21 atmosphere_file21
  --atmosphere_file00tomorrow atmosphere_file00tomorrow]
  --ice_file ice_file
  --doy doy
  [--first_bin first_bin]
  [--last_bin last_bin]
  --outfile outfile
 where
 rrs_type is the sensor for the Rrs, S for SeaWiFS and A for MODISA.
 rrs_file is the Rrs file.
 atmosphere_type is the dataset for the atmospheric products, I for ISCCP and
  M for MODISA.
 atmosphere_file is the file for the atmospheric products for one day.
  This argument is only valid with atmosphere_type = M.
 atmosphere_file00 is the file for the atmospheric products at 00 UTC.
  This argument is only valid with atmosphere_type = I.
 atmosphere_file03 is the file for the atmospheric products at 03 UTC.
  This argument is only valid with atmosphere_type = I.
 atmosphere_file06 is the file for the atmospheric products at 06 UTC.
  This argument is only valid with atmosphere_type = I.
 atmosphere_file09 is the file for the atmospheric products at 09 UTC.
  This argument is only valid with atmosphere_type = I.
 atmosphere_file12 is the file for the atmospheric products at 12 UTC.
  This argument is only valid with atmosphere_type = I.
 atmosphere_file15 is the file for the atmospheric products at 15 UTC.
  This argument is only valid with atmosphere_type = I.
 atmosphere_file18 is the file for the atmospheric products at 18 UTC.
  This argument is only valid with atmosphere_type = I.
 atmosphere_file21 is the file for the atmospheric products at 21 UTC.
  This argument is only valid with atmosphere_type = I.
 atmosphere_file00tomorrow is the file for the atmospheric products at 00 UTC
  tomorro. This argument is only valid with atmosphere_type = I.
 ice_file is the sea-ice concentration file.
 doy is the day of year.
 first_bin is the first bin in the in the L3BIN grid to be treated (inclusive).
  Default is 0.
 last_bin is the last bin in the in the L3BIN grid to be treated (inclusive).
  Default is the last bin of the L3BIN grid.
 outfile is the NetCDF outfile.
 
 description:
 Compute primary productivity and other on the L3BIN grid.
 
 The algorithm is based on Belanger, S., and M. Babin (March 2011), Algorithm
 Theoretical Basis Document (ATBD) for satellite-based Arctic Ocean Primary
 Productivity assessment, version 1.0.
 
 We keep only the pixels from the L3BIN images verifying the following
 conditions:
 - Latitude > 45 N.
 - Ice concentration is documented.
 - Rrs > 0.
 - 0 < Chlorophyll-a concentration < 100 mg Chl-a m^-3.
 - Inherent Optical Properties (IOPs) > 0.
 - If documented, ice concentration is < 10%.
 
 Write the results in a NetCDF file.
 The primary production units are (mg C m^-2 d^-1).
 
 uses: TODO
 
 keywords: primary, productivity
 
 ============================================================================ */

#include <errno.h>
#include <getopt.h>
#include <math.h>
#include <memory.h>
#include <netcdf.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <hdf.h>
#include <mfhdf.h>


#include "Arctic.h"
#include "chl_cota.h"
#include "chl_oc.h"
#include "chl_gsm.h"
#include "Color.h"
#include "fonction_read_i2f_FD_srf.h"
#include "fonction_read_SSMI_TNB.h"
#include "function_read_L3BIN_MODIS_Atmosphere.h"
#include "get_args.h"
#include "qaa4ppMA.h"
#include "qaa4ppSW.h"
#include "rscalc.h"
#include "takuvik.h"

/* ====================== TYPEDEF ====================== */
/* Type of output. */
typedef enum {hdf4_out_false, hdf4_out_true} hdf4_out_type;
typedef enum {netcdf_out_false, netcdf_out_true} netcdf_out_type;

/* ====================== CONSTANTES A MODIFIER ====================== */

/* Affiche une ligne de quatre colonnes. */
#define DEBUG_0 0
/* Affiche un output minimal. */
#define DEBUG_1 1
/* Affiche un output maximal. */
#define DEBUG_2 2
/* Type d'output souhaite. */
#define DEBUG DEBUG_0
/* Type of output. */
#define HDF4_OUT hdf4_out_true
/* Pixel used for debugging. */
//#define IPIX 822     // SeaWiFS, ISCCP, L3BIN
#define IPIX 804650    // SeaWiFS, MODIS-Atmosphere, L3BIN
//#define IPIX 0       // MODISA, ISCCP, L3BIN
//#define IPIX 0       // MODISA, MODIS-Atmosphere, L3BIN
/* Type of output. */
#define NETCDF_OUT netcdf_out_true

/* ====================== CONSTANTES ====================== */
/*
 * Spectral slope of absorption due to colored dissolved and detrital 
 * organic matters in the GSM IOP model.
 * The symbol is S in Maritorena et al. 2002. The symbol is S_CDM in 
 * Ben Mustapha et al. 2012.
 * The units are nm^-1.
 * See Maritorena et al. 2002, table 2 for a description of the model.
 * The value is from SeaDAS 7.0.1 in ocssw/build/src/l2gen/giop.c.
 */
#define ADG_S_CHL_GSM 0.02061
/*
 * Spectral slope of absorption due to colored dissolved and detrital 
 * organic matters in the GSM IOP model.
 * The symbol is S in Maritorena et al. 2002. The symbol is S_CDM in 
 * Ben Mustapha et al. 2012.
 * The units are nm^-1.
 * The value is from Ben Mustapha et al. 2012, p. 547.
 */
#define ADG_S_CHL_GSM_MUSTAPHA 0.018
/*
 * Empirical coefficient used in the computing of aCDOM(412).
 * From Belanger et al. 2008 table 2 for Beaufort Sea.
 */
#define ALPHA -0.514
#define ANNEE 1998 /* Annee */
/* Coefficient A_phi(lambda) derive empiriquement pour la loi de puissance
 * servant a calculer a_phi(lambda), le coefficient d'absorption du 
 * phytoplancton selon Matsuoka et al. 2007 tableau 2.
*/
#define A_PHI_443 0.0288
/*
 * Power law exponent for particulate backscattering coefficient in the 
 * GSM IOP model.
 * The symbol is the greek letter eta.
 * Unitless.
 * See Maritorena et al. 2002, table 2 for a description of the model.
 * The value is from SeaDAS 7.0.1 in ocssw/build/src/l2gen/giop.c.
 */
#define BBP_S_CHL_GSM 1.03373
/*
 * Power law exponent for particulate backscattering coefficient in the 
 * GSM IOP model.
 * The symbol is the greek letter eta.
 * Unitless.
 * The value is from Ben Mustapha et al. 2012, p. 547.
 */
#define BBP_S_CHL_GSM_MUSTAPHA 1.4
/*
 * Empirical coefficient used in the computing of aCDOM(412).
 * From Belanger et al. 2008 table 2 for Beaufort Sea.
 */
#define BETA -0.546
/*
 * Empirical coefficient used in the computing of aCDOM(412).
 * From Belanger et al. 2008 table 2 for Beaufort Sea.
 */
#define CHI 0.480
#define CLASS_NAME_ACDOM_412 "aCDOM_412"
#define CLASS_NAME_CF00 "CF_00UTC"
#define CLASS_NAME_O300 "O3_00UTC"
#define CLASS_NAME_PAR_CLEAR "PAR_clear"
#define CLASS_NAME_PAR_CLOUD "PAR_cloud"
#define CLASS_NAME_TAUCLD00 "TauCld_00UTC"
#define CLASS_NAME_THETAS00 "thetas_00UTC"
#define CLASS_NAME_ZEU "Zeu"
/*
 * Empirical coefficient used in the computing of aCDOM(412).
 * From Belanger et al. 2008 table 2 for Beaufort Sea.
 */
#define DELTA -0.454
/* Coefficient E_phi(lambda) derive empiriquement pour la loi de puissance
 * servant a calculer a_phi(lambda), le coefficient d'absorption du 
 * phytoplancton selon Matsuoka et al. 2007 tableau 2.
*/
#define E_PHI_443 0.82
/* Handle errors by printing an error message and exiting with a
 * non-zero status. */
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}
#define FIELD_NAME_ACDOM_412 "aCDOM_412"
#define FIELD_NAME_CF00 "CF_00UTC"
#define FIELD_NAME_O300 "O3_00UTC"
#define FIELD_NAME_PAR_CLEAR "PAR_clear"
#define FIELD_NAME_PAR_CLOUD "PAR_cloud"
#define FIELD_NAME_TAUCLD00 "TauCld_00UTC"
#define FIELD_NAME_THETAS00 "thetas_00UTC"
#define FIELD_NAME_ZEU "Zeu"
#define FORMAT_A "a%d"
#define FORMAT_BB "bb%d"
#define FORMAT_KD "Kd%d"
#define FORMAT_RRS "Rrs%d"
/* Index de l'heure a 00h UTC. */
#define I00UTC 0
/* emplacement de lambda 490 dans les longueurs d'onde 400-700 */
#define I490_FROM400 18
/* emplacement de lambda 443 dans les longueurs d'onde de SeaWiFS */
#define I443SW 1
/* emplacement de lambda 490 dans les longueurs d'onde de SeaWiFS */
#define I490SW 2
/* emplacement de lambda 550 dans les longueurs d'onde 400-700 */
#define I550_FROM400 30
#define ID400 22 /* emplacement de lambda 400 dans les longeurs d'onde 290-700*/
#define ID440 30 
/* Intervalle de temps entre deux fichiers du ISCCP en heures. */
#define INTERVALLE_TEMPS_ISCCP 3.
#define KD700_PURE_WATER 0.624 /* Kd(700)=0.62 from Morel for pure water */
#define ITIME_TOMORROW_MIDNIGHT 8
/* Une valeur de 168 dans un fichier SSM/I represente de la terre*/
#define LAND 168
#define LENGTH_DAY 2
#define LENGTH_MONTH 2
#define LENGTH_YEAR 4
#define LUT "/Volumes/taku-njall/LUTS/Ed0moins_LUT.dat"
/* Une valeur de 157 dans un fichier SSM/I represente une donnee manquante.*/
#define MISSING 157
#define MODISA 'A'
/* Names of the products in the output files. */
#define NAME_LAT "lat"
#define NAME_LON "lon"
#define NAME_INDLAT "indlat"
#define NAME_INDLON "indlon"
#define NAME_PP "PP"
#define NAME_CHL_GSM "chl_gsm"
/*
 * The chlorophyll-a concentration from Ben Mustapha et al. 2012.
 * Units:  (mg Chl-a m^-3).
 */
#define NAME_CHL_GSM_MUSTAPHA "chl_gsm_mustapha"
#define NAME_CHL_COTA "chl_cota"
#define NAME_CHL_OC "chl_oc"
#define NAME_APHY443 "aphy443"
#define NAME_ICE "Ice"
#define NAME_ZEU "Zeu"
#define NAME_POC "poc"
#define NAME_ACDOM_412 "aCDOM_412"
#define NAME_PAR_CLOUD "PAR_cloud"
#define NAME_PAR_CLEAR "PAR_clear"
#define NAME_CF_00UTC "CF_00UTC"
#define NAME_O3_00UTC "O3_00UTC"
#define NAME_TAUCLD_00UTC "TauCld_00UTC"
#define NAME_CF_MEAN "CF_mean"
#define NAME_O3_MEAN "O3_mean"
#define NAME_TAUCLD_MEAN "TauCld_mean"
#define NAME_KPUR "KPUR"
#define NAME_PP_COTA "PP_cota"
#define NAME_PP_GSM_MUSTAPHA "PP_gsm_mustapha"
#define NAME_PP_OC "PP_oc"
#define NBANDS 6 /* nb lambda SeaWiFS */
#define NB_COL_FILE_IOP 20 /* Nombre de colonnes dans le fichier d'entree */
#define NBDEPTHS 12 /* nombre de profondeurs */
/* Nombre de longitudes a la latitude la plus basse de notre sous-grille EQ
 * du ISCCP au nord de 45 degres Nord. Il s'agit de 100 parce qu'il y a 100
 * longitudes a la latitude 46.25 degres Nord de la grille EQ du ISCCP. */
#define NB_LONG_A_LAT_MIN_ISCCP_EQ 100
#define NBHOURS 9
#define NBTHETAS 3
#define NBWL 83 /* nb longeurs d'onde 290-700nm step 5 */
#define NVIS 61 /* nb longeurs d'onde 400-700nm step 5 */
#define PI 3.14159265
#define POS_412 0 //Position de la longueur d'onde 412 dans SeaWiFS et MODISA.
#define POS_488 2 //Position de la longueur d'onde 488 dans MODISA.
#define POS_490 2 //Position de la longueur d'onde 490 dans SeaWiFS.
#define POS_555 4 //Position de la longueur d'onde 555 dans SeaWiFS et MODISA.
#define POS_A 3 /* Position de la premiere absorption totale dans le fichier d'entree */
#define POS_BB 9 /* Position de la premiere retrodiffusion totale dans le fichier d'entree */
#define POS_CHL 16 /* Position de la chlorophylle dans le fichier d'entree */
#define POS_LAT 1 /* Position de la latitude dans le fichier d'entree */
#define POS_LON 2 /* Position de la longitude dans le fichier d'entree */
#define POS_YEAR_FROM_END 22
#define POS_ZEU 10 /* Position de la profondeur optique 0.01 (la profondeur
		    * euphotique) dans le tableau des profondeurs. */
#define PRODUIT_CF 1
#define PRODUIT_O3 2
#define PRODUIT_TAUCLD 3
#define SEAWIFS 'S'
#define VDATA_NAME_ACDOM_412 "aCDOM_412"
#define VDATA_NAME_CF00 "CF_00UTC"
#define VDATA_NAME_O300 "O3_00UTC"
#define VDATA_NAME_PAR_CLEAR "PAR_clear"
#define VDATA_NAME_PAR_CLOUD "PAR_cloud"
#define VDATA_NAME_TAUCLD00 "TauCld_00UTC"
#define VDATA_NAME_THETAS00 "thetas_00UTC"
#define VDATA_NAME_ZEU "Zeu"

/* ================= FONCTIONS EXTERNES ================= */

/* fonction de : subroutine_locate.f */
extern void locate_(int *gtype,
		    int *ihem,
		    int *itrans,
		    int *i,
		    int *j,
		    float *lat,
		    float *lon);
/* fonction de : subroutine_Ed0moins_at_pixel.f */
extern void ed0moins_(int *jday,
		      float *rtime,
		      float *lat,
		      float *lon,
		      float *o3,
		      float *tcl,
		      float *cf,
		      float Ed_pixel[NBWL],
		      float *ed_lut,
		      float *thetas);
extern void read_ed0moins_lut_(char *lut_fic, float *ed_lut);
extern void sunpos_(int *jday,
		    float *rtime,
		    float *lat,
		    float *lon,
		    float *thetas,
		    float *azim);

/* ====================== PROTOTYPES ====================== */

float calc_aCDOM_412(float at_412,
		     float Rrs_412,
		     float Rrs_490,
		     float Rrs_555);
void calc_aphy(float aphy443, float vis_aphy[NVIS]);
void calc_aphy443(float chla, float* aphy443);
void calc_array1d_iprof_PUR(float array2d_ivis_iprof_E0[NVIS][NBDEPTHS],
                            float vis_aphy[NVIS],
                            float aphy443,
                            float array1d_ivis_lambda[NVIS],
                            float array1d_iprof_PUR[NBDEPTHS]);
void calc_array1d_iprof_Z(float kd550, float array1d_iprof_Z[NBDEPTHS]);
void calc_array2d_h_ivis_Ed_pixel(float Ed_pixel[NBWL],
				  float array2d_h_ivis_Ed_pixel[NBHOURS][NVIS],
				  int h);
void calc_array2d_ivis_iprof_E0(float Kd[NVIS],
				float array1d_iband_a[NBANDS],
				float bb[NBANDS],
				float array1d_iprof_Z[NBDEPTHS],
				float Ed_pixel[NBWL],
				float E0_pixel_z_t[NVIS][NBDEPTHS]);
void calc_E0_pixel_z_t(float Kd[NVIS],
		       float array1d_iband_a[NBANDS],
		       float bb[NBANDS],
		       float Z[NBDEPTHS][NBHOURS],
		       float Ed_pixel[NBWL],
		       int h,
		       float E0_pixel_z_t[NVIS][NBDEPTHS]);
void calc_Ek(float Ek[NBDEPTHS], float meanPUR[NBDEPTHS]);
void calc_isccp_eq(float O3[NBHOURS][NB_PIXELS_ISCCP],
		   float CF[NBHOURS][NB_PIXELS_ISCCP],
		   float TauCld[NBHOURS][NB_PIXELS_ISCCP],
		   int itime,
		   float lat,
		   float lon180,
		   int isccp_grid_lat[NB_PIXELS_ISCCP],
		   float isccp_grid_lon[NB_PIXELS_ISCCP],
		   float* oneO3,
		   float* oneCF,
		   float* oneTauCld);
float calc_isccp_eq_one_product(float product[NBHOURS][NB_PIXELS_ISCCP],
				int itime,
				int iNW,
				int iSW,
				int iNE,
				int iSE,
				float rlat,
				float rlon);
void calc_kdsimon(float array1d_iband_a[NBANDS],
		  float bb[NBANDS],
		  float Kd[NBANDS],
		  float thetas);
void calc_meanPUR(float PUR[NBDEPTHS][NBHOURS], float meanPUR[NBDEPTHS], float photoperiode);
float calc_PAR(float array2d_h_ivis_Ed_pixel[NBHOURS][NVIS]);
float calc_poc(float  Rrs_488, float Rrs_555);
float calc_PP(float chl,
	      float PUR[NBDEPTHS][NBHOURS],
	      float Ek[NBDEPTHS],
	      float Z[NBDEPTHS][NBHOURS]);
void calc_PUR(float array2d_iprof_h_PUR[NBDEPTHS][NBHOURS],
              float array2d_ivis_iprof_E0[NVIS][NBDEPTHS],
              int h,
              float vis_aphy[NVIS],
              float aphy443,
              float array1d_ivis_lambda[NVIS]);
void calc_Z(float Z[NBDEPTHS][NBHOURS], int h, float kd550);
int get_are_pixels_on_current_line(char rrs_type,
                                   short i,
                                   int array1d_irow_begin_MODISA
                                   [SZLAT_MODISA],
                                   int array1d_irow_begin_SeaWiFS
                                   [SZLAT_SEAWIFS]);
void get_array1d_iprof_val_from_array2d_iprof_itime_val(float array2d_iprof_itime_val
							[NBDEPTHS][NBHOURS],
							int itime,
							float array1d_iprof_val
							[NBDEPTHS]);
void get_array1d_itime_KPUR_from_PUR_and_Z(float array2d_iprof_itime_PUR
					   [NBDEPTHS][NBHOURS],
					   float array2d_iprof_itime_Z
					   [NBDEPTHS][NBHOURS],
					   float array1d_itime_KPUR[NBHOURS]);
int get_array2d_iband_itime_Kd(int doy,
			       float lat,
			       float lon,
			       float array1d_iband_a[NBANDS],
			       float array1d_iband_bb[NBANDS],
                               float array2d_iband_itime_Kd[NBANDS][NBHOURS]);
void get_day(char *in_ptr_tnbfilename,
	     int *out_ptr_y,
	     int *out_ptr_m,
	     int * out_ptr_day);
void get_isccp_mean(float CF[NBHOURS][NB_PIXELS_ISCCP],
		    float O3[NBHOURS][NB_PIXELS_ISCCP],
		    float TauCld[NBHOURS][NB_PIXELS_ISCCP],
		    float lat,
		    float lon180,
		    int isccp_grid_lat[NB_PIXELS_ISCCP],
		    float isccp_grid_lon[NB_PIXELS_ISCCP],
		    float* ptr_CFmean,
		    float* ptr_O3mean,
		    float* ptr_TauCldmean);
int get_is_ice_10(float glace);
int get_is_ice_val(float glace);
int get_is_modis_atm_val(float oneCFday,
			 float oneO3day,
			 float oneTauCldday);
int get_is_PUR_val(float array1d_iprof_PUR[NBDEPTHS]);
int get_is_Rrs_val_L3BIN(int ibinimage,
                         float* array1d_ptr_array1d_iband_ibinimage_Rrs
                         [NBANDS]);
float get_itime_at_local_zenith_from_lon(float lonMinus180_180);
float get_KPUR_from_PUR_and_Z(float array1d_iprof_PUR[NBDEPTHS],
			      float array1d_iprof_Z[NBDEPTHS]);
float get_KPUR_from_PUR_and_Z_and_itime(float array2d_iprof_itime_PUR
				    [NBDEPTHS][NBHOURS],
				    float array2d_iprof_itime_Z
				    [NBDEPTHS][NBHOURS],
                                    int itime);
int get_ibinimage(int idx,
		  int ibinimage_firstdocumentedbin_currentrow,
		  int ibinimage_firstdocumentedbin_nextrow,
		  unsigned int* ptr_bin_num);
float get_modis_atm(float array2d_ilat_ilon_atm
		    [NBLAT_MODIS_ATMOSPHERE][NBLON_MODIS_ATMOSPHERE],
		    float lat,
		    float lon);
void interp_Kd(float Kd[NBANDS],
	       float vis_Kd[NVIS],
	       float lo[NBANDS],
	       char rrs_type);
float interp_line(float x1, float x2, float x, float y1, float y2);
int is_array1d_full(char*[], int);
int is_array1d_null(char*[], int);
int is_valide(float chla, float a[NBANDS], float bb[NBANDS]);
float lon180_2_lon360(float lon180);
void parse(char input[500], int* ptr_first_number, float arr[NB_COL_FILE_IOP]);
void print_aphy(float in_tab_aphy[NVIS],
                float array1d_ivis_lambda[NVIS]);
void print_aphy443(float aphy443);
void print_array2d_h_ivis_Ed_pixel(float array2d_h_ivis_Ed_pixel
				   [NBHOURS][NVIS]);
void print_calc_isccp_eq(int itime,
			 float lat,
			 float lon180,
			 float lon, // Longitude de 0 a 360.
			 int latN,
			 int latS,
			 int iNW,
			 int iNE,
			 int iSW,
			 int iSE,
			 int iLATN,
			 int iLATS,
			 int last,
			 float rlat,
			 float rlonS,
			 float rlonN,
			 float rlon,
			 float oneO3,
			 float oneCF,
			 float oneTauCld);
void print_calc_isccp_eq_one_product(float productNW,
				     float productSW,
				     float xlatW,
				     float productNE,
				     float productSE,
				     float xlatE);
void print_eclairement(float Ed_pixel[NBWL],
		       float E0_pixel_z_t[NVIS][NBDEPTHS]);
void print_final(int is_clim,
		 int ice_val,
		 int ice_not_val,
		 int l3b_val,
		 int Rrs_val,
		 int oc_val,
		 int ice_10,
		 int clim_val,
		 int pp_val);
void print_final_noclim(int l3b_val,
			int Rrs_val,
			int oc_val,
			int ice_10,
			int pp_val);
void print_glace(int gtype,
		 int ihem,
		 int itrans,
		 int iTNB,
		 int jTNB,
		 float xlat,
		 float xlon,
		 float glace);
void print_glace_exhaustive(float array2d_jTNB_iTNB_ice[NBLINE][NBCOL]);
void print_iops(float in_tab_aphy[NVIS],
                float in_tab_atotal[NVIS],
                float in_tab_bbtotal[NVIS],
                float array1d_ivis_lambda[NVIS]);
void print_isccp_interpolation(float hh, float o3, float tauCld, float cf);
void print_kd(float a[NBANDS],
	      float bb[NBANDS],
	      float thetas,
	      float Kd[NBANDS],
	      float lo[NBANDS]);
void print_KPUR(float array1d_itime_KPUR[NBHOURS]);
void print_meanpur_ek(float meanPUR[NBDEPTHS], float Ek[NBDEPTHS]);
void print_photoperiode(float photoperiode);
void print_pp(float pp);
void print_pp_temp(float pp);
void print_pur(int itime, float PUR[NBDEPTHS][NBHOURS]);
void print_result(int no_pixel_grille_isin, float lat, float lon, float pp);
void print_thetas(int doy, float hh, float xlat, float xlon, float thetas);
void print_thetas45();
void print_vis_Kd(float vis_Kd[NVIS],
                  float array1d_ivis_lambda[NVIS]);
void print_Z(float Z[NBDEPTHS][NBHOURS], int h);
void read_l3b(char *filename,
              char rrs_type,
              int array1d_irow_start_num_MODISA[SZLAT_MODISA],
              int array1d_irow_start_num_SeaWiFS[SZLAT_SEAWIFS],
              int array1d_irow_begin_MODISA[SZLAT_MODISA],
              int array1d_irow_begin_SeaWiFS[SZLAT_SEAWIFS],
              int array1d_irow_extent_MODISA[SZLAT_MODISA],
              int array1d_irow_extent_SeaWiFS[SZLAT_SEAWIFS],
              int array1d_irow_max_MODISA[SZLAT_MODISA],
              int array1d_irow_max_SeaWiFS[SZLAT_SEAWIFS],
              unsigned int** ptr_ptr_bin_num,
              float* array1d_ptr_array1d_iband_ibinimage_Rrs[NBANDS]);
void reproject_isccp_grid_to_sat_grid(
				      float O3[NBHOURS][NB_PIXELS_ISCCP],
				      float CF[NBHOURS][NB_PIXELS_ISCCP],
				      float TauCld[NBHOURS][NB_PIXELS_ISCCP],
				      float *array1d_ptr_array1d_itime_ipix_O3sat[NBHOURS],
				      float *array1d_ptr_array1d_itime_ipix_CFsat[NBHOURS],
				      float *array1d_ptr_array1d_itime_ipix_TauCldsat[NBHOURS],
				      int szlat,
                                      int nbpix);
void write_hdf(char *fic_out,
               int nb_valid,
               float lo[NBANDS],
               float *ptr_array1d_ipix_lat,
               float *ptr_array1d_ipix_lon,
               short *ptr_array1d_ipix_indlat,
               short *ptr_array1d_ipix_indlon,
               float *array1d_ptr_array1d_iband_ipix_Rrs[NBANDS],
               float *array1d_ptr_array1d_iband_ipix_a[NBANDS],
               float *array1d_ptr_array1d_iband_ipix_bb[NBANDS],
               float *array1d_ptr_array1d_iband_ipix_Kd[NBANDS],
               float *ptr_array1d_ipix_Zeu_bar,
               float *ptr_array1d_ipix_poc,
               float *ptr_array1d_ipix_aCDOM_412,
               float *ptr_array1d_ipix_PAR_cloud,
               float *ptr_array1d_ipix_PAR_clear,
               float *ptr_array1d_ipix_CF00,
               float *ptr_array1d_ipix_O300,
               float *ptr_array1d_ipix_tauCld00,
               float *ptr_array1d_ipix_thetas00,
               float *ptr_array1d_ipix_CFmean,
               float *ptr_array1d_ipix_O3mean,
               float *ptr_array1d_ipix_tauCldmean,
               float *ptr_array1d_ipix_KPUR,
               float *ptr_array1d_ipix_chl_cota,
               float *ptr_array1d_ipix_chl_oc,
               float *ptr_array1d_ipix_chl_gsm_mustapha,
               float *ptr_array1d_ipix_PP_cota,
               float *ptr_array1d_ipix_PP_gsm_mustapha,
               float *ptr_array1d_ipix_PP_oc,
               float *ptr_array1d_ipix_aphy443,
               float *ptr_array1d_ipix_chl_gsm,
               float *ptr_array1d_ipix_PP_chl_gsm,
               float *ptr_array1d_ipix_ice);
void write_netcdf(char *fic_out,
                  int nb_valid,
                  float lo[NBANDS],
                  float *ptr_array1d_ipix_lat,
                  float *ptr_array1d_ipix_lon,
                  short *ptr_array1d_ipix_indlat,
                  short *ptr_array1d_ipix_indlon,
                  float *array1d_ptr_array1d_iband_ipix_Rrs[NBANDS],
                  float *array1d_ptr_array1d_iband_ipix_a[NBANDS],
                  float *array1d_ptr_array1d_iband_ipix_bb[NBANDS],
                  float *array1d_ptr_array1d_iband_ipix_Kd00[NBANDS],
                  float *ptr_array1d_ipix_Zeu_bar,
                  float *ptr_array1d_ipix_poc,
                  float *ptr_array1d_ipix_aCDOM_412,
                  float *ptr_array1d_ipix_PAR_cloud,
                  float *ptr_array1d_ipix_PAR_clear,
                  float *ptr_array1d_ipix_CF00,
                  float *ptr_array1d_ipix_O300,
                  float *ptr_array1d_ipix_tauCld00,
                  float *ptr_array1d_ipix_thetas00,
                  float *ptr_array1d_ipix_CFmean,
                  float *ptr_array1d_ipix_O3mean,
                  float *ptr_array1d_ipix_tauCldmean,
                  float *ptr_array1d_ipix_KPUR,
                  float *ptr_array1d_ipix_chl_cota,
                  float *ptr_array1d_ipix_chl_oc,
                  float *ptr_array1d_ipix_chl_gsm_mustapha,
                  float *ptr_array1d_ipix_PP_cota,
                  float *ptr_array1d_ipix_PP_gsm_mustapha,
                  float *ptr_array1d_ipix_PP_oc,
                  float *ptr_array1d_ipix_aphy443,
                  float *ptr_array1d_ipix_chl_gsm,
                  float *ptr_array1d_ipix_PP_chl_gsm,
                  float *ptr_array1d_ipix_ice);

/* ============== GLOBAL VARIABLES (ONLY ARRAYS OF CONSTANTS) ============== */

/*
 * Array of dimension 9.
 * The first dimension is the hours UTC from 0 to 24 by step of 3 h.
 * The values are the hours UTC.
 * The units are hours UTC.
 */
float ARRAY1D_ITIME_HOUR[NBHOURS]={0., 3., 6., 9., 12., 15., 18., 21., 0.};
/*
 * Array of dimension NBANDS = 6.					       
 * The first dimension are the wavelengths of the bands of the satellite.    
 * The values are the band of MODIS-AQUA.
 * The units are nm.                                                      
 */
float ARRAY1D_IBAND_BANDMODISA[NBANDS] = {412., 443., 488., 531., 555., 667.};
/*
 * Array of dimension NBANDS = 6.					       
 * The first dimension are the wavelengths of the bands of the satellite.    
 * The values are the band of SeaWiFS.
 * The units are nm.                                                      
 */
float ARRAY1D_IBAND_BANDSEAWIFS[NBANDS] = {412., 443., 490., 510., 555., 670.};
/*
 * Array of dimensions 12.
 * The first dimension are the optical depths. 			     
 * The values are the optical depths.
 * Unitless.
 */
float ARRAY1D_IDEPTH_OPTICALDEPTHS[NBDEPTHS]
= {1.,.9, .8, .7, .6, .5, .4, .3, .2, .1, .01, .001};

/* ====================== FONCTIONS ====================== */

/* --------------------------------------------------------------
calcul de aphy a partir de la chloro (Matsuoka et al 2007)
*/
void calc_aphy(float aphy443, float vis_aphy[NVIS])
{
  float A[NVIS] = {0.0209, 0.0232, 0.0252, 0.0266, 0.0275, 0.0281, 0.0291, 0.0304, 0.0306, 0.0291, 0.0271, 0.0257, 0.0253, 0.0249, 0.0242, 0.0228, 0.0241, 0.0202, 0.0189, 0.0172, 0.0156, 0.0141, 0.0126, 0.0113, 0.0103, 0.0093, 0.0085, 0.0077, 0.0070, 0.0064, 0.0057, 0.0049, 0.0043, 0.0039, 0.0036, 0.0036, 0.0035, 0.0036, 0.0037, 0.0036, 0.0034, 0.0033, 0.0035, 0.0038, 0.0041, 0.0042, 0.0045, 0.0048, 0.0050, 0.0052, 0.0054, 0.0058, 0.0072, 0.0100, 0.0127, 0.0140, 0.0128, 0.0093, 0.0054, 0.0029, 0.0018};
  float B[NVIS] = {0.881, 0.898, 0.902, 0.891, 0.881, 0.878, 0.866, 0.863, 0.860, 0.857, 0.856, 0.855, 0.858, 0.854, 0.854, 0.847, 0.843, 0.845, 0.853, 0.858, 0.866, 0.874, 0.877, 0.877, 0.888, 0.894, 0.909, 0.924, 0.950, 0.965, 0.975, 0.970, 0.972, 0.947, 0.932, 0.946, 0.894, 0.904, 0.903, 0.890, 0.880, 0.867, 0.891, 0.895, 0.912, 0.907, 0.909, 0.898, 0.898, 0.870, 0.883, 0.903, 0.920, 0.930, 0.922, 0.921, 0.932, 0.951, 0.956, 0.969, 0.993};
  int i;  
  for(i=0; i<NVIS; i++)
    vis_aphy[i] 
      = A[i] * (float) pow((double) aphy443 / A_PHI_443, B[i] / E_PHI_443);
}

/*
 * Calcul des profondeurs.
 * array1d_iprof_Z: Les profondeurs geometriques (en m).
 * kd550: Le coefficient d'attenuation a 550 nm (en m^-1).
 */
void calc_array1d_iprof_Z(float kd550, float array1d_iprof_Z[NBDEPTHS]){
  int j;
  for(j = 0; j < NBDEPTHS; j++){
    array1d_iprof_Z[j]
      = (- (float) log((double) ARRAY1D_IDEPTH_OPTICALDEPTHS[j]) / kd550);
  }
}

/*
 * IN
 * array2d_ivis_iprof_E0: 
 *  Array of dimension 61 * 12.
 *  The first dimension is the wavelenght from 400 nm to 700 nm by step of 5
 *  nm.
 *  The second dimension is the optical depths:
 *  {1.,.9, .8, .7, .6, .5, .4, .3, .2, .1, .01, .001}
 *  The values are the scalar irradiances for a pixel on a given day.
 *  The units are uEinstein*m^-2*s^-1*nm^-1.
 * vis_aphy:
 *  Array of dimension 61.
 *  The first dimension is the wavelenght from 400 nm to 700 nm by step of 5
 *  nm.
 *  The values are the phytoplancton absorption coefficients.
 *  The units are m^-1.
 * aphy443:
 *  A scalar.
 *  The phytoplancton absorption coefficient at 443 nm.
 *  The units are m^-1.
 * array1d_ivis_lambda:
 *  Array of dimension 61.
 *  The first dimension is the index of the wavelength.
 *  The values are the wavelengths from 400 to 700 by step 5.
 *  The units are nm.
 * OUT
 * array1d_iprof_PUR:
 *  The values are the photosynthetically usable radiation (PUR) for a pixel
 *  on a given day.
 *  The units are uEinstein*m^-2*s^-1.
 * Compute the PUR (equation 6).
 * Integration with the method of the trapezes.
 */
void calc_array1d_iprof_PUR(float array2d_ivis_iprof_E0[NVIS][NBDEPTHS],
                            float vis_aphy[NVIS],
                            float aphy443,
                            float array1d_ivis_lambda[NVIS],
                            float array1d_iprof_PUR[NBDEPTHS]){
  int iprof;
  int ivis;
  float a1;
  float a2;
  for(iprof = 0; iprof < NBDEPTHS; iprof++){
    array1d_iprof_PUR[iprof] = 0.;
    for(ivis = 0; ivis < (NVIS - 1); ivis++){
      a1 = 
	vis_aphy[ivis] / aphy443 * array2d_ivis_iprof_E0[ivis][iprof];
      a2 =
	vis_aphy[ivis + 1] / aphy443 * array2d_ivis_iprof_E0[ivis + 1][iprof];
      array1d_iprof_PUR[iprof] +=
	( (array1d_ivis_lambda[ivis + 1] - array1d_ivis_lambda[ivis])
   * (a1 + a2) / 2.);
    }
  }
}

/*
 * Calcul de aphy443.
 * chla: La chlorophylle a.
 * aphy443: Le coefficient d'absorption du phytoplancton a 443 nm.
 */
void calc_aphy443(float chla, float* aphy443){
  *aphy443 = A_PHI_443 * (float)pow((double)chla, E_PHI_443);
}

/*
 * Ed_pixel:
 *  Array of dimension 83.						 
 *  The first dimension is the index of the wavelength from 290 nm to 700 nm 
 *  by step of 5 nm.
 *  The values are the downward irradiances just below the surface water 
 *  (Ed0moins) for the pixel.						 
 *  The units are uEinstein*m^-2*s^-1*nm^-1.                             
 * array2d_h_ivis_Ed_pixel:
 *  Array of dimensions 9 * 61.						  
 *  The first dimension is the hours from 0 to 24 by step of 3 h.
 *  The second dimension is the wavelenght from 400 nm to 700 nm by step of 5
 *  nm.
 *  The values are the downward irradiances just below the surface water 
 *  (Ed0moins) for the pixel.
 *  The units are uEinstein*m^-2*s^-1*nm^-1.                             
 * h:
 *  The index of the hour for the hours from 0 to 24 by step of 3 h.
 */
void calc_array2d_h_ivis_Ed_pixel(float Ed_pixel[NBWL],
				  float array2d_h_ivis_Ed_pixel[NBHOURS][NVIS],
				  int h){
  // Index on the wavelengths from 400 nm to 700 nm by step of 5 nm.
  int ivis;
  for(ivis = 0; ivis < NVIS; ivis++){
    array2d_h_ivis_Ed_pixel[h][ivis] = Ed_pixel[ID400 + ivis];
  }
}

/*
 * Calcul de E0_pixel_z_t (equation 7)
 *
 * IN
 * array1d_iband_a:
 *  Array of dimension NBANDS = 6.					       
 *  The first dimension are the wavelengths of the bands of the satellite.    
 *  The values are the total absorption coefficients (a or a_t).	       
 *  The units are m^-1.                                                       
 * OUT
 */

void calc_array2d_ivis_iprof_E0(float Kd[NVIS],
				float array1d_iband_a[NBANDS],
				float bb[NBANDS],
				float array1d_iprof_Z[NBDEPTHS],
				float Ed_pixel[NBWL],
				float E0_pixel_z_t[NVIS][NBDEPTHS]){
  int i,j;
  float one_z;
  
  for(i=0; i<NVIS; i++){
    for(j=0; j<NBDEPTHS; j++){
      one_z = array1d_iprof_Z[j];
      E0_pixel_z_t[i][j]
	= (Ed_pixel[i+ID400] * (float)exp((double) -Kd[i]*one_z))
	* (Kd[I490_FROM400]/(array1d_iband_a[I490SW] + bb[I490SW]));
    }
  }
}

/* --------------------------------------------------------------
 Calcul de E0_pixel_z_t (equation 7)
 * IN
 * array1d_iband_a:
 *  Array of dimension NBANDS = 6.					       
 *  The first dimension are the wavelengths of the bands of the satellite.    
 *  The values are the total absorption coefficients (a or a_t).	       
 *  The units are m^-1.                                                       
 * OUT
*/
void calc_E0_pixel_z_t(float Kd[NVIS],
		       float array1d_iband_a[NBANDS],
		       float bb[NBANDS],
		       float Z[NBDEPTHS][NBHOURS],
		       float Ed_pixel[NBWL],
		       int h,
		       float E0_pixel_z_t[NVIS][NBDEPTHS])
{
  int iprof;
  float array1d_iprof_Z[NBDEPTHS];
  for(iprof = 0; iprof < NBDEPTHS; iprof++){
    array1d_iprof_Z[iprof] = Z[iprof][h];
  }
  calc_array2d_ivis_iprof_E0(Kd,
			     array1d_iband_a,
			     bb,
			     array1d_iprof_Z,
			     Ed_pixel,
			     E0_pixel_z_t);
}

/* --------------------------------------------------------------
 Calcul de Ek (equation 8 & 9)
 */
void calc_Ek(float Ek[NBDEPTHS], float meanPUR[NBDEPTHS])
{
  int i;
  double Ek_max=80., B;
  
  B = exp( 1.089 - 2.12*log10(Ek_max));
  for(i=0; i<NBDEPTHS; i++){
    Ek[i] = (float)( Ek_max / (1+2*exp(-B * meanPUR[i])) );
  }
}

/*
 * O3: Les valeurs d'ozone lues dans les fichiers du ISCCP sur notre 
 *     sous-grille de la grille EQ du ISCCP au nord de 45 degres Nord.
 *     Il y a 966 valeurs pour chaque heure.
 * CF: Les valeurs de fraction nuageuse lues dans les fichiers du ISCCP sur 
 *     notre sous-grille de la grille EQ du ISCCP au nord de 45 degres Nord.
 *     Il y a 966 valeurs pour chaque heure.
 * TauCld: Les valeurs de l'epaisseur optique lues dans les fichiers du ISCCP 
 *         sur notre sous-grille de la grille EQ du ISCCP au nord de 45 degres 
 *         Nord.
 *         Il y a 966 valeurs pour chaque heure.
 * itime: Indice de l'heure dans ARRAY1D_ITIME_HOUR.
 * lat: Latitude (de 45 a 90).
 * lon180: Longitude (de -180 a 180).
 * isccp_grid_lat: Les latitudes pour chaque pixel de notre sous-grille de la 
 *                 grille EQ du ISCCP au nord de 45 degres Nord.
 *      	   Les latitudes sont multipliees par un facteur de 100.
 *      	   Les latitudes vont de 4500 a 9000.
 * isccp_grid_lon: Les longitudes pour chaque pixel de notre sous-grille de 
 *                 la grille EQ du ISCCP au nord de 45 degres Nord.
 *                 Les longitudes vont de 0 a 360.
 * oneO3: Une valeur d'ozone.
 * oneCF: Une valeur de fraction nuageuse.
 * oneTauCld: Une valeur d'epaisseur optique.
 */
void calc_isccp_eq(float O3[NBHOURS][NB_PIXELS_ISCCP],
		   float CF[NBHOURS][NB_PIXELS_ISCCP],
		   float TauCld[NBHOURS][NB_PIXELS_ISCCP],
		   int itime,
		   float lat,
		   float lon180,
		   int isccp_grid_lat[NB_PIXELS_ISCCP],
		   float isccp_grid_lon[NB_PIXELS_ISCCP],
		   float* oneO3,
		   float* oneCF,
		   float* oneTauCld){
  float lon; // Longitude de 0 a 360.
  int i; // Indice sur notre sous-grille de la grille EQ du ISCCP.
  // Soient les boites NW, NE, SW et SE au nord-ouest, au nord-est, 
  // au sud-ouest et au sud-est du pixel sur la grille EQ du ISCCP.
  int latN; // Latitude des boites NW et NE.
            // Les latitudes sont multipliees par 100.
  int latS; // Latitude des boites SW et SE.
            // Les latitudes sont multipliees par 100.
  int iNW = 0.; // Indice de la boite NW.
  int iNE = 0.; // Indice de la boite NE.
  int iSW; // Indice de la boite SW.
  int iSE; // Indice de la boite SE.
  int iLATN; // Indice du premier pixel de la rangee nord (contenant NW et NE).
  int iLATS; // Indice du deuxieme pixel de la rangee sud (contenant SW et SE).
  int last; // Indice du dernier pixel de la rangee nord (contenant NW et NE).
  float rlat; // Proportion indiquant la distance du pixel d'avec la latitude 
              // de la rangee nord.
              // 0 signifie sur la rangee nord.
              // 1 signifie sur la rangee sud.
  float rlonS; // Proportion indiquant la distance du pixel d'avec la
               // longitude de la boite SW.
               // 0 signifie sur la longitude de la boite SW.
               // 1 signifie sur la longitude de la boite SE.
  float rlonN = 0.; // Proportion indiquant la distance du pixel d'avec la
               // longitude de la boite NW.
               // 0 signifie sur la longitude de la boite NW.
               // 1 signifie sur la longitude de la boite NE.
  float rlon; // Moyenne ponderee par rlat de rlonS et rlonN.
  lon = lon180_2_lon360(lon180);
  if(lat > 88.75){
    lat = 88.75;
  }
  // Boucle pour calculer rlat.
  // Cas general.
  // lat >= 46.25 degres Nord.
  if(lat > isccp_grid_lat[0] / 100.){
    for(i = 0; i < NB_PIXELS_ISCCP; i++){
      if(lat <= isccp_grid_lat[i] / 100.){
	latN = isccp_grid_lat[i];
	latS = latN - 250;
	rlat = (latN / 100. - lat) / 2.5;
	iLATN = i;
	break;
      }
    }
  }
  // Cas special.
  // Le pixel est sur le bord de la limite sud de la grille.
  else{
    latN = isccp_grid_lat[NB_LONG_A_LAT_MIN_ISCCP_EQ]; // 4875
    latS = isccp_grid_lat[0]; // 4625
    rlat = (latN / 100. - lat) / 2.5;
    iLATN = NB_LONG_A_LAT_MIN_ISCCP_EQ; // 100
  }
  iLATS = 0;
  // Boucle pour calculer rlonS.
  for(i = 0; i < NB_PIXELS_ISCCP - 1; i++){
    // Pour trouver iLATS.
    if(latS == isccp_grid_lat[i]
       && iLATS == 0){
      iLATS = i;
    }
    // Cas general.
    if(latS == isccp_grid_lat[i]
       && lon > isccp_grid_lon[i]
       && lon <= isccp_grid_lon[i + 1]){
      iSW = i;
      iSE = i + 1;
      rlonS = (lon - isccp_grid_lon[i])
	/ (isccp_grid_lon[i + 1] - isccp_grid_lon[i]);
      break;
    }
    // Cas special.
    // Le pixel est sur le bord de la limite ouest de la grille.
    if(latS == isccp_grid_lat[i]
       && lon <= isccp_grid_lon[i]){
      iSE = i; // C'est le premier pixel de la rangee sud.
      iSW = iLATN - 1; // C'est le dernier pixel de la rangee sud.
      rlonS = (lon + (360. - isccp_grid_lon[iSW]))
	/ (isccp_grid_lon[i] + (360. - isccp_grid_lon[iSW]));
      break;
    }
    // Cas special.
    // Le pixel est sur le bord de la limite est de la grille.
    if(latS == isccp_grid_lat[i]
       && lon > isccp_grid_lon[iLATN - 1]){
      iSW = iLATN - 1;
      iSE = iLATS;
      rlonS = (lon - isccp_grid_lon[iSW])
	/ (isccp_grid_lon[iSE] + (360. - isccp_grid_lon[iSW]));
      break;
    }
  }
  // Pour trouver last.
  for(i = 0; i < NB_PIXELS_ISCCP; i++){
    if(latN + 250 == isccp_grid_lat[i]){
      last = i - 1;
      break;
    }
  }
  // Boucle pour calculer rlonN.
  for(i = 0; i < NB_PIXELS_ISCCP -1; i++){
    // Cas general.
    if(latN == isccp_grid_lat[i]
       && lon > isccp_grid_lon[i]
       && lon <= isccp_grid_lon[i + 1]){
      iNW = i;
      iNE = i + 1;
      rlonN = (lon - isccp_grid_lon[i])
	/ (isccp_grid_lon[i + 1] - isccp_grid_lon[i]);
      break;
    }
    // Cas special.
    // Le pixel est sur le bord de la limite ouest de la grille.
    if(latN == isccp_grid_lat[i]
       && lon <= isccp_grid_lon[i]){
      iNE = i; // C'est le premier pixel de la rangee nord.
      iNW = last; // C'est le dernier pixel de la rangee nord.
      rlonN = (lon + (360. - isccp_grid_lon[iNW]))
	/ (isccp_grid_lon[i] + (360. - isccp_grid_lon[iNW]));
      break;
    }
    // Cas special.
    // Le pixel est sur le bord de la limite est de la grille.
    if(latN == isccp_grid_lat[i]
       && lon > isccp_grid_lon[last]){
      iNW = last;
      iNE = iLATN;
      rlonN = (lon - isccp_grid_lon[iNW])
	/ (isccp_grid_lon[iNE] + (360. - isccp_grid_lon[iNW]));
      break;
    }
  }
  rlon = rlonN * (1. - rlat) + rlonS * rlat;
  *oneO3 = calc_isccp_eq_one_product(O3, itime, iNW, iSW, iNE, iSE, rlat, rlon);
  *oneCF = calc_isccp_eq_one_product(CF, itime, iNW, iSW, iNE, iSE, rlat, rlon);
  *oneTauCld = calc_isccp_eq_one_product(TauCld, itime, iNW, iSW, iNE, iSE,
					 rlat, rlon);
  
  /*print_calc_isccp_eq(itime,
		      lat,
		      lon180,
		      lon,
		      latN,
		      latS,
		      iNW,
		      iNE,
		      iSW,
		      iSE,
		      iLATN,
		      iLATS,
		      last,
		      rlat,
		      rlonS,
		      rlonN,
		      rlon,
		      *oneO3,
		      *oneCF,
		      *oneTauCld);*/
		      
}

/*
 * product: Les valeurs de produit atmospherique lues dans les fichiers du 
 *          ISCCP sur notre sous-grille de la grille EQ du ISCCP au nord de 
 *          45 degres Nord.
 * itime: Indice de l'heure dans ARRAY1D_ITIME_HOUR.
 * Soient les boites NW, NE, SW et SE au nord-ouest, au nord-est, 
 * au sud-ouest et au sud-est du pixel sur la grille EQ du ISCCP.
 * iNW: Indice de la boite NW.
 * iNE: Indice de la boite NE.
 * iSW: Indice de la boite SW.
 * iSE: Indice de la boite SE.
 * rlat: Proportion indiquant la distance du pixel d'avec la latitude 
 *       de la rangee contenant NW et NE.
 *       0 signifie sur la rangee nord.
 *       1 signifie sur la rangee sud.
 * rlon: Proportion indiquant la distance du pixel d'avec la
 *       longitude de la boite NW.
 *       0 signifie sur la longitude de la boite NW.
 *       1 signifie sur la longitude de la boite NE.
 */
float calc_isccp_eq_one_product(float product[NBHOURS][NB_PIXELS_ISCCP],
				int itime,
				int iNW,
				int iSW,
				int iNE,
				int iSE,
				float rlat,
				float rlon){
  float xlatW; // Moyenne ponderee du produit dans les boites NW et SW.
  float xlatE; // Moyenne ponderee du produit dans les boites NE et SE.
  float one_product; // Le produit calculee suite a l'interpolation lineaire
                     // lineaire en deux dimensions.
  xlatW = product[itime][iNW] * (1. - rlat) + product[itime][iSW] * rlat;
  xlatE = product[itime][iNE] * (1. - rlat) + product[itime][iSE] * rlat;
  one_product = xlatW * (1. - rlon) + xlatE * rlon;
  /*print_calc_isccp_eq_one_product(product[itime][iNW],
				  product[itime][iSW],
				  xlatW,
				  product[itime][iNE],
				  product[itime][iSE],
				  xlatE);*/
  // Cela peut arriver lorsque le pixel est sur le bord de la limite sud de la 
  // grille parce que l'on fait alors une extrapolation lineaire.
  if(one_product < 0){
    one_product = 0.;
  }
  return one_product;
}
/* -------------------------------------------------------------- */

/* Computing of Kd (equation 4 of the ATBD).
 *
 * IN
 * array1d_iband_a:
 *  Array of dimension NBANDS = 6.					       
 *  The first dimension are the wavelengths of the bands of the satellite.    
 *  The values are the total absorption coefficients (a or a_t).	       
 *  The units are m^-1.                                                       
 * OUT
 */
void calc_kdsimon(float array1d_iband_a[NBANDS],
		  float bb[NBANDS],
		  float Kd[NBANDS],
		  float thetas)
{
  short j,ithetas;
  float m[4][NBTHETAS] = { {1.044, 1.108, 1.32},
			   {4.173, 4.245, 4.120},
			   {0.530, 0.526, 0.504},
			   {11.157, 10.942, 10.304}
  };
  // m3(10) de Simon
  float thetas_lut[3] = {10.,30.,60.};
  float f,minter[4];
  int iband;
  if(thetas <= thetas_lut[0]){ 
    for(iband = 0; iband < NBANDS; iband++){
      Kd[iband] = m[0][0] * array1d_iband_a[iband]
	+ m[1][0]
	* ( 1. - m[2][0]
	    * (float)exp( (double)-m[3][0] * array1d_iband_a[iband]) )
	* bb[iband];
    }
  }
  if(thetas >= thetas_lut[2]){
    for(iband = 0; iband < NBANDS; iband++){
      Kd[iband] = m[0][2] * array1d_iband_a[iband]
	+ m[1][2]
	* ( 1. - m[2][2]
	    * (float)exp( (double)-m[3][2] * array1d_iband_a[iband]) )
	* bb[iband];
    }
  }
  if(thetas > 10. && thetas<60){
    /* Interpolation des coefficients */
    for (j=0; j<=1; j++){
      if(thetas >= thetas_lut[j] && thetas < thetas_lut[j+1]){
	ithetas = j;
      }
    }
    f = (thetas - thetas_lut[ithetas])
      / (thetas_lut[ithetas+1] - thetas_lut[ithetas]);
    for(j=0; j<=3; j++){
      minter[j] = (1. - f) * m[j][ithetas] + f * m[j][ithetas+1];
    }
    for (iband = 0; iband < NBANDS; iband++){
      Kd[iband] = minter[0] * array1d_iband_a[iband]
	+ minter[1]
	* ( 1. - minter[2]
	    * (float)exp( (double)-minter[3] * array1d_iband_a[iband]) )
	* bb[iband];
    }
  }
}

/* --------------------------------------------------------------
 Calcul du PUR moyen (equation 10)
 integration par la methode des trapezes
*/
void calc_meanPUR(float PUR[NBDEPTHS][NBHOURS], float meanPUR[NBDEPTHS], float photoperiode)
{
  int i, j;
  
  for(i=0; i<NBDEPTHS; i++){
    meanPUR[i] = 0.;
    for(j=0; j<NBHOURS-1; j++){
      meanPUR[i] 
	+= ( INTERVALLE_TEMPS_ISCCP * 3600. * (PUR[i][j+1] + PUR[i][j]) / 2. );
    }
    meanPUR[i] /= photoperiode;
  }
}

/* ------------------------------------------------------------------
 * array2d_h_ivis_Ed_pixel:
 *  Array of dimensions 9 * 61.						  
 *  The first dimension is the hours from 0 to 24 by step of 3 h.
 *  The second dimension is the wavelenght from 400 nm to 700 nm by step of 5
 *  nm.
 *  The values are the downward irradiances just below the surface water 
 *  (Ed0moins) for the pixel.
 *  The units are uEinstein*m^-2*s^-1*nm^-1.                             
 * Return PAR. PAR is the photosynthetically available irradiance (PAR) for a 
 * pixel on a given day.
 * The units are uEinstein*m^-2.
 * The method used is the trapeze method.
 */
float calc_PAR(float array2d_h_ivis_Ed_pixel[NBHOURS][NVIS]){
  /*
   * Array of dimension 9.
   * The first dimension is the hours.
   * The values are the integration of the different downward irradiances 
   * just below the surface water on the wavelengths.
   * The units are uEinstein*m^-2*s^-1.
   * The method used for the integration is the trapeze method.
   */
  float array1d_h_PAR_s[NBHOURS];
  int h; // The index of the hour.
  /*
   * Work variable when computing PAR_s.
   * The units of hauteur_Ed_pixel are uEinstein*m^-2*s^-1nm^-1.
   */
  float hauteur_Ed_pixel;
  /*
   * Work variable when computing PAR.
   * The units of hauteur_PAR_s are uEinstein*m^-1*s^-1.
   */
  float hauteur_PAR_s;
  int ivis; // The index of the wavelength.
  /*
   * PAR is the photosynthetically available irradiance (PAR) for a pixel on a 
   * given day.
   * The units are uEinstein*m^-2.
   */
  float PAR;
  PAR = 0.;
  for(h = 0; h < NBHOURS; h++){
    array1d_h_PAR_s[h] = 0.;
    // Integration on the wavelengths.
    for(ivis = 0; ivis <= NVIS - 2; ivis++){
      hauteur_Ed_pixel
	= (array2d_h_ivis_Ed_pixel[h][ivis]
	   + array2d_h_ivis_Ed_pixel[h][ivis + 1]
	   ) / 2.;
      array1d_h_PAR_s[h] += hauteur_Ed_pixel * 5.;
    } // End of the integration on the wavelengths.
  }
  // Integration on the hours.
  for(h = 0; h <= NBHOURS - 2; h++){
    hauteur_PAR_s = (array1d_h_PAR_s[h]
		     + array1d_h_PAR_s[h + 1]
		     ) / 2.;
    PAR += (hauteur_PAR_s * 3. * 3600.);
  }// End of the integration on the hours.
  return PAR;
}

/* ------------------------------------------------------------------
 * Rrs_488: Rrs at 488 nm.
 * Rrs_555: Rrs at 555 nm.
 * Return the Particulate organic carbon (poc).
 */
float calc_poc(float  Rrs_488, float Rrs_555){
  float poc = (float)pow(196.98 * (Rrs_488 / Rrs_555), -1.312);
  return poc;
}

/* --------------------------------------------------------------
 Calcul de la production primaire (equation 5)
*/
float calc_PP(float chl, float PUR[NBDEPTHS][NBHOURS], float Ek[NBDEPTHS], float Z[NBDEPTHS][NBHOURS])
{
  int i, j;
  /* P_Bmax avec les unites mgC (mg Chl-a)^-1 h^-1 */
  float P_Bmax_unit_h = 2.0;
  /* P_Bmax avec les unites mgC (mg Chl-a)^-1 d^-1 */
  float P_Bmax_unit_d;
  float a1, a2, tmp[NBHOURS], dintegrale=0.;

  //P_Bmax_unit_d = P_Bmax_unit_h * 24;
  P_Bmax_unit_d = P_Bmax_unit_h;
  for(i=0; i<NBHOURS-1; i++){
    if(i==0){
      tmp[i] = 0.;
      for(j=0; j<NBDEPTHS-1; j++){
        a1 = 1. - exp((-PUR[j][i]) / Ek[j] );
        a2 = 1. - exp((-PUR[j+1][i]) / Ek[j+1] );
        tmp[i] += ( (Z[j+1][i] - Z[j][i]) * (a1 + a2) / 2. );
      }
    }
    tmp[i+1] = 0.;
    for(j=0; j<NBDEPTHS-1; j++){
      a1 = 1. - exp((-PUR[j][i+1]) / Ek[j] );
      a2 = 1. - exp((-PUR[j+1][i+1]) / Ek[j+1] );
      tmp[i+1] += ( (Z[j+1][i+1] - Z[j][i+1]) * (a1 + a2) / 2. );
    }
    dintegrale += ( INTERVALLE_TEMPS_ISCCP * (tmp[i] + tmp[i+1]) / 2. );
  }
  dintegrale *= (chl*P_Bmax_unit_d);
  
  return(dintegrale);
}

/* 
 * OUT
 * array2d_iprof_h_PUR:
 *  Array of dimensions 12 * 9.
 *  The first dimension is the optical depths:
 *  {1.,.9, .8, .7, .6, .5, .4, .3, .2, .1, .01, .001}.
 *  The second dimension is the hours from 0 to 24 by step of 3 h.
 * IN
 * array2d_ivis_iprof_E0: 
 *  Array of dimension 61 * 12.
 *  The first dimension is the wavelenght from 400 nm to 700 nm by step of 5
 *  nm.
 *  The second dimension is the optical depths:
 *  {1.,.9, .8, .7, .6, .5, .4, .3, .2, .1, .01, .001}
 *  The values are the scalar irradiances for a pixel on a given day.
 *  The units are uEinstein*m^-2*s^-1*nm^-1.
 * h:
 *  The index of the hour for the hours from 0 to 24 by step of 3 h.
 * vis_aphy:
 *  Array of dimension 61.
 *  The first dimension is the wavelenght from 400 nm to 700 nm by step of 5
 *  nm.
 *  The values are the phytoplancton absorption coefficients.
 *  The units are m^-1.
 * aphy443:
 *  A scalar.
 *  The phytoplancton absorption coefficient at 443 nm.
 *  The units are m^-1.
 * array1d_ivis_lambda:
 *  Array of dimension 61.
 *  The first dimension is the index of the wavelength.
 *  The values are the wavelengths from 400 to 700 by step 5.
 *  The units are nm.
 * Compute the PUR for a specific time and write it in a 2d array of the PUR
 * for each depth and each time.
 */
void calc_PUR(float array2d_iprof_h_PUR[NBDEPTHS][NBHOURS],
              float array2d_ivis_iprof_E0[NVIS][NBDEPTHS],
              int h,
              float vis_aphy[NVIS],
              float aphy443,
              float array1d_ivis_lambda[NVIS]){
  /*
   *  The values are the photosynthetically usable radiation (PUR) for a pixel
   *  on a given day.
   *  The units are uEinstein*m^-2*s^-1.
   */
  float array1d_iprof_PUR[NBDEPTHS];
  int iprof;
  calc_array1d_iprof_PUR(array2d_ivis_iprof_E0,
                         vis_aphy,
                         aphy443,
                         array1d_ivis_lambda,
                         array1d_iprof_PUR);
  for(iprof = 0; iprof < NBDEPTHS; iprof++){
    array2d_iprof_h_PUR[iprof][h] = array1d_iprof_PUR[iprof];
  }
}

/*
 * Calcul des profondeurs.
 * Z: Les profondeurs geometriques (en m).
 * h: Index de l'heure.
 * kd550: Le coefficient d'attenuation a 550 nm (en m^1).
 */
void calc_Z(float Z[NBDEPTHS][NBHOURS], int h, float kd550){
  int iprof;
  float array1d_iprof_Z[NBDEPTHS];
  calc_array1d_iprof_Z(kd550, array1d_iprof_Z);
  for(iprof = 0; iprof < NBDEPTHS; iprof++){
    Z[iprof][h] = array1d_iprof_Z[iprof];
  }
}

/*
 * at_412: Total absorption coefficient at 412 nm (m^-1).
 * Rrs_412: Remote sensing reflectance above the sea surface at 412 nm
 *          (nm^-1 sr^-1).
 * Rrs_490: Remote sensing reflectance above the sea surface at 490 nm
 *          (nm^-1 sr^-1).
 * Rrs_555: Remote sensing reflectance above the sea surface at 555 nm
 *          (nm^-1 sr^-1).
 * Return the absorption coefficient for CDOM at 412 nm (m^-1).
 * Reference : Belanger et al. 2008 equation 5.
 */
float calc_aCDOM_412(float at_412,
		     float Rrs_412,
		     float Rrs_490,
		     float Rrs_555){
  /*
   * Absorption coefficient for CDOM at 412 nm (m^-1).
   */
  float aCDOM_412;
  aCDOM_412 = at_412 * ( ALPHA
			 + BETA * log10(Rrs_412 / Rrs_555)
			 + CHI  * log10(Rrs_490 / Rrs_555)
			 + DELTA * log10(Rrs_555)
			 );
  if(aCDOM_412 > 4.){
    aCDOM_412 = -999.;
  }
  return aCDOM_412;
}

/*
 * rrs_type:
 *  S pour SeaWiFS et A pour MODISA
 * i:
 *  Row number on the L3BIN grid.					       
 *  Corresponds to ROW in						       
 *  http://oceancolor.gsfc.nasa.gov/SeaWiFS/TECH_REPORTS/PreLPDF/PreLVol32.pdf
 *  appendix A.							       
 *  From 1 to number of rows on the grid.                                     
 * array1d_irow_begin_MODISA:
 *  Array of dimension SZLAT_MODISA (the number of rows in the MODISA grid).
 *  The first dimension is the index (0-based) of the row in the MODISA grid.
 *  The value is the index (1-based) of the bin number (bin_num) of the first
 *  documented bin in the current row.
 *  From 1 to 23761676.
 *  Fill value: 0.
 *  Unitless.
 * array1d_irow_begin_SeaWiFS:
 *  Array of dimension SZLAT_SEAWIFS (the number of rows in the SeaWiFS grid).
 *  The first dimension is the index (0-based) of the row in the SeaWiFS grid.
 *  The value is the index (1-based) of the bin number (bin_num) of the first
 *  documented bin in the current row.
 *  From 1 to 5940422.
 *  Fill value: 0.
 *  Unitless.
 * Retourne 1 s'il y a des pixels sur la ligne.
 * Retourne 0 sinon.
 */
int get_are_pixels_on_current_line(char rrs_type,
                                   short i,
                                   int array1d_irow_begin_MODISA
                                   [SZLAT_MODISA],
                                   int array1d_irow_begin_SeaWiFS
                                   [SZLAT_SEAWIFS]){
  int out = 0;
  i--; // 1-based to 0-based.
  if((rrs_type == SEAWIFS && array1d_irow_begin_SeaWiFS[i])
     || (rrs_type == MODISA && array1d_irow_begin_MODISA[i])){
    out = 1;
  }
  return out;
}
/* ------------------------------------------------------------------ */

/*
 * IN
 * array2d_iprof_itime_val:
 *  Array of dimensions 12 * 9.
 *  The first dimension is the optical depths:
 *  {1.,.9, .8, .7, .6, .5, .4, .3, .2, .1, .01, .001}.
 *  The second dimension is the hours from 0 to 24 by step of 3 h.
 *  The values are any values.
 * itime:
 *  The index of the hour for the hours from 0 to 24 by step of 3 h.
 * OUT
 *  array1d_iprof_val:
 *  Array of dimensions 12.
 *  The first dimension is the optical depths:
 *  {1.,.9, .8, .7, .6, .5, .4, .3, .2, .1, .01, .001}.
 *  The values are the values from array2d_iprof_itime_val for the time at index
 *  itime.
 */
void get_array1d_iprof_val_from_array2d_iprof_itime_val(float array2d_iprof_itime_val
						    [NBDEPTHS][NBHOURS],
						    int itime,
						    float array1d_iprof_val
						    [NBDEPTHS]){
  int iprof;
  for(iprof = 0; iprof < NBDEPTHS; iprof++){
    array1d_iprof_val[iprof] = array2d_iprof_itime_val[iprof][itime];
  }
}
/* ------------------------------------------------------------------ */

/*
 * IN
 * array2d_iprof_itime_PUR:
 *  Array of dimensions 12 * 9.
 *  The first dimension is the optical depths:
 *  {1.,.9, .8, .7, .6, .5, .4, .3, .2, .1, .01, .001}.
 *  The second dimension is the hours from 0 to 24 by step of 3 h.
 *  The values are the photosynthetically usable radiation.
 *  The units are uE*s^-1*m^-2.
 * array2d_iprof_itime_Z:
 *  Array of dimensions 12 * 9.
 *  The first dimension is the optical depths:
 *  {1.,.9, .8, .7, .6, .5, .4, .3, .2, .1, .01, .001}.
 *  The second dimension is the geometrical depths.
 *  The units are m.
 * OUT
 *  Array of dimensions 12.					      
 *  The first dimension is the hours UTC from 0 to 24 by step of 3 h.
 *  The values are vertical attenuation of photosynthetically usable 
 *  radiation.							      
 *  The units are m^-1.                                              
 */
void get_array1d_itime_KPUR_from_PUR_and_Z(float array2d_iprof_itime_PUR
				       [NBDEPTHS][NBHOURS],
				       float array2d_iprof_itime_Z
				       [NBDEPTHS][NBHOURS],
				       float array1d_itime_KPUR[NBHOURS]){
  int itime;
  for(itime = 0; itime < NBHOURS; itime++){
    array1d_itime_KPUR[itime]
      = get_KPUR_from_PUR_and_Z_and_itime(array2d_iprof_itime_PUR,
				      array2d_iprof_itime_Z,
				      itime);
  }
}
/*---------------------------------------------------------------*/

/*
 * IN
 * doy:
 *  Day of year.
 * lat:
 *  Latitudes from -90 to 90.
 *  The units are degrees North.
 * lon:
 *  Longitudes from -180 to 180.   
 *  The units are degrees East.
 * array1d_iband_a:
 *  Array of dimension NBANDS = 6.					       
 *  The first dimension are the wavelengths of the bands of the satellite.    
 *  The values are the total absorption coefficients (a or a_t).	       
 *  The units are m^-1.                                                       
 * array1d_iband_bb:
 *  Array of dimension NBANDS = 6.					       
 *  The first dimension are the wavelengths of the bands of the satellite.    
 *  The values are the total backscattering coefficients (bb or b_bt) at      
 *  the surface (z = 0).						       
 *  The units are m^-1.                                                       
 * OUT
 * array2d_iband_itime_Kd:
 *  Array of dimension 6 * 9.					       	       
 *  The first dimension are the wavelengths of the bands of the satellite.    
 *  The second dimension is the hours UTC from 0 to 24 by step of 3 h.	       
 *  The values are the diffuse attenuation coefficient (K_d) at 	       
 *  the surface (z = 0).						       
 *  The units are m^-1.                                                       
 */
int get_array2d_iband_itime_Kd(int doy,
			       float lat,
			       float lon,
			       float array1d_iband_a[NBANDS],
			       float array1d_iband_bb[NBANDS],
			       float array2d_iband_itime_Kd[NBANDS][NBHOURS]){
  /*
   * Array of dimension 6.					       	      
   * The first dimension are the wavelengths of the bands of the satellite.    
   * The values are the diffuse attenuation coefficient (K_d) at 	      
   * the surface (z = 0).						      
   * The units are m^-1.                                                       
   */
  float array1d_iband_Kd[NBANDS];
  /*
   * Index of the band of the sensor.
   */
  int iband;
  /* 
   * The index of the hour for the hours from 0 to 24 by step of 3 h.
   */
  int itime;
  /*
   * Solar azimuth angle.
   * Units: Degrees.
   */
  float phi;
  /*
   * 0 if normal execution.
   * 1 if not.
   */
  int ret = 0;
  /*
   * Solar zenith angle.
   * Units: Degrees.
   */
  float thetas;
  /*
   * doy or doy + 1
   */
  int tmpdoy;
  if(!is_array1df_valid(array1d_iband_a, NBANDS)
     ||
     !is_array1df_valid(array1d_iband_bb, NBANDS)
     ){
    ret = -1;
  }
  if(!ret){
    for(itime = 0; itime < NBHOURS; itime++){
      ///// Compute doy. /////
      if(itime == ITIME_TOMORROW_MIDNIGHT){
	tmpdoy = doy + 1;
      }else{
	tmpdoy = doy;
      }
      for(iband = 0; iband < NBANDS; iband++){
	array1d_iband_Kd[iband] = -999.;
      }
      sunpos_(&tmpdoy,
	      &ARRAY1D_ITIME_HOUR[itime],
	      &lat,
	      &lon,
	      &thetas,
	      &phi);
      calc_kdsimon(array1d_iband_a,
		   array1d_iband_bb,
		   array1d_iband_Kd,
		   thetas);
      for(iband = 0; iband < NBANDS; iband++){
	array2d_iband_itime_Kd[iband][itime] = array1d_iband_Kd[iband];
      }
    }
  }
  return ret;
}
/*---------------------------------------------------------------*/

/*
  ENTREES:
  in_ptr_tnbfilename : Pointeur sur le nom du fichier TNB
  out_ptr_y : Pointeur sur l'annee
  out_ptr_m : Pointeur sur le mois
  out_ptr_day : Pointeur sur le jour
 */
void get_day(char *in_ptr_tnbfilename, int *out_ptr_y, int *out_ptr_m, int * out_ptr_day){
  int day, len_tnbfilename, m, pos_year, y;
  char *ptr_y_string = malloc(sizeof(char) * (LENGTH_YEAR + 1));
  char *ptr_m_string = malloc(sizeof(char) * (LENGTH_MONTH + 1));
  char *ptr_day_string = malloc(sizeof(char) * (LENGTH_DAY + 1));
  len_tnbfilename = strlen(in_ptr_tnbfilename);
  pos_year = len_tnbfilename - POS_YEAR_FROM_END;
  ptr_y_string = strncpy(ptr_y_string,
			 in_ptr_tnbfilename + pos_year,
			 LENGTH_YEAR);
  y = atoi(ptr_y_string);
  *out_ptr_y = y;
  ptr_m_string = strncpy(ptr_m_string,
			 in_ptr_tnbfilename + pos_year + LENGTH_YEAR,
			 LENGTH_MONTH);
  m = atoi(ptr_m_string);
  *out_ptr_m = m;
  ptr_day_string= strncpy(ptr_day_string,
			  in_ptr_tnbfilename + pos_year + LENGTH_YEAR + LENGTH_MONTH,
			  LENGTH_DAY);
  day = atoi(ptr_day_string);
  *out_ptr_day = day;
}
/* ------------------------------------------------------------------ */

/*
 * CF: Les valeurs de fraction nuageuse lues dans les fichiers du ISCCP sur 
 *     notre sous-grille de la grille EQ du ISCCP au nord de 45 degres Nord.
 *     Il y a 966 valeurs pour chaque heure.
 * O3: Les valeurs d'ozone lues dans les fichiers du ISCCP sur notre 
 *     sous-grille de la grille EQ du ISCCP au nord de 45 degres Nord.
 *     Il y a 966 valeurs pour chaque heure.
 * TauCld: Les valeurs de l'epaisseur optique lues dans les fichiers du ISCCP 
 *         sur notre sous-grille de la grille EQ du ISCCP au nord de 45 degres 
 *         Nord.
 *         Il y a 966 valeurs pour chaque heure.
 * lat: Latitude (de 45 a 90).
 * lon180: Longitude (de -180 a 180).
 * isccp_grid_lat: Les latitudes pour chaque pixel de notre sous-grille de la 
 *                 grille EQ du ISCCP au nord de 45 degres Nord.
 *      	   Les latitudes sont multipliees par un facteur de 100.
 *      	   Les latitudes vont de 4500 a 9000.
 * isccp_grid_lon: Les longitudes pour chaque pixel de notre sous-grille de 
 *                 la grille EQ du ISCCP au nord de 45 degres Nord.
 *                 Les longitudes vont de 0 a 360.
 * CFmean: The mean cloud fraction for a time interval of one day.
 *	   It is the mean of the cloud fraction at 00UTC, 03UTC, 06UTC, 09UTC, 
 *         12UTC, 15UTC, 18UTC and 21UTC.
 *	   Unitless.
 * O3mean: The mean ozone for a time interval of one day.
 *	   It is the mean of the ozone at 00UTC, 03UTC, 06UTC, 09UTC, 
 *         12UTC, 15UTC, 18UTC and 21UTC.
 *	   Dobson units.
 * TauCldmean: The mean optical thickness for a time interval of one day.
 *	       It is the mean of the optical thickness at 00UTC, 03UTC, 06UTC,
 *             09UTC, 12UTC, 15UTC, 18UTC and 21UTC.
 *	       Unitless.
 */
void get_isccp_mean(float CF[NBHOURS][NB_PIXELS_ISCCP],
		    float O3[NBHOURS][NB_PIXELS_ISCCP],
		    float TauCld[NBHOURS][NB_PIXELS_ISCCP],
		    float lat,
		    float lon180,
		    int isccp_grid_lat[NB_PIXELS_ISCCP],
		    float isccp_grid_lon[NB_PIXELS_ISCCP],
		    float* ptr_CFmean,
		    float* ptr_O3mean,
		    float* ptr_TauCldmean){
  int itime; // Indice de l'heure dans ARRAY1D_ITIME_HOUR.
  /*
   * Array of dimension 8.
   * The first dimension is the index of the hour in ARRAY1D_ITIME_HOUR.
   * The values are the cloud fractions (0 to 1).
   * The values are unitless.
   */
  float array1d_itime_CF[NBHOURS - 1];
  /*
   * Array of dimension 8.
   * The first dimension is the index of the hour in ARRAY1D_ITIME_HOUR.
   * The values are the total ozone column.
   * The values are Dobson units.
   */
  float array1d_itime_O3[NBHOURS - 1];
  /*
   * Array of dimension 8.
   * The first dimension is the index of the hour in ARRAY1D_ITIME_HOUR.
   * The values are the cloud optical thicknesses.
   * Unitless.
   */
  float array1d_itime_tauCld[NBHOURS - 1];
  for(itime = 0; itime < NBHOURS - 1; itime++){
    calc_isccp_eq(O3,
		  CF,
		  TauCld,
		  itime,
		  lat,
		  lon180,
		  isccp_grid_lat,
		  isccp_grid_lon,
		  &(array1d_itime_O3[itime]),
		  &(array1d_itime_CF[itime]),
		  &(array1d_itime_tauCld[itime])
		  );
  }
  *ptr_CFmean = mean(array1d_itime_CF, NBHOURS - 1);
  *ptr_O3mean = mean(array1d_itime_O3, NBHOURS - 1);
  *ptr_TauCldmean = mean(array1d_itime_tauCld, NBHOURS - 1);
}
/* ------------------------------------------------------------------ */

/*
 * glace: Concentration de glace ou valeur speciale (de 0 a 1.02).
 * Retourne 1 si la concentration de glace est < 10%.
 * Retourne 0 sinon. Retourne 0 si la concentration de glace est >= 10% ou 
 * une valeur speciale.
 */
int get_is_ice_10(float glace){
  int out = 0;
  if(glace < 0.102){
    out = 1;
  }
  return out;
}

/*
 * La concentration de glace.
 * Retourne 1 si [glace] est valide.
 * Retourne 0 sinon.
 */
int get_is_ice_val(float glace){
  int out = 0;
  if(glace < 1.002){
    out = 1;
  }
  return out;
}
/* ------------------------------------------------------------------ */

/*
 * oneCFday:
 *  The cloud fraction for a time interval of one day. Used with 
 *  MODIS-Atmosphere.
 * oneO3day:
 *  The ozone for a time interval of one day. Used with MODIS-Atmosphere.
 * oneTauCldday:
 *  The optical thickness for a time interval of one day. Used with 
 *  MODIS-Atmosphere.
 * Return 1 if the atmospheric products from MODIS-Atmosphere are valid.
 * Return 0 if not.
 */
int get_is_modis_atm_val(float oneCFday,
			 float oneO3day,
			 float oneTauCldday){
  int is_atm_val = 1;
  if(oneCFday < -0.5 || oneO3day < -0.5 || oneTauCldday < -0.5){
    is_atm_val = 0;
  }
  return is_atm_val;
}
/* ------------------------------------------------------------------ */

/*
 * IN
 * array1d_iprof_PUR:
 *  Array of dimensions 12.
 *  The first dimension is the optical depths:
 *  {1.,.9, .8, .7, .6, .5, .4, .3, .2, .1, .01, .001}.
 *  The values are the photosynthetically usable radiation.
 *  The units are uE*s^-1*m^-2.
 * Return 1 if each PUR is > 0 and each PUR is > the PUR of the depth below.
 */
int get_is_PUR_val(float array1d_iprof_PUR[NBDEPTHS]){
  int is_PUR_val = 1;
  int iprof = 0;
  while(is_PUR_val && iprof <= NBDEPTHS - 2){
    if(!gt(array1d_iprof_PUR[iprof], 0.) || !gt(array1d_iprof_PUR[iprof], array1d_iprof_PUR[iprof + 1])){
      is_PUR_val = 0;
    }
    iprof++;
  }
  return is_PUR_val;
}

/* ------------------------------------------------------------------ */

/*
 * ibinimage:
 *  Index (0-based) of the bin in the image.	  
 *  -1 if the bin is not in the image.           
 * array1d_ptr_array1d_iband_ibinimage_Rrs:
 *  Array of dimensions 6 * nbinimage
 *  The first dimension is band of the sensor.
 *  The second dimension is the bins in the image. Note that not all bins of
 *  the grid are in the image.
 *  The values are the remote-sensing reflectances (Rrs).
 *  The units are sr^-1.
 * Retourne 1 si les Rrs sont valides.
 * Retourne 0 sinon.
 */
int get_is_Rrs_val_L3BIN(int ibinimage,
                         float* array1d_ptr_array1d_iband_ibinimage_Rrs
                         [NBANDS]){
  int out = 1;
  int iband;
  for(iband = 0; iband < NBANDS; iband++){
    if( *(array1d_ptr_array1d_iband_ibinimage_Rrs[iband] + ibinimage) < 0.){
      out = 0;
    }
  }
  return out;
}
/* ------------------------------------------------------------------ */

/*
 * IN
 * lonMinus180_180: The longitude between -180 degrees and 180 degrees.
 * Return the index of the hour in ARRAY1D_ITIME_HOUR 
 * that is nearer to the local zenith (local noon).
 */
float get_itime_at_local_zenith_from_lon(float lonMinus180_180){
  float timeUTC_at_local_zenith;
  int itime;
  timeUTC_at_local_zenith = 12. - 24. * lonMinus180_180 / 360.;
  itime = lroundf(timeUTC_at_local_zenith / 3.);
  return itime;
}

/*
 * IN
 * array1d_iprof_PUR:
 *  Array of dimensions 12.
 *  The first dimension is the optical depths:
 *  {1.,.9, .8, .7, .6, .5, .4, .3, .2, .1, .01, .001}.
 *  The values are the photosynthetically usable radiation.
 *  The units are uE*s^-1*m^-2.
 * array1d_iprof_Z:
 *  Array of dimensions 12.					   
 *  The first dimension is the optical depths:			   
 *  {1.,.9, .8, .7, .6, .5, .4, .3, .2, .1, .01, .001}.	   
 *  The second dimension is the hours from 0 to 24 by step of 3 h.
 *  The values are the geometrical depths.			   
 *  The units are m.                                              
 * RETURN
 * The vertical attenuation of photosynthetically usable radiation.
 * The units are m^-1.
 */
float get_KPUR_from_PUR_and_Z(float array1d_iprof_PUR[NBDEPTHS],
			      float array1d_iprof_Z[NBDEPTHS]){
  float KPUR;
  if(!get_is_PUR_val(array1d_iprof_PUR)){
    KPUR = -999.;
  }else{
    /*
     * Array of dimensions 12.
     * The first dimension is the optical depths:
     * {1.,.9, .8, .7, .6, .5, .4, .3, .2, .1, .01, .001}.
     * The values are the proportion (between 0 and 1) of the PUR at depth
     * iprof on the PUR just below the surface.
     * The are no units.
     */
    float array1d_iprof_fPUR[NBDEPTHS];
    int iprof;
    int iz = -1;
    float rz;
    float z1p100;
    for(iprof = 0; iprof < NBDEPTHS; iprof++){
      array1d_iprof_fPUR[iprof]
        = array1d_iprof_PUR[iprof] / array1d_iprof_PUR[0];
    }
    for(iprof = 0; iprof <= NBDEPTHS - 2; iprof++){
      if(array1d_iprof_fPUR[iprof] > 0.01
         && array1d_iprof_fPUR[iprof + 1] < 0.01){
        iz = iprof;
      }
    }
    /*
     * The PUR at an optical depth of .001 is more than 1% of the PUR
     * just below the surface!
     */
    if(iz == -1){
      KPUR = -999.;
    }else{
      rz = (0.01 - array1d_iprof_fPUR[iz])
	/ (array1d_iprof_fPUR[iz + 1] - array1d_iprof_fPUR[iz]);
      z1p100 = array1d_iprof_Z[iz] * (1 - rz) + array1d_iprof_Z[iz + 1] * rz;
      KPUR = 4.6 / z1p100;
    }
  }
  return KPUR;
}
/* ------------------------------------------------------------------ */

/*
 * IN
 * array2d_iprof_itime_PUR:
 *  Array of dimensions 12 * 9.
 *  The first dimension is the optical depths:
 *  {1.,.9, .8, .7, .6, .5, .4, .3, .2, .1, .01, .001}.
 *  The second dimension is the hours from 0 to 24 by step of 3 h.
 *  The values are the photosynthetically usable radiation.
 *  The units are uE*s^-1*m^-2.
 * array2d_iprof_itime_Z:
 *  Array of dimensions 12 * 9.
 *  The first dimension is the optical depths:
 *  {1.,.9, .8, .7, .6, .5, .4, .3, .2, .1, .01, .001}.
 *  The second dimension is the geometrical depths.
 *  The units are m.
 * itime:
 *  The index of the hour for the hours from 0 to 24 by step of 3 h.
 * RETURN
 * The vertical attenuation of photosynthetically usable radiation.
 * The units are m^-1.
 */
float get_KPUR_from_PUR_and_Z_and_itime(float array2d_iprof_itime_PUR
				    [NBDEPTHS][NBHOURS],
				    float array2d_iprof_itime_Z
				    [NBDEPTHS][NBHOURS],
				    int itime){
  /*
   * Array of dimensions 12.
   * The first dimension is the optical depths:
   * {1.,.9, .8, .7, .6, .5, .4, .3, .2, .1, .01, .001}.
   * The values are the photosynthetically usable radiation.
   * The units are uE*s^-1*m^-2.
   */
  float array1d_iprof_PUR[NBDEPTHS];
  /*
   * Array of dimensions 12.					   
   * The first dimension is the optical depths:			   
   * {1.,.9, .8, .7, .6, .5, .4, .3, .2, .1, .01, .001}.	   
   * The second dimension is the hours from 0 to 24 by step of 3 h.
   * The values are the geometrical depths.			   
   * The units are m.                                              
   */
  float array1d_iprof_Z[NBDEPTHS];
  /*
   * Vertical attenuation of photosynthetically usable radiation.
   * Units: m^-1.
   */
  float KPUR;
  get_array1d_iprof_val_from_array2d_iprof_itime_val(array2d_iprof_itime_PUR,
						 itime,
						 array1d_iprof_PUR);
  get_array1d_iprof_val_from_array2d_iprof_itime_val(array2d_iprof_itime_Z,
						 itime,
						 array1d_iprof_Z);
  KPUR = get_KPUR_from_PUR_and_Z(array1d_iprof_PUR,
				 array1d_iprof_Z);
  return KPUR;
}
/* ------------------------------------------------------------------ */

/*
 * idx: 
 *  The bin index number on the (full) L3BIN grid.			       
 *  Corresponds to IDX in						       
 *  http://oceancolor.gsfc.nasa.gov/SeaWiFS/TECH_REPORTS/PreLPDF/PreLVol32.pdf
 *  appendix A.							       
 *  From 1 to the total number of possible bins in the grid.                  
 * ibinimage_firstdocumentedbin_currentrow:
 *  ibinimage (0-based) of the first pixel of the current row of the grid.
 *  Number of (documented) bins in the image in the latitudes south
 *  of the current latitude.
 * ibinimage_firstdocumentedbin_nextrow:
 *  ibinimage (0-based) of the first pixel of the next row of the grid.      
 *  Number of (documented) bins in the image in the latitudes south or equal 
 *  to the current latitude						      
 *  (ibinimage_firstdocumentedbin_currentrow + extent[i]).                   
 * ptr_bin_num:
 *  Array of dimension nbinimage.
 *  The first dimension is the index of the pixel in the image.
 *  The value is the index of the pixel in the grid.
 *  From 1 to nbinimage.
 *  Unitless.
 * Return the ndex (0-based) of the bin idx in the image.
 * -1 if the bin idx is not in the image.
 */
int get_ibinimage(int idx,
	    int ibinimage_firstdocumentedbin_currentrow,
	    int ibinimage_firstdocumentedbin_nextrow,
	    unsigned int* ptr_bin_num){
  int trouve = 0;
  int ibinimage;
  for(ibinimage = ibinimage_firstdocumentedbin_currentrow;
      !trouve && ibinimage < ibinimage_firstdocumentedbin_nextrow;
      ibinimage++){
    if(*(ptr_bin_num + ibinimage) == idx){
      trouve = 1;
    }
  }
  if(trouve){
    ibinimage--;
  }else{
    ibinimage = -1;
  }
  return ibinimage;
}
/* ------------------------------------------------------------------ */

/*
 * array2d_ilat_ilon_atm:
 *  Array of dimensions 46 x 360.				    
 *  The first dimension is the latitude from 89.5 N to 44.5 N.	    
 *  The second dimension is the longitude from -179.5 E to 179.5 E.
 *  The values are the atmospheric products read on the 	    
 *  MODIS-Atmosphere file north of 45 degrees North.               
 * lat:
 *  Latitude from 45 N to 90 N.
 * lon:
 *  Longitude from -180 E to 180 E.
 * Get the value of the atmospheric product read on the MODIS-Atmosphere file
 * for one latitude and one longitude.
 */
float get_modis_atm(float array2d_ilat_ilon_atm
		    [NBLAT_MODIS_ATMOSPHERE][NBLON_MODIS_ATMOSPHERE],
		    float lat,
		    float lon){
  float atm; // The atmospheric product.
  int ilat = (int)lroundf(89.5 - lat);
  int ilon = (int)lroundf(179.5 + lon);
  atm = array2d_ilat_ilon_atm[ilat][ilon];
  return atm;
}
/* ------------------------------------------------------------------ */

/*
 * Interpolation et extrapolation du Kd des bandes du satellite a tout le 
 * spectre visible.
 * Kd: Un tableau de taille 6 contenant les Kd.
 * vis_Kd: Un tableau sur tout le spectre visible contenant les Kd.
 * lo: Un tableau de taille 6 contenant les longueurs d'onde.
 * rrs_type: Le satellite, S pour SeaWiFS et A pour MODISA.
 */
void interp_Kd(float Kd[NBANDS],
	       float vis_Kd[NVIS],
	       float lo[NBANDS],
	       char rrs_type){
  float S;
  int i_vis;
  /*
   * Longueur d'onde maximale qui est interpolee en utilisant les longueurs 
   * d'onde aux positions 2 et 3 du satellite.
   * Les positions commencent a 0.
   * Vaut 505 nm dans le cas de SeaWiFS et 525 nm dans le cas de MODISA.
   */
  int borne_sup_3;
  /*
   * Longueur d'onde maximale qui est interpolee en utilisant les longueurs 
   * d'onde aux positions 4 et 5 du satellite.
   * Les positions commencent a 0.
   * Vaut 670 nm dans le cas de SeaWiFS et 660 nm dans le cas de MODISA.
   */
  int borne_sup_5;
  if(rrs_type == SEAWIFS){
    borne_sup_3 = 505;
    borne_sup_5 = 670;
  }else if(rrs_type == MODISA){
    borne_sup_3 = 525;
    borne_sup_5 = 660;
  }
  i_vis = 0;
  S = (log(Kd[1] / Kd[0])) / (lo[0] - lo[1]);
  while(400 + 5 * i_vis <= 410){
    vis_Kd[i_vis] = Kd[0] * exp(S * (lo[0] - (400. + 5 * i_vis)));
    i_vis++;
  }
  while(400 + 5 * i_vis <= 440){
    vis_Kd[i_vis] 
      = interp_line(lo[0], lo[1], 400. + 5 * i_vis, Kd[0], Kd[1]);
    i_vis++;
  }
  while(400 + 5 * i_vis <= 485){
    vis_Kd[i_vis]
      = interp_line(lo[1], lo[2], 400. + 5 * i_vis, Kd[1], Kd[2]);
    i_vis++;
  }
  while(400 + 5 * i_vis <= borne_sup_3){
    vis_Kd[i_vis]
      = interp_line(lo[2], lo[3], 400. + 5 * i_vis, Kd[2], Kd[3]);
    i_vis++;
  }
  while(400 + 5 * i_vis <= 550){
    vis_Kd[i_vis]
      = interp_line(lo[3], lo[4], 400. + 5 * i_vis, Kd[3], Kd[4]);
    i_vis++;
  }
  while(400 + 5 * i_vis <= borne_sup_5){
    vis_Kd[i_vis]
      = interp_line(lo[4], lo[5], 400. + 5 * i_vis, Kd[4], Kd[5]);
    i_vis++;
  }
  while(400 + 5 * i_vis <= 700){
    vis_Kd[i_vis]
      = interp_line(lo[5], 700., 400. + 5 * i_vis, Kd[5], KD700_PURE_WATER);
    i_vis++;
  }
}

/*-----------------------------------------------------------------------------
Linear interpolation
-----------------------------------------------------------------------------*/
float interp_line(float x1, float x2, float x, float y1, float y2){
  float a,b,y;
  a=(y2-y1)/(x2-x1);
  b=y1-(a*x1);
  y=a*x+b;
  return y;
}

/*
 * a: Array of strings in one dimension.
 * length: Length of a.
 * Return 1 if and only if all the strings in a are not null.
 * Return 0 if not.
 */
int is_array1d_full(char* a[], int length){
  int i;
  int ret = 1;
  for(i = 0; i < length && ret; i++){
    if(a[i] == NULL){
      ret = 0;
    }
  }
  return ret;
}

/*
 * a: Array of strings in one dimension.
 * length: Length of a.
 * Return 1 if and only if all the strings in a are null.
 * Return 0 if not.
 */
int is_array1d_null(char* a[], int length){
  int i;
  int ret = 1;
  for(i = 0; i < length && ret; i++){
    if(a[i] != NULL){
      ret = 0;
    }
  }
  return ret;
}

/*
 * chla: La concentration de chlorophylle-a (mg Chl-a m^-3).
 * a: Tableau des coefficients d'absorption (m^-1).
 * bb: Tableau des coefficients de retrodiffusion (m^-1).
 *
 * Retourne 1 si et seulement si la chlorophylle-a et les IOPs sont tous 
 * valides.
 * Retourne 0 sinon.
 */
int is_valide(float chla, float a[NBANDS], float bb[NBANDS]){
  int out;
  int iband;
  out = 1; // Valide
  if(chla < 0. || chla  > 100.){
    out = 0; // Invalide
  }
  for(iband = 0; iband < NBANDS; iband++){
    if(a[iband] < 0. || bb[iband] < 0.){
      out = 0; // Invalide
    }
  }
  return out;
}

/*
 * lon180: Longitude de -180 a 180.
 * Retourne la longitude de 0 a 360.
 */
float lon180_2_lon360(float lon180){
  float lon360;
  if(lon180 < 0){
    lon360 = lon180 + 360;
  }else{
    lon360 = lon180;
  }
  return lon360;
}

/*
 * Parse input in arr.
 * input: A string of one int and many floats. The separators are one or many 
 * space(s) between each number.
 * ptr_first_number: Pointeur vers le int d'input.
 * arr: An array of the floats.
 */
void parse(char input[500], int* ptr_first_number, float arr[NB_COL_FILE_IOP]){
  char temp[500];
  char c_input;
  int i_arr = 0;
  int i_input = 0;
  int i_temp = 0;
  int first_number = 1;
  while(i_arr < NB_COL_FILE_IOP){
    c_input = input[i_input];
    while(c_input != ' ' && c_input != '\n'){
      temp[i_temp] = input[i_input];
      i_input++;
      i_temp++;
      c_input = input[i_input];
    }
    temp[i_temp] = '\0';
    i_input++;
    i_temp = 0;
    if(temp[0] != '\0'){
      if(first_number){
	*ptr_first_number = atoi(temp);
	first_number = 0;
	arr[i_arr] = 0.;
	i_arr++;
      }else{
	arr[i_arr] = atof(temp);
	i_arr++;
      }
    }
  }
}

/*
 * Affiche a_phy pour les longueurs d'onde du visible par pas de 5
 * nm.
 * array1d_ivis_lambda:
 *  Array of dimension 61.
 *  The first dimension is the index of the wavelength.
 *  The values are the wavelengths from 400 to 700 by step 5.
 *  The units are nm.
 * in_tab_aphy: Les aphy pour les longueurs d'onde du visible par pas de 5 nm.
 */
void print_aphy(float in_tab_aphy[NVIS],
                float array1d_ivis_lambda[NVIS]){
  int ivis;
  float wl, aphy;
  printf(" Longueur d'onde aphy\n");
  printf("            (nm) (m^-1)\n");
  for(ivis = 0; ivis < NVIS; ivis++){
    wl = array1d_ivis_lambda[ivis];
    aphy = in_tab_aphy[ivis];
    printf("             %.0f %f\n", wl, aphy);
  }
}

/*
 * Affiche a_phy_443.
 * aphy443: Le coefficient d'absorption du phytoplancton a 443 nm.
 */
void print_aphy443(float aphy443){
  printf(" aphy443 (m^-1): %f\n", aphy443);
}

/*
 * array2d_h_ivis_Ed_pixel:
 *  Array of dimensions 9 * 61.						  
 *  The first dimension is the hours from 0 to 24 by step of 3 h.
 *  The second dimension is the wavelenght from 400 nm to 700 nm by step of 5
 *  nm.
 *  The values are the downward irradiances just below the surface water 
 *  (Ed0moins) for the pixel.
 *  The units are uEinstein*m^-2*s^-1*nm^-1.                             
 * Print the content of the array.
 */
void print_array2d_h_ivis_Ed_pixel(float array2d_h_ivis_Ed_pixel[NBHOURS][NVIS]){
  int h;
  int ivis;
  int wl;
  printf("Interpolation of Ed_cloud 0UTC 3UTC 6UTC 9UTC 12UTC 15UTC 18UTC 21UTC 0UTC\n");
  printf("(uEinstein*m^-2*s^-1*nm^-1)\n");
  for(ivis = 0; ivis < NVIS; ivis++){
    wl = 400 + ivis * 5;
    printf(" %d", wl);
    for(h = 0; h < NBHOURS; h++){
      printf(" %.3f",array2d_h_ivis_Ed_pixel[h][ivis]);
    }
    printf("\n");
  }
}

/*
 * Affiche les donnees liees au calcul des produits ISCCP sur la grille EQ.
 * Voir la description de chaque donnee dans les commentaire de calc_isccp_eq.
 */
void print_calc_isccp_eq(int itime,
			 float lat,
			 float lon180,
			 float lon, // Longitude de 0 a 360.
			 int latN,
			 int latS,
			 int iNW,
			 int iNE,
			 int iSW,
			 int iSE,
			 int iLATN,
			 int iLATS,
			 int last,
			 float rlat,
			 float rlonS,
			 float rlonN,
			 float rlon,
			 float oneO3,
			 float oneCF,
			 float oneTauCld){
  printf("  Interpolation 2D du ISCCP\n");
  printf("   IN\n");
  printf("    itime:      %3d\n", itime);
  printf("    lat:    %7.2f\n", lat);
  printf("    lon180: %7.2f\n", lon180);
  printf("    lon360: %7.2f\n", lon);
  printf("   CALCUL\n");
  printf("    latN:   %3d\n", latN);
  printf("    latS:   %3d\n", latS);
  printf("    iNW:    %3d\n", iNW);
  printf("    iNE:    %3d\n", iNE);
  printf("    iSW:    %3d\n", iSW);
  printf("    iSE:    %3d\n", iSE);
  printf("    iLATN:  %3d\n", iLATN);
  printf("    iLATS:  %3d\n", iLATS);
  printf("    last:   %3d\n", last);
  printf("    rlat:   %7.2f\n", rlat);
  printf("    rlonS:  %7.2f\n", rlonS);
  printf("    rlonN:  %7.2f\n", rlonN);
  printf("    rlon:   %7.2f\n", rlon);
  printf("   OUT\n");
  printf("    oneO3:     %7.2f\n", oneO3);
  printf("    oneCF:     %7.2f\n", oneCF);
  printf("    oneTauCld: %7.2f\n", oneTauCld);

}

/*
 * Affiche les donnees liees au calcul d'un produitdu ISCCP sur la grille EQ.
 * Voir la description de chaque donnee dans les commentaire de 
 * calc_isccp_eq_one_product.
 */
void print_calc_isccp_eq_one_product(float productNW,
				     float productSW,
				     float xlatW,
				     float productNE,
				     float productSE,
				     float xlatE){
  printf("  Trace de calc_isccp_eq_one_product\n");
  printf("   productNW: %f\n", productNW);
  printf("   productSW: %f\n", productSW);
  printf("   xlatW:     %f\n", xlatW);
  printf("   productNE: %f\n", productNE);
  printf("   productSE: %f\n", productSE);
  printf("   xlatE:     %f\n", xlatE);
}

/*
 * Affiche les eclairements vers le bas juste sous la surface de l'eau et a 
 * differentes longueurs d'onde.
 * Ed_pixel: Les eclairements vers le bas juste sous la surface de l'eau.
 * E0_pixel_z_t: Les eclairements vers le bas a differentes longueurs d'onde.
 */
void print_eclairement(float Ed_pixel[NBWL], float E0_pixel_z_t[NVIS][NBDEPTHS]){
  int i,iprof, ivis, iwl;
  float ed0moins ,tmp, wl;
  float LWALL[NBWL];
  tmp=290;
  for(i=0; i<NBWL; i++){
    LWALL[i] = tmp;
    tmp+=5;
  }
  printf("  Longueur d'onde Ed0-    E0_1     E0_.9    E0_.8    E0_.7    E0_.6    E0_.5    E0_.4    E0_.3    E0_.2    E0_.1    E0_.01   E0_.001\n");
  printf("             (nm) Unites de Ed0- et E0 (uE*s^-1*m^-2):\n");
  for(iwl = 0; iwl < NBWL; iwl++){
    wl = LWALL[iwl];
    ed0moins = Ed_pixel[iwl];
    printf("             %.0f %f", wl, ed0moins);
    if(iwl >= ID400){
      ivis = iwl - ID400;
      for(iprof = 0; iprof < NBDEPTHS; iprof++){
	printf(" %f", E0_pixel_z_t[ivis][iprof]);
      }
    }
    printf("\n");
  }  
}

/*
 * is_clim:
 *  Vaut 1 si et seulement si la climatologie doit etre utilisee.
 *  Vaut 0 sinon.                                                
 * ice_val: 
 *  Nombre de pixels pour lesquels la concentration de glace est valide. 
 *  (On compte tous les pixels de la grille du satellite et pas seulement
 *  ceux de l'image l3BIN.)           
 * ice_not_val:
 *  Nombre de pixels pour lesquels la concentration de glace n'est pas valide.
 *  (On compte tous les pixels de la grille du satellite et pas seulement     
 *  ceux de l'image l3BIN.)                                                   
 * l3b_val:
 *  Nombre de pixels dans l'image qui sont au nord de 45N.
 * Rrs_val:
 *  Nombre de pixels dans l'image qui sont au nord de 45N et qui ont des Rrs 
 *  valides.                    
 * oc_val:
 *  Nombre de pixels dans l'image qui sont au nord de 45N, qui ont des Rrs  
 *  valides, qui ont une chlorophylle calculee par GSM valide et qui ont des
 *  IOPs calcules par QAA valides.                                          
 * ice_10:
 *  Dans le cas sans la climatologie :					       
 *   Nombre de pixels dans l'image qui sont au nord de 45N, qui ont des Rrs   
 *   valides, qui ont une chlorophylle calculee par GSM valide, qui ont des   
 *   IOPs calcules par QAA valides et qui ont une concentration de glace      
 *   < 10%.                            				       
 *  Dans le cas avec la climatologie :  				       
 *   Nombre de pixels qui ont une concentration de glace < 10%.	       
 *   (On compte tous les pixels de la grille du satellite et pas seulement    
 *   ceux de l'image l3BIN.)                                                  
 * clim_val:
 *  Nombre de pixels pour lesquels l'algorithme utilise la climatologie.
 * pp_val:
 *  Nombre de pixels dans l'image qui sont au nord de 45 N, qui ont des Rrs 
 *  valides, qui ont une concentration de chlorophylle et des IOPs valides  
 *  tout en ayant une concentration de glace < 10% ou bien qui ont une      
 *  climatologie valide et dont la production primaire est valide.          
 */
void print_final(
		 int is_clim,
		 int ice_val,
		 int ice_not_val,
		 int l3b_val,
		 int Rrs_val,
		 int oc_val,
		 int ice_10,
		 int clim_val,
		 int pp_val){
  if(!is_clim){
    print_final_noclim(
		       l3b_val,
		       Rrs_val,
		       oc_val,
		       ice_10,
		       pp_val);
  }
}

/*
 * l3b_val:
 *  Nombre de pixels dans l'image qui sont au nord de 45N.
 * Rrs_val:
 *  Nombre de pixels dans l'image qui sont au nord de 45N et qui ont des Rrs 
 *  valides.                    
 * oc_val:
 *  Nombre de pixels dans l'image qui sont au nord de 45N, qui ont des Rrs  
 *  valides, qui ont une chlorophylle calculee par GSM valide et qui ont des
 *  IOPs calcules par QAA valides.                                          
 * ice_10:
 *  Nombre de pixels dans l'image qui sont au nord de 45N, qui ont des Rrs 
 *  valides, qui ont une chlorophylle calculee par GSM valide, qui ont des 
 *  IOPs calcules par QAA valides et qui ont une concentration de glace    
 *  < 10%.                                                                 
 * pp_val:
 *  Nombre de pixels dans l'image qui sont au nord de 45 N, qui ont des Rrs 
 *  valides, qui ont une concentration de chlorophylle et des IOPs valides  
 *  tout en ayant une concentration de glace < 10% et dont la production 
 * primaire est valide.          
 */
void print_final_noclim(
			int l3b_val,
			int Rrs_val,
			int oc_val,
			int ice_10,
			int pp_val){
  printf("L3BIN: %d\n", l3b_val);
  printf("Rrs: %d\n", Rrs_val);
  printf("chla et les IOPs: %d\n", oc_val);
  printf("glace<10%%: %d\n", ice_10);
  printf("pp: %d\n", pp_val);
}

/*
 * Affiche la glace.
 * gtype: gtype
 * ihem: ihem
 * itrans: itrans
 * iTNB: iTNB
 * jTNB: jTNB
 * xlat: Latitude (en degres)
 * xlon: Longitude (en degres)
 * Glace: La glace (sans unites).
 */
void print_glace(int gtype, int ihem, int itrans, int iTNB, int jTNB, float xlat, float xlon, float glace){
  printf(" Glace\n");
  printf("  gtype: %d\n", gtype);
  printf("  ihem: %d\n", ihem);
  printf("  itrans: %d\n", itrans);
  printf("  iTNB: %d\n", iTNB);
  printf("  jTNB: %d\n", jTNB);
  printf("  Latitude (en degres): %f\n", xlat);
  printf("  Longitude (en degres): %f\n", xlon);
  printf("  Glace (sans unites): %f", glace);
  printf("\n");
}

/*
 * array2d_jTNB_iTNB_ice:
 *  Array of dimensions 448 * 304.
 *  The first dimension is the line.
 *  The second dimension is the column.
 *  The values are the ice concentrations.
 *  Concentration de glace de 0 a 1.02 (sans unites).
 *  Les valeurs de 0 a 1 sont des concentrations de glace.
 *  1.004 (251/250.0) est le masque circulaire utilise en Arctique pour
 *  couvrir le manque de donnees de forme irreguliere autour du pole (cause
 *  par l'inclinaison de l'orbite et le swath de l'instrument).
 *  1.008 (252/250.0) est une valeur inutilisee.
 *  1.012 (253/250.0) est une ligne de cote.
 *  1.016 (254/250.0) est un masque de terre superpose.
 *  1.02 (255/250.0) est une donnee manquante.
 *  Affiche la glace de la grille SSM/I TNB du NSIDC.
 */
void print_glace_exhaustive(float array2d_jTNB_iTNB_ice[NBLINE][NBCOL]){
  int iTNB, jTNB;
  printf("Glace (%%) sur la grille TNB (NSIDC)\n");
  for(jTNB = 0; jTNB < NBLINE; jTNB++){
    for(iTNB = 0; iTNB < NBCOL; iTNB++){
      printf(" %f", array2d_jTNB_iTNB_ice[jTNB][iTNB]);
    }
    printf("\n");
  }
}

/*
 * Affiche a_phy et les IOPs pour les longueurs d'onde du visible par pas de 5
 * nm.
 * in_tab_aphy: Les aphy pour les longueurs d'onde du visible par pas de 5 nm.
 * in_tab_atotal: Les a_total pour les longueurs d'onde du visible par pas de 5 
 *                nm.
 * in_tab_bbtotal: Les bb_total pour les longueurs d'onde du visible par pas 
 *                 de 5 nm.
 * array1d_ivis_lambda:
 *  Array of dimension 61.
 *  The first dimension is the index of the wavelength.
 *  The values are the wavelengths from 400 to 700 by step 5.
 *  The units are nm.
 */
void print_iops(float in_tab_aphy[NVIS],
                float in_tab_atotal[NVIS],
                float in_tab_bbtotal[NVIS],
                float array1d_ivis_lambda[NVIS]){
  int ivis;
  float wl, aphy, atotal, bbtotal;
  printf(" Longueur d'onde aphy     atotal   bbtotal\n");
  printf("            (nm) (m^-1)   (m^-1)   (m^-1)\n");
  for(ivis = 0; ivis < NVIS; ivis++){
    wl = array1d_ivis_lambda[ivis];
    aphy = in_tab_aphy[ivis];
    atotal = in_tab_atotal[ivis];
    bbtotal = in_tab_bbtotal[ivis];
    printf("             %.0f %f %f %f \n", wl, aphy, atotal, bbtotal);
  }
}
/*
 * Affiche les donnees lues a partir des fichiers du ISCCP.
 * hh: Heure UTC
 * o3: Ozone
 * tauCld: Densite optique
 * cf: Fraction nuageuse
 */
void print_isccp_interpolation(float hh, float o3, float tauCld, float cf){
  printf("  Donnees du ISCCP\n");
  printf("   heure UTC (heures): %0.f\n", hh);
  printf("   O3 (DU): %f\n", o3);
  printf("   CF (sans unites): %f\n", cf);
  printf("   TauCld (sans unites): %f\n", tauCld);
}

/*
 * Affiche Kd pour les longueurs d'onde par pas de 5 nm.
 * vis_a: Les a_total pour les longueurs d'onde par pas de 5 nm.
 * vis_bb: Les bb_total pour les longueurs d'onde par pas de 5 nm.
 * thetas
 * Kd: Les Kd pour les longueurs d'onde par pas de 5 nm.
 * lo: Un tableau de taille 6 contenant les longueurs d'onde.
 */
void print_kd(float a[NBANDS],
	      float bb[NBANDS],
	      float thetas,
	      float Kd[NBANDS],
	      float lo[NBANDS]){
  int iband;
  float wl, atotal, bbtotal, kd;
  printf("  Longueur d'onde atotal  bbtotal  thetas    Kd\n");
  printf("             (nm) (m^-1)  (m^-1)   (degres)  (m^1)\n");
  for(iband = 0; iband < NBANDS; iband++){
    wl = lo[iband];
    atotal = a[iband];
    bbtotal = bb[iband];
    kd = Kd[iband];
    printf("             %.0f %f %f %f %f \n", wl, atotal, bbtotal, thetas, kd);
  }
}
/* ------------------------------------------------------------------ */

/*
 * IN
 *  Array of dimensions 12.					      
 *  The first dimension is the hours UTC from 0 to 24 by step of 3 h.
 *  The values are vertical attenuation of photosynthetically usable 
 *  radiation.							      
 *  The units are m^-1.                                              
 */
void print_KPUR(float array1d_itime_KPUR[NBHOURS]){
  int itime;
  float KPUR;
  float hourUTC;
  printf(" Time  KPUR\n");
  printf(" (UTC) (m^-1)\n");
  for(itime = 0; itime < NBHOURS; itime++){
    hourUTC = ARRAY1D_ITIME_HOUR[itime];
    KPUR = array1d_itime_KPUR[itime];
    printf("    %.0f %f\n", hourUTC, KPUR);
  }
}
/*
void print_KPUR(int itime, float KPUR){
  printf(" TIME (UTC) %.0f\n", ARRAY1D_ITIME_HOUR[itime]);
  printf("  KPUR (m^-1): %f\n", KPUR);
}
*/

/* ------------------------------------------------------------------ */
/*
 * Affiche le PUR moyen et le Ek pour chaque profondeur.
 * meanPUR: Le PUR moyen
 * Ek: Le Ek
 */
void print_meanpur_ek(float meanPUR[NBDEPTHS], float Ek[NBDEPTHS]){
  int iprof;
  printf(" Unites de Profondeur: sans unites\n");
  printf(" Unites de PUR moyen:  uE*s^-1*m^-2\n");
  printf(" Unites de Ek:         uE*s^-1*m^-2\n");
  printf(" Profondeur PUR moyen  Ek\n");
  for(iprof = 0; iprof < NBDEPTHS; iprof++){
    printf("  E0_%0.3f  %f %f\n",
	   ARRAY1D_IDEPTH_OPTICALDEPTHS[iprof],
	   meanPUR[iprof],
	   Ek[iprof]);
  }
}

/*
 * Affiche la photoperiode.
 * photoperiode: La photoperiode.
 */
void print_photoperiode(float photoperiode){
  printf(" Photoperiode (s): %f\n", photoperiode);
}

/*
 * Affiche la production primaire.
 * pp: La production primaire.
 */
void print_pp(float pp){
  printf(" Production primaire (mgC*m^-2*d^-1): %f\n", pp);
  //printf("%f\n", pp);
}

/*
 * Affiche la production primaire sans tenir compte de la glace.
 * pp: La production primaire sans tenir compte de la glace.
 */
void print_pp_temp(float pp){
  printf(" Production primaire sans la glace (mgC*m^-2*d^-1): %f\n", pp);
}

/*
 * Affiche le PUR pour chaque profondeur.
 * itime: L'indice de l'heure
 * PUR: Le PUR
 */
void print_pur(int itime, float PUR[NBDEPTHS][NBHOURS]){
  int iprof;
  printf("  PUR\n");
  for(iprof = 0; iprof < NBDEPTHS; iprof++){
    printf("   E0_%0.3f (uE*s^-1*m^-2): %f\n",
	   ARRAY1D_IDEPTH_OPTICALDEPTHS[iprof],
	   PUR[iprof][itime]);
  }
}

/*
 * Affiche le numero du pixel sur la grille ISIN, la latitude, la longitude
 * et la production primaire.
 * no_pixel_grille_isin: Le numero du pixel sur la grille ISIN.
 * lat: La latitude.
 * lon: La longitude.
 * pp: La production primaire.
 */
void print_result(int no_pixel_grille_isin, float lat, float lon, float pp){
  printf("%d %f %f %f\n", no_pixel_grille_isin, lat, lon, pp);
}

/*
 * Affiche les parametres d'entree pour le calcul de thetas et affiche thetas.
 * doy: Jour julien
 * hh: Heure UTC
 * xlat: Latitude
 * xlon: Longitude
 * thetas
 */
void print_thetas(int doy, float hh, float xlat, float xlon, float thetas){
  printf("  Calcul de thetas pour Kd\n");
  printf("   ENTREES\n");
  printf("    doy %d\n", doy);
  printf("    hh (heure UTC en heures): %.0f\n", hh);
  printf("    xlat (degres): %f\n", xlat);
  printf("    xlon (degres): %f\n", xlon);
  printf("   SORTIE\n");
  printf("    thetas pour le calcul de Kd (degres): %f\n", thetas);
}

/*
 * Affiche que thetas vaut 45 parce que CF > 0.3.
 */
void print_thetas45(){
  printf("  CF > 0.3 donc thetas pour le calcul de Kd (en degres): 45.\n");
}

/*
 * Affiche Kd pour les longueurs d'onde du visible par pas de 5 nm.
 * vis_Kd: les Kd pour les longueurs d'onde du visible par pas de 5 nm.
 * array1d_ivis_lambda:
 *  Array of dimension 61.
 *  The first dimension is the index of the wavelength.
 *  The values are the wavelengths from 400 to 700 by step 5.
 *  The units are nm.
 */
void print_vis_Kd(float vis_Kd[NVIS],
                  float array1d_ivis_lambda[NVIS]){
  int i_vis;
  float wl, kd;
  printf("  Longueur d'onde Kd\n");
  printf("            (nm) (m^-1)\n");
  for(i_vis = 0; i_vis < NVIS; i_vis++){
    wl = array1d_ivis_lambda[i_vis];
    kd = vis_Kd[i_vis];
    printf("             %.0f %f\n", wl, kd);
  }
}

/*
 * Affiche les profondeurs geometriques.
 * Z: Les profondeurs geometriques (en m).
 * h: L'index des heures.
 */
void print_Z(float Z[NBDEPTHS][NBHOURS], int h){
  int i;
  printf("  Profondeurs optiques (sans unites) Profondeurs geometriques (m)\n");
  for(i = 0; i < NBDEPTHS; i++){
    printf("                               %.3f",
	   ARRAY1D_IDEPTH_OPTICALDEPTHS[i]);
    printf(" %f\n", Z[i][h]);
  }
}
/* ------------------------------------------------------------------ */

/*
 * IN
 * filename: Name of the ocean color file name.
 * rrs_type: S for SeaWiFS and A for MODISA.
 * OUT
 * array1d_irow_start_num_MODISA:
 *  Array of dimension SZLAT_MODISA (the number of rows in the MODISA grid).
 *  The first dimension is the index (0-based) of the row in the MODISA grid.
 *  The value is the index (1-based) of the bin number (bin_num) of the first
 *  bin in the current row.
 *  From 1 to 23761674.
 *  Unitless.
 * array1d_irow_start_num_SeaWiFS:
 *  Array of dimension SZLAT_SEAWIFS (the number of rows in the SeaWiFS grid).
 *  The first dimension is the index (0-based) of the row in the SeaWiFS grid.
 *  The value is the index (1-based) of the bin number (bin_num) of the first
 *  bin in the current row.
 *  From 1 to 5940420.
 *  Unitless.
 * array1d_irow_begin_MODISA:
 *  Array of dimension SZLAT_MODISA (the number of rows in the MODISA grid).
 *  The first dimension is the index (0-based) of the row in the MODISA grid.
 *  The value is the index (1-based) of the bin number (bin_num) of the first
 *  documented bin in the current row.
 *  From 1 to 23761676.
 *  Fill value: 0.
 *  Unitless.
 * array1d_irow_begin_SeaWiFS:
 *  Array of dimension SZLAT_SEAWIFS (the number of rows in the SeaWiFS grid).
 *  The first dimension is the index (0-based) of the row in the SeaWiFS grid.
 *  The value is the index (1-based) of the bin number (bin_num) of the first
 *  documented bin in the current row.
 *  From 1 to 5940422.
 *  Fill value: 0.
 *  Unitless.
 * array1d_irow_extent_MODISA:
 *  Array of dimension SZLAT_MODISA (the number of rows in the MODISA grid).
 *  The first dimension is the index (0-based) of the row in the MODISA grid.
 *  The value is the number of documented bins in the current row.
 *  From 0 to the total number of possible bins in each row.
 *  The units are bins.
 * array1d_irow_extent_SeaWiFS:
 *  Array of dimension SZLAT_SEAWIFS (the number of rows in the SeaWiFS grid).
 *  The first dimension is the index (0-based) of the row in the SeaWiFS grid.
 *  The value is the number of documented bins in the current row.
 *  From 0 to the total number of possible bins in each row.
 *  The units are bins.
 * array1d_irow_max_MODISA:
 *  Array of dimension SZLAT_MODISA (the number of rows in the MODISA grid).
 *  The first dimension is the index (0-based) of the row in the MODISA grid.
 *  The value is the number of possible bins in the current row.
 *  From 3 to 8640.
 *  The units are bins.
 * array1d_irow_max_SeaWiFS:
 *  Array of dimension SZLAT_SEAWIFS (the number of rows in the SeaWiFS grid).
 *  The first dimension is the index (0-based) of the row in the SeaWiFS grid.
 *  The value is the number of possible bins in the current row.
 *  From 3 to 8640.
 *  The units are bins.
 * ptr_ptr_bin_num:
 *  Array of dimension nbinimage.
 *  The first dimension is the index of the pixel in the image.
 *  The value is the index of the pixel in the grid.
 *  From 1 to nbinimage.
 *  Unitless.
 * array1d_ptr_array1d_iband_ibinimage_Rrs:
 *  Array of dimensions 6 * nbinimage
 *  The first dimension is band of the sensor.
 *  The second dimension is the bins in the image. Note that not all bins of 
 *  the grid are in the image.
 *  The values are the remote-sensing reflectances (Rrs).
 *  The units are sr^-1.
 * Read the Rrs in the ocean color file.
 */
void read_l3b(char *filename,
              char rrs_type,
              int array1d_irow_start_num_MODISA[SZLAT_MODISA],
              int array1d_irow_start_num_SeaWiFS[SZLAT_SEAWIFS],
              int array1d_irow_begin_MODISA[SZLAT_MODISA],
              int array1d_irow_begin_SeaWiFS[SZLAT_SEAWIFS],
              int array1d_irow_extent_MODISA[SZLAT_MODISA],
              int array1d_irow_extent_SeaWiFS[SZLAT_SEAWIFS],
              int array1d_irow_max_MODISA[SZLAT_MODISA],
              int array1d_irow_max_SeaWiFS[SZLAT_SEAWIFS],
              unsigned int** ptr_ptr_bin_num,
              float* array1d_ptr_array1d_iband_ibinimage_Rrs[NBANDS]){
  
  /////////// declaration of the variables ///////////
  
  /* The index (1-based) of the bin number (bin_num) of the first
   * documented bin in the current row. */
  int begin;
  /*
   * The bin index number of the global L3BIN grid.
   * From 1 to TOTBINS.
   */
  unsigned int bin_num = 0;
  /* The number of bins in the image. */
  size_t data_bins;
  /* Dimension ID of BinIndex */
  int dimid_binindex = -1;
  /* Dimension ID of BinList */
  int dimid_binlist = -1;
  /* The number of documented bins in the current row. */
  int extent;
  // The band of the sensor.
  int iband;
  /*
   * Index of the data bin.
   * From 0 to data_bins - 1.
   */
  int idata_bin = -1;
  /*
   * Index (0-based) of the row in the MODIS or SeaWiFS grid.
   * 0 to 2159 for SeaWiFS.
   * 0 to 4319 for MODIS.
   */
  int irow;
  /* Group ID for level-3_binned_data */
  int grpid = -1;
  /* The number of possible bins in the current row. */
  int max;
  /* NetCDF ID for the file. */
  int ncid = -1;
  // The fieldname.
  char* prodname;
  /*
   * The instrument (aka Sensor Name).
   * "MODIS" or "SeaWiFS".
   */
  char* sensor;
  /*
   * The (integer) number of rows in the grid (equal to 2160 for SeaWiFS and to
   * 4320 for MODIS).
   */
  size_t numrows = -1;
  /* Return value */
  float retval = -1;
  char* array1d_iband_RrsMODIS[NBANDS]={
    "Rrs_412",
    "Rrs_443",
    "Rrs_488",
    "Rrs_531",
    "Rrs_555",
    "Rrs_667"
  };
  char* array1d_iband_RrsSeaWiFS[NBANDS]={
    "Rrs_412",
    "Rrs_443",
    "Rrs_490",
    "Rrs_510",
    "Rrs_555",
    "Rrs_670"
  };
  /*
   * The lenght of the sensor variable.
   */
  size_t s_lenght;
  /*
   * The index (1-based) of the bin number (bin_num) of the first bin in 
   * the current row.
   */
  int start_num;
  // The sum of the data.
  float sum;
  // Value of the data.
  float val;
  /* Compound data variable ID for BinIndex. */
  int varid_binindex = -1;
  /* Compound data variable ID for BinList. */
  int varid_binlist = -1;
  // Compound data variable ID for product.
  int varid_prodname = -1;
  // The weight of a data.
  float weights;
  
  /////////// read BinIndex ///////////
  
  /* Open the file. */
  if ((retval = nc_open(filename, NC_NOWRITE, &ncid))){ERR(retval);}
  
  /* get sensor name */
  if( (retval = nc_inq_attlen(ncid, NC_GLOBAL, "instrument", &s_lenght)) ){
    ERR(retval);
  }
  sensor = (char*) malloc(s_lenght + 1);
  if ( (retval=nc_get_att_text(ncid, NC_GLOBAL, "instrument", sensor) )){
    ERR(retval);
  }

  /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   Important !: read the number of rows before declaring the structures to read
   the data
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
  
  /* Get the group id of level-3_binned_data. */
  if ((retval = nc_inq_ncid(ncid, "level-3_binned_data", &grpid))){
    ERR(retval);
  }
  /* Get the varid of the compound data variable, based on its name,
   * in grpid. */
  if ((retval = nc_inq_varid(grpid, "BinIndex", &varid_binindex))){
    ERR(retval);
  }
  /* Get the dimid of the dimension of the compound data variable. */
  if ((retval = nc_inq_vardimid(grpid, varid_binindex, &dimid_binindex))){
    ERR(retval);
  }
  /* Get the length of the dimension of the compound data variable. */
  if ((retval = nc_inq_dimlen(grpid, dimid_binindex, &numrows))){
    ERR(retval);
  }
  // Struct for BinIndex
  struct sBinIndex{
    int start_num;
    int begin;
    int extent;
    int max;
  };
  // struct sBinList binIndex[numrows]
  struct sBinIndex* binIndex = malloc(numrows * sizeof(struct sBinIndex));
  /* Read the data. */
  if ((retval = nc_get_var(grpid, varid_binindex, &binIndex[0]))){
    ERR(retval);
  }
  for(irow = 0; irow < numrows; irow++){
    start_num = binIndex[irow].start_num;
    begin     = binIndex[irow].begin;
    extent    = binIndex[irow].extent;
    max       = binIndex[irow].max;
    if(!strcmp(sensor, "MODIS")){
      array1d_irow_start_num_MODISA[irow] = start_num;
      array1d_irow_begin_MODISA[irow]     = begin;
      array1d_irow_extent_MODISA[irow]    = extent;
      array1d_irow_max_MODISA[irow]       = max;
    }else if(!strcmp(sensor, "SeaWiFS")){
      array1d_irow_start_num_SeaWiFS[irow] = start_num;
      array1d_irow_begin_SeaWiFS[irow]     = begin;
      array1d_irow_extent_SeaWiFS[irow]    = extent;
      array1d_irow_max_SeaWiFS[irow]       = max;
    }
  }

	/////////// read BinList ///////////

	/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   Important !: read the number of data bins before declaring the structures to
	read the data
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */

  /* Get the varid of the compound data variable, based on its name,
   * in grpid. */
  if ((retval = nc_inq_varid(grpid, "BinList", &varid_binlist))){ERR(retval);}
  /* Get the dimid of the dimension of the compound data variable. */
  if ((retval = nc_inq_vardimid(grpid, varid_binlist, &dimid_binlist))){
    ERR(retval);
  }
  /* Get the length of the dimension of the compound data variable. */
  if ((retval = nc_inq_dimlen(grpid, dimid_binlist, &data_bins))){
    ERR(retval);
  }
  // Struct of BinList.
  struct sBinList{
    unsigned int bin_num ;
    short nobs ;
    short nscenes ;
    float weights ;
    float time_rec ;
  };
  // struct sBinList binList[data_bins]
  struct sBinList* binList = malloc(data_bins * sizeof(struct sBinList));
  /* Read the data. */
  if ((retval = nc_get_var(grpid, varid_binlist, &binList[0]))){
    ERR(retval);
  }
  // Allocate memory for *ptr_ptr_bin_num.
  *ptr_ptr_bin_num = (unsigned int*)malloc(data_bins * sizeof(unsigned int));
  for(idata_bin = 0; idata_bin < data_bins; idata_bin++){
    bin_num = binList[idata_bin].bin_num;
    (*ptr_ptr_bin_num)[idata_bin] = bin_num;
  }
  
	/////////// read Rrs ///////////

	/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   Important !: read the number of data bins before declaring the structures to
   read the data
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */

  // Struct of product.
  struct sProduct{
    float sum;
    float sum_squared;
  };
  // struct sProduct product[data_bins].
  struct sProduct* product = malloc(data_bins * sizeof(struct sProduct));
  for(iband = 0; iband < NBANDS; iband++){
    // Allocate memory for the Rrs.
    array1d_ptr_array1d_iband_ibinimage_Rrs[iband]
    	= malloc(data_bins * sizeof(float));
    if(!strcmp(sensor, "MODIS")){
      prodname = array1d_iband_RrsMODIS[iband];
    }else if(!strcmp(sensor, "SeaWiFS")){
      prodname = array1d_iband_RrsSeaWiFS[iband];
    }
    /* Get the varid of the compound data variable, based on its name,
     * in grpid. */
    if ((retval = nc_inq_varid(grpid, prodname, &varid_prodname))){
      ERR(retval);
    }
    /* Read the data. */
    if ((retval = nc_get_var(grpid, varid_prodname, &product[0]))){
      ERR(retval);
    }
    for(idata_bin = 0; idata_bin < data_bins; idata_bin++){
      sum = product[idata_bin].sum;
      weights = binList[idata_bin].weights;
      val = sum / weights;
      array1d_ptr_array1d_iband_ibinimage_Rrs[iband][idata_bin] = val;
    }
  }
}
/* ------------------------------------------------------------------ */

/*
 * O3:
 *  Tableau de dimensions 9 x 966.					   
 *  Ce tableau contient les valeurs des produits atmospheriques lus dans  
 *  les fichiers du ISCCP sur notre sous-grille de de la grille EQ du 	   
 *  ISCCP au nord de 45 degres Nord.					   
 *  Il y a 966 valeurs pour chaque heure.				   
 *  O3 est l'ozone.							   
 *  CF est la fraction nuageuse.					   
 *  TauCld est l'epaisseur optique.                                       
 * CF:
 *  Tableau de dimensions 9 x 966.					   
 *  Ce tableau contient les valeurs des produits atmospheriques lus dans  
 *  les fichiers du ISCCP sur notre sous-grille de de la grille EQ du 	   
 *  ISCCP au nord de 45 degres Nord.					   
 *  Il y a 966 valeurs pour chaque heure.				   
 *  O3 est l'ozone.							   
 *  CF est la fraction nuageuse.					   
 *  TauCld est l'epaisseur optique.                                       
 * TauCld:
 *  Tableau de dimensions 9 x 966.					   
 *  Ce tableau contient les valeurs des produits atmospheriques lus dans  
 *  les fichiers du ISCCP sur notre sous-grille de de la grille EQ du 	   
 *  ISCCP au nord de 45 degres Nord.					   
 *  Il y a 966 valeurs pour chaque heure.				   
 *  O3 est l'ozone.							   
 *  CF est la fraction nuageuse.					   
 *  TauCld est l'epaisseur optique.                                       
 * array1d_ptr_array1d_itime_ipix_O3sat:
 *  Tableau de pointeurs de dimensions 9 x nbpix.			     
 *  La premiere dimension est les heures.				     
 *  La deuxieme dimension est le pixel sur la grille du satellite.	     
 *  Les valeurs sont les produits atmospheriques lus sur la grille du ISCCP 
 *  reprojetees sur la grille du satellite au nord de 45 degres Nord.	     
 *  O3 est l'ozone.							     
 *  CF est la fraction nuageuse.					     
 *  TauCld est l'epaisseur optique.                                         
 * array1d_ptr_array1d_itime_ipix_CFsat:
 *  Tableau de pointeurs de dimensions 9 x nbpix.			     
 *  La premiere dimension est les heures.				     
 *  La deuxieme dimension est le pixel sur la grille du satellite.	     
 *  Les valeurs sont les produits atmospheriques lus sur la grille du ISCCP 
 *  reprojetees sur la grille du satellite au nord de 45 degres Nord.	     
 *  O3 est l'ozone.							     
 *  CF est la fraction nuageuse.					     
 *  TauCld est l'epaisseur optique.                                         
 * array1d_ptr_array1d_itime_ipix_TauCldsat:
 *  Tableau de pointeurs de dimensions 9 x nbpix.			     
 *  La premiere dimension est les heures.				     
 *  La deuxieme dimension est le pixel sur la grille du satellite.	     
 *  Les valeurs sont les produits atmospheriques lus sur la grille du ISCCP 
 *  reprojetees sur la grille du satellite au nord de 45 degres Nord.	     
 *  O3 est l'ozone.							     
 *  CF est la fraction nuageuse.					     
 *  TauCld est l'epaisseur optique.                                         
 * szlat:
 *  Nombre de lignes (latitudes) de 90S a 90N sur la grille du satellite.
 * nbpix:
 *  Nombre de pixels au nord de 45N sur la grille du satellite.
 */
void reproject_isccp_grid_to_sat_grid(
				      float O3[NBHOURS][NB_PIXELS_ISCCP],
				      float CF[NBHOURS][NB_PIXELS_ISCCP],
				      float TauCld[NBHOURS][NB_PIXELS_ISCCP],
				      float *array1d_ptr_array1d_itime_ipix_O3sat[NBHOURS],
				      float *array1d_ptr_array1d_itime_ipix_CFsat[NBHOURS],
				      float *array1d_ptr_array1d_itime_ipix_TauCldsat[NBHOURS],
				      int szlat,
				      int nbpix){
  int i; //Indice de la ligne sur la grille SeaWiFS ou MODISA et sur l'image.
  int ipix; //Index du pixel courant sur la grille du satellite.
  /*
   * Les latitudes pour chaque pixel de notre sous-grille de la grille EQ du 
   * ISCCP au nord de 45 degres Nord.					    
   * Les latitudes sont multiplies par un facteur de 100.		    
   * Les latitudes vont de 4500 a 9000.                                     
   */
  int isccp_grid_lat[NB_PIXELS_ISCCP];
  /*
   * Les longitudes pour chaque pixel de notre sous-grille de la grille EQ du 
   * ISCCP au nord de 45 degres Nord.
   * Les longitudes vont de 0 a 360.
   */
  float isccp_grid_lon[NB_PIXELS_ISCCP];
  int j; // Indice des pixels sur la grille du satellite.
  int itime; // Indice des heures.
  /* Nombre de bins sur la ligne courante de la grille du satellite. */
  int numbin_local;
  float oneCF; // Une valeur de fraction nuageuse.
  float oneO3; // Une valeur d'ozone.
  float oneTauCld; // Une valeur d'epaisseur optique.
  float xlat; // La latitude de 0N a 90N.
  float xlon; // La longitude de -180E a 180E.
  /* Allocation de l'espace memoire pour les tableaux de pointeurs. */
  for(itime = 0; itime < NBHOURS; itime++){
    array1d_ptr_array1d_itime_ipix_O3sat[itime]
      = (float*)malloc(sizeof(float) * nbpix);
    array1d_ptr_array1d_itime_ipix_CFsat[itime]
      = (float*)malloc(sizeof(float) * nbpix);
    array1d_ptr_array1d_itime_ipix_TauCldsat[itime] =
      (float*)malloc(sizeof(float) * nbpix);
  }
  /* initialisation des Latitudes/Longitudes des donnees SRF sur la grille EQ */
  calc_isccp_grid(isccp_grid_lat, isccp_grid_lon);
  if(DEBUG >= DEBUG_2){
    print_isccp_grid_eq(isccp_grid_lat, isccp_grid_lon);
    print_isccp_grid_eq_parabole(isccp_grid_lat, isccp_grid_lon);
  }
  ipix = 0;
  for (i=0; i < szlat; i++){
    xlat = ((i + .5) * 180. / szlat) - 90.;
    numbin_local = (int) (2*szlat*cos((xlat)*PI/180.)+.5);
    // Le pixel est-il au nord de 45N?
    if (xlat > 45.){
      // Boucle sur les colonnes (ie: les longitudes) de la grille SeaWiFS.
      for(j = 1; j <= numbin_local; j++){
	xlon = (360. * (j - .5) / numbin_local) - 180.;
	// Boucle sur les heures.
	for(itime = 0; itime < NBHOURS; itime++){
	  calc_isccp_eq(O3,
			CF,
			TauCld,
			itime,
			xlat,
			xlon,
			isccp_grid_lat,
			isccp_grid_lon,
			&oneO3,
			&oneCF,
			&oneTauCld);
	  *(array1d_ptr_array1d_itime_ipix_O3sat[itime] + ipix) = oneO3;
	  *(array1d_ptr_array1d_itime_ipix_CFsat[itime] + ipix) = oneCF;
	  *(array1d_ptr_array1d_itime_ipix_TauCldsat[itime] + ipix) = oneTauCld;
	  if(DEBUG >= DEBUG_2){
	    print_isccp_products_eq_parabole(O3,
					     CF,
					     TauCld,
					     ARRAY1D_ITIME_HOUR,
					     itime);
	    print_isccp_interpolation(ARRAY1D_ITIME_HOUR[itime],
				      oneO3,
				      oneTauCld,
				      oneCF);
	  }
	} // Fin de la boucle sur les heures.
	ipix++;
      } /* Fin de la boucle sur les colonnes (ie: les longitudes) de la 
	 * grille SeaWiFS. */
    } // Fin du si (Le pixel est-il au nord de 45N?)
  }/* Fin de la boucle sur les lignes (ie: les latitudes) de la grille du
   * satellite.*/
  // Pour s'assurer que le nombre de pixels traites correspond bien au 
  // nombre de pixels sur la grille du satellite au nord de 45N.
  if(ipix != nbpix){
    printf("prodprim_l3b.c: reproject_isccp_grid_to_sat_grid(...) ipix: %d nbpix: %d\n",
	   ipix,
	   nbpix);
    exit(-1);
  }
}

/* -------------------------------------------------------
 * fic_out:
 *  Nom du fichier de sortie.
 * nb_valid:
 *  Nombre de pixels a ecrire.
 * lo: 
 *  Tableau de dimension 6.			      
 *  Contient les longueurs d'onde du satellite en nm.
 * ptr_array1d_ipix_lat:
 *  Array of dimension nbpix.
 *  The first dimension is the pixel.
 *  The values are the latitudes from -90 to 90.
 *  The units are degrees N.
 * ptr_array1d_ipix_lon:
 *  Array of dimension nbpix.
 *  The first dimension is the pixel.
 *  The values are the latitudes from -180 to 180.
 *  The units are degrees E.
 * ptr_array1d_ipix_indlat:
 *  Array of dimension nbpix.
 *  The first dimension is the pixel.
 *  The values are the row number on the L3BIN grid.
 *  Corresponds to ROW - 1 in
 *  http://oceancolor.gsfc.nasa.gov/SeaWiFS/TECH_REPORTS/PreLPDF/PreLVol32.pdf
 *  appendix A.
 * Array of dimension nbpix.
 *  The first dimension is the pixel.
 *  The values are the bin number on this row.
 *  Corresponds to COL in
 *  http://oceancolor.gsfc.nasa.gov/SeaWiFS/TECH_REPORTS/PreLPDF/PreLVol32.pdf
 *  appendix A.
 *  From 1 to the total number of possible bins in this row.
 * array1d_ptr_array1d_iband_ipix_Rrs:
 *  Array of dimensions 6 X nbpix.
 *  (Array of pointers.)
 *  The first dimension is the band of the satellite.
 *  The second dimension is the pixel.
 *  The values are the remote-sensing reflectances (Rrs).
 *  The units are sr^-1.
 * array1d_ptr_array1d_iband_ipix_a:
 *  Array of dimensions 6 X nbpix.
 *  (Array of pointers.)
 *  The first dimension is the band of the satellite.
 *  The second dimension is the pixel.
 *  The values are the absorption coefficients.
 *  The units are m^-1.
 * array1d_ptr_array1d_iband_ipix_bb:
 *  Array of dimensions 6 X nbpix.
 *  (Array of pointers.)
 *  The first dimension is the band of the satellite.
 *  The second dimension is the pixel.
 *  The values are the backscattering coefficients (bb).
 *  The units are m^-1.
 * array1d_ptr_array1d_iband_ipix_Kd:
 *  Array of dimensions 6 X nbpix.				       
 *  (Array of pointers.)					       
 *  The first dimension is the band of the satellite.		       
 *  The second dimension is the pixel.				       
 *  The values are the vertical attenuation coefficients (Kd) for one 
 *  wavelength and one pixel at the local noon.				  
 *  The units are m^-1.                                                   
 * ptr_array1d_ipix_Zeu_bar:
 *  Array of dimension nbpix.						     
 *  The first dimension is the pixel.
 *  The values are the mean euphotic depths.
 *  The units are m.
 * ptr_array1d_ipix_poc:
 *  Array of dimension nbpix.			   
 *  The first dimension is the pixel.		   
 *  The values are the particulate organic carbon.                          
 * ptr_array1d_ipix_aCDOM_412:
 *  Array of dimension nbpix.						        
 *  The columns contain the absorption coefficient for CDOM at 412 nm (m^-1).
 * ptr_array1d_ipix_PAR_cloud:
 *  Array of dimension nbpix.						     
 *  The first dimension is the pixel.					     
 *  The values are the photosynthetically available irradiances (PAR) for a 
 *  pixel on a given day.						     
 *  The units are uEinstent*m^-1.                                           
 * float ptr_array1d_ipix_CF00:
 *  Array of dimension nbpix.				      
 *  The first dimension is the pixel.			      
 *  The values are the cloud fractions (0 to 1) at 00h UTC.  
 *  The values are unitless.                                 
 * ptr_array1d_ipix_O300:
 *  Array of dimension nbpix.			      
 *  The first dimension is the pixel.		      
 *  The values are the total ozone column at 00h UTC.
 *  The units are Dobson units.                      
 * ptr_array1d_ipix_tauCld00:
 *  Array of dimension nbpix.				      
 *  The first dimension is the pixel.			      
 *  The values are the cloud optical thickness at 00h UTC.   
 *  The values are unitless.                                 
 * ptr_array1d_ipix_thetas00:
 *  Array of dimension nbpix.				
 *  The first dimension is the pixel.			
 *  The values are the solar zenith angles at 00h UTC.	
 *  The values are in degrees.                         
 * ptr_array1d_ipix_CFmean:
 *  Array of dimension nbpix.				      
 *  The first dimension is the pixel.			      
 *  The values are the mean cloud fractions (0 to 1).        
 *  The values are unitless.                                 
 * ptr_array1d_ipix_O300:
 *  Array of dimension nbpix.			      
 *  The first dimension is the pixel.		      
 *  The values are the total ozone column at 00h UTC.
 *  The units are Dobson units.                      
 * ptr_array1d_ipix_tauCldmean:
 *  Array of dimension nbpix.				      
 *  The first dimension is the pixel.			      
 *  The values are the mean cloud optical thickness 	      
 *  The values are unitless.                                 
 * ptr_array1d_ipix_KPUR:
 *  Array of dimension nbpix.			      			   
 *  The first dimension is the pixel.		      			   
 *  The values are the vertical attenuation of photosynthetically usable  
 *  radiation.								   
 *  The units are m^-1.                                                   
 * ptr_array1d_ipix_chl_cota:
 *  Array of dimension nbpix.			      			     
 *  The first dimension is the pixel.		      			     
 *  The values are the chlorophyll-a concentration from Cota et al. 2004.
 *  The units are (mg Chl-a m^-3).
 * ptr_array1d_ipix_chl_oc:
 *  Array of dimension nbpix.			      			     
 *  The first dimension is the pixel.		      			     
 *  The values are the chlorophyll-a concentration from OC3 or OC4.	   
 *  The units are (mg Chl-a m^-3).
 * ptr_array1d_ipix_chl_gsm_mustapha:
 *  Array of dimension nbpix.			      			       
 *  The first dimension is the pixel.		      			       
 *  The values are the chlorophyll-a concentration from Ben Mustapha et       
 *  al. 2012.								       
 *  The units are (mg Chl-a m^-3).                                            
 * ptr_array1d_ipix_PP_cota:
 *  Array of dimension nbpix.						     
 *  The first dimension is the pixel.					     
 *  The values are the primary productivity computed with Cota et al. 2004. 
 *  The units are mgC*m^-2*d^-1.                                            
 * ptr_array1d_ipix_PP_gsm_mustapha:
 *  Array of dimension nbpix.						       
 *  The first dimension is the pixel.					       
 *  The values are the primary productivity computed with Ben Mustapha et al. 
 *  2012.								       
 *  The units are mgC*m^-2*d^-1.                                              
 * ptr_array1d_ipix_PP_oc:
 *  Array of dimension nbpix.						     
 *  The first dimension is the pixel.					     
 *  The values are the primary productivity computed with OC3 or OC4.	     
 *  The units are mgC*m^-2*d^-1.                                            
 * ptr_array1d_ipix_aphy443:
 *  Array of dimension nbpix.
 *  The first dimension is the pixel.
 *  The values are phytoplankton absorption coefficient at 443 nm (a_phy at
 *  443 nm).
 *  The units are m^-1.
 * ptr_array1d_ipix_chl_gsm:
 *  Array of dimension nbpix.
 *  The first dimension is the pixel.
 *  The values are the chlorophyll-a concentration from GSM.
 *  The units are (mg Chl-a m^-3).
 * ptr_array1d_ipix_PP_chl_gsm:
 *  Array of dimension nbpix.
 *  The first dimension is the pixel.
 *  The values are the primary productivity computed with GSM.
 *  The units are mgC*m^-2*d^-1.
 * ptr_array1d_ipix_ice:
 *  Array of dimension nbpix.
 *  The first dimension is the pixel.
 *  The values are the sea-ice concentration.
 *  Normal values from 0 to 1.
 *  Special values > 1.
 *  Unitless.
 * Mise en forme et ecriture du fichier HDF4.
 */
void write_hdf(char *fic_out,
               int nb_valid,
               float lo[NBANDS],
               float *ptr_array1d_ipix_lat,
               float *ptr_array1d_ipix_lon,
               short *ptr_array1d_ipix_indlat,
               short *ptr_array1d_ipix_indlon,
               float *array1d_ptr_array1d_iband_ipix_Rrs[NBANDS],
               float *array1d_ptr_array1d_iband_ipix_a[NBANDS],
               float *array1d_ptr_array1d_iband_ipix_bb[NBANDS],
               float *array1d_ptr_array1d_iband_ipix_Kd[NBANDS],
               float *ptr_array1d_ipix_Zeu_bar,
               float *ptr_array1d_ipix_poc,
               float *ptr_array1d_ipix_aCDOM_412,
               float *ptr_array1d_ipix_PAR_cloud,
               float *ptr_array1d_ipix_PAR_clear,
               float *ptr_array1d_ipix_CF00,
               float *ptr_array1d_ipix_O300,
               float *ptr_array1d_ipix_tauCld00,
               float *ptr_array1d_ipix_thetas00,
               float *ptr_array1d_ipix_CFmean,
               float *ptr_array1d_ipix_O3mean,
               float *ptr_array1d_ipix_tauCldmean,
               float *ptr_array1d_ipix_KPUR,
               float *ptr_array1d_ipix_chl_cota,
               float *ptr_array1d_ipix_chl_oc,
               float *ptr_array1d_ipix_chl_gsm_mustapha,
               float *ptr_array1d_ipix_PP_cota,
               float *ptr_array1d_ipix_PP_gsm_mustapha,
               float *ptr_array1d_ipix_PP_oc,
               float *ptr_array1d_ipix_aphy443,
               float *ptr_array1d_ipix_chl_gsm,
               float *ptr_array1d_ipix_PP_chl_gsm,
               float *ptr_array1d_ipix_ice){
  intn status;
  int32 status_32;
  int32 file_id, vdata_ref;
  int iband;
  char nomvar[250];
  if((file_id = Hopen(fic_out, DFACC_CREATE, 0)) == FAIL) { printf("ouverture\n"); HEprint(stdout,2); exit(-1);}
  status = Vstart(file_id);  
  vdata_ref = VHstoredata(file_id,
                          NAME_LAT,
                          (uint8 *)ptr_array1d_ipix_lat,
                          nb_valid,
                          DFNT_FLOAT32,
                          NAME_LAT,
                          NAME_LAT);
  vdata_ref = VHstoredata(file_id,
                          NAME_LON,
                          (uint8 *)ptr_array1d_ipix_lon,
                          nb_valid,
                          DFNT_FLOAT32,
                          NAME_LON,
                          NAME_LON);
  vdata_ref = VHstoredata(file_id,
                          NAME_INDLAT,
                          (uint8 *)ptr_array1d_ipix_indlat,
                          nb_valid,
                          DFNT_INT16,
                          NAME_INDLAT,
                          NAME_INDLAT);
  
  vdata_ref = VHstoredata(file_id,
                          NAME_INDLON,
                          (uint8 *)ptr_array1d_ipix_indlon,
                          nb_valid,
                          DFNT_INT16,
                          NAME_INDLON,
                          NAME_INDLON);
  vdata_ref = VHstoredata(file_id,
                          NAME_PP,
                          (uint8 *)ptr_array1d_ipix_PP_chl_gsm,
                          nb_valid,
                          DFNT_FLOAT32,
                          NAME_PP,
                          NAME_PP);
  vdata_ref = VHstoredata(file_id,
                          NAME_CHL_GSM,
                          (uint8 *)ptr_array1d_ipix_chl_gsm,
                          nb_valid,
                          DFNT_FLOAT32,
                          NAME_CHL_GSM,
                          NAME_CHL_GSM);
  vdata_ref = VHstoredata(file_id,
                          NAME_APHY443,
                          (uint8 *)ptr_array1d_ipix_aphy443,
                          nb_valid,
                          DFNT_FLOAT32,
                          NAME_APHY443,
                          NAME_APHY443);
  vdata_ref = VHstoredata(file_id,
                          NAME_ICE,
                          (uint8 *)ptr_array1d_ipix_ice,
                          nb_valid,
                          DFNT_FLOAT32,
                          NAME_ICE,
                          NAME_ICE);
  vdata_ref = VHstoredata(file_id, NAME_ZEU, (uint8*)ptr_array1d_ipix_Zeu_bar, nb_valid, DFNT_FLOAT32, NAME_ZEU, NAME_ZEU);
  vdata_ref = VHstoredata(file_id, NAME_POC, (uint8*)ptr_array1d_ipix_poc, nb_valid, DFNT_FLOAT32, NAME_POC, NAME_POC);
  vdata_ref = VHstoredata(file_id, NAME_ACDOM_412, (uint8*)ptr_array1d_ipix_aCDOM_412, nb_valid, DFNT_FLOAT32, NAME_ACDOM_412, NAME_ACDOM_412);
  vdata_ref = VHstoredata(file_id, NAME_PAR_CLOUD, (uint8*)ptr_array1d_ipix_PAR_cloud, nb_valid, DFNT_FLOAT32, NAME_PAR_CLOUD, NAME_PAR_CLOUD);
  vdata_ref = VHstoredata(file_id, NAME_PAR_CLEAR, (uint8*)ptr_array1d_ipix_PAR_clear, nb_valid, DFNT_FLOAT32, NAME_PAR_CLEAR, NAME_PAR_CLEAR);
  vdata_ref = VHstoredata(file_id, NAME_CF_00UTC, (uint8 *)ptr_array1d_ipix_CF00, nb_valid, DFNT_FLOAT32, NAME_CF_00UTC, NAME_CF_00UTC);
  vdata_ref = VHstoredata(file_id, NAME_O3_00UTC, (uint8 *)ptr_array1d_ipix_O300, nb_valid, DFNT_FLOAT32, NAME_O3_00UTC, NAME_O3_00UTC);
  vdata_ref = VHstoredata(file_id, NAME_TAUCLD_00UTC, (uint8 *)ptr_array1d_ipix_tauCld00, nb_valid, DFNT_FLOAT32, NAME_TAUCLD_00UTC, NAME_TAUCLD_00UTC);
  vdata_ref = VHstoredata(file_id, NAME_CF_MEAN, (uint8 *)ptr_array1d_ipix_CFmean, nb_valid, DFNT_FLOAT32, NAME_CF_MEAN, NAME_CF_MEAN);
  vdata_ref = VHstoredata(file_id, NAME_O3_MEAN, (uint8 *)ptr_array1d_ipix_O3mean, nb_valid, DFNT_FLOAT32, NAME_O3_MEAN, NAME_O3_MEAN);
  vdata_ref = VHstoredata(file_id, NAME_TAUCLD_MEAN, (uint8 *)ptr_array1d_ipix_tauCldmean, nb_valid, DFNT_FLOAT32, NAME_TAUCLD_MEAN, NAME_TAUCLD_MEAN);
  vdata_ref = VHstoredata(file_id, NAME_KPUR, (uint8 *)ptr_array1d_ipix_KPUR, nb_valid, DFNT_FLOAT32, NAME_KPUR, NAME_KPUR);
  vdata_ref = VHstoredata(file_id, NAME_CHL_COTA, (uint8 *)ptr_array1d_ipix_chl_cota, nb_valid, DFNT_FLOAT32, NAME_CHL_COTA, NAME_CHL_COTA);
  vdata_ref = VHstoredata(file_id, NAME_CHL_OC, (uint8 *)ptr_array1d_ipix_chl_oc, nb_valid, DFNT_FLOAT32, NAME_CHL_OC, NAME_CHL_OC);
  vdata_ref = VHstoredata(file_id, NAME_CHL_GSM_MUSTAPHA, (uint8 *)ptr_array1d_ipix_chl_gsm_mustapha, nb_valid, DFNT_FLOAT32, NAME_CHL_GSM_MUSTAPHA, NAME_CHL_GSM_MUSTAPHA);
  vdata_ref = VHstoredata(file_id, NAME_PP_COTA, (uint8 *)ptr_array1d_ipix_PP_cota, nb_valid, DFNT_FLOAT32, NAME_PP_COTA, NAME_PP_COTA);
  vdata_ref = VHstoredata(file_id, NAME_PP_GSM_MUSTAPHA, (uint8 *)ptr_array1d_ipix_PP_gsm_mustapha, nb_valid, DFNT_FLOAT32, NAME_PP_GSM_MUSTAPHA, NAME_PP_GSM_MUSTAPHA);
  vdata_ref = VHstoredata(file_id, NAME_PP_OC, (uint8 *)ptr_array1d_ipix_PP_oc, nb_valid, DFNT_FLOAT32, NAME_PP_OC, NAME_PP_OC);
  for(iband = 0; iband < NBANDS; iband++){
    sprintf(nomvar,FORMAT_A,(int)lo[iband]);
    vdata_ref = VHstoredata(file_id,
                            nomvar,
                            (uint8 *)array1d_ptr_array1d_iband_ipix_a[iband],
                            nb_valid,
                            DFNT_FLOAT32,
                            nomvar,
                            nomvar);
    sprintf(nomvar,FORMAT_BB,(int)lo[iband]);
    vdata_ref = VHstoredata(file_id,
                            nomvar,
                            (uint8 *)array1d_ptr_array1d_iband_ipix_bb[iband],
                            nb_valid,
                            DFNT_FLOAT32,
                            nomvar,
                            nomvar);
    sprintf(nomvar,FORMAT_KD,(int)lo[iband]);
    vdata_ref = VHstoredata(file_id,
                            nomvar,
                            (uint8 *)array1d_ptr_array1d_iband_ipix_Kd[iband],
                            nb_valid,
                            DFNT_FLOAT32,
                            nomvar,
                            nomvar);
    sprintf(nomvar,FORMAT_RRS,(int)lo[iband]);
    vdata_ref = VHstoredata(file_id,
                            nomvar,
                            (uint8 *)array1d_ptr_array1d_iband_ipix_Rrs[iband],
                            nb_valid,
                            DFNT_FLOAT32,
                            nomvar,
                            nomvar);
  }
  status = Vend(file_id);
  status_32 = Hclose(file_id);
}

/* -------------------------------------------------------
 * fic_out:
 *  Nom du fichier de sortie.
 * nb_valid:
 *  Nombre de pixels a ecrire.
 * lo: 
 *  Tableau de dimension 6.			      
 *  Contient les longueurs d'onde du satellite en nm.
 * ptr_array1d_ipix_lat:
 *  Array of dimension nbpix.
 *  The first dimension is the pixel.
 *  The values are the latitudes from -90 to 90.
 *  The units are degrees N.
 * ptr_array1d_ipix_lon:
 *  Array of dimension nbpix.
 *  The first dimension is the pixel.
 *  The values are the latitudes from -180 to 180.
 *  The units are degrees E.
 * ptr_array1d_ipix_indlat:
 *  Array of dimension nbpix.
 *  The first dimension is the pixel.
 *  The values are the row number on the L3BIN grid.
 *  Corresponds to ROW - 1 in
 *  http://oceancolor.gsfc.nasa.gov/SeaWiFS/TECH_REPORTS/PreLPDF/PreLVol32.pdf
 *  appendix A.
 * ptr_array1d_ipix_indlon:
 *  Array of dimension nbpix.
 *  The first dimension is the pixel.
 *  The values are the bin number on this row.
 *  Corresponds to COL in
 *  http://oceancolor.gsfc.nasa.gov/SeaWiFS/TECH_REPORTS/PreLPDF/PreLVol32.pdf
 *  appendix A.
 *  From 1 to the total number of possible bins in this row.
 * array1d_ptr_array1d_iband_ipix_Rrs:
 *  Array of dimensions 6 X nbpix.
 *  (Array of pointers.)
 *  The first dimension is the band of the satellite.
 *  The second dimension is the pixel.
 *  The values are the remote-sensing reflectances (Rrs).
 *  The units are sr^-1.
 * array1d_ptr_array1d_iband_ipix_a:
 *  Array of dimensions 6 X nbpix.
 *  (Array of pointers.)
 *  The first dimension is the band of the satellite.
 *  The second dimension is the pixel.
 *  The values are the absorption coefficients.
 *  The units are m^-1.
 * array1d_ptr_array1d_iband_ipix_bb:
 *  Array of dimensions 6 X nbpix.
 *  (Array of pointers.)
 *  The first dimension is the band of the satellite.
 *  The second dimension is the pixel.
 *  The values are the backscattering coefficients (bb).
 *  The units are m^-1.
 * array1d_ptr_array1d_iband_ipix_Kd:
 *  Array of dimensions 6 X nbpix.				       
 *  (Array of pointers.)					       
 *  The first dimension is the band of the satellite.		       
 *  The second dimension is the pixel.				       
 *  The values are the vertical attenuation coefficients (Kd) for one 
 *  wavelength and one pixel at the local noon.
 *  The units are m^-1.                                                   
 * ptr_array1d_ipix_Zeu_bar:
 *  Array of dimension nbpix.						     
 *  The first dimension is the pixel.
 *  The values are the mean euphotic depths.
 *  The units are m.
 * ptr_array1d_ipix_poc:
 *  Array of dimension nbpix.			   
 *  The first dimension is the pixel.		   
 *  The values are the particulate organic carbon.                          
 * ptr_array1d_ipix_aCDOM_412
 *  Array of dimension nbpix.						        
 *  The columns contain the absorption coefficient for CDOM at 412 nm (m^-1).
 * ptr_array1d_ipix_PAR_cloud:
 *  Array of dimension nbpix.						     
 *  The first dimension is the pixel.					     
 *  The values are the photosynthetically available irradiances (PAR) for a 
 *  pixel on a given day.						     
 *  The units are uEinstent*m^-1.                                           
 * float ptr_array1d_ipix_CF00:
 *  Array of dimension nbpix.				      
 *  The first dimension is the pixel.			      
 *  The values are the cloud fractions (0 to 1) at 00h UTC.  
 *  The values are unitless.                                 
 * ptr_array1d_ipix_O300:
 *  Array of dimension nbpix.			      
 *  The first dimension is the pixel.		      
 *  The values are the total ozone column at 00h UTC.
 *  The units are Dobson units.                      
 * ptr_array1d_ipix_tauCld00:
 *  Array of dimension nbpix.				      
 *  The first dimension is the pixel.			      
 *  The values are the cloud optical thickness at 00h UTC.   
 *  The values are unitless.                                 
 * ptr_array1d_ipix_thetas00:
 *  Array of dimension nbpix.				
 *  The first dimension is the pixel.			
 *  The values are the solar zenith angles at 00h UTC.	
 *  The values are in degrees.                               
 * ptr_array1d_ipix_CFmean:
 *  Array of dimension nbpix.				      
 *  The first dimension is the pixel.			      
 *  The values are the mean cloud fractions (0 to 1).        
 *  The values are unitless.                                 
 * ptr_array1d_ipix_O300:
 *  Array of dimension nbpix.			      
 *  The first dimension is the pixel.		      
 *  The values are the total ozone column at 00h UTC.
 *  The units are Dobson units.                      
 * ptr_array1d_ipix_tauCldmean:
 *  Array of dimension nbpix.				      
 *  The first dimension is the pixel.			      
 *  The values are the mean cloud optical thickness 	      
 *  The values are unitless.                                 
 * ptr_array1d_ipix_KPUR:
 *  Array of dimension nbpix.			      			   
 *  The first dimension is the pixel.		      			   
 *  The values are the vertical attenuation of photosynthetically usable  
 *  radiation.								   
 *  The units are m^-1.                                                   
 * ptr_array1d_ipix_chl_cota:
 *  Array of dimension nbpix.			      			     
 *  The first dimension is the pixel.		      			     
 *  The values are the chlorophyll-a concentration from Cota et al. 2004.
 *  The units are (mg Chl-a m^-3).
 * ptr_array1d_ipix_chl_oc:
 *  Array of dimension nbpix.			      			     
 *  The first dimension is the pixel.		      			     
 *  The values are the chlorophyll-a concentration from OC3 or OC4.	   
 *  The units are (mg Chl-a m^-3).
 * ptr_array1d_ipix_chl_gsm_mustapha:
 *  Array of dimension nbpix.			      			       
 *  The first dimension is the pixel.		      			       
 *  The values are the chlorophyll-a concentration from Ben Mustapha et       
 *  al. 2012.								       
 *  The units are (mg Chl-a m^-3).
 * ptr_array1d_ipix_PP_cota:
 *  Array of dimension nbpix.						     
 *  The first dimension is the pixel.					     
 *  The values are the primary productivity computed with Cota et al. 2004. 
 *  The units are mgC*m^-2*d^-1.                                            
 * ptr_array1d_ipix_PP_gsm_mustapha:
 *  Array of dimension nbpix.						       
 *  The first dimension is the pixel.					       
 *  The values are the primary productivity computed with Ben Mustapha et al. 
 *  2012.								       
 *  The units are mgC*m^-2*d^-1.                                              
 * ptr_array1d_ipix_PP_oc:
 *  Array of dimension nbpix.						     
 *  The first dimension is the pixel.					     
 *  The values are the primary productivity computed with OC3 or OC4.
 * ptr_array1d_ipix_aphy443:
 *  Array of dimension nbpix.
 *  The first dimension is the pixel.
 *  The values are phytoplankton absorption coefficient at 443 nm (a_phy at
 *  443 nm).
 *  The units are m^-1.
 * ptr_array1d_ipix_chl_gsm:
 *  Array of dimension nbpix.
 *  The first dimension is the pixel.
 *  The values are the chlorophyll-a concentration from GSM.
 *  The units are (mg Chl-a m^-3).
 * ptr_array1d_ipix_PP_chl_gsm:
 *  Array of dimension nbpix.
 *  The first dimension is the pixel.
 *  The values are the primary productivity computed with GSM.
 *  The units are mgC*m^-2*d^-1.
 * ptr_array1d_ipix_ice:
 *  Array of dimension nbpix.
 *  The first dimension is the pixel.
 *  The values are the sea-ice concentration.
 *  Normal values from 0 to 1.
 *  Special values > 1.
 *  Unitless.
 * Prepare the format and write the netCDF file.
 */
void write_netcdf(char *fic_out,
                  int nb_valid,
                  float lo[NBANDS],
                  float *ptr_array1d_ipix_lat,
                  float *ptr_array1d_ipix_lon,
                  short *ptr_array1d_ipix_indlat,
                  short *ptr_array1d_ipix_indlon,
                  float *array1d_ptr_array1d_iband_ipix_Rrs[NBANDS],
                  float *array1d_ptr_array1d_iband_ipix_a[NBANDS],
                  float *array1d_ptr_array1d_iband_ipix_bb[NBANDS],
                  float *array1d_ptr_array1d_iband_ipix_Kd[NBANDS],
                  float *ptr_array1d_ipix_Zeu_bar,
                  float *ptr_array1d_ipix_poc,
                  float *ptr_array1d_ipix_aCDOM_412,
                  float *ptr_array1d_ipix_PAR_cloud,
                  float *ptr_array1d_ipix_PAR_clear,
                  float *ptr_array1d_ipix_CF00,
                  float *ptr_array1d_ipix_O300,
                  float *ptr_array1d_ipix_tauCld00,
                  float *ptr_array1d_ipix_thetas00,
                  float *ptr_array1d_ipix_CFmean,
                  float *ptr_array1d_ipix_O3mean,
                  float *ptr_array1d_ipix_tauCldmean,
                  float *ptr_array1d_ipix_KPUR,
                  float *ptr_array1d_ipix_chl_cota,
                  float *ptr_array1d_ipix_chl_oc,
                  float *ptr_array1d_ipix_chl_gsm_mustapha,
                  float *ptr_array1d_ipix_PP_cota,
                  float *ptr_array1d_ipix_PP_gsm_mustapha,
                  float *ptr_array1d_ipix_PP_oc,
                  float *ptr_array1d_ipix_aphy443,
                  float *ptr_array1d_ipix_chl_gsm,
                  float *ptr_array1d_ipix_PP_chl_gsm,
                  float *ptr_array1d_ipix_ice){
  /* identifiers */
  int file_id, ipixdim_id;
  int lat_id, lon_id, indlat_id, indlon_id, PP_id, chl_gsm_id;
  int aphy443_id, Ice_id;
  int Zeu_id, poc_id, aCDOM_412_id, PAR_clear_id, PAR_cloud_id;
  int CF_00UTC_id, O3_00UTC_id, TauCld_00UTC_id;
  int CF_mean_id, O3_mean_id, TauCld_mean_id;
  int KPUR_id, chl_cota_id, chl_oc_id, chl_gsm_mustapha_id;
  int PP_cota_id, PP_gsm_mustapha_id, PP_oc_id;
  int a_id[NBANDS], bb_id[NBANDS], Kd_id[NBANDS], Rrs_id[NBANDS];
  /* sizes */
  int ipixdims[1];
  /* status */
  int status;
  /* Other variables. */
  int iband; // The index of the wavelength.
  char nomvar[250]; // The name of the variable.
  /* Create a new file - clobber anything existing. */
  /*
   * Signification of |NC_NETCDF4|NC_CLASSIC_MODEL
   * It will be in the HDF5 format, but will not allow the use of any 
   * netCDF-4 advanced features. That is, it will conform to the classic 
   * netCDF-3 data model.
   */
  status = nc_create(fic_out,
		     NC_CLOBBER|NC_NETCDF4|NC_CLASSIC_MODEL,
		     &file_id);
  if(status != NC_NOERR){
    fprintf(stderr, "Could not open file %s\n", fic_out);
    exit(-1);
  }
  /* Define the dimensions. */
  nc_def_dim(file_id, "ipix", nb_valid, &ipixdim_id);
  /* Now that we've defined the dimensions, we can define variables on them. */
  ipixdims[0] = ipixdim_id;
  nc_def_var(file_id, NAME_LAT, NC_FLOAT, 1, ipixdims, &lat_id);
  nc_def_var(file_id, NAME_LON, NC_FLOAT, 1, ipixdims, &lon_id);
  nc_def_var(file_id, NAME_INDLAT, NC_SHORT, 1, ipixdims, &indlat_id);
  nc_def_var(file_id, NAME_INDLON, NC_SHORT, 1, ipixdims, &indlon_id);
  nc_def_var(file_id, NAME_PP, NC_FLOAT, 1, ipixdims, &PP_id);
  nc_def_var(file_id, NAME_CHL_GSM, NC_FLOAT, 1, ipixdims, &chl_gsm_id);
  nc_def_var(file_id, NAME_APHY443, NC_FLOAT, 1, ipixdims, &aphy443_id);
  nc_def_var(file_id, NAME_ICE, NC_FLOAT, 1, ipixdims, &Ice_id);
  nc_def_var(file_id, NAME_ZEU, NC_FLOAT, 1, ipixdims, &Zeu_id);
  nc_def_var(file_id, NAME_POC, NC_FLOAT, 1, ipixdims, &poc_id);
  nc_def_var(file_id, NAME_ACDOM_412, NC_FLOAT, 1, ipixdims, &aCDOM_412_id);
  nc_def_var(file_id, NAME_PAR_CLOUD, NC_FLOAT, 1, ipixdims, &PAR_cloud_id);
  nc_def_var(file_id, NAME_PAR_CLEAR, NC_FLOAT, 1, ipixdims, &PAR_clear_id);
  nc_def_var(file_id, NAME_CF_00UTC, NC_FLOAT, 1, ipixdims, &CF_00UTC_id);
  nc_def_var(file_id, NAME_O3_00UTC, NC_FLOAT, 1, ipixdims, &O3_00UTC_id);
  nc_def_var(file_id, NAME_TAUCLD_00UTC, NC_FLOAT, 1, ipixdims,
	     &TauCld_00UTC_id);
  nc_def_var(file_id, NAME_KPUR, NC_FLOAT, 1, ipixdims, &KPUR_id);
  nc_def_var(file_id, NAME_CHL_COTA, NC_FLOAT, 1, ipixdims, &chl_cota_id);
  nc_def_var(file_id, NAME_CHL_OC, NC_FLOAT, 1, ipixdims, &chl_oc_id);
  nc_def_var(file_id, NAME_CHL_GSM_MUSTAPHA, NC_FLOAT, 1, ipixdims, 
	     &chl_gsm_mustapha_id);
  nc_def_var(file_id, NAME_CF_MEAN, NC_FLOAT, 1, ipixdims, &CF_mean_id);
  nc_def_var(file_id, NAME_O3_MEAN, NC_FLOAT, 1, ipixdims, &O3_mean_id);
  nc_def_var(file_id, NAME_TAUCLD_MEAN, NC_FLOAT, 1, ipixdims, &TauCld_mean_id);
  nc_def_var(file_id, NAME_PP_COTA, NC_FLOAT, 1, ipixdims, &PP_cota_id);
  nc_def_var(file_id, NAME_PP_GSM_MUSTAPHA, NC_FLOAT, 1, ipixdims,
	     &PP_gsm_mustapha_id);
  nc_def_var(file_id, NAME_PP_OC, NC_FLOAT, 1, ipixdims, &PP_oc_id);
  for(iband = 0; iband < NBANDS; iband++){
    sprintf(nomvar,FORMAT_A,(int)lo[iband]);
    nc_def_var(file_id, nomvar, NC_FLOAT, 1, ipixdims, &(a_id[iband]));
    sprintf(nomvar,FORMAT_BB,(int)lo[iband]);
    nc_def_var(file_id, nomvar, NC_FLOAT, 1, ipixdims, &(bb_id[iband]));
    sprintf(nomvar,FORMAT_KD,(int)lo[iband]);
    nc_def_var(file_id, nomvar, NC_FLOAT, 1, ipixdims, &(Kd_id[iband]));
    sprintf(nomvar,FORMAT_RRS,(int)lo[iband]);
    nc_def_var(file_id, nomvar, NC_FLOAT, 1, ipixdims, &(Rrs_id[iband]));
  }
  /* We are now done defining variables and their attributes. */
  nc_enddef(file_id);
  /* Write out the data in the variables we've defined. */
  nc_put_var_float(file_id, lat_id, ptr_array1d_ipix_lat);
  nc_put_var_float(file_id, lon_id, ptr_array1d_ipix_lon);
  nc_put_var_short(file_id, indlat_id, ptr_array1d_ipix_indlat);
  nc_put_var_short(file_id, indlon_id, ptr_array1d_ipix_indlon);
  nc_put_var_float(file_id, PP_id, ptr_array1d_ipix_PP_chl_gsm);
  nc_put_var_float(file_id, chl_gsm_id, ptr_array1d_ipix_chl_gsm);
  nc_put_var_float(file_id, aphy443_id, ptr_array1d_ipix_aphy443);
  nc_put_var_float(file_id, Ice_id, ptr_array1d_ipix_ice);
  nc_put_var_float(file_id, Zeu_id, ptr_array1d_ipix_Zeu_bar);
  nc_put_var_float(file_id, poc_id, ptr_array1d_ipix_poc);
  nc_put_var_float(file_id, aCDOM_412_id, ptr_array1d_ipix_aCDOM_412);
  nc_put_var_float(file_id, PAR_cloud_id, ptr_array1d_ipix_PAR_cloud);
  nc_put_var_float(file_id, PAR_clear_id, ptr_array1d_ipix_PAR_clear);
  nc_put_var_float(file_id, CF_00UTC_id, ptr_array1d_ipix_CF00);
  nc_put_var_float(file_id, O3_00UTC_id, ptr_array1d_ipix_O300);
  nc_put_var_float(file_id, TauCld_00UTC_id, ptr_array1d_ipix_tauCld00);
  nc_put_var_float(file_id, CF_mean_id, ptr_array1d_ipix_CFmean);
  nc_put_var_float(file_id, O3_mean_id, ptr_array1d_ipix_O3mean);
  nc_put_var_float(file_id, TauCld_mean_id, ptr_array1d_ipix_tauCldmean);
  nc_put_var_float(file_id, KPUR_id, ptr_array1d_ipix_KPUR);
  nc_put_var_float(file_id, chl_cota_id, ptr_array1d_ipix_chl_cota);
  nc_put_var_float(file_id, chl_oc_id, ptr_array1d_ipix_chl_oc);
  nc_put_var_float(file_id, chl_gsm_mustapha_id,
		   ptr_array1d_ipix_chl_gsm_mustapha);
  nc_put_var_float(file_id, PP_cota_id, ptr_array1d_ipix_PP_cota);
  nc_put_var_float(file_id, PP_gsm_mustapha_id,
		   ptr_array1d_ipix_PP_gsm_mustapha);
  nc_put_var_float(file_id, PP_oc_id, ptr_array1d_ipix_PP_oc);
  for(iband = 0; iband < NBANDS; iband++){
    sprintf(nomvar, FORMAT_A, (int)lo[iband]);
    nc_put_var_float(file_id,
		     a_id[iband],
		     array1d_ptr_array1d_iband_ipix_a[iband]);
    sprintf(nomvar, FORMAT_BB, (int)lo[iband]);
    nc_put_var_float(file_id,
		     bb_id[iband],
		     array1d_ptr_array1d_iband_ipix_bb[iband]);
    sprintf(nomvar, FORMAT_KD, (int)lo[iband]);
    nc_put_var_float(file_id,
		     Kd_id[iband],
		     array1d_ptr_array1d_iband_ipix_Kd[iband]);
    sprintf(nomvar, FORMAT_RRS, (int)lo[iband]);
    nc_put_var_float(file_id,
		     Rrs_id[iband],
		     array1d_ptr_array1d_iband_ipix_Rrs[iband]);
  }
  nc_close(file_id);
  return;
}

/* ================================= MAIN ================================= */
int main(int argc, char *argv[]){
  
  /////////// Declaration of the variables. ///////////
  
  /*
   * Absorption coefficient for CDOM at 412 nm (m^-1).
   */
  float aCDOM_412;
  /* Le coefficient d'absorption du phytoplancton a 443 nm. */
  float aphy443;
  /*
   * Vaut 1 s'il y a des pixels sur la ligne.
   * Vaut 0 sinon.
   */
  int are_pixels_on_current_line = 0;
  /*
   * Array of dimension NBANDS = 6.					       
   * The first dimension are the wavelengths of the bands of the satellite.    
   * The values are the total absorption coefficients (a or a_t).	       
   * The units are m^-1.                                                       
   */
  float array1d_iband_a[NBANDS];
  /*
   * Array of dimension NBANDS = 6.					       
   * The first dimension are the wavelengths of the bands of the satellite.    
   * The values are the total backscattering coefficients (bb or b_bt) at      
   * the surface (z = 0).						       
   * The units are m^-1.                                                       
   */
  float array1d_iband_bb[NBANDS];
  /*
   * Array of dimension NBANDS = 6.					       
   * The first dimension are the wavelengths of the bands of the satellite.    
   * The values are the diffuse attenuation coefficient (K_d) at 
   * the surface (z = 0).
   * The units are m^-1.                                                       
   */
  float array1d_iband_Kd[NBANDS];
  /*
   * Array of dimension NBANDS = 6.					       
   * The first dimension are the wavelengths of the bands of the satellite.    
   * The values are the remote-sensing reflectances (Rrs).
   * The units are sr^-1.                                                      
   */
  float array1d_iband_Rrs[NBANDS];
  /*
   * Array of dimensions 12.
   * The first dimension are the optical depths. 			     
   * The values are the saturation irradiances.
   * The units are umol photons*m^-2*s^-1.
   */
  float array1d_idepth_Ek[NBDEPTHS];
  /*
   * Array of dimensions 12.
   * The first dimension are the optical depths. 			     
   * The values are the daily mean of the photosynthetically usable radiation.
   * The units are umol photons*m^-2*s^-1.
   */
  float array1d_idepth_meanPUR[NBDEPTHS];
  /*
   * Array of dimensions 9.
   * The first dimension is the hours UTC from 0 to 24 by step of 3 h.
   * The values are the names of the atmosphere files.
   */
  char* array1d_ifile_atmosphere_file[NBHOURS];
  /*
   * Array of dimension SZLAT_MODISA (the number of rows in the MODISA grid).
   * The first dimension is the index (0-based) of the row in the MODISA grid.
   * The value is the index (1-based) of the bin number (bin_num) of the first
   * documented bin in the current row.
   * From 1 to 23761676.
   * Fill value: 0.
   * Unitless.
   */
  int array1d_irow_begin_MODISA[SZLAT_MODISA];
  /*
   * Array of dimension SZLAT_SEAWIFS (the number of rows in the SeaWiFS grid).
   * The first dimension is the index (0-based) of the row in the SeaWiFS grid.
   * The value is the index (1-based) of the bin number (bin_num) of the first
   * documented bin in the current row.
   * From 1 to 5940422.
   * Fill value: 0.
   * Unitless.
   */
  int array1d_irow_begin_SeaWiFS[SZLAT_SEAWIFS];
  /*
   * Array of dimension SZLAT_MODISA (the number of rows in the MODISA grid).
   * The first dimension is the index (0-based) of the row in the MODISA grid.
   * The value is the number of documented bins in the current row.
   * From 0 to the total number of possible bins in each row.
   * The units are bins.
   */
  int array1d_irow_extent_MODISA[SZLAT_MODISA];
  /*
   * Array of dimension SZLAT_SEAWIFS (the number of rows in the SeaWiFS grid).
   * The first dimension is the index (0-based) of the row in the SeaWiFS grid.
   * The value is the number of documented bins in the current row.
   * From 0 to the total number of possible bins in each row.
   * The units are bins.
   */
  int array1d_irow_extent_SeaWiFS[SZLAT_SEAWIFS];
  /*
   * Array of dimension SZLAT_MODISA (the number of rows in the MODISA grid).
   * The first dimension is the index (0-based) of the row in the MODISA grid.
   * The value is the number of possible bins in the current row.
   * From 3 to 8640.
   * The units are bins.
   */
  int array1d_irow_max_MODISA[SZLAT_MODISA];
  /*
   * Array of dimension SZLAT_SEAWIFS (the number of rows in the SeaWiFS grid).
   * The first dimension is the index (0-based) of the row in the SeaWiFS grid.
   * The value is the number of possible bins in the current row.
   * From 3 to 8640.
   * The units are bins.
   */
  int array1d_irow_max_SeaWiFS[SZLAT_SEAWIFS];
  /*
   * Array of dimension SZLAT_MODISA (the number of rows in the MODISA grid).
   * The first dimension is the index (0-based) of the row in the MODISA grid.
   * The value is the index (1-based) of the bin number (bin_num) of the first
   * bin in the current row.
   * From 1 to 4320.
   * Unitless.
   */
  int array1d_irow_start_num_MODISA[SZLAT_MODISA];
  /*
   * Array of dimension SZLAT_SEAWIFS (the number of rows in the SeaWiFS grid).
   * The first dimension is the index (0-based) of the row in the SeaWiFS grid.
   * The value is the index (1-based) of the bin number (bin_num) of the first
   * bin in the current row.
   * From 1 to 5940420.
   * Unitless.
   */
  int array1d_irow_start_num_SeaWiFS[SZLAT_SEAWIFS];
  /*
   * Array of dimensions 12.
   * The first dimension is the hours UTC from 0 to 24 by step of 3 h.
   * The values are vertical attenuation of photosynthetically usable
   * radiation at local zenith (local noon).
   * The units are m^-1.
   */
  float array1d_itime_KPUR[NBHOURS];
  /*
   * Array of dimension NBANDS = 6.					       
   * The first dimension are the wavelength from 400 to 700 nm by step 5.
   * The values are the diffuse attenuation coefficient (K_d) at 
   * the surface (z = 0).
   * The units are m^-1.                                                       
   */
  float array1d_ivis_Kd[NVIS];
  /*
   * Array of dimension 61.
   * The first dimension is the index of the wavelength.
   * The values are the wavelengths from 400 to 700 by step 5.
   * The units are nm.
   */
  float array1d_ivis_lambda[NVIS];
  /*
   * Array of dimensions NBWL = 83.
   * The first dimension are the wavelengths from 290 to 700 by step 5 nm.
   * The values are the downward irradiance just below the water surface.
   * The units are umol photons*m^-2*s^-1*nm^-1. 
   */
  float array1d_iwl_Ed0minus[NBWL];
  /*
   * Array of dimensions NBWL = 83.
   * The first dimension are the wavelengths from 290 to 700 by step 5 nm.
   * The values are the downward irradiance just below the water surface
   * with the cloud fraction forced to 0.
   * The units are umol photons*m^-2*s^-1*nm^-1. 
   */
  float array1d_iwl_Ed0minus_clear[NBWL];
  /*
   * Array of dimensions 6 * nbinimage
   * The first dimension is band of the sensor.
   * The second dimension is the bins in the image. Note that not all bins of
   * the grid are in the image.
   * The values are the remote-sensing reflectances (Rrs).
   * The units are sr^-1.
   */
  float* array1d_ptr_array1d_iband_ibinimage_Rrs[NBANDS];
  /*
   * Array of dimensions 6 X nbpix.
   * (Array of pointers.)
   * The first dimension is the band of the satellite.
   * The second dimension is the pixel.
   * The values are the absorption coefficients.
   * The units are m^-1.
   */
  float *array1d_ptr_array1d_iband_ipix_a[NBANDS];
  /*
   * Array of dimensions 6 X nbpix.
   * (Array of pointers.)
   * The first dimension is the band of the satellite.
   * The second dimension is the pixel.
   * The values are the backscattering coefficients (bb).
   * The units are m^-1.
   */
  float *array1d_ptr_array1d_iband_ipix_bb[NBANDS];
  /*
   * Array of dimensions 6 X nbpix.
   * (Array of pointers.)
   * The first dimension is the band of the satellite.
   * The second dimension is the pixel.
   * The values are the vertical attenuation coefficients (Kd) for one
   * wavelength and one pixel at the local noon.
   * The units are m^-1.
   */
  float *array1d_ptr_array1d_iband_ipix_Kd[NBANDS];
  /*
   * Array of dimensions 6 X nbpix.
   * (Array of pointers.)
   * The first dimension is the band of the satellite.
   * The second dimension is the pixel.
   * The values are the remote-sensing reflectances (Rrs).
   * The units are sr^-1.
   */
  float *array1d_ptr_array1d_iband_ipix_Rrs[NBANDS];
  /*
   * Tableau de dimensions 9.
   * La premiere dimension est les heures.
   * Les valeurs sont les profondeurs de la zone euphotique (en m).
   */
  float array1d_itime_Zeu[NBHOURS];
  /*
   * Array of dimension 6 * 9.					       	       
   * The first dimension are the wavelengths of the bands of the satellite.    
   * The second dimension is the hours UTC from 0 to 24 by step of 3 h.	       
   * The values are the diffuse attenuation coefficient (K_d) at 	       
   * the surface (z = 0).						       
   * The units are m^-1.                                                       
   */
  float array2d_iband_itime_Kd[NBANDS][NBHOURS];
  /*
   * Array of dimensions 9 * 61.
   * The first dimension is the hours from 0 to 24 by step of 3 h.
   * The second dimension is the wavelenght from 400 nm to 700 nm by step of 5
   * nm.
   * The values are the downward irradiances just below the surface water
   * (Ed0moins) for one pixel with a cloud fraction forced to 0.
   * The units are uEinstein*m^-2*s^-1*nm^-1.
   */
  float array2d_h_ivis_Ed_pixel_clear[NBHOURS][NVIS];
  /*
   * Array of dimensions 9 * 61.
   * The first dimension is the hours from 0 to 24 by step of 3 h.
   * The second dimension is the wavelenght from 400 nm to 700 nm by step of 5
   * nm.
   * The values are the downward irradiances just below the surface water
   * (Ed0moins) for one pixel with the real cloud fraction.
   * The units are uEinstein*m^-2*s^-1*nm^-1.
   */
  float array2d_h_ivis_Ed_pixel_cloud[NBHOURS][NVIS];
  /*
   * Array of dimensions 46 x 360.
   * The first dimension is the latitude from 89.5 N to 44.5 N.
   * The second dimension is the longitude from -179.5 E to 179.5 E.
   * The values are the cloud fraction day mean read on the
   * MODIS-Atmosphere file north of 45 degrees North.
   * No units.
   * Valid values: 0 to 1.
   * Fill value: -999.
   */
  float array2d_ilat_ilon_cfdm[NBLAT_MODIS_ATMOSPHERE][NBLON_MODIS_ATMOSPHERE];
  /*
   * Array of dimensions 46 x 360.
   * The first dimension is the latitude from 89.5 N to 44.5 N.
   * The second dimension is the longitude from -179.5 E to 179.5 E.
   * The values are the cloud optical thickness combined mean read on the
   * MODIS-Atmosphere file north of 45 degrees North.
   * No units.
   * Fill value: -999.
   */
  float array2d_ilat_ilon_cotcm[NBLAT_MODIS_ATMOSPHERE][NBLON_MODIS_ATMOSPHERE];
  /* Tableau des coefficients de retrodiffusion (m^-1). */
  /*
   * Array of dimensions 46 x 360.
   * The first dimension is the latitude from 89.5 N to 44.5 N.
   * The second dimension is the longitude from -179.5 E to 179.5 E.
   * The values are the total ozone mean read on the
   * MODIS-Atmosphere file north of 45 degrees North.
   * Units: Dobson units.
   * Fill value: -999.
   */
  float array2d_ilat_ilon_to3m[NBLAT_MODIS_ATMOSPHERE][NBLON_MODIS_ATMOSPHERE];
  /*
   * Array of dimensions 12 * 9.
   * The first dimension is the index of the optical depths.
   * The second dimension is the hours from 0 to 24 by step of 3 h.
   * The values are the photosynthetically usable radiation.
   * The units are umol photons*m^-2*s^-1.
   */
  float array2d_idepth_itime_PUR[NBDEPTHS][NBHOURS];
  /*
   * Array of dimensions 12 * 9.
   * The first dimension is the index of the optical depths.
   * The second dimension is the hours from 0 to 24 by step of 3 h.
   * The values are the geometrical depths.
   * The units are m.
   */
  float array2d_idepth_itime_Z[NBDEPTHS][NBHOURS];
  /*
   * Array of dimensions NBDEPTHS * NVIS = 12 * 61.			     
   * The first dimension are the optical depths. 			     
   * The second dimension are the wavelengths from 400 nm to 700 nm by step 
   * of 5. Units: nm.							    
   * The values are the scalar irradiances.				    
   * The units are umol photons.m^-2.s^-1.nm^-1.                            
   */
  float array2d_ivis_idepth_E0[NVIS][NBDEPTHS];
  /*
   * Array of dimensions 448 * 304.
   * The first dimension is the line.
   * The second dimension is the column.
   * The values are the ice concentrations.
   * Concentration de glace de 0 a 1.02 (sans unites).
   * Les valeurs de 0 a 1 sont des concentrations de glace.
   * 1.004 (251/250.0) est le masque circulaire utilise en Arctique pour
   * couvrir le manque de donnees de forme irreguliere autour du pole (cause
   * par l'inclinaison de l'orbite et le swath de l'instrument).
   * 1.008 (252/250.0) est une valeur inutilisee.
   * 1.012 (253/250.0) est une ligne de cote.
   * 1.016 (254/250.0) est un masque de terre superpose.
   * 1.02 (255/250.0) est une donnee manquante.
   */
  float array2d_jTNB_iTNB_ice[NBLINE][NBCOL];
  /*
   * File for the atmospheric products for one day.
   * This argument is only valid with atmosphere_type = M.
   */
  char* atmosphere_file = NULL;
  /*
   * The type for the atmospheric products.
   * isccp or modis_atmosphere.
   * See get_args.h.
   */
  atmosphere_type_type atmosphere_type;
  /*
   * Index number of the first bin in the current row.
   * From 1 to the number of possible bins in the grid.
   */
  int basebin;
  /*
   * Index number of the first bin in the previous row.
   * From 1 to the number of possible bins in the grid.
   */
  int basebinM1;
  // Return of getopt_long.
  int c;
  /*
   * Tableau de dimensions 9 x 966.
   * Ce tableau contient les valeurs des produits atmospheriques lus dans
   * les fichiers du ISCCP sur notre sous-grille de de la grille EQ du
   * ISCCP au nord de 45 degres Nord.
   * Il y a 966 valeurs pour chaque heure.
   * O3 est l'ozone.
   * CF est la fraction nuageuse.
   * TauCld est l'epaisseur optique.
   */
  float CF[NBHOURS][NB_PIXELS_ISCCP];
  /*
   * The mean cloud fraction for a time interval of one day.
   * With ISCCP, it is the mean of the cloud fraction at 00UTC, 03UTC, 06UTC,
   * 09UTC, 12UTC, 15UTC, 18UTC and 21UTC.
   * With MODIS-Atmosphere, it is directly the cloud fraction for that day.
   * There is no mean to do in that case.
   * Unitless.
   */
  float CFmean;
  /*
   * The chlorophyll-a concentration from Cota et al. 2004.
   * Units:  (mg Chl-a m^-3).
   */
  float chl_cota;
  /* 
   * The chlorophyll-a concentration from GSM (Garver and Siegel 1997. 
   * Maritorena et al. 2002.)
   * Units:  (mg Chl-a m^-3).
   */
  float chl_gsm;
  /*
   * The chlorophyll-a concentration from Ben Mustapha et al. 2012.
   * Units:  (mg Chl-a m^-3).
   */
  float chl_gsm_mustapha;
  /*
   * Chlorophyll-a concentration from OC3 or OC4.
   * Units:  (mg Chl-a m^-3).
   */
  float chl_oc;
  /*
   * Nombre de pixels pour lesquels l'algorithme utilise la climatologie.
   */
  int clim_val = 0;
  int day;
  // Day of year.
  int doy;
  char errmsg[] = "Usage:\n"
  "./prodprim\n"
 	"--rrs_type rrs_type\n"
  "--rrs_file rrs_file\n"
  "--atmosphere_type atmosphere_type\n"
  "[--atmosphere_file atmosphere_file]\n"
  "[--atmosphere_file00 atmosphere_file00\n"
  "--atmosphere_file03 atmosphere_file03\n"
  "--atmosphere_file06 atmosphere_file06\n"
  "--atmosphere_file09 atmosphere_file09\n"
  "--atmosphere_file12 atmosphere_file12\n"
  "--atmosphere_file15 atmosphere_file15\n"
  "--atmosphere_file18 atmosphere_file18\n"
  "--atmosphere_file21 atmosphere_file21\n"
  "--atmosphere_file00tomorrow atmosphere_file00tomorrow]\n"
  "--ice_file ice_file\n"
  "--doy doy\n"
  "[--first_bin first_bin]\n"
  "[--last_bin last_bin]\n"
  "--outfile outfile\n"
  "where\n"
  "rrs_type is the sensor for the Rrs, S for SeaWiFS and A for MODISA.\n"
  "rrs_file is the Rrs file.\n"
  "atmosphere_type is the dataset for the atmospheric products, I for ISCCP\n"
  " andM for MODISA.\n"
  "atmosphere_file is the file for the atmospheric products for one day.\n"
  " This argument is only valid with atmosphere_type = M.\n"
  "atmosphere_file00 is the file for the atmospheric products at 00 UTC.\n"
  " This argument is only valid with atmosphere_type = I.\n"
  "atmosphere_file03 is the file for the atmospheric products at 03 UTC.\n"
  " This argument is only valid with atmosphere_type = I.\n"
  "atmosphere_file06 is the file for the atmospheric products at 06 UTC.\n"
  " This argument is only valid with atmosphere_type = I.\n"
  "atmosphere_file09 is the file for the atmospheric products at 09 UTC.\n"
  " This argument is only valid with atmosphere_type = I.\n"
  "atmosphere_file12 is the file for the atmospheric products at 12 UTC.\n"
  " This argument is only valid with atmosphere_type = I.\n"
  "atmosphere_file15 is the file for the atmospheric products at 15 UTC.\n"
  " This argument is only valid with atmosphere_type = I.\n"
  "atmosphere_file18 is the file for the atmospheric products at 18 UTC.\n"
  " This argument is only valid with atmosphere_type = I.\n"
  "atmosphere_file21 is the file for the atmospheric products at 21 UTC.\n"
  " This argument is only valid with atmosphere_type = I.\n"
  "atmosphere_file00tomorrow is the file for the atmospheric products at\n"
  " 00 UTC tomorrow. This argument is only valid with atmosphere_type = I.\n"
  "ice_file is the sea-ice concentration file.\n"
  "doy is the day of year.\n"
  "outfile is the NetCDF outfile.\n"
  "first_bin is the first bin in the in the L3BIN grid to be treated\n"
  " (inclusive).\n"
  " Default is 0.\n"
  "last_bin is the last bin in the in the L3BIN grid to be treated\n"
  " (inclusive).\n"
  " Default is the last bin of the L3BIN grid.\n";
  // First bin in the L3BIN grid to be treated north of 45N.
  int first_bin = 0;
  /*
   * Vaut 1 si l'algorithme rejette le pixel.
   * Vaut 0 sinon.
   */
  int get_out = 0;
  /*
   * La concentration de glace lue dans le fichier de glace si la elle est
   * valide,
   * 50% si elle est invalide.
   */
  float glace;
  /*
   * La concentration de glace lue dans le fichier de glace.
   */
  float glace_lue;
  /*
   * Entier decrivant une des deux dimensions de cellule de grille (1=12.5 km,
   * 2=25.0 km).
   */
  int gtype = 2;
  /* Type of output. */
  hdf4_out_type hdf4_out;
  /* 
   * Row number on the L3BIN grid.					       
   * Corresponds to ROW in						       
   * http://oceancolor.gsfc.nasa.gov/SeaWiFS/TECH_REPORTS/PreLPDF/PreLVol32.pdf
   * appendix A.							       
   * From 1 to number of rows on the grid.                                     
   */
  short i;
  /*
   * Index of the band of the sensor.
   */
  int iband;
  /*
   * Index (0-based) of the bin in the image.	  
   * -1 if the bin is not in the image.           
   */
  int ibinimage;
  /*
   *  Nombre de pixels dans l'image qui sont au nord de 45N, qui ont des Rrs
   *  valides, qui ont une chlorophylle calculee par GSM valide, qui ont des
   *  IOPs calcules par QAA valides et qui ont une concentration de glace
   *  < 10%.
   */
  int ice_10 = 0;
  /*
   * The sea-ice concentration file.
   */
  char* ice_file = NULL;
  /*
   * Nombre de pixels pour lesquels la concentration de glace n'est pas valide.
   * (On compte tous les pixels de la grille du satellite et pas seulement
   * ceux de l'image l3BIN.)
   */
  int ice_not_val = 0;
  /*
   * Nombre de pixels pour lesquels la concentration de glace est valide.
   * (On compte tous les pixels de la grille du satellite et pas seulement
   * ceux de l'image l3BIN.)
   */
  int ice_val = 0;
  /*
   * The bin index number on the (full) L3BIN grid.			       
   * Corresponds to IDX in						       
   * http://oceancolor.gsfc.nasa.gov/SeaWiFS/TECH_REPORTS/PreLPDF/PreLVol32.pdf
   * appendix A.							       
   * From 1 to the total number of possible bins in the grid.                  
   */
  int idx;
  /*
   * The index of the file in array1d_ifile_atmosphere_file.
   */
  int ifile;
  /*
   * Entier decrivant une des deux regions polaires ( 1=Nord, 2=Sud).
   */
  int ihem = 1;
  /*
   * Index of the current pixel on the grid north of 45N.
   */
  int ipix;
  /*
   * Is 1 if and only if the atmospheric products files are passed in the
   * arguments.
   * Is 0 if not.
   */
  int is_atm = 0;
  /*
   * Les latitudes pour chaque pixel de notre sous-grille de la grille EQ du
   * ISCCP au nord de 45 degres Nord.
   * Les latitudes sont multiplies par un facteur de 100.
   * Les latitudes vont de 4500 a 9000.
   */
  int isccp_grid_lat[NB_PIXELS_ISCCP];
  /*
   * Les longitudes pour chaque pixel de notre sous-grille de la grille EQ du
   * ISCCP au nord de 45 degres Nord.
   * Les longitudes vont de 0 a 360.
   */
  float isccp_grid_lon[NB_PIXELS_ISCCP];
  /*
   * Vaut 1 si et seulement si la climatologie doit etre utilisee.
   * Vaut 0 sinon.
   */
  int is_clim = 0;
  /*
   * Vaut 1 si la concentration de glace est < 10%.
   * Vaut 0 sinon. Retourne 0 si la concentration de glace est >= 10% ou
   * une valeur speciale.
   */
  int is_ice_10 = 0;
  /*
   * Vaut 1 si et seulement si la concentration de glace est valide.
   * Vaut 0 sinon.
   */
  int is_ice_val = 0;
  /*
   * Vaut 1 si la concentration de chlorophylle-a et les IOPs sont valides.
   * Vaut 0 sinon.
   */
  int is_oc_val = 0;
  /*
   * Vaut 1 si les Rrs sont valides.
   * Vaut 0 sinon.
   */
  int is_Rrs_val = 0;
  /*
   * Index of the time at local noon.
   * Unitless.
   */
  int itime_local_noon;
  /* 
   * The index of the hour for the hours from 0 to 24 by step of 3 h.
   */
  int itime;
  /*
   * Index commencant a 1 de la colonne de la grille du NSIDC.
   */
  int iTNB;
  /*
   * Entier decrivant le type de transformation LOCATE fera
   * (1=I,J-vers-Lat,Lon; 2=Lat,Lon-vers-I,J).
   */
  int itrans = 2;
  /*
   * Index of the wavelength from 400 to 700 nm by steps of 5 nm.
   */
  int ivis;
  /*
   * Bin number on this row.
   * Corresponds to COL in
   * http://oceancolor.gsfc.nasa.gov/SeaWiFS/TECH_REPORTS/PreLPDF/PreLVol32.pdf
   * appendix A.
   * From 1 to the total number of possible bins in this row.
   */
  short j;
  /*
   * Index commencant a 1 de la rangee de la grille du NSIDC.
   */
  int jTNB;
  /*
   * Vertical attenuation of photosynthetically usable radiation at local
   * zenith (local noon).
   * Units: m^-1.
   */
  float KPUR;
  // Last bin in the L3BIN grid to be treated (inclusive) north of 45N.
  int last_bin = 0;
  /*
   * Nombre de pixels dans l'image qui sont au nord de 45N.
   */
  int l_val = 0;
  int month;
  /* Number of pixels in the grid north of 45N. */
  int nbpix = 0;
  /* Type of output. */
  netcdf_out_type netcdf_out;
  /*
   * ibinimage (0-based) of the first pixel of the next row of the grid.      
   * Number of (documented) bins in the image in the latitudes south or equal 
   * to the current latitude						      
   * (ibinimage_firstdocumentedbin_currentrow + extent[i]).                   
   */
  int ibinimage_firstdocumentedbin_nextrow;
  /*
   * Total number of possible bins in the current row.
   */
  int numbin;
  /*
   * Total number of possible bins in the previous row.
   */
  int numbinM1;
  /*
   * ibinimage (0-based) of the first pixel of the current row of the grid.
   * Number of (documented) bins in the image in the latitudes south
   * of the current latitude.
   */
  int ibinimage_firstdocumentedbin_currentrow;
  /*
   * Tableau de dimensions 9 x 966.
   * Ce tableau contient les valeurs des produits atmospheriques lus dans
   * les fichiers du ISCCP sur notre sous-grille de de la grille EQ du
   * ISCCP au nord de 45 degres Nord.
   * Il y a 966 valeurs pour chaque heure.
   * O3 est l'ozone.
   * CF est la fraction nuageuse.
   * TauCld est l'epaisseur optique.
   */
  float O3[NBHOURS][NB_PIXELS_ISCCP];
  /*
   * The mean ozone for a time interval of one day.
   * With ISCCP, it is the mean of the ozone at 00UTC, 03UTC, 06UTC,
   * 09UTC, 12UTC, 15UTC, 18UTC and 21UTC.
   * With MODIS-Atmosphere, it is directly the ozone for that day.
   * There is no mean to do in that case.
   * Unitless.
   */
  float O3mean;
  /*
   * Nombre de pixels dans l'image qui sont au nord de 45N, qui ont des Rrs
   * valides, qui ont une chlorophylle calculee par GSM valide et qui ont des
   * IOPs calcules par QAA valides.
   */
  int oc_val = 0;
  // A cloud fraction forced to 0.
  float oneCF_clear = 0.;
  float oneCF;
  /*
   * The cloud fraction for a time interval of one day. Used with
   * MODIS-Atmosphere.
   */
  float oneCFday;
  float oneO3;
  /*
   * The ozone for a time interval of one day. Used with MODIS-Atmosphere.
   */
  float oneO3day;
  float oneTauCld;
  /*
   * The optical thickness for a time interval of one day. Used with
   * MODIS-Atmosphere.
   */
  float oneTauCldday;
  /* Name of the output file without extension. */
  char *outfile;
  /* Name of the HDF4 output file */
  char outfile_hdf4[500];
  /*
   * The photosynthetically available irradiance (PAR) for a pixel on a given
   * day with the cloud fraction forced to 0.
   * The units are uEinstent*m^-1.
   */
  float PAR_clear;
  /*
   * The photosynthetically available irradiance (PAR) for a pixel on a given
   * day.
   * The units are uEinstent*m^-1.
   */
  float PAR_cloud;
  float phi;
  float photoperiode;
  float poc; // Particulate organic carbon.
  /*
   * The primary productivity computed with GSM.
   * The units are mgC*m^-2*d^-1.
   */
  float pp;
  /*
   * The primary productivity computed with GSM regardless of the ice.
   * The units are mgC*m^-2*d^-1.
   */
  float pp_temp;
  /*
   * The primary productivity computed with Cota et al. 2004.
   * The units are mgC*m^-2*d^-1.
   */
  float PP_cota;
  /*
   * The primary productivity computed with Ben Mustapha et al. 2012.
   * The units are mgC*m^-2*d^-1.
   */
  float PP_gsm_mustapha;
  /*
   * The primary productivity computed with OC3 or OC4.
   * The units are mgC*m^-2*d^-1.
   */
  float PP_oc;
  /*
   * Nombre de pixels dans l'image qui sont au nord de 45 N, qui ont des Rrs
   * valides, qui ont une concentration de chlorophylle et des IOPs valides
   * tout en ayant une concentration de glace < 10% et dont la production
   * primaire est valide.
   */
  int pp_val = 0;
  /*
   * Array of dimension NBANDS = 6.					       
   * The first dimension are the wavelengths of the bands of the satellite.    
   * The values are the bands of the satellite.
   * The units are nm.                                                       
   */
  float *ptr_array1d_iband_band;
  /*
   * Array of dimension nbpix.
   * The columns contain the absorption coefficient for CDOM at 412 nm (m^-1).
   */
  float *ptr_array1d_ipix_aCDOM_412;
  /*
   * Array of dimension nbpix.
   * The first dimension is the pixel.
   * The values are phytoplankton absorption coefficient at 443 nm (a_phy at
   * 443 nm).
   * The units are m^-1.
   */
  float *ptr_array1d_ipix_aphy443;
  /*
   * Array of dimension nbpix.
   * The first dimension is the pixel.
   * The values are the cloud fractions (0 to 1) at 00h UTC.
   * The values are unitless.
   */
  float *ptr_array1d_ipix_CF00;
  /*
   * Array of dimension nbpix.
   * The first dimension is the pixel.
   * The values are the chlorophyll-a concentration from GSM.
   * The units are (mg Chl-a m^-3).
   */
  float *ptr_array1d_ipix_chl_gsm;
  /*
   * Array of dimension nbpix.
   * The first dimension is the pixel.
   * The values are the chlorophyll-a concentration from Ben Mustapha et
   * al. 2012.
   * The units are (mg Chl-a m^-3).
   */
  float *ptr_array1d_ipix_chl_gsm_mustapha;
  /*
   * Array of dimension nbpix.
   * The first dimension is the pixel.
   * The values are the chlorophyll-a concentration from Cota et al. 2004.
   * The units are (mg Chl-a m^-3).
   */
  float *ptr_array1d_ipix_chl_cota;
  /*
   * Array of dimension nbpix.
   * The first dimension is the pixel.
   * The values are the chlorophyll-a concentration from OC3 or OC4.
   * The units are (mg Chl-a m^-3).
   */
  float *ptr_array1d_ipix_chl_oc;
  /*
   * Array of dimension nbpix.
   * The first dimension is the pixel.
   * The values are the mean cloud fractions (0 to 1).
   * The values are unitless.
   */
  float *ptr_array1d_ipix_CFmean;
  /*
   * Array of dimension nbpix.
   * The first dimension is the pixel.
   * The values are the sea-ice concentration.
   * Normal values from 0 to 1.
   * Special values > 1.
   * Unitless.
   */
  float *ptr_array1d_ipix_ice;
  /*
   * Array of dimension nbpix.
   * The first dimension is the pixel.
   * The values are the row number on the L3BIN grid.
   * Corresponds to ROW in
   * http://oceancolor.gsfc.nasa.gov/SeaWiFS/TECH_REPORTS/PreLPDF/PreLVol32.pdf
   * appendix A.
   * From 1 to number of rows on the grid.
   */
  short *ptr_array1d_ipix_indlat;
  /*
   * Array of dimension nbpix.
   * The first dimension is the pixel.
   * The values are the bin number on this row.
   * Corresponds to COL in
   * http://oceancolor.gsfc.nasa.gov/SeaWiFS/TECH_REPORTS/PreLPDF/PreLVol32.pdf
   * appendix A.
   * From 1 to the total number of possible bins in this row.
   */
  short *ptr_array1d_ipix_indlon;
  /*
   * Array of dimension nbpix.
   * The first dimension is the pixel.
   * The values are the vertical attenuation of photosynthetically usable
   * radiation.
   * The units are m^-1.
   */
  float *ptr_array1d_ipix_KPUR;
  /*
   * Array of dimension nbpix.
   * The first dimension is the pixel.
   * The values are the total ozone column at 00h UTC.
   * The units are Dobson units.
   */
  float *ptr_array1d_ipix_O300;
  /*
   * Array of dimension nbpix.
   * The first dimension is the pixel.
   * The values are the mean total ozone column.
   * The units are Dobson units.
   */
  float *ptr_array1d_ipix_O3mean;
  /*
   * Array of dimension nbpix.
   * The first dimension is the pixel.
   * The values are the latitudes from -90 to 90.
   * The units are degrees N.
   */
  float *ptr_array1d_ipix_lat;
  /*
   * Array of dimension nbpix.
   * The first dimension is the pixel.
   * The values are the latitudes from -180 to 180.
   * The units are degrees E.
   */
  float *ptr_array1d_ipix_lon;
  /*
   * Array of dimension nbpix.
   * The first dimension is the pixel.
   * The values are the photosynthetically available irradiances (PAR) for a
   * pixel on a given day with the cloud fraction forced to 0.
   * The units are uEinstent*m^-1.
   */
  float *ptr_array1d_ipix_PAR_clear;
  /*
   * Array of dimension nbpix.
   * The first dimension is the pixel.
   * The values are the photosynthetically available irradiances (PAR) for a
   * pixel on a given day.
   * The units are uEinstent*m^-1.
   */
  float *ptr_array1d_ipix_PAR_cloud;
  /*
   * Array of dimension nbpix.
   * The first dimension is the pixel.
   * The values are the particulate organic carbon.
   */
  float *ptr_array1d_ipix_poc;
  /*
   * Array of dimension nbpix.
   * The first dimension is the pixel.
   * The values are the primary productivity computed with GSM.
   * The units are mgC*m^-2*d^-1.
   */
  float *ptr_array1d_ipix_PP_chl_gsm;
  /*
   * Array of dimension nbpix.
   * The first dimension is the pixel.
   * The values are the primary productivity computed with Cota et al. 2004.
   * The units are mgC*m^-2*d^-1.
   */
  float *ptr_array1d_ipix_PP_cota;
  /*
   * Array of dimension nbpix.
   * The first dimension is the pixel.
   * The values are the primary productivity computed with Ben Mustapha et al.
   * 2012.
   * The units are mgC*m^-2*d^-1.
   */
  float *ptr_array1d_ipix_PP_gsm_mustapha;
  /*
   * Array of dimension nbpix.
   * The first dimension is the pixel.
   * The values are the primary productivity computed with OC3 or OC4.
   * The units are mgC*m^-2*d^-1.
   */
  float *ptr_array1d_ipix_PP_oc;
  /*
   * Array of dimension nbpix.
   * The first dimension is the pixel.
   * The values are the cloud optical thickness at 00h UTC.
   * The values are unitless.
   */
  float *ptr_array1d_ipix_tauCld00;
  /*
   * Array of dimension nbpix.
   * The first dimension is the pixel.
   * The values are the mean cloud optical thickness
   * The values are unitless.
   */
  float *ptr_array1d_ipix_tauCldmean;
  /*
   * Array of dimension nbpix.
   * The first dimension is the pixel.
   * The values are the solar zenith angles at 00h UTC.
   * The values are in degrees.
   */
  float *ptr_array1d_ipix_thetas00;
  /*
   * Tableau de floats de dimensions nbpix.
   * Les colonnes representent les profondeurs euphotiques moyennes (en m).
   */
  float *ptr_array1d_ipix_Zeu_bar;
  /*
   * Pointer to an array of dimension 
   * ntaucl=8 * nozone=8 * nthetas=19 * nwl=83.
   * The first dimension is the mean cloud optical thickness.
   * The second dimension is the total ozone column.
   * The third dimension is the solar zenith angle.
   * The fourth dimension is the wavelengths from 290 to 700 by step 5 nm.
   * The values are the downward irradiance just below the water surface
   * from the lookup table.
   * The units are umol photons*m^-2*s^-1*nm^-1. 
   */
  float* ptr_array4d_itaucl_iozone_ithetas_iwl_Ed0minus;
  /*
   * The type for the atmospheric products.
   * I for ISCCP and M for MODIS.
   */
  char* ptr_atmosphere_type;
  /*
   * Array of dimension nbinimage.
   * The first dimension is the index of the pixel in the image.
   * The value is the index of the pixel in the grid.
   * From 1 to nbinimage.
   * Unitless.
   */
  unsigned int* ptr_bin_num = NULL;
  /*
   * Day of year.
   */
  char* ptr_doy;
  // First bin in the L3BIN grid to be treated.
  char* ptr_first_bin = NULL;
  // Last bin in the L3BIN grid to be treated (inclusive).
  char* ptr_last_bin = NULL;
  /*
   * S for SeaWiFS and A for MODIS-AQUA.
   */
  char* ptr_rrs_type;
  /*
   * The Rrs file.
   */
  char* rrs_file;
  /*
   * S for SeaWiFS and A for MODIS-AQUA.
   */
  char rrs_type;
  /*
   * Nombre de pixels dans l'image qui sont au nord de 45N et qui ont des Rrs
   * valides.
   */
  int Rrs_val = 0;
  /* Nombre de lignes (latitudes) de 90S a 90N sur la grille du satellite. */
  int szlat;
  /* Index de la ligne de la grille SeaWiFS qui est a la latitude de 45N. */
  int szlat_45;
  /*
   * Tableau de dimensions 9 x 966.
   * Ce tableau contient les valeurs des produits atmospheriques lus dans
   * les fichiers du ISCCP sur notre sous-grille de de la grille EQ du
   * ISCCP au nord de 45 degres Nord.
   * Il y a 966 valeurs pour chaque heure.
   * O3 est l'ozone.
   * CF est la fraction nuageuse.
   * TauCld est l'epaisseur optique.
   */
  float TauCld[NBHOURS][NB_PIXELS_ISCCP];
  /*
   * The mean optical thickness for a time interval of one day.
   * With ISCCP, it is the mean of the optical thickness at 00UTC, 03UTC, 06UTC,
   * 09UTC, 12UTC, 15UTC, 18UTC and 21UTC.
   * With MODIS-Atmosphere, it is directly the optical thickness for that day.
   * There is no mean to do in that case.
   * Unitless.
   */
  float TauCldmean;
  float thetas;
  /*
   * Time at local noon.
   * The units are hours UTC.
   */
  float time_local_noon;
  float tmp;
  /*
   * doy or doy + 1
   */
  int tmpdoy;
  /*
   * Tableaux de dimensions 966.
   * Ces tableaux contiennent les valeurs des produits atmospheriques lus dans
   * les fichiers du ISCCP sur notre sous-grille de de la grille EQ du
   * ISCCP au nord de 45 degres Nord.
   * Ces tableaux ne contiennent les valeurs que pour une seule heure precise.
   */
  float tmpO3[NB_PIXELS_ISCCP],
  tmpCF[NB_PIXELS_ISCCP],
  tmpTauCld[NB_PIXELS_ISCCP];
  /*
   * Vaut 1 si l'algorithme utilise le pixel.
   * Vaut 0 sinon.
   */
  int use_l3bin = 0;
  //float vis_a[NVIS];
  float vis_aphy[NVIS];
  //float vis_bb[NVIS];
  //  float vis_Kd[NVIS];
  float xlat;
  float xlon;
  int year;
  /* La profondeur de la zone euphotique (en m). */
  float Zeu;
  /* La moyenne de la profondeur de la zone euphotique pour les differentes
   * heures (0, 3, 6, 9, 12, 15, 18, 21) pour un pixel.
   * L'unite est le m. */
  float Zeu_bar;
  
	/////////// Get the arguments. ///////////
  // See http://www.gnu.org/software/libc/manual/html_node/Getopt.html#Getopt.
  
  for(ifile = 0; ifile < NBHOURS; ifile++){
    array1d_ifile_atmosphere_file[ifile] = NULL;
  }
  
  while(1){
    static struct option long_options[] ={
      {"rrs_type", required_argument, 0, 's'},
      {"rrs_file", required_argument, 0, 'r'},
      {"atmosphere_type", required_argument, 0, 'b'},
      {"atmosphere_file", required_argument, 0, 'a'},
      {"atmosphere_file00", required_argument, 0, '0'},
      {"atmosphere_file03", required_argument, 0, '1'},
      {"atmosphere_file06", required_argument, 0, '2'},
      {"atmosphere_file09", required_argument, 0, '3'},
      {"atmosphere_file12", required_argument, 0, '4'},
      {"atmosphere_file15", required_argument, 0, '5'},
      {"atmosphere_file18", required_argument, 0, '6'},
      {"atmosphere_file21", required_argument, 0, '7'},
      {"atmosphere_file00tomorrow", required_argument, 0, '8'},
      {"ice_file", required_argument, 0, 'i'},
      {"doy", required_argument, 0, 'd'},
      {"first_bin", required_argument, 0, 'f'},
      {"last_bin", required_argument, 0, 'l'},
      {"outfile", required_argument, 0, 'o'},
      {0, 0, 0, 0}
    };
    /* getopt_long stores the option index here. */
    int option_index = 0;
    c = getopt_long (argc,
                     argv,
                     "a:d:f:i:o:r:s:",
                     long_options,
                     &option_index);
    /* Detect the end of the options. */
    if (c == -1){
      break;
    }
    switch (c){
      case 0:
        /* If this option set a flag, do nothing else now. */
        if(long_options[option_index].flag != 0){
          break;
        }
      case '0':
        printf ("option --atmosphere_file00 with value `%s'\n", optarg);
        array1d_ifile_atmosphere_file[0] = optarg;
        break;
      case '1':
        printf ("option --atmosphere_file03 with value `%s'\n", optarg);
        array1d_ifile_atmosphere_file[1] = optarg;
        break;
      case '2':
        printf ("option --atmosphere_file06 with value `%s'\n", optarg);
        array1d_ifile_atmosphere_file[2] = optarg;
        break;
      case '3':
        printf ("option --atmosphere_file09 with value `%s'\n", optarg);
        array1d_ifile_atmosphere_file[3] = optarg;
        break;
      case '4':
        printf ("option --atmosphere_file12 with value `%s'\n", optarg);
        array1d_ifile_atmosphere_file[4] = optarg;
        break;
      case '5':
        printf ("option --atmosphere_file15 with value `%s'\n", optarg);
        array1d_ifile_atmosphere_file[5] = optarg;
        break;
      case '6':
        printf ("option --atmosphere_file18 with value `%s'\n", optarg);
        array1d_ifile_atmosphere_file[6] = optarg;
        break;
      case '7':
        printf ("option --atmosphere_file21 with value `%s'\n", optarg);
        array1d_ifile_atmosphere_file[7] = optarg;
        break;
      case '8':
        printf ("option --atmosphere_file00tomorrow with value `%s'\n", optarg);
        array1d_ifile_atmosphere_file[8] = optarg;
        break;
      case 'a':
        printf ("option --atmosphere_file with value `%s'\n", optarg);
        atmosphere_file = optarg;
        break;
      case 'b':
        printf ("option --atmosphere_type with value `%s'\n", optarg);
        ptr_atmosphere_type = optarg;
        break;
      case 'd':
        printf ("option --doy with value `%s'\n", optarg);
        ptr_doy = optarg;
        break;
      case 'f':
        printf ("option --first_bin with value `%s'\n", optarg);
        ptr_first_bin = optarg;
        break;
      case 'i':
        printf ("option --ice_file with value `%s'\n", optarg);
        ice_file = optarg;
        break;
      case 'l':
        printf ("option --last_bin with value `%s'\n", optarg);
        ptr_last_bin = optarg;
        break;
      case 'o':
        printf ("option --outfile with value `%s'\n", optarg);
        outfile = optarg;
        break;
      case 'r':
        printf ("option -rrs_file with value `%s'\n", optarg);
        rrs_file = optarg;
        break;
      case 's':
        printf ("option --rrs_type with value `%s'\n", optarg);
        ptr_rrs_type = optarg;
        break;
      case '?':
        fprintf(stderr,
                errmsg,
                NULL);
        exit(-1);
        break;
      default:
        exit(-1);
    }
  }
  
  /////////// Test the arguments. ///////////
  
  rrs_type = get_args_sensor(ptr_rrs_type, errmsg);
  printf("rrs_type = %c\n", rrs_type);
  
  if(rrs_type == SEAWIFS){
    nbpix = NBPIX_SEAWIFS;
    szlat = SZLAT_SEAWIFS;
    szlat_45 = SZLAT_45_SEAWIFS;
    ptr_array1d_iband_band = ARRAY1D_IBAND_BANDSEAWIFS;
  }else if(rrs_type == MODISA){
    nbpix = NBPIX_MODISA;
    szlat = SZLAT_MODISA;
    szlat_45 = SZLAT_45_MODISA;
    ptr_array1d_iband_band = ARRAY1D_IBAND_BANDMODISA;
  }
  
  rrs_file = get_args_file(rrs_file, errmsg);
  printf("rrs_file = %s\n", rrs_file);
  
  atmosphere_type = get_args_atmosphere_type(ptr_atmosphere_type, errmsg);
  
  if(!
     (atmosphere_type == isccp
      && atmosphere_file == NULL
      && is_array1d_full(array1d_ifile_atmosphere_file, NBHOURS)
      )
     && !
     (atmosphere_type == modis_atmosphere
      && atmosphere_file != NULL
      && is_array1d_null(array1d_ifile_atmosphere_file, NBHOURS)
      )
     ){
    fprintf(stderr,
            "Either use\n"
            "--atmosphere_type M\n"
            "with\n"
            "--atmosphere_file\n"
            "but without\n"
            "--atmosphere_file00,\n"
            "--atmosphere_file03,\n"
            "--atmosphere_file06,\n"
            "--atmosphere_file09,\n"
            "--atmosphere_file12,\n"
            "--atmosphere_file15,\n"
            "--atmosphere_file18,\n"
            "--atmosphere_file21,\n"
            "--atmosphere_file00tomorrow,\n"
            "OR\n"
            "--atmosphere_type I\n"
            "with\n"
            "--atmosphere_file00,\n"
            "--atmosphere_file03,\n"
            "--atmosphere_file06,\n"
            "--atmosphere_file09,\n"
            "--atmosphere_file12,\n"
            "--atmosphere_file15,\n"
            "--atmosphere_file18,\n"
            "--atmosphere_file21,\n"
            "--atmosphere_file00tomorrow,\n"
            "but without\n"
            "--atmosphere_file.\n",
            NULL);
    exit(-1);
  }
  
  printf("atmosphere_type = %c\n", atmosphere_type);
  
  printf("atmosphere_file = %s\n", atmosphere_file);
  
  printf("atmosphere_file00 = %s\n", array1d_ifile_atmosphere_file[0]);
  printf("atmosphere_file03 = %s\n", array1d_ifile_atmosphere_file[1]);
  printf("atmosphere_file06 = %s\n", array1d_ifile_atmosphere_file[2]);
  printf("atmosphere_file09 = %s\n", array1d_ifile_atmosphere_file[3]);
  printf("atmosphere_file12 = %s\n", array1d_ifile_atmosphere_file[4]);
  printf("atmosphere_file15 = %s\n", array1d_ifile_atmosphere_file[5]);
  printf("atmosphere_file18 = %s\n", array1d_ifile_atmosphere_file[6]);
  printf("atmosphere_file21 = %s\n", array1d_ifile_atmosphere_file[7]);
  printf("atmosphere_file00tomorrow = %s\n", array1d_ifile_atmosphere_file[8]);
  
  ice_file = get_args_file(ice_file, errmsg);
  printf("ice_file = %s\n", ice_file);
  
  doy = (short)get_args_doy(ptr_doy, errmsg);
  printf("doy = %d\n", doy);
  
  first_bin = get_args_bin(ptr_first_bin,
                           0,
                           nbpix - 1,
                           0,
                           "first_bin");
  printf("first_bin = %d\n", first_bin);
  
  last_bin = get_args_bin(ptr_last_bin,
                          0,
                          nbpix - 1,
                          nbpix - 1,
                          "last_bin");
  printf("last_bin = %d\n", last_bin);
  
  
  
  outfile = get_args_outfile(outfile, errmsg);
  printf("outfile = %s\n", outfile);
  
  /////////// Read constant parameters. ///////////
  
  /* Type of output. */
  hdf4_out = HDF4_OUT;
  netcdf_out = NETCDF_OUT;
  
  /////////// Prepare Atm documented? ///////////
  
  if(atmosphere_type == isccp){
    is_atm = 1;
  }
  
  /////////// Compute wavelengths. ///////////
  
  tmp=400;
  for(ivis = 0; ivis < NVIS; ivis++){
    array1d_ivis_lambda[ivis] = tmp;
    tmp += 5;
  }
  
  /////////// Test the arguments. ///////////
  
  if(doy >= 365){
    printf("Cette version de valide_pp.c ne supporte pas les jours juliens 365 et 366.\n");
    exit(-1);
  }
  
  /////////// Read lookup table of Ed0-. ///////////
  
  // Memory allocation.
  // nwl=83, nthetas=19, nozone=8, ntaucl=8
  ptr_array4d_itaucl_iozone_ithetas_iwl_Ed0minus
    = malloc(sizeof(float)*83*19*8*8);
  read_ed0moins_lut_(LUT, ptr_array4d_itaucl_iozone_ithetas_iwl_Ed0minus);
  
  /////////// Read ice. ///////////
  
  read_SSMI_TNB(ice_file, array2d_jTNB_iTNB_ice);
  
  /////////// Read atmospheric products. ///////////
  
  if(atmosphere_type == isccp){
    //Initialize Latitudes / Longitudes ISCCP data on the EQ grid.
    calc_isccp_grid(isccp_grid_lat, isccp_grid_lon);
    /* Trace */
    if(DEBUG >= DEBUG_2){
      print_isccp_grid_eq(isccp_grid_lat, isccp_grid_lon);
      print_isccp_grid_eq_parabole(isccp_grid_lat, isccp_grid_lon);
    } /* End of trace */
    // Read the atmospheric products from the ISCCP files.
    for(itime = 0; itime < NBHOURS; itime++){
      lecture_un_fic_srf(array1d_ifile_atmosphere_file[itime],
			 tmpO3,
			 tmpCF,
			 tmpTauCld);
      memcpy(O3[itime],
	     tmpO3,
	     NB_PIXELS_ISCCP * sizeof(float));
      memcpy(CF[itime],
	     tmpCF,
	     NB_PIXELS_ISCCP * sizeof(float));
      memcpy(TauCld[itime],
	     tmpTauCld,
	     NB_PIXELS_ISCCP * sizeof(float));
    }
  }else if(atmosphere_type == modis_atmosphere){
    // Read the atmospheric products from the MODIS (atmospheric) files.
    read_L3BIN_MODIS_Atmosphere(atmosphere_file,
                                array2d_ilat_ilon_cfdm,
                                array2d_ilat_ilon_to3m,
                                array2d_ilat_ilon_cotcm);
  }
  
  /////////// Read Rrs(lambda_6). ///////////
  
  read_l3b(rrs_file,
           rrs_type,
           array1d_irow_start_num_MODISA,
           array1d_irow_start_num_SeaWiFS,
           array1d_irow_begin_MODISA,
           array1d_irow_begin_SeaWiFS,
           array1d_irow_extent_MODISA,
           array1d_irow_extent_SeaWiFS,
           array1d_irow_max_MODISA,
           array1d_irow_max_SeaWiFS,
           &ptr_bin_num,
           array1d_ptr_array1d_iband_ibinimage_Rrs);
  
  /////////// Memory allocatin for the output variables. ///////////

  ptr_array1d_ipix_indlat           = (short*)malloc(sizeof(short) * nbpix);
  ptr_array1d_ipix_indlon           = (short*)malloc(sizeof(short) * nbpix);
  ptr_array1d_ipix_lat              = (float*)malloc(sizeof(float) * nbpix);
  ptr_array1d_ipix_lon              = (float*)malloc(sizeof(float) * nbpix);
  ptr_array1d_ipix_chl_gsm          = (float*)malloc(sizeof(float) * nbpix);
  ptr_array1d_ipix_aphy443          = (float*)malloc(sizeof(float) * nbpix);
  ptr_array1d_ipix_ice              = (float*)malloc(sizeof(float) * nbpix);
  ptr_array1d_ipix_PP_chl_gsm       = (float*)malloc(sizeof(float) * nbpix);
  ptr_array1d_ipix_Zeu_bar          = (float*)malloc(sizeof(float) * nbpix);
  ptr_array1d_ipix_poc              = (float*)malloc(sizeof(float) * nbpix);
  ptr_array1d_ipix_aCDOM_412        = (float*)malloc(sizeof(float) * nbpix);
  ptr_array1d_ipix_PAR_cloud        = (float*)malloc(sizeof(float) * nbpix);
  ptr_array1d_ipix_PAR_clear        = (float*)malloc(sizeof(float) * nbpix);
  ptr_array1d_ipix_CF00             = (float*)malloc(sizeof(float) * nbpix);
  ptr_array1d_ipix_O300             = (float*)malloc(sizeof(float) * nbpix);
  ptr_array1d_ipix_tauCld00         = (float*)malloc(sizeof(float) * nbpix);
  ptr_array1d_ipix_thetas00         = (float*)malloc(sizeof(float) * nbpix);
  ptr_array1d_ipix_CFmean           = (float*)malloc(sizeof(float) * nbpix);
  ptr_array1d_ipix_O3mean           = (float*)malloc(sizeof(float) * nbpix);
  ptr_array1d_ipix_tauCldmean       = (float*)malloc(sizeof(float) * nbpix);
  ptr_array1d_ipix_KPUR             = (float*)malloc(sizeof(float) * nbpix);
  ptr_array1d_ipix_chl_cota         = (float*)malloc(sizeof(float) * nbpix);
  ptr_array1d_ipix_chl_gsm_mustapha = (float*)malloc(sizeof(float) * nbpix);
  ptr_array1d_ipix_chl_oc           = (float*)malloc(sizeof(float) * nbpix);
  ptr_array1d_ipix_PP_cota          = (float*)malloc(sizeof(float) * nbpix);
  ptr_array1d_ipix_PP_gsm_mustapha  = (float*)malloc(sizeof(float) * nbpix);
  ptr_array1d_ipix_PP_oc            = (float*)malloc(sizeof(float) * nbpix);
  for(iband = 0; iband < NBANDS; iband++){
    array1d_ptr_array1d_iband_ipix_Rrs[iband]
      = (float*)malloc(sizeof(float)*nbpix);
    array1d_ptr_array1d_iband_ipix_a[iband]
      = (float*)malloc(sizeof(float) * nbpix);
    array1d_ptr_array1d_iband_ipix_bb[iband]
      = (float*)malloc(sizeof(float) * nbpix);
    array1d_ptr_array1d_iband_ipix_Kd[iband]
    = (float*)malloc(sizeof(float) * nbpix);
  }
  
  /* Trace */
  if(DEBUG >= DEBUG_1){
    printf("*******************************************************\n");
    printf("CALCUL DE LA PRODUCTION PRIMAIRE\n");
    //    printf("Fichier lu %s\n", argv[4]);
    printf("Day of year: %d\n", doy);
  } /* End of trace */
  /* Trace */
  if(DEBUG >= DEBUG_2){
    print_glace_exhaustive(array2d_jTNB_iTNB_ice);
  } /* End of trace */
  
  /////////////////////////////////////////////////////////////////////////////
  //                     Computing of a L3BIN image                          //
  /////////////////////////////////////////////////////////////////////////////
  
  ipix=0;
  basebin = 1;
  ibinimage_firstdocumentedbin_currentrow = 0;
  /* Initial count of the calculation for latitude 45N. */
  for (i = 1; i <= szlat_45; i++){
    if(rrs_type == SEAWIFS){
      ibinimage_firstdocumentedbin_nextrow
	= ibinimage_firstdocumentedbin_currentrow
	+ array1d_irow_extent_SeaWiFS[i - 1];
    }else if(rrs_type == MODISA){
      ibinimage_firstdocumentedbin_nextrow
	= ibinimage_firstdocumentedbin_currentrow
	+ array1d_irow_extent_MODISA[i - 1];
    }
    ibinimage_firstdocumentedbin_currentrow
      = ibinimage_firstdocumentedbin_nextrow;
  } /* End of the initial count of the calculation for latitude 45N. */
  /* Loop on lines (ie: latitudes) of the L3BIN grid. */
  for(i = 1; i <= szlat; i++){
    if(rrs_type == SEAWIFS){
      ibinimage_firstdocumentedbin_nextrow
	= ibinimage_firstdocumentedbin_currentrow
	+ array1d_irow_extent_SeaWiFS[i - 1];
    }else if(rrs_type == MODISA){
      ibinimage_firstdocumentedbin_nextrow
	= ibinimage_firstdocumentedbin_currentrow
	+ array1d_irow_extent_MODISA[i - 1];
    }
    xlat = ( (i - .5) * 180. / szlat) - 90.;
    numbin = (int) (2 * szlat * cos( (xlat) * PI / 180.) + .5);
    if(i > 1) basebin = basebinM1 + numbinM1;
    numbinM1  = numbin;
    basebinM1 = basebin;
    /* Is the pixel north of 45N? */
    if (xlat > 45.){
      /* Loop on the columns (ie: longitudes) of the L3BIN grid. */
      for(j=1;j<=numbin;j++){
	/* Initialize variables with foo values. */
        get_out                    = 0;
        use_l3bin                  = 0;
        is_ice_val                 = 0;
        are_pixels_on_current_line = 0;
        is_Rrs_val                 = 0;
        is_oc_val                  = 0;
        is_ice_10                  = 0;
	itime_local_noon           = -1;
	time_local_noon            = -999.;
        chl_gsm             = -999.;
        aphy443          = -999.;
        pp               = -999.;
        Zeu_bar          = -999.;
        poc              = -999.;
        aCDOM_412        = -999.;
        PAR_cloud        = -999.;
        PAR_clear        = -999.;
        *(ptr_array1d_ipix_CF00 + ipix)     = -999.;
        *(ptr_array1d_ipix_O300 + ipix)     = -999.;
        *(ptr_array1d_ipix_tauCld00 + ipix) = -999.;
        *(ptr_array1d_ipix_thetas00 + ipix) = -999.;
        CFmean           = -999.;
        O3mean           = -999.;
        TauCldmean       = -999.;
        KPUR             = -999.;
        chl_cota         = -999.;
        chl_gsm_mustapha = -999.;
        chl_oc           = -999.;
        PP_cota          = -999.;
        PP_gsm_mustapha  = -999.;
        PP_oc            = -999.;
        for(iband = 0; iband < NBANDS; iband++){
          array1d_iband_a[iband]      = -999.;
          array1d_iband_bb[iband]     = -999.;
          array1d_iband_Kd[iband]     = -999.;
          array1d_iband_Rrs[iband] = -999.;
          *(array1d_ptr_array1d_iband_ipix_Kd[iband] + ipix) = -999.;
	  for(itime = 0; itime < NBHOURS; itime++){
	    array2d_iband_itime_Kd[iband][itime] = -999.;
	  }
        }
        /* Is ipix in the interval selected from the arguments? */
        if(first_bin <= ipix && ipix <= last_bin){
          
          xlon = (360.*(j-.5)/numbin) - 180.;
          idx=basebin+j-1;

	  /////////// Local noon. ///////////
	  time_local_noon = get_timeGMT_local_noon(xlon);
	  itime_local_noon = get_itime(time_local_noon);
          
          /////////// Matchup ice. ///////////
          locate_(&gtype, &ihem, &itrans, &iTNB, &jTNB, &xlat, &xlon);
          /* Le pixel est en-dehors de la grille du NSIDC. */
          if(iTNB < 1 || iTNB > NBCOL || jTNB < 1 || jTNB > NBLINE){
            glace_lue = 1.02;
            /* Le pixel est dans la grille du NSIDC. */
          }else{
            glace_lue = array2d_jTNB_iTNB_ice[jTNB - 1][iTNB - 1];
          }
          /////////// Matchup ice. ///////////
          
          /// Matchup of the daily means of the atmospheric products. ///
          if(atmosphere_type == isccp){
            get_isccp_mean(CF,
                           O3,
                           TauCld,
                           xlat,
                           xlon,
                           isccp_grid_lat,
                           isccp_grid_lon,
                           &CFmean,
                           &O3mean,
                           &TauCldmean);
          }else if(atmosphere_type == modis_atmosphere){
            CFmean = get_modis_atm(array2d_ilat_ilon_cfdm, xlat, xlon);
            O3mean = get_modis_atm(array2d_ilat_ilon_to3m, xlat, xlon);
            TauCldmean = get_modis_atm(array2d_ilat_ilon_cotcm, xlat, xlon);
          }
          /// End of Matchup of the daily means of the atmospheric products. ///
          
          /////////// In L3BIN? ///////////
          if(!get_out){
            are_pixels_on_current_line
            = get_are_pixels_on_current_line(rrs_type,
                                             i,
                                             array1d_irow_begin_MODISA,
                                             array1d_irow_begin_SeaWiFS);
            if(!are_pixels_on_current_line){
              get_out = 1;
            }
          }
          if(!get_out){
            ibinimage = get_ibinimage(idx,
				      ibinimage_firstdocumentedbin_currentrow,
				      ibinimage_firstdocumentedbin_nextrow,
				      ptr_bin_num);
            if(ibinimage >= 0){
              l_val++;
            }else{
              get_out = 1;
            }
          }
          /////////// End of In L3BIN? ///////////
          
          /////////// Matchup of Rrs(lambda_6). ///////////
          if(!get_out){
            for(iband = 0; iband < NBANDS; iband++){
              array1d_iband_Rrs[iband]
		= *(array1d_ptr_array1d_iband_ibinimage_Rrs[iband] + ibinimage);
            }
          }
          /////////// Matchup of Rrs(lambda_6). ///////////
          
          /////////// Computing of chl_oc. ///////////
          if(!get_out){
            chl_oc = get_chl_oc(array1d_iband_Rrs, rrs_type);
          }
          /////////// End of the computing of chl_oc. ///////////
          
          /////////// Computing of chl_cota. ///////////
          if(!get_out){
            chl_cota = get_chl_cota(array1d_iband_Rrs, rrs_type);
          }
          /////////// End of the computing of chl_cota. ///////////
          
          /////////// Rrs documented? ///////////
          if(!get_out){
            is_Rrs_val
            = get_is_Rrs_val_L3BIN(ibinimage,
                                   array1d_ptr_array1d_iband_ibinimage_Rrs);
            if(is_Rrs_val){
              // Rrs valides.
              Rrs_val++;
            }else{
              get_out = 1;
            }
          }
          /////////// End of Rrs documented? ///////////
          
          /////////// Computing of poc. ///////////
          // Used only to write the poc product in the output file.
          // Not implemented for SeaWiFS.
          if(rrs_type == MODISA && !get_out){
            poc = calc_poc(array1d_iband_Rrs[POS_488],
			   array1d_iband_Rrs[POS_555]);
          }
          /////////// End of Computing of poc. ///////////
          
          /////////// Computing of chl_gsm_mustapha. ///////////
          if(!get_out){
            chl_gsm_mustapha = get_chl_gsm(array1d_iband_Rrs,
                                           rrs_type,
                                           ADG_S_CHL_GSM_MUSTAPHA,
                                           BBP_S_CHL_GSM_MUSTAPHA);
          }
          /////////// End of the computing of chl_mustapha. ///////////
          
          /////////// [chl] and all IOPs documented? ///////////
          if(!get_out){
            chl_gsm = get_chl_gsm(array1d_iband_Rrs,
                               rrs_type,
                               ADG_S_CHL_GSM,
                               BBP_S_CHL_GSM);
            if(rrs_type == SEAWIFS){
	      /* QAA pour a et bb */
              qaaSW(array1d_iband_Rrs, array1d_iband_a, array1d_iband_bb);
            }else if(rrs_type == MODISA){
	      /* QAA pour a et bb */
              qaaMA(array1d_iband_Rrs, array1d_iband_a, array1d_iband_bb);
            }
            is_oc_val = is_valide(chl_gsm, array1d_iband_a, array1d_iband_bb);
            if(is_oc_val){
              // Chlorophylle et IOPs valides.
              use_l3bin = 1;
              oc_val++;
            }else{
              get_out = 1;
            }
          }
          /////////// End of [chl] and all IOPs documented? ///////////

	  /////////// Compute Kd(lambda_6, t) ///////////
	  get_array2d_iband_itime_Kd(doy,
				     xlat,
				     xlon,
				     array1d_iband_a,
				     array1d_iband_bb,
				     array2d_iband_itime_Kd);
	  get_array1df_i_from_array2df_i_j(NBANDS,
					   NBHOURS,
					   array2d_iband_itime_Kd,
					   itime_local_noon,
					   array1d_iband_Kd);
	  for(iband = 0; iband < NBANDS; iband++){
	    *(array1d_ptr_array1d_iband_ipix_Kd[iband] + ipix)
	      = array1d_iband_Kd[iband];
	  }
	  /////////// End of compute Kd(lambda_6, t) ///////////
          
          /////////// Computing of aCDOM_412. ///////////
          // Used only to write the aCDOM_412 product in the output file.
          if(is_oc_val){
            aCDOM_412 = calc_aCDOM_412(array1d_iband_a[POS_412],
                                       array1d_iband_Rrs[POS_412],
                                       array1d_iband_Rrs[POS_490],
                                       array1d_iband_Rrs[POS_555]);
          }
          /////////// End of Computing of aCDOM_412. ///////////
          
          /////////// [ice] documented? ///////////
          glace = glace_lue;
          is_ice_val = get_is_ice_val(glace_lue);
          if(is_ice_val){
            // Ice valid.
            ice_val++;
          }else{
            // Ice not valid.
            ice_not_val++;
            get_out = 1;
          }
          /////////// End of [ice] documented? ///////////
          
          /////////// [ice] < 10%. ///////////
          if(!get_out){
            is_ice_10 = get_is_ice_10(glace);
            if(is_ice_10){
              // Glace < 10%
              ice_10++;
            }else{
              // glace >= 10% ou valeur speciale.
              get_out = 1;
            }
          }
          /////////// End of [ice] < 10%. ///////////
          
          /////////// Matchup and Atm documented? ///////////
          if(atmosphere_type == modis_atmosphere){
            // Reading of the atmospheric products from the MODIS-Atmosphere
            // file.
            oneCFday = get_modis_atm(array2d_ilat_ilon_cfdm, xlat, xlon);
            oneO3day = get_modis_atm(array2d_ilat_ilon_to3m, xlat, xlon);
            oneTauCldday = get_modis_atm(array2d_ilat_ilon_cotcm, xlat, xlon);
            is_atm = get_is_modis_atm_val(oneCFday, oneO3day, oneTauCldday);
            
            // Atmospheric product not valid.
            if(!is_atm){
              get_out = 1;
            }
          }
          /////////// End of Matchup and Atm documented? ///////////
          
          ///// Compute daily P rate for open water (P_temp). /////
          /* Will primary productivity be computed for this pixel? */
          if(!get_out && (use_l3bin)){
            calc_aphy443(chl_gsm, &aphy443);
            ///// Compute aphy(lambda_61) from chl_gsm  /////
            // From chl_gsm (Matsuoka et al. 2007).
            calc_aphy(aphy443, vis_aphy);
            ///// End of Compute aphy. /////
            /* Trace. */
            if(DEBUG >= DEBUG_1 && ipix == IPIX){
              printf("ipix: %d\n", ipix);
              printf(" Latitude (degres): %f\n", xlat);
              printf(" Longitude (degres): %f\n", xlon);
              printf(" chl_gsm (mgChla*m^-3): %f\n", chl_gsm);
              if(atmosphere_type == modis_atmosphere){
                printf(" CF (unitless): %0.f\n", oneCFday);
                printf(" O3 (DU): %0.f\n", oneO3day);
                printf(" TauCld (unitless): %0.f\n", oneTauCldday);
              }
              print_aphy443(aphy443);
              //print_iops(vis_aphy, vis_a, vis_bb);
              print_aphy(vis_aphy,
                         array1d_ivis_lambda);
            } /* End of trace. */
            /* Loop on time steps of 3 hours. */
            for(itime=0; itime<NBHOURS; itime++){
              ///// Compute time. /////
              if(itime == ITIME_TOMORROW_MIDNIGHT){
                tmpdoy = doy + 1;
              }else{
                tmpdoy = doy;
              }
              ///// End of Compute time. /////
              /* Trace. */
              if(DEBUG >= DEBUG_1 && ipix == IPIX){
                printf(" HEURE UTC %.0f\n",
		       ARRAY1D_ITIME_HOUR[itime]);
              } /* End of trace. */
	      /////////// Get Kd(ivis) ///////////
	      get_array1df_i_from_array2df_i_j(NBANDS,
					       NBHOURS,
					       array2d_iband_itime_Kd,
					       itime,
					       array1d_iband_Kd);
              interp_Kd(array1d_iband_Kd,
			array1d_ivis_Kd,
			ptr_array1d_iband_band,
			rrs_type);

	      /////////// End of get Kd(ivis) ///////////
	      
              ////////// Matchup of the atmospheric products. //////////
              if(atmosphere_type == isccp){
                calc_isccp_eq(O3,
                              CF,
                              TauCld,
                              itime,
                              xlat,
                              xlon,
                              isccp_grid_lat,
                              isccp_grid_lon,
                              &oneO3,
                              &oneCF,
                              &oneTauCld);
                /* Trace. */
                if(DEBUG >= DEBUG_1 && ipix == IPIX){
                  print_isccp_interpolation(ARRAY1D_ITIME_HOUR[itime],
                                            oneO3,
                                            oneTauCld,
                                            oneCF);
                } /* End of trace. */
                /* Trace. */
                if(DEBUG >= DEBUG_2 && ipix == IPIX){
                  print_isccp_products_eq_parabole(O3,
						   CF,
						   TauCld,
						   ARRAY1D_ITIME_HOUR,
						   itime);
                } /* End of trace. */
              }else if(atmosphere_type == modis_atmosphere){
                oneCF = oneCFday;
                oneO3 = oneO3day;
                oneTauCld = oneTauCldday;
              }
              ////////// End of Matchup of the atmospheric products. //////////
              ////////// Compute thetas. //////////
	      /* recuperation pour le calcul de Kd */
	      sunpos_(&tmpdoy,
		      &ARRAY1D_ITIME_HOUR[itime],
		      &xlat,
		      &xlon,
		      &thetas,
		      &phi);
	      /* Trace. */
	      if(DEBUG >= DEBUG_1 && ipix == IPIX){
		print_thetas(tmpdoy,
			     ARRAY1D_ITIME_HOUR[itime],
			     xlat,
			     xlon,
			     thetas);
	      } /* End of trace. */
              ////////// End of compute thetas. //////////
              ////////// Compute Ed0-(lambda_61). //////////
              ed0moins_(&tmpdoy,
                        &ARRAY1D_ITIME_HOUR[itime],
                        &xlat,
                        &xlon,
                        &oneO3,
                        &oneTauCld,
                        &oneCF,
                        array1d_iwl_Ed0minus,
                        ptr_array4d_itaucl_iozone_ithetas_iwl_Ed0minus,
                        &thetas);
              ////////// End of Compute Ed0-(lambda_61). //////////
              ////////// Compute Z(z, t). //////////
              calc_Z(array2d_idepth_itime_Z,
		     itime,
		     array1d_ivis_Kd[I550_FROM400]);
              ////////// End of Compute Z(z, t). //////////
              ////////// Compute E0(lambda_61, z). //////////
              /* equation 7 de l'ATBD */
              calc_E0_pixel_z_t(array1d_ivis_Kd,
                                array1d_iband_a,
                                array1d_iband_bb,
                                array2d_idepth_itime_Z,
                                array1d_iwl_Ed0minus,
                                itime,
                                array2d_ivis_idepth_E0);
              ////////// End of Compute E0(lambda_61, z). //////////
              ////////// Compute PUR(z, t). //////////
              /*equation 6 */
              calc_PUR(array2d_idepth_itime_PUR,
                       array2d_ivis_idepth_E0,
                       itime,
                       vis_aphy,
                       aphy443,
                       array1d_ivis_lambda);
              ////////// End of Compute PUR(z, t). //////////
              /* Store */
                *(ptr_array1d_ipix_CF00 + ipix)     = oneCF;
                *(ptr_array1d_ipix_O300 + ipix)     = oneO3;
                *(ptr_array1d_ipix_tauCld00 + ipix) = oneTauCld;
                *(ptr_array1d_ipix_thetas00 + ipix) = thetas;
              /* End of Store */
              /* Trace */
              if(DEBUG >= DEBUG_1 && ipix == IPIX){
                print_kd(array1d_iband_a,
			 array1d_iband_bb,
			 thetas,
			 array1d_iband_Kd,
			 ptr_array1d_iband_band);
                print_vis_Kd(array1d_ivis_Kd,
                             array1d_ivis_lambda);
                print_Z(array2d_idepth_itime_Z, itime);
                print_eclairement(array1d_iwl_Ed0minus,
				  array2d_ivis_idepth_E0);
                print_pur(itime, array2d_idepth_itime_PUR);
              } /* End of trace */
              ////////// Preparation for Compute Zeu_bar. //////////
              Zeu = array2d_idepth_itime_Z[POS_ZEU][itime];
              array1d_itime_Zeu[itime] = Zeu;
              ////////// End of Preparation for Compute Zeu_bar. //////////
              //// Preparation for Compute PAR_cloud and PAR_clear. ////
              calc_array2d_h_ivis_Ed_pixel(array1d_iwl_Ed0minus,
                                           array2d_h_ivis_Ed_pixel_cloud,
                                           itime);
              ed0moins_(&tmpdoy,
                        &ARRAY1D_ITIME_HOUR[itime],
                        &xlat,
                        &xlon,
                        &oneO3,
                        &oneTauCld,
                        &oneCF_clear,
                        array1d_iwl_Ed0minus_clear,
                        ptr_array4d_itaucl_iozone_ithetas_iwl_Ed0minus,
                        &thetas);
              calc_array2d_h_ivis_Ed_pixel(array1d_iwl_Ed0minus_clear,
                                           array2d_h_ivis_Ed_pixel_clear,
                                           itime);
               //// End of Preparation for Compute PAR_cloud and PAR_clear. ////
            } /* End of the loop on time steps of 3 hours. */
            ///// Compute day of month /////
            get_day(ice_file, &year, &month, &day);
            ///// End of Compute day of month /////
            ///// Compute daylength /////
            photoperiode = daylength(year, month, day, xlat);
            ///// End of Compute daylength /////
            ///// Compute meanPUR(z) /////
            calc_meanPUR(array2d_idepth_itime_PUR,
			 array1d_idepth_meanPUR,
			 photoperiode); /* equation 10 */
            ///// End of Compute meanPUR(z) /////
            ///// Compute Ek(z) /////
            calc_Ek(array1d_idepth_Ek, array1d_idepth_meanPUR); /* equation 8 */
            ///// End of Compute Ek(z) /////
            ///// Compute P_temp /////
            /* equation 5 */
            pp              = calc_PP(chl_gsm,
				      array2d_idepth_itime_PUR,
				      array1d_idepth_Ek,
				      array2d_idepth_itime_Z);
            PP_cota         = calc_PP(chl_cota,
				      array2d_idepth_itime_PUR,
				      array1d_idepth_Ek,
				      array2d_idepth_itime_Z);
            PP_gsm_mustapha = calc_PP(chl_gsm_mustapha,
				      array2d_idepth_itime_PUR,
				      array1d_idepth_Ek,
				      array2d_idepth_itime_Z);
            PP_oc           = calc_PP(chl_oc,
				      array2d_idepth_itime_PUR,
				      array1d_idepth_Ek,
				      array2d_idepth_itime_Z);
            ///// End of Compute P_temp /////
            ///// End of Compute daily P rate for open water (P_temp). /////
            
            /////////// P <- P_temp * (1 - [ice]) ///////////
            pp_temp         = pp;
            pp              = pp * (1. - glace);
            PP_cota         = PP_cota * (1. - glace);
            PP_gsm_mustapha = PP_gsm_mustapha * (1. - glace);
            PP_oc           = PP_oc * (1. - glace);
            /////////// End of P <- P_temp * (1 - [ice]) ///////////
            /* Is primary productivity valid? */
            if(!isnan(pp) && !isinf(pp)){
              pp_val++;
            }else{
              pp = -999.;
            }
            if(isnan(PP_cota) || isinf(PP_cota)){
              PP_cota = -999.;
            }
            if(isnan(PP_gsm_mustapha) || isinf(PP_gsm_mustapha)){
              PP_gsm_mustapha = -999.;
            }
            if(isnan(PP_oc) || isinf(PP_oc)){
              PP_oc = -999.;
            }
            /* End of if (Is primary productivity valid?). */
            ////////// Computing of Zeu_bar. //////////
            /*
             * It does not consider the values at midnight the next day for not 
             * giving more weight to midnight values.
             */
            Zeu_bar = mean(array1d_itime_Zeu, NBHOURS - 1);
            ////////// End of the computing of Zeu_bar. //////////
            ////////// Computing of PAR_cloud and PAR_clear. //////////
            PAR_cloud = calc_PAR(array2d_h_ivis_Ed_pixel_cloud);
            PAR_clear = calc_PAR(array2d_h_ivis_Ed_pixel_clear);
            //////// End of the computing of PAR_cloud and PAR_clear. ////////
            ////////// Computing of KPUR at local noon. //////////
            get_array1d_itime_KPUR_from_PUR_and_Z(array2d_idepth_itime_PUR,
						  array2d_idepth_itime_Z,
						  array1d_itime_KPUR);
            itime = get_itime_at_local_zenith_from_lon(xlon);
            KPUR = array1d_itime_KPUR[itime];
            ////////// End of Computing of KPUR at local noon. //////////
            /* Trace */
            if(DEBUG >= DEBUG_1 && ipix == IPIX){
              //printf("cloud\n");
              //print_array2d_h_ivis_Ed_pixel(array2d_h_ivis_Ed_pixel_cloud);
              //printf("PAR_cloud: %.3f\n", PAR_cloud);
              print_KPUR(array1d_itime_KPUR);
              print_photoperiode(photoperiode);
              print_meanpur_ek(array1d_idepth_meanPUR, array1d_idepth_Ek);
              print_pp_temp(pp_temp);
              print_glace(gtype, ihem, itrans, iTNB, jTNB, xlat, xlon, glace);
              print_pp(pp);
            } /* End of trace */
          } /* 
             * End of if (Will primary productivity be computed for this
             * pixel?).
             */
        } /* End of if (Is ipix in the interval selected from the arguments?).*/
        
        /* Store */
        *(ptr_array1d_ipix_indlat+ipix) = i;
        *(ptr_array1d_ipix_indlon+ipix) = j;
        *(ptr_array1d_ipix_lon+ipix) = xlon;
        *(ptr_array1d_ipix_lat+ipix) = xlat;
        *(ptr_array1d_ipix_chl_gsm+ipix) = chl_gsm;
        *(ptr_array1d_ipix_aphy443 + ipix) = aphy443;
        *(ptr_array1d_ipix_ice + ipix) = glace_lue;
        *(ptr_array1d_ipix_PP_chl_gsm + ipix) = pp;
        *(ptr_array1d_ipix_Zeu_bar + ipix) = Zeu_bar;
        *(ptr_array1d_ipix_poc + ipix) = poc;
        *(ptr_array1d_ipix_aCDOM_412 + ipix) = aCDOM_412;
        *(ptr_array1d_ipix_PAR_cloud + ipix) = PAR_cloud;
        *(ptr_array1d_ipix_PAR_clear + ipix) = PAR_clear;
        *(ptr_array1d_ipix_CFmean + ipix) = CFmean;
        *(ptr_array1d_ipix_O3mean + ipix) = O3mean;
        *(ptr_array1d_ipix_tauCldmean + ipix) = TauCldmean;
        *(ptr_array1d_ipix_KPUR + ipix) = KPUR;
        *(ptr_array1d_ipix_chl_cota + ipix) = chl_cota;
        *(ptr_array1d_ipix_chl_gsm_mustapha + ipix) = chl_gsm_mustapha;
        *(ptr_array1d_ipix_chl_oc + ipix) = chl_oc;
        *(ptr_array1d_ipix_PP_cota + ipix) = PP_cota;
        *(ptr_array1d_ipix_PP_gsm_mustapha + ipix) = PP_gsm_mustapha;
        *(ptr_array1d_ipix_PP_oc + ipix) = PP_oc;
        for(iband=0; iband<NBANDS; iband++){
          *(array1d_ptr_array1d_iband_ipix_a[iband]+ipix)
	    = array1d_iband_a[iband];
          *(array1d_ptr_array1d_iband_ipix_bb[iband]+ipix)
	    = array1d_iband_bb[iband];
          *(array1d_ptr_array1d_iband_ipix_Rrs[iband]+ipix)
	    = array1d_iband_Rrs[iband];
        }
        /* End of Store */
        ipix++; // Increment ipix.
      } /* End of the loop on the columns (ie: longitudes) of the L3BIN grid. */
      ibinimage_firstdocumentedbin_currentrow
	= ibinimage_firstdocumentedbin_nextrow;
    } /* End of if (Is the pixel north of 45N?). */
  } /* End of the loop on lines (ie: latitudes) of the L3BIN grid. */
  /* Trace */
  print_final(is_clim,
              ice_val,
              ice_not_val,
              l_val,
              Rrs_val,
              oc_val,
              ice_10,
              clim_val,
              pp_val);
  /* End of trace */
  
  /////////////////////////////////////////////////////////////////////////////
  //                  End of Computing of a L3BIN image                      //
  /////////////////////////////////////////////////////////////////////////////
  
  /////////// Write. ///////////
  
  if(hdf4_out == hdf4_out_true){
    // Remove ".nc".
    strncat(outfile_hdf4,
            outfile,
            strlen(outfile) - 3);
    strcat(outfile_hdf4, ".hdf");
    write_hdf(outfile_hdf4,
              ipix,
              ptr_array1d_iband_band,
              ptr_array1d_ipix_lat,
              ptr_array1d_ipix_lon,
              ptr_array1d_ipix_indlat,
              ptr_array1d_ipix_indlon,
              array1d_ptr_array1d_iband_ipix_Rrs,
              array1d_ptr_array1d_iband_ipix_a,
              array1d_ptr_array1d_iband_ipix_bb,
              array1d_ptr_array1d_iband_ipix_Kd,
              ptr_array1d_ipix_Zeu_bar,
              ptr_array1d_ipix_poc,
              ptr_array1d_ipix_aCDOM_412,
              ptr_array1d_ipix_PAR_cloud,
              ptr_array1d_ipix_PAR_clear,
              ptr_array1d_ipix_CF00,
              ptr_array1d_ipix_O300,
              ptr_array1d_ipix_tauCld00,
              ptr_array1d_ipix_thetas00,
              ptr_array1d_ipix_CFmean,
              ptr_array1d_ipix_O3mean,
              ptr_array1d_ipix_tauCldmean,
              ptr_array1d_ipix_KPUR,
              ptr_array1d_ipix_chl_cota,
              ptr_array1d_ipix_chl_oc,
              ptr_array1d_ipix_chl_gsm_mustapha,
              ptr_array1d_ipix_PP_cota,
              ptr_array1d_ipix_PP_gsm_mustapha,
              ptr_array1d_ipix_PP_oc,
              ptr_array1d_ipix_aphy443,
              ptr_array1d_ipix_chl_gsm,
              ptr_array1d_ipix_PP_chl_gsm,
              ptr_array1d_ipix_ice
              );
  }
  if(netcdf_out == netcdf_out_true){
    write_netcdf(outfile,
                 ipix,
                 ptr_array1d_iband_band,
                 ptr_array1d_ipix_lat,
                 ptr_array1d_ipix_lon,
                 ptr_array1d_ipix_indlat,
                 ptr_array1d_ipix_indlon,
                 array1d_ptr_array1d_iband_ipix_Rrs,
                 array1d_ptr_array1d_iband_ipix_a,
                 array1d_ptr_array1d_iband_ipix_bb,
                 array1d_ptr_array1d_iband_ipix_Kd,
                 ptr_array1d_ipix_Zeu_bar,
                 ptr_array1d_ipix_poc,
                 ptr_array1d_ipix_aCDOM_412,
                 ptr_array1d_ipix_PAR_cloud,
                 ptr_array1d_ipix_PAR_clear,
                 ptr_array1d_ipix_CF00,
                 ptr_array1d_ipix_O300,
                 ptr_array1d_ipix_tauCld00,
                 ptr_array1d_ipix_thetas00,
                 ptr_array1d_ipix_CFmean,
                 ptr_array1d_ipix_O3mean,
                 ptr_array1d_ipix_tauCldmean,
                 ptr_array1d_ipix_KPUR,
                 ptr_array1d_ipix_chl_cota,
                 ptr_array1d_ipix_chl_oc,
                 ptr_array1d_ipix_chl_gsm_mustapha,
                 ptr_array1d_ipix_PP_cota,
                 ptr_array1d_ipix_PP_gsm_mustapha,
                 ptr_array1d_ipix_PP_oc,
                 ptr_array1d_ipix_aphy443,
                 ptr_array1d_ipix_chl_gsm,
                 ptr_array1d_ipix_PP_chl_gsm,
                 ptr_array1d_ipix_ice
                 );
  }
  
  /////////// End. ///////////
  
  exit(0);
}
 
