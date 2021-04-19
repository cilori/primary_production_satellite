#include <Rcpp.h>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
using namespace std;
using namespace Rcpp;

// [[Rcpp::export]]
double GetAverage(NumericVector numbers) {
    double sum;
    for (int n=0; n < numbers.size(); n++) {
        sum = sum + numbers[n];
    }
    return(sum / numbers.size());
}

// Source: http://www.cplusplus.com/forum/general/216928/
// Returns interpolated value at x from parallel arrays ( xData, yData )
// Assumes that xData has at least two elements, is sorted and is strictly monotonic increasing
// boolean argument extrapolate determines behaviour beyond ends of array (if needed)
// [[Rcpp::export]]
double interpolate( vector<double> &xData, vector<double> &yData, double x, bool extrapolate )
{
    int size = xData.size();
    
    int i = 0;                                                                  // find left end of interval for interpolation
    if ( x >= xData[size - 2] )                                                 // special case: beyond right end
    {
        i = size - 2;
    }
    else
    {
        while ( x > xData[i+1] ) i++;
    }
    double xL = xData[i], yL = yData[i], xR = xData[i+1], yR = yData[i+1];      // points on either side (unless beyond ends)
    if ( !extrapolate )                                                         // if beyond ends of array and not extrapolating
    {
        if ( x < xL ) yR = yL;
        if ( x > xR ) yL = yR;
    }
    
    double dydx = ( yR - yL ) / ( xR - xL );                                    // gradient
    
    return yL + dydx * ( x - xL );                                              // linear interpolation
}

// [[Rcpp::export]]
List bird_C(double ZEN, NumericVector wl) {
    
    // ZEN   = solar zenith angle, radians
    // wl    = wavelengths where direct and diffuse should be computed, micrometres
    //
    // Compute direct and diffuse components of spectral irradiance at sea level,
    // given solar zenith angle in air (< 80 degrees). The direct component is
    // calculated normal to the sun's direction.
    // UNITS: watts * metre^-2 * micrometre^-1
    // SOURCE: Bird 1984, Bird and Riordan 1986
    //
    // AIRMAS = AIRMASS
    // C1     = C, INTERPOLATED FOR THE ZENITH ANGLE OF INTEREST
    // C2     = C1, INTERPOLATED FOR THE 24 LAMBDAS.
    // DIF    = DIFFUSE SPECTRAL IRRADIANCE AT 24 WAVELENGTHS
    // DIFUSE = DIFFUSE SPECTRAL IRRADIANCE AT 61 WAVELENGTHS
    // DIR    = DIRECT  SPECTRAL IRRADIANCE AT 24 WAVELENGTHS
    // DIRECT = DIRECT  SPECTRAL IRRADIANCE AT 61 WAVELENGTHS
    // EM0    = M(0), THE AIRMASS EXPRESSION FOR OZONE
    // ROS    = RO (S), THE ALBEDO OF AIR.
    // TA     = AEROSOL TRANSMITTANCE
    // TO     = OZONE TRANSMITTANCE
    // TR     = RAYLEIGH TRANSMITTANCE
    // TW     = WATER VAPOUR TRANSMITTANCE
    // ZEND   = ZENITH ANGLE IN DEGREES

    //==========================================================================
    // VARIABLES

    double Rsize = 24;

    // Value for alpha2 corrected March 1993, changed again Feb 1998
    double ALPHA1 = 1.0274;
    double BETA1 = 0.13238;
    double ALPHA2 = 1.206;
    double BETA2 = 0.116981;

    // LAMBDA = WAVELENGTHS AT WHICH TRANSMITTANCE IS CALCULATED
    static const double LAMBDA[] = {400,410,420,430,440,450,460,470,480,490,500,510,520,530,540,550,570,593,610,630,656,667.6,690,710};
    vector<double> LAMPTH = {0.4,0.41,0.42,0.43,0.44,0.45,0.46,0.47,0.48,0.49,0.5,0.51,0.52,0.53,0.54,0.55,0.57,0.593,0.61,0.63,0.656,0.6676,0.69,0.71};
    
    // EXTER = EXTRA TERRESTRIAL SPECTRAL IRRADIANCE
    static const double EXTER[] = {1479.1,1701.3,1740.4,1587.2,1837.0,2005.0,2043.0,1987.0,2027.0,1896.0,1909.0,1927.0,1831.0,1891.0,1898.0,1892.0,1840.0,1768.0,1728.0,1658.0,1524.0,1531.0,1420.0,1399.0};

    // AV = WATER VAPOUR ABSORPTION COEFFICIENT
    // Feb 1998 - corrected value at wavelength 710 from 0.05 to 0.0125
    static const double AV[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.075,0,0,0,0,0.016,0.0125};

    // AO = OZONE ABSORPTION COEFFICIENT
    // Feb 1998 - corrected value at wavelengths 550, 610, 630, 6676 from 0.095 to 0.085,
    // from 0.132 to 0.120, from 0.120 to 0.090, and from 0.060 to 0.051 respectively
    static const double AO[] = {0,0,0,0,0,0.003,0.006,0.009,0.014,0.021,0.030,0.040,0.048,0.063,0.075,0.085,0.120,0.119,0.120,0.090,0.065,0.051,0.028,0.018};

    // AU = ABSORPTION COEFFICIENT AND GASEOUS AMOUNT
    static const double AU[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.15,0};
    
    
    //==========================================================================
    // DIRECT NORMAL IRRADIANCE FOR 24 WAVELENGTHS

    // Compute direct irradiance transmittance components using AIRMAS=1.9,
    // then use them to compute ROS (air albedo).

    // Bird and Riordan, 1986: M0 = (1 + ho/6370) / sqrt( (cos(ZEN)^2) + 2*ho/6370),
    // where ho=22km, approx max height of ozone
    // Below: TO = exp(-AO * 0.344 * EM0)
    double EM0 = 35. / sqrt(1224. * pow(cos(ZEN),2) + 1.);

    double AIRMAS = 1.9;
    double W = 2; // precipitable water vapor in a vertical path (cm), value from George's script

    NumericVector TR(Rsize);
    NumericVector TW(Rsize);
    NumericVector TO(Rsize);
    NumericVector TA(Rsize);
    NumericVector TU(Rsize);
    NumericVector ROS(Rsize);
    vector<double> DIR(Rsize);
    
    for (int i=0; i < Rsize; i++) {
        // RAYLEIGH SCATTERING
        TR[i] = exp(-AIRMAS / (pow(LAMPTH[i],4) * ( 115.6406 - 1.335/pow(LAMPTH[i],2) )));
        // WATER VAPOUR ABSORPTION
        // Feb 1998 - corrected TW equation as per Bird and Riordan 1986
        TW[i] = exp( (-0.2385 * AV[i] * W *  AIRMAS) / pow((1. + 20.07 * AV[i] * W * AIRMAS),0.45));
        // OZONE
        // From Campbell and Aarup 1989, assume O3=0.344cm (corrected Feb 1998),
        // where O3 = ozone amount in atmosphere
        TO[i] = exp (-AO[i]*0.344*EM0);
        // AEROSOL SCATTERING AND ABSORPTION (TURBIDITY ASSUMED TO BE 0.27)
        if (i < 10) {
            TA[i] = exp(-BETA1 * pow(LAMPTH[i],-ALPHA1) * AIRMAS);
        } else {
            TA[i] = exp(-BETA2 * pow(LAMPTH[i],-ALPHA2) * AIRMAS);
        }
        // UNIFORMLY MIXED GAS ABSORPTION
        // Note: AU = 0 for all but lambda=690nm, so TU = 1 for those values
        // Feb 1998 - Corrected Leckner's value from 118.93 to 118.3 as described in
        // Bird and Riordan, 1986
        // Also, AU(690) changed from .30 to .15
        TU[i] = exp(-1.41 * AU[i] * AIRMAS / pow((1. + 118.3 * AU[i] * AIRMAS),0.45) );
        // ALBEDO OF AIR
        ROS[i] = TO[i] * TW[i] * (TA[i] * (1. - TR[i]) *0.5 + TR[i] * (1. - TA[i]) * 0.22 * 0.928) * TU[i];
    }
    
    // Get ZEND (zenith angle in degrees), and recompute AIRMAS and all direct
    // irradiance transmittance components that use the AIRMAS variable.
    // ***NOTE: do NOT recompute ROS, just use the airmas=1.9 ROS after this
    double ZEND = ZEN * 180/M_PI;
    AIRMAS = 1 / (cos(ZEN) + 0.15 * pow(93.885 - ZEND,-1.253));
    if (AIRMAS < 1) {AIRMAS = 1;}


    
    for (int i=0; i < Rsize; i++) {
        LAMPTH[i] = LAMBDA[i]/1000;
        TR[i] = exp(-AIRMAS / (pow(LAMPTH[i],4) * ( 115.6406 -1.335/pow(LAMPTH[i],2) )));
        TW[i] = exp( (-0.2385 * AV[i] * W *  AIRMAS) / pow((1. + 20.07 * AV[i] * W * AIRMAS),0.45));
        if (i < 10) {
            TA[i] = exp(-BETA1 * pow(LAMPTH[i],-ALPHA1) * AIRMAS);
        } else {
            TA[i] = exp(-BETA2 * pow(LAMPTH[i],-ALPHA2) * AIRMAS);
        }
        TU[i] = exp (-1.41 * AU[i] * AIRMAS / pow((1. + 118.3 * AU[i] * AIRMAS),0.45) );
        // TOTAL DIRECT IRRADIANCE
        DIR[i] = EXTER[i] * TR[i] * TA[i] * TW[i] * TO[i] * TU[i];
    }
    
    
    //==========================================================================
    // DIFFUSE HORIZONTAL IRRADIANCE FOR 24 WAVELENGTHS

    NumericVector XX(Rsize);
    NumericVector R(Rsize);
    NumericVector A(Rsize);
    NumericVector G(Rsize);
    vector<double> DIF(Rsize);

    for (int i=0; i < Rsize; i++) {
        XX[i] = EXTER[i] * cos(ZEN) * TO[i] * TU[i] * TW[i];
        // Rayleigh scattering component
        R[i] = XX[i] * TA[i]*(1-TR[i])*0.5;
        // Aerosol scattering component
        A[i] = XX[i] * TR[i]*(1-TA[i])*0.928*0.82;
        // Component accounting for multiple reflection of irradiance between ground and air
        G[i] = (DIR[i] * cos(ZEN) + (R[i]+A[i])) * 0.05 * ROS[i] / (1 - 0.05 * ROS[i]);
        // Total, including correction factor for sum of R and A
        DIF[i] = (R[i]+A[i]) + G[i];
    }

    
    //==========================================================================
    // INTERPOLATION TO OTHER WAVELENGTHS wl

    vector<double> DIRECT;
    vector<double> DIFFUSE;

    for ( double x : wl ) {
        double ydir = interpolate( LAMPTH, DIR, x, true );
        DIRECT.push_back( ydir );
        double ydif = interpolate( LAMPTH, DIF, x, true );
        DIFFUSE.push_back( ydif );
    }

    // A correction factor (CORF) is applied here to LAMBDA in the original Fortran script,
    // but looks like it might be a Fortran quirk rather than something in the model.

    List Output;
    Output["DIRECT"] = DIRECT;
    Output["DIFFUSE"] = DIFFUSE;
    return(Output);
    
}

// [[Rcpp::export]]
List sam_penguin_C(NumericVector chl, NumericVector I0_direct, NumericVector I0_diffuse, double alphaB, NumericVector Z, double ZENW, NumericVector LAMBDA) {
    
    double nstep = Z.size();
    double dz = Z[1]-Z[0];
    double Rsize = LAMBDA.size();
    
    // https://stackoverflow.com/questions/2236197/what-is-the-easiest-way-to-initialize-a-stdvector-with-hardcoded-elements
    
    // ATLANTIC ZONE
    static const double pc1[] = {0.0744, 0.0773, 0.0778, 0.0776, 0.0773, 0.0745, 0.0746, 0.0748, 0.0723, 0.0668, 0.0612, 0.0571, 0.0546, 0.0516, 0.0488, 0.0446, 0.0401, 0.036, 0.0322, 0.0283, 0.0247, 0.0217, 0.019, 0.017, 0.0151, 0.0134, 0.012, 0.0106, 0.0093, 0.0081, 0.007, 0.0061, 0.0056, 0.0056, 0.0064, 0.007, 0.0089, 0.0104, 0.0103, 0.0118, 0.0105, 0.0102, 0.0113, 0.0112, 0.0116, 0.0119, 0.0128, 0.0135, 0.0116, 0.0103, 0.0104, 0.012, 0.0168, 0.0224, 0.0277, 0.0296, 0.0269, 0.0207, 0.0135, 0.0082, 0.0072};
    static const double pc2[] = {0.0104, 0.0105, 0.0112, 0.0115, 0.0114, 0.0117, 0.0119, 0.0124, 0.0125, 0.0121, 0.0114, 0.0109, 0.0106, 0.0103, 0.0096, 0.009, 0.0084, 0.0081, 0.0079, 0.0078, 0.0077, 0.0073, 0.007, 0.0066, 0.0063, 0.0061, 0.0058, 0.0056, 0.0053, 0.005, 0.0046, 0.0042, 0.0037, 0.0034, 0.0031, 0.003, 0.0029, 0.0028, 0.0028, 0.0025, 0.0024, 0.0025, 0.0027, 0.003, 0.0033, 0.0034, 0.0035, 0.0035, 0.0036, 0.0036, 0.0037, 0.0043, 0.0056, 0.0074, 0.0089, 0.0092, 0.0081, 0.0061, 0.0039, 0.0023, 0.0011};
    static const double rate[] = {0.5582, 0.5706, 0.611, 0.6465, 0.6668, 0.7038, 0.7335, 0.7657, 0.7975, 0.8253, 0.8571, 0.8993, 0.9557, 1.0389, 1.0929, 1.1765, 1.259, 1.3455, 1.3954, 1.4338, 1.4485, 1.4396, 1.434, 1.415, 1.4161, 1.399, 1.365, 1.3307, 1.3254, 1.293, 1.2783, 1.2116, 1.0704, 0.8557, 0.6616, 0.5854, 0.47, 0.429, 0.4673, 0.4162, 0.4494, 0.4422, 0.3926, 0.4014, 0.4103, 0.4147, 0.4125, 0.4179, 0.5266, 0.6462, 0.6428, 0.5651, 0.4828, 0.5287, 0.5984, 0.6433, 0.6482, 0.5811, 0.4876, 0.4428, 0.3376};
    
    static const double AW[] = {0.00663, 0.00530, 0.00473, 0.00444, 0.00454, 0.00478, 0.00495, 0.00530, 0.00635, 0.00751, 0.00922, 0.00962, 0.00979, 0.01011, 0.0106, 0.0114, 0.0127, 0.0136, 0.0150, 0.0173, 0.0204, 0.0256, 0.0325, 0.0396, 0.0409, 0.0417, 0.0434, 0.0452, 0.0474, 0.0511, 0.0565, 0.0596, 0.0619, 0.0642, 0.0695, 0.0772, 0.0896, 0.1100, 0.1351, 0.1672, 0.2224, 0.2577, 0.2644, 0.2678, 0.2755, 0.2834, 0.2916, 0.3012, 0.3108, 0.325, 0.340, 0.371, 0.410, 0.429, 0.439, 0.448, 0.465, 0.486, 0.516, 0.559, 0.624};
    
    double BW500 = 0.00288;
    double BR488 = 0.00027;
    
    NumericVector BW(Rsize);
    NumericVector BBR(Rsize);
    NumericVector AY(Rsize);
    NumericVector Iz_scalar(Rsize);
    NumericVector Iz_planar(Rsize);
    NumericVector mu_d(Rsize);
    
    // Fill the vectors
    for (int i=0; i < Rsize; i++) {
        
        BW[i] = BW500 * pow(LAMBDA[i] / 500, -4.3);
        BBR[i] = 0.5 * BR488 * pow(LAMBDA[i] / 488, -5.3);
        
        double x = -0.014 * (LAMBDA[i] - 440);
        if (x < -600) {
            AY[i] = 0;
        } else {
            AY[i] = exp(x);
        }
        
        // Get total scalar and planar irradiance Iz
        Iz_scalar[i] = I0_direct[i] + I0_diffuse[i];
        Iz_planar[i] = I0_direct[i] * cos(ZENW) + I0_diffuse[i];//should I0_diffuse be multiplied by 0.831 here?
        mu_d[i] = (I0_direct[i] * cos(ZENW) + I0_diffuse[i] * 0.831) / Iz_scalar[i];
        
    }
    
    // PAR and irradiance at each depth
    NumericVector PARz(nstep);
    NumericVector XXz(nstep);
    NumericVector yelsubz(nstep);
    
    for (int i=0; i < nstep; i++) {
        
        // Chlorophyll at depth i
        double chlz = chl[i];
        
        // Absorption and backscattering for depth i, Devred 2006
        NumericVector AC(Rsize);
        for (int j=0; j < Rsize; j++) {
            
            double x = (-1) * rate[j] * chlz;
            if (x < -600) {
                AC[j] = pc1[j] + pc2[j] * chlz;
            } else {
                AC[j] = pc1[j] * (1 - exp(x)) + pc2[j] * chlz;
            }
            
        }
        
        // Remember that indexing starts at 0 in C++, so AC[8] is the 9th value
        double AC440 = AC[8];
        double yelsub;
        
        if (chlz > 0) {
            
            yelsub = (log10(chlz) + 2) * 0.3;
            
            if (chlz > 1) {
                yelsub = yelsub * (log10(chlz) + 1);
            }
            if (yelsub < 0.0001) {
                yelsub = 0.0001;
            }
            
        } else if (chlz == 0) {
            
            yelsub =  0.0001;
            
        }
        
        yelsubz[i] = yelsub;
        
        double AY440 = yelsub * AC440;
        double ACmean = GetAverage(AC);
        
        NumericVector A(Rsize);
        for (int j=0; j < Rsize; j++) {
            A[j] = AW[j] + AC[j] + AY440 * AY[j] + 2 * BBR[j];
        }
        
        // Loisel and Morel 1998
        double bbtilda = (0.78 + 0.42 * (-1) * log10(chlz)) * 0.01;
        if (bbtilda < 0.0005) {
            bbtilda = 0.0005;
        } else if (bbtilda > 0.01) {
            bbtilda = 0.01;
        }
        double BC660 = 0.407 * pow(chlz, 0.795);
        
        double tmpPARz = 0;
        double tmpXXz = 0;
        
        for (int j=0; j < Rsize; j++) {
            
            double BC;
            double BB;
            
            double x = BC660 * pow(660 / LAMBDA[j], (-1) * log10(chlz));
            if (x < 0) {
                BC = 0;
            } else {
                BC = x;
            }
            
            BB = BC * bbtilda + BW[j] * 0.5;
            
            // Compute PAR for this depth (integration over wavelength, using Riemann sum)
            tmpPARz = tmpPARz + (Iz_scalar[j] * 5);
            
            tmpXXz = tmpXXz + (alphaB * Iz_planar[j] * AC[j] * 6022 / ACmean / mu_d[j] / 2.77 / 36) * 5;
            
            // Find K, diffuse attenuation coefficient, for this layer at each wavelength, and use it to update Iz at this depth for each wavelength.
            double K = (A[j] + BB) / mu_d[j];
            
            double y = (-1) * K * dz;
            if (y < -600) {
                Iz_planar[j] = 0;
                Iz_scalar[j] = 0;
            } else {
                Iz_planar[j] = Iz_planar[j] * exp(y);
                Iz_scalar[j] = Iz_scalar[j] * exp(y);
            }
            
        }
        
        XXz[i] = tmpXXz;
        
        // Convert BIO units (E/m^2/hr) to Laval units (uE/m^2/s) for easier comparison
        PARz[i] = tmpPARz * pow(10, 6) / 3600;
        
        // End calculations if we've reached euphotic depth (when PAR <= 1% surface PAR)
        if ((i > 0) && (PARz[i] < (0.01 * PARz[0]))) {
            
            // Return PAR, XX, euphotic depth, index of euphotic depth
            
            // https://stackoverflow.com/questions/52187622/retrieve-object-from-a-list-rcpp
            List Output;
            Output["PARz"] = PARz;  
            Output["XXz"] = XXz;
            Output["EuphoD"] = Z[i]; // Z is zero-indexed, but i also starts at 0 in this loop, so this is right
            Output["id"] = i + 1; // this is the index returned to the R section of the script, so for example if i=0 here, you must add 1 since R indexing starts at 1
            Output["yelsubz"] = yelsubz;
            return(Output);
            
        }
        
    }
    
    // Remember C++ is zero-indexed (and R starts at 1), so adjust indices here to match with the original R code
    List Output;
    Output["PARz"] = PARz;  
    Output["XXz"] = XXz;
    Output["EuphoD"] = Z[nstep - 1]; // last index in C++
    Output["id"] = nstep; // last index in R (this value is returned and used in the R section of the script)
    Output["yelsubz"] = yelsubz;
    return(Output);
    
}
