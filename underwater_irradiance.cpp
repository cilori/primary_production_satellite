#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
#include <cmath>
#include <vector>
using namespace std;
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix underwater_irradiance_c(NumericMatrix Eqdirw, NumericMatrix Eqdifw, NumericVector chlpix,
                                      NumericVector lambda, NumericVector pc1f, NumericVector pc2f,
                                      NumericVector ratef, double fracys, double Sy, NumericVector asw,
                                      NumericVector zendR, NumericVector zendw, double alphaB, double PBm,
                                      NumericVector BW, int hr1, int hr2) {
    
    double nstep = chlpix.size();
    double Rsize = lambda.size();
    NumericMatrix PP(23, nstep);
    NumericVector Kddir(Rsize);
    NumericVector Kddif(Rsize);
    NumericVector tmp_Eqdirw(Rsize);
    NumericVector tmp_Eqdifw(Rsize);
    double PIz;
    double apbar;
    double chlz = chlpix[0];
    
    // To compute the underwater light field, we need to know total absorption and total backscattering
    // a = aw + aphy + acdom, bb = bbw + bbp, this as to be define for each wavelength
    
    // backscattering is modelled according to Loisel and Morel 1998
    double bbtilda = (0.78 + 0.42 * (-1) * log10(chlz)) * 0.01; // bbtilda et le rapport entre bb et b.
    if (bbtilda < 0.0005) {
        bbtilda = 0.0005;
    } else if (bbtilda > 0.01) {
        bbtilda = 0.01;
    }
    double BC660 = 0.407 * pow(chlz, 0.795); // ici, bc corresponds to b_chl a 660 nm (scattering by chlorophyll-a)
    
    NumericVector BC = BC660*pow((660/lambda),(-log10(chlz))); // BC spectral dependence depends on chlz,
    // this means that at chl = 1 tere is no spectral dependence.This is weird and has to be checked
    BC[BC < 0] = 0;
    NumericVector bbt = BC*bbtilda + BW*0.5; // Here, we get total backscattering.
    
    //phytoplankton absorption //## !!!!! THIS HAS TO BE UPDATED.
    NumericVector aphy = pc1f*(1-exp(-ratef*chlz)) + pc2f*chlz;
    
    // Yellow substances and detritus absorption is derived from chl
    NumericVector ays = fracys * aphy[43] * exp(-Sy * (lambda -443));
    
    // Now we can compute total absorption
    NumericVector at = asw + ays + aphy;
    
    // Here we compute the photosynthetic action spectrum according to Kyewangala et al. 1997
    apbar = mean(aphy);
    NumericVector alphaBc = alphaB * aphy / apbar;
    
    
    for (int it=hr1; it <= hr2; it++) { // loop on the time of day (for daylight zendR)
        
        // get E0 at surface, to be reduced at each depth based on Kd
        tmp_Eqdirw = Eqdirw(it, _);
        tmp_Eqdifw = Eqdifw(it, _);
        
        // calculate PP for the first layer
        PIz = 1/cos(zendw[it]) * sum(alphaBc * tmp_Eqdirw) + 1.20 * sum(alphaBc * tmp_Eqdifw);
        PP(it, 0) = chlz * PIz/sqrt(1 + pow((PIz/PBm),2));
        
        // calculate attenuation coefficient for each wavelength (for uniform profile, it's the same at every depth)
        Kddir = exp(-0.5 * (at + bbt) / cos(zendw[it]));
        Kddif = exp(-0.5 * (at + bbt) / 0.83);
        
        for (int iz=1; iz < nstep; iz++) { // loop on depth, we start at the second step as the first one is Eqdw (note: C++ uses zero-indexing)
            
            tmp_Eqdirw = tmp_Eqdirw * Kddir;
            tmp_Eqdifw = tmp_Eqdifw * Kddif;

            PIz = 1/cos(zendw[it]) * sum(alphaBc * tmp_Eqdirw) + 1.20 * sum(alphaBc * tmp_Eqdifw);
            PP(it, iz) = chlz * PIz/sqrt(1 + pow((PIz/PBm),2));
            
        } // end loop on depth
        
    } // end loop on time of day (xhr)
    
    return(PP);
    
}

// [[Rcpp::export]]
NumericMatrix underwater_irradiance_NU_c(NumericMatrix E0dir, NumericMatrix E0dif, NumericVector chlpix,
                                         NumericVector lambda, NumericVector pc1f, NumericVector pc2f,
                                         NumericVector ratef, double fracys, double Sy, NumericVector asw,
                                         NumericVector zendR, NumericVector zendw, double alphaB, double PBm,
                                         NumericVector BW, int hr1, int hr2) {
    
    double nstep = chlpix.size();
    double Rsize = lambda.size();
    NumericMatrix PP(23, nstep);
    NumericVector BC(Rsize);
    NumericVector bbt(Rsize);
    NumericVector aphy(Rsize);
    NumericVector ays(Rsize);
    NumericVector at(Rsize);
    NumericVector Kddif(Rsize);
    NumericVector alphaBc(Rsize);
    NumericVector E0dir_num(Rsize);
    double chlz;
    double bbtilda;
    double BC660;
    double apbar;
    double PIz;
    NumericMatrix tmp_E0dir(23, Rsize);
    NumericMatrix tmp_E0dif(23, Rsize);
    
    for (int iz=1; iz < nstep; iz++) { // loop on depth, we start at the second step as the first one is Eqdw (note: C++ uses zero-indexing)
        
        chlz = chlpix[iz];

        // To compute the underwater light field, we need to know total absorption and total backscattering
        // a = aw + aphy + acdom, bb = bbw + bbp, this as to be define for each wavelength

        // backscattering is modelled according to Loisel and Morel 1998
        bbtilda = (0.78 + 0.42 * (-1) * log10(chlz)) * 0.01; // bbtilda et le rapport entre bb et b.
        if (bbtilda < 0.0005) {
            bbtilda = 0.0005;
        } else if (bbtilda > 0.01) {
            bbtilda = 0.01;
        }
        BC660 = 0.407 * pow(chlz, 0.795); // ici, bc corresponds to b_chl a 660 nm (scattering by chlorophyll-a)

        BC = BC660*pow((660/lambda),(-log10(chlz))); // BC spectral dependence depends on chlz,
        // this means that at chl = 1 tere is no spectral dependence.This is weird and has to be checked
        BC[BC < 0] = 0;
        bbt = BC*bbtilda + BW*0.5; // Here, we get total backscattering.

        //phytoplankton absorption //## !!!!! THIS HAS TO BE UPDATED.
        aphy = pc1f*(1-exp(-ratef*chlz)) + pc2f*chlz;

        // Yellow substances and detritus absorption is derived from chl
        ays = fracys * aphy[43] * exp(-Sy * (lambda -443));

        // Now we can compute total absorption
        at = asw + ays + aphy;

        // Here we compute the photosynthetic action spectrum according to Kyewangala et al. 1997
        apbar = mean(aphy);
        alphaBc = alphaB * aphy / apbar;
        
        Kddif = exp(-0.5 * (at + bbt) / 0.83);
        E0dir_num = -0.5 * (at + bbt);

        if (iz==1) {
            for (int it=hr1; it <= hr2; it++) {
                PIz = 1/cos(zendw[it]) * sum(alphaBc * E0dir(it, _)) + 1.20 * sum(alphaBc * E0dif(it, _));
                PP(it, 0) = chlz * PIz/sqrt(1 + pow((PIz/PBm),2));
                tmp_E0dir(it, _) = E0dir(it, _) * exp(E0dir_num / cos(zendw[it]));
                tmp_E0dif(it, _) = E0dif(it, _) * Kddif;
                PIz = 1/cos(zendw[it]) * sum(alphaBc * tmp_E0dir(it, _)) + 1.20 * sum(alphaBc * tmp_E0dif(it, _));
                PP(it, 1) = chlz * PIz/sqrt(1 + pow((PIz/PBm),2));
            }
        } else {
            for (int it=hr1; it <= hr2; it++) {
                tmp_E0dir(it, _) = tmp_E0dir(it, _) * exp(E0dir_num / cos(zendw[it]));
                tmp_E0dif(it, _) = tmp_E0dif(it, _) * Kddif;
                PIz = 1/cos(zendw[it]) * sum(alphaBc * tmp_E0dir(it, _)) + 1.20 * sum(alphaBc * tmp_E0dif(it, _));
                PP(it, iz) = chlz * PIz/sqrt(1 + pow((PIz/PBm),2));
            }
        }
        
    } // end loop on depth
    
    return(PP);
    
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
vector<double> higher_res(vector<double> xold, vector<double> xnew, vector<double> yold) {
    vector<double> ynew;
    for ( double x : xnew ) {
        double ynew_tmp = interpolate( xold, yold, x, true );
        ynew.push_back( ynew_tmp );
    }
    return(ynew);
}


