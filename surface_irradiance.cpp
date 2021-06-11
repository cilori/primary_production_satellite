#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace std;
using namespace Rcpp;

// // [[Rcpp::export]]
// NumericVector Watt_perm_2_to_microMol_per_m2_per_s(NumericVector Ed, NumericVector lambda){
//     double h = 6.626e-34; // in J.s
//     double c = 2.998e8; // in m.s-1
//     double Na = 6.022e23; // Avogadro number
//     NumericVector Np = Ed * lambda * 1e-9/(h*c);
//     NumericVector new_Ed = Np / (Na * 1e-6);
//     return(new_Ed);
// }


// [[Rcpp::export]]
double comp_sco3_c(double ylat,double xlon,int jday,double ro3,double rad) {
    double rlat = ylat/rad;
    double rlon = xlon/rad;
    double to3 = 235.0 + (150.0+40.0*sin(0.9865*(jday-30.0)) + 20.0*sin(3.0*(rlon)))*sin(1.28*rlat)*sin(1.28*rlat);
    if (ro3 > 0) {
        to3 = ro3;
    } else {
        to3 = 235.0 + (150.0+40.0*sin(0.9865*(jday-30.0)) + 20.0*sin(3.0*(rlon)))*sin(1.28*rlat)*sin(1.28*rlat);
    }
    double sco3 = to3*1e-3;
    return(sco3);
}


// [[Rcpp::export]]
NumericVector solar_irr_c(NumericVector Fobar, int jday, double pi2) {
    //  Correct for Earth-Sun distance
    //  modify date if it is last day of a leap year
    if(jday >365) {jday=365.0;}
    NumericVector Fo = Fobar * pow((1.0+1.67E-2*cos(pi2*(jday-3)/365.0)),2);
    return(Fo);
}


// [[Rcpp::export]]
List sfcrfl_c(double rad, double theta, double ws) {

    //  Computes surface reflectance for direct (rod) and diffuse (ros)
    //  components separately, as a function of theta, wind speed or stress.
    
    double rn = 1.341;      //index of refraction of pure seawater
    double roair = 1.2E3;   //density of air g/m3
    double cn;
    double rof;
    double rosps;
    double rospd;
    double rtheta;
    double sintr;
    double rthetar;
    double rmin;
    double rpls;
    double sinp;
    double tanp;
    double a;
    double b;

    //  Foam and diffuse reflectance
    if (ws > 4.0) {
        if (ws < 7.0) {
            cn = 6.2E-4 + 1.56E-3/ws;
            rof = roair*cn*2.2E-5*ws*ws - 4.0E-4;
        } else {
            cn = 0.49E-3 + 0.065E-3*ws;
            rof = (roair*cn*4.5E-5 - 4.0E-5)*ws*ws;
        }
        rosps = 0.057;
    } else {
        rof = 0.0;
        rosps = 0.066;
    }

    //  Direct
    //  Fresnel reflectance for theta < 40, ws < 2 m/s
    if (theta < 40.0 || ws < 2.0) {
        if (theta == 0.0) {
            rospd = 0.0211;
        } else {
            rtheta = theta/rad;
            sintr = sin(rtheta)/rn;
            rthetar = asin(sintr);
            rmin = rtheta - rthetar;
            rpls = rtheta + rthetar;
            sinp = (sin(rmin)*sin(rmin))/(sin(rpls)*sin(rpls));
            tanp = (tan(rmin)*tan(rmin))/(tan(rpls)*tan(rpls));
            rospd = 0.5*(sinp + tanp);
        }
    } else { // Empirical fit otherwise
        a = 0.0253;
        b = -7.14E-4*ws + 0.0618;
        rospd = a*exp(b*(theta-40.0));
    }

    //  Reflectance totals
    double rodaux = rospd + rof;
    double rosaux = rosps + rof;
    
    List Output;
    Output["rod"] = rodaux;
    Output["ros"] = rosaux;
    return(Output);
    
}





// [[Rcpp::export]]
List navaer_c(double rh,double am,double wsm,double ws,double vis) {
    
    //  Computes aerosol parameters according to a simplified version
    //  of the Navy marine aerosol model.
    
    NumericVector ro = {0.03,0.24,2.0};
    NumericVector r = {0.1,1.0,10.0};
    
    double rlam = 0.55;
    
    //  Relative humidity factor
    if (rh >= 100.0) rh = 99.9;
    double rnum = 2.0 - rh/100.0;
    double rden = 6.0*(1.0-rh/100.0);
    double frh = pow((rnum/rden),0.333);

    //  Size distribution amplitude components
    NumericVector a = {2000.0*am*am,5.866*(wsm-2.2),0.01527*(ws-2.2)*0.05};
    if (a[1] < 0.5) a[1] = 0.5;
    if (a[2] < 1.4E-5) a[2] = 1.4E-5;

    //  Compute size distribution at three selected radii according to Navy method
    NumericVector dndr(3);
    for (int n=0; n < 3; n++) {
        for (int i=0; i < 3; i++) {
            dndr[n] = dndr[n] + a[i]*exp((-1)*pow(log(r[n]/frh/ro[i]),2))/frh;
        }
    }
    
    //  Least squares approximation
    double sumxy = sum(log10(r)*log10(dndr));
    double sumx2 = sum(pow(log10(r),2));
    double gama = sumxy/sumx2;
    double alphaux = -(gama+3.0);

    //  Compute beta
    double cext = 3.91/vis;
    double betaux = cext*pow(rlam,alphaux);

    //  Compute asymmetry parameter -- a function of alpha
    double asympaux;
    if (alphaux > 1.2) {
        asympaux = 0.65;
    } else if (alphaux > 0.0) {
        asympaux = 0.82;
    } else {
        asympaux = -0.14167*alphaux + 0.82;
    }
    
    //  Single scattering albedo at 550; function of RH
    double waux = (-0.0032*am + 0.972)*exp(3.06E-4*rh);
    
    List Output;
    Output["alpha"] = alphaux;
    Output["beta"] = betaux;
    Output["asymp"] = asympaux;
    Output["wa"] = waux;
    return(Output);

}



// [[Rcpp::export]]
List atmodd_v2_c(double rod, double ros, double rad, NumericVector lam, double theta,
                 NumericVector oza, NumericVector ag, NumericVector aw,
                 double sco3, double p, double p0, double wv, double rh, double am,
                 double wsm, double ws, double vis, NumericVector Fo) {
    
    // rod and ros calculated before calling this version of atmodd
    
    //  Model for atmospheric transmittance of solar irradiance through
    //  a maritime atmosphere.  Computes direct and diffuse separately.
    //  Includes water vapor and oxygen absorption.
            
    //  Compute atmospheric path lengths (air mass); pressure-corrected
    double cosunz = cos(theta/rad);
    //  Modified March, 1994 according to Kasten and Young 1989.
    
    double rex = -1.6364;
    double rtmp = pow((96.07995-theta),rex);
    double rm = 1.0/(cosunz+0.50572*rtmp);
    double rmp = p/p0*rm;
    double otmp = pow((cosunz*cosunz+44.0/6370.0),0.5);
    double rmo = (1.0+22.0/6370.0)/otmp;

    //  Obtain aerosol parameters; simplified Navy aerosol model
    List res_aer = navaer_c(rh,am,wsm,ws,vis);
    double alpha = res_aer["alpha"];
    double beta = res_aer["beta"];
    double asymp = res_aer["asymp"];
    double wa = res_aer["wa"];
    double eta = -alpha;
    
    //   Forward scattering probability
    double alg = log(1.0-asymp);
    double afs = alg*(1.459+alg*(.1595+alg*.4129));
    double bfs = alg*(.0783+alg*(-.3824-alg*.5874));
    double Fa = 1.0 - 0.5*exp((afs+bfs*cosunz)*cosunz);

    // Compute spectral irradiance
    //   Rayleigh, by Bird's method
    NumericVector rlam = lam*1.0E-3;
    NumericVector tr = 1.0/(115.6406*pow(rlam,4) - 1.335*pow(rlam,2));
    NumericVector rtra = exp(-tr*rmp);
    //    Ozone
    NumericVector to = oza*sco3;   //optical thickness
    NumericVector otra = exp(-to*rmo);   //transmittance
    //   Aerosols
    NumericVector ta = beta*pow(rlam,eta);
    NumericVector atra = exp(-ta*rm);
    NumericVector taa = exp(-(1.0-wa)*ta*rm);
    NumericVector tas = exp(-wa*ta*rm);
    //   Oxygen/gases
    NumericVector gtmp = pow((1.0 + 118.3*ag*rmp),0.45);
    NumericVector gtmp2 = -1.41*ag*rmp;
    NumericVector gtra = exp(gtmp2/gtmp);
    //   Water Vapor
    NumericVector  wtmp = pow((1.0+20.07*aw*wv*rm),0.45);
    NumericVector wtmp2 = -0.2385*aw*wv*rm;
    NumericVector wtra = exp(wtmp2/wtmp);

    //  Direct irradiance
    NumericVector Ediraux = Fo*cosunz*rtra*otra*atra*gtra*wtra*(1.0-rod);

    //   Diffuse irradiance
    NumericVector dray = Fo*cosunz*gtra*wtra*otra*taa*0.5*(1.0-pow(rtra,.95));
    NumericVector daer = Fo*cosunz*gtra*wtra*otra*pow(rtra,1.5)*taa*Fa*(1.0-tas);

    //  Total diffuse
    NumericVector Edifaux = (dray + daer)*(1.0-ros);
    
    NumericVector Edaux = Ediraux + Edifaux;
    
    List Output;
    Output["Edir"] = Ediraux;
    Output["Edif"] = Edifaux;
    Output["Ed"] = Edaux;
    return(Output);

}


// [[Rcpp::export]]
List surface_irr_c(NumericVector zendR, double rad, NumericVector lam,
                   NumericVector oza, NumericVector ag, NumericVector aw,
                   double sco3, double p0, double wv, double rh, double am,
                   double wsm, double ws, double pxlvis, NumericVector Fo) {

    // Ed(0+,lambda)
    NumericMatrix Ed(23,301); // total in mW.m-2.nm-1
    
    // Eq(0-)
    NumericMatrix Eqdirw(23,301); // direct
    NumericMatrix Eqdifw(23,301); // diffus
    // NumericMatrix Eqdw(23,301); // total in Einstein.m-2
    
    List irr_gc;
    List irrm_gc;
    List res_sfc;
    double rod;
    double ros;
    
    NumericVector Edir(23);
    NumericVector Edif(23);
    NumericVector Ed_row(301);
    
    for (int i=0; i < 23; i++) {
        
        if (zendR[i] < 90.) {
            
            // above water spectral irradiance
            rod = 0.0;
            ros = 0.0;
            irr_gc = atmodd_v2_c(rod, ros, rad, lam, zendR[i], oza, ag, aw, sco3,
                                 p0, p0, wv, rh, am, wsm, ws, pxlvis, Fo);
            Ed_row = irr_gc["Ed"];
            Ed(i, _) = Ed_row;
                
            // Under water spectral irradiance.
            res_sfc = sfcrfl_c(rad,zendR[i],ws);
            rod = res_sfc["rod"];
            ros = res_sfc["ros"];
            irrm_gc = atmodd_v2_c(rod, ros, rad, lam, zendR[i], oza, ag, aw, sco3,
                                  p0, p0, wv, rh, am, wsm, ws, pxlvis, Fo);
            
            Edir = irrm_gc["Edir"];
            Eqdirw(i, _) = Edir;
            Edif = irrm_gc["Edif"];
            Eqdifw(i, _) = Edif;

        }
    }
    
    List Output;
    Output["Ed"] = Ed;
    Output["Eqdirw"] = Eqdirw;
    Output["Eqdifw"] = Eqdifw;
    return(Output);
    
}

