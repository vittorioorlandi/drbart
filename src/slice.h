#ifndef slice_h
#define slice_h

#include "funs.h"

class logdensity {
  public:
  virtual double val(double y) = 0;
};

class ld_norm: public logdensity {
  public:
  double mu;
  double sigma;
  double val(double y) { return(R::dnorm(y, mu, sigma, 1)); }
  ld_norm(double mu_, double sigma_) { mu=mu_; sigma=sigma_; }
};

class ld_bartU: public logdensity {
  public:
  double f; //fit that doesn't depend on u
  double sigma;
  bool scalemix;
  
  size_t i; //observation index
  std::vector<tree> using_u; //set of trees that split on u
  xinfo xi; 
  dinfo di; //pointers to xi, di
  
  std::vector<tree> using_uprec; //set of trees that split on u
  xinfo xiprec; 
  dinfo diprec; //pointers to xi, di
  
  double yobs;
  int j, p;
  double val(double y) { 
    double oldx = *(di.x + i*di.p);
    *(di.x + i*di.p) = y;
    if(scalemix) *(diprec.x + i*diprec.p) = y;
    double mm = f + fit_i(i, using_u, xi, di);
    double pp = sigma;
    if(scalemix) {
      pp /= sqrt(fit_i_mult(i, using_uprec, xiprec, diprec));
    }
    *(di.x + i*di.p) = oldx;
    if(scalemix) *(diprec.x + i*diprec.p) = oldx;
    return(R::dnorm(yobs, mm, pp, 1)); 
  }
  ld_bartU(double f_, double sigma_) { f=f_; sigma=sigma_; scalemix=false;}
};

double slice(double x0, logdensity* g, double w=1., double m=INFINITY, 
             double lower=-INFINITY, double upper=INFINITY);
#endif slice_h
