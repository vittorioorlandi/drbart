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

// not even sure what the point of this being a class is... 
// like yes the inheritance hierarchy is nice but we really don't use it 
// and we only end up making one copy of this so it seems like it would 
// be more natural for it to be a function. for now keeping it a class for 
// practice's sake. 
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
  //   // temporary method to agree with virtual method in base class...? 
  //   return(0); 
  // }
  // 
  // double val(double y, dinfo& di, dinfo& diprec, 
  //            std::vector<tree>& using_u, std::vector<tree>& using_uprec) { 
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
  // ld_bartU(double f, double sigma, bool scalemix) {
  //   this->f = f; 
  //   this->sigma = sigma;
  //   this->scalemix = scalemix;
  // }
  // 
  // // ld_bartU(double f_, double sigma_, dinfo& di, dinfo& diprec) {
  // //   f = f_; 
  // //   sigma = sigma_;
  // //   scalemix = false;
  // // }
};

double slice(double x0, logdensity* g, double w=1., double m=INFINITY, 
             double lower=-INFINITY, double upper=INFINITY);
#endif
