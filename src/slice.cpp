#include "slice.h"

double slice(double x0, logdensity* g, double w, double m, 
             double lower, double upper) {
  double x1;             
  
  double gx0 = g->val(x0);
  double logy = gx0 - R::rexp(1.);
  double u = R::runif(0., w);
  double L = x0 - u;
  double R = x0 + (w-u);
  
  while(true) {
    R_CheckUserInterrupt();
    if(L<=lower) { break; }
    if(g->val(L)<=logy) { break; }
    L -= w;
  }
  while(true) {
    R_CheckUserInterrupt();
    if(R>=upper) { break; }
    if(g->val(R)<=logy) { break; }
    R += w;
  }
  
  if(L<lower) {L=lower;}
  if(R>upper) {R=upper;}
  
  while(true) {
    R_CheckUserInterrupt();
    x1 = R::runif(L, R);
    double gx1 = g->val(x1);
    if(gx1>=logy) { break; }
    if(x1>x0) {
      R = x1;
    } else {
      L = x1;
    }
  }
  
  return(x1);
}
