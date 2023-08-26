#include "slice.h"

// typically called with w = 1, m = INFINITY, lower = 0, upper = 1
double slice(double x0, logdensity* g, double w, double m, 
             double lower, double upper) { // , 
             // dinfo& di, dinfo& diprec, 
             // std::vector<tree>& using_u, std::vector<tree>& using_uprec) {
  double x1; // new sample 

	// std::ofstream treef;
	// treef.open ("slice.txt");
  
  // double gx0 = g->val(x0, di, diprec, using_u, using_uprec); // current loglik
  double gx0 = g->val(x0); // current loglik
 	// treef << "basic comps" << std::endl; 
  double logy = gx0 - R::rexp(1.);
  double u = R::runif(0., w); 
  double L = x0 - u;
  double R = x0 + (w - u);
	// MAYBE CAN AUTOMATICALLY GET A LARGE ENOUGH INTERVAL 
	// DIRECTLY FROM THE CUTPOINTS 
  while(true) {
    R_CheckUserInterrupt();
    if(L<=lower) { break; }
    // if(g->val(L, di, diprec, using_u, using_uprec) <= logy) { break; }
    if(g->val(L) <= logy) { break; }
    L -= w;
  }
  while(true) {
    R_CheckUserInterrupt();
    if(R>=upper) { break; }
    // if(g->val(R, di, diprec, using_u, using_uprec) <= logy) { break; }
    if(g->val(R) <= logy) { break; }
    R += w;
  }
//  	treef << "found interval" << std::endl; 
// 	treef.close();
  if(L<lower) {L=lower;}
  if(R>upper) {R=upper;}

	// [L, R] is our interval to sample x1 uniformly from 
  
  while(true) {
    R_CheckUserInterrupt();
    x1 = R::runif(L, R);
    // double gx1 = g->val(x1, di, diprec, using_u, using_uprec);
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
