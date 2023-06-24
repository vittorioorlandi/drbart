#ifndef GUARD_info_h
#define GUARD_info_h
#include <cmath>

//data
class dinfo {
public:
   dinfo() {p=0;n=0;x=0;y=0;}
   size_t p;  //number of vars
   size_t n;  //number of observations
   double *x; // jth var of ith obs is *(x + p*i+j)
   double *y; // ith y is *(y+i) or y[i]
};

//prior and mcmc
struct pinfo
{
   //mcmc info
   double pbd = 1.0; // prob of birth / death
   double pb = 0.5;  // prob of birth given birth / death
   
   //prior info
   // prior prob a bot node splits is alpha / (1 + depth) ^ beta
   double alpha = 0.95;
   double beta = 2.0;
   double tau = 1.0;
   
   //sigma
   double sigma = 1.0;
   
   pinfo() = default; 
   
   pinfo(double pbd, double pb, 
         double alpha, double beta, double tau, 
         double sigma) {
      this->pbd = pbd;
      this->pb = pb;
      this->alpha = alpha;
      this->beta = beta;
      this->tau = tau;
      this->sigma = sigma;
   }
   
   pinfo(double pbd, double pb, double alpha, double beta, 
         double miny, double maxy, double kfac, int m, 
         double sigma) {
      this->pbd = pbd;
      this->pb = pb;
      this->alpha = alpha;
      this->beta = beta;
      this->tau = (maxy - miny) / (2 * kfac * sqrt((double) m));
      this->sigma = sigma;
   }

};

//sufficient statistics for 1 node
class sinfo
{
public:
   sinfo() {n0=0.0;n=0;sy=0.0;sy2=0.0;}
   double n0; //unweighted sample counts
   double n;
   double sy;
   double sy2;
};

#endif
