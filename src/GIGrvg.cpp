
/*---------------------------------------------------------------------------*/
/* header files */

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>

#include "GIGrvg.h"


/*---------------------------------------------------------------------------*/
/* Prototypes of private functions                                           */
/*---------------------------------------------------------------------------*/

static double _gig_mode(double lambda, double omega);

/* Type 1 */
static void _rgig_ROU_noshift (double *res, int n, double lambda, double lambda_old, double omega, double alpha);

/* Type 4 */
static void _rgig_newapproach1 (double *res, int n, double lambda, double lambda_old, double omega, double alpha);

/* Type 8 */
static void _rgig_ROU_shift_alt (double *res, int n, double lambda, double lambda_old,  double omega, double alpha);

static double _unur_bessel_k_nuasympt (double x, double nu, int islog, int expon_scaled);


/*****************************************************************************/
/* API                                                                       */
/*****************************************************************************/

SEXP dgig(SEXP sexp_x, SEXP sexp_lambda, SEXP sexp_chi, SEXP sexp_psi, SEXP sexp_logvalue)
/*---------------------------------------------------------------------------*/
/* Evaluate density of GIG distribution.                                     */
/*                                                                           */
/* Parameters:                                                               */
/*   x ....... argument(s) for density                                       */
/*   lambda .. parameter for distribution                                    */
/*   chi   ... parameter for distribution                                    */
/*   psi   ... parameter for distribution                                    */
/*   logvalue  if TRUE the logarithm of the density will be returned         */
/*                                                                           */
/* Return:                                                                   */
/*   density at x                                                            */
/*---------------------------------------------------------------------------*/
{ 
  /* Arguments */
  double *x;                 /* array of points */
  int nx;                    /* length of x */
  double lambda, chi, psi;   /* parameters of distribution */
  int logvalue;
  double LOGNORMCONSTANT;    /* log of normalization constant for density */
  int i;
  SEXP sexp_res;             /* results of computations */
  double *res;
  double res_err;            /* result in case error */
  int err;                   /* indicate invalid input */

  /* get array of points for which the density has to be evaluated */
  nx = length(sexp_x);
  x = REAL(AS_NUMERIC(sexp_x));

  /* compute logarithm? */
  logvalue = *(LOGICAL( AS_LOGICAL(sexp_logvalue) ));

  /* extract parameters of distribution */
  lambda = *(REAL(AS_NUMERIC(sexp_lambda)));
  chi = *(REAL(AS_NUMERIC(sexp_chi)));
  psi = *(REAL(AS_NUMERIC(sexp_psi)));

  /* allocate array for results of density evaluations */
  PROTECT(sexp_res = NEW_NUMERIC(nx));
  res = REAL(sexp_res);

  /* check GIG parameters: */
  err = FALSE;
  /* we handle invalid input as in the R core function dnorm() */
  if (! (R_FINITE(lambda) && R_FINITE(chi) && R_FINITE(psi)) ) {
    if (ISNAN(lambda) || ISNAN(chi)  || ISNAN(psi)) {
      /* NA or NaN */
      err = TRUE;
      res_err = (ISNA(lambda) || ISNA(chi)  || ISNA(psi)) ? NA_REAL : R_NaN;
    }
    else {
      /* Inf */
      err = TRUE;
      res_err = R_NaN;
      warning("NaNs produced");
    }
  } 
  else {
    /* check GIG parameters */
    if ( (chi <  0. || psi < 0)      || 
	 (chi == 0. && lambda <= 0.) ||
	 (psi == 0. && lambda >= 0.) ) {
      err = TRUE;
      res_err = R_NaN;
      warning("NaNs produced");
    }
  }

  /* invalid arguments */
  if (err) {
    warning("invalid parameters for GIG distribution: lambda=%g, chi=%g, psi=%g",
	    lambda, chi, psi);
    for(i=0; i<nx; i++)
      res[i] = ISNA(x[i]) ? NA_REAL : res_err;
    UNPROTECT(1);
    return sexp_res;
  }

  /* compute normalization constant */
  if (psi == 0.) {
    /* case: Inverse Gamma distribution */
    LOGNORMCONSTANT = -lambda * log(0.5*chi) - lgammafn(-lambda);
  }
  else if (chi == 0.) {
    /* case: Gamma distribution */
    LOGNORMCONSTANT = lambda * log(0.5*psi) - lgammafn(lambda);
  }
  else {
    /* general case: */
    /*   pow(psi/chi, lambda/2.) / (2. * bessel_k(sqrt(psi*chi),lambda,1)); */
    double alambda = fabs(lambda);
    double beta = sqrt(psi*chi);
    LOGNORMCONSTANT = 0.5*lambda*log(psi/chi) - M_LN2;
    if (alambda < 50.) {
      /* threshold value 50 is selected by experiments */
      LOGNORMCONSTANT -= log(bessel_k(beta, alambda, 2)) - beta;
    }
    else {
      LOGNORMCONSTANT -= _unur_bessel_k_nuasympt(beta, alambda, TRUE, FALSE);
    }
  }

  /* evaluate density */
  for(i=0; i<nx; i++) {
    if (ISNAN(x[i])) {
      res[i] = x[i];
    }
    else if (!R_FINITE(x[i]) || x[i] <= 0.) {
      res[i] =  (logvalue) ? R_NegInf : 0.;
    }
    else {
      res[i] = LOGNORMCONSTANT + ((lambda-1.)*log(x[i]) - 0.5*(chi/x[i]+psi*x[i]));
      if (!logvalue) res[i] = exp(res[i]);
    }
  }

  /* return result to R */
  UNPROTECT(1);
  return sexp_res;
  
} /* end of dgig() */

/*---------------------------------------------------------------------------*/

#define ZTOL (DBL_EPSILON*10.0)

SEXP rgig(SEXP sexp_n, SEXP sexp_lambda, SEXP sexp_chi, SEXP sexp_psi)
/*---------------------------------------------------------------------------*/
/* Draw sample from GIG distribution.                                        */
/* Wrapper for do_rgig() with GetRNGstate() ... PutRNGstate()                */
/*                                                                           */
/* Parameters:                                                               */
/*   n ....... sample size (positive integer)                                */
/*   lambda .. parameter for distribution                                    */
/*   chi   ... parameter for distribution                                    */
/*   psi   ... parameter for distribution                                    */
/*                                                                           */
/* Return:                                                                   */
/*   random sample of size 'n'                                               */
/*---------------------------------------------------------------------------*/
{
  int n;                   /* sample size */
  double lambda, chi, psi; /* parameters of distribution */
  SEXP sexp_res;           /* results */

  /* extract sample size */
  n = *(INTEGER (AS_INTEGER (sexp_n)));

  /* extract parameters of distribution */
  lambda = *(REAL(AS_NUMERIC(sexp_lambda)));
  chi = *(REAL(AS_NUMERIC(sexp_chi)));
  psi = *(REAL(AS_NUMERIC(sexp_psi)));

  /* Get state of R uniform PRNG */
  GetRNGstate();

  /* run generator */
  sexp_res = do_rgig(n, lambda, chi, psi);

  /* Return state of PRNG to R */
  PutRNGstate();

  /* return result to R */
  return (sexp_res);

} /* end of rgig() */

/*---------------------------------------------------------------------------*/

SEXP do_rgig(int n, double lambda, double chi, double psi)
/*---------------------------------------------------------------------------*/
/* Draw sample from GIG distribution.                                        */
/* without calling GetRNGstate() ... PutRNGstate()                           */
/*                                                                           */
/* Parameters:                                                               */
/*   n ....... sample size (positive integer)                                */
/*   lambda .. parameter for distribution                                    */
/*   chi   ... parameter for distribution                                    */
/*   psi   ... parameter for distribution                                    */
/*                                                                           */
/* Return:                                                                   */
/*   random sample of size 'n'                                               */
/*---------------------------------------------------------------------------*/
{
  double omega, alpha;     /* parameters of standard distribution */
  SEXP sexp_res;           /* results */
  double *res;
  int i;

  /* check sample size */
  if (n<=0) {
    error("sample size 'n' must be positive integer.");
  }

  /* check GIG parameters: */
  if ( !(R_FINITE(lambda) && R_FINITE(chi) && R_FINITE(psi)) ||
       (chi <  0. || psi < 0)      || 
       (chi == 0. && lambda <= 0.) ||
       (psi == 0. && lambda >= 0.) ) {
    error("invalid parameters for GIG distribution: lambda=%g, chi=%g, psi=%g",
	  lambda, chi, psi);
  }

  /* allocate array for random sample */
  PROTECT(sexp_res = NEW_NUMERIC(n));
  res = REAL(sexp_res);

  if (chi < ZTOL) { 
    /* special cases which are basically Gamma and Inverse Gamma distribution */
    if (lambda > 0.0) {
      for (i=0; i<n; i++) res[i] = rgamma(lambda, 2.0/psi); 
    }
    else {
      for (i=0; i<n; i++) res[i] = 1.0/rgamma(-lambda, 2.0/psi); 
    }    
  }

  else if (psi < ZTOL) {
    /* special cases which are basically Gamma and Inverse Gamma distribution */
    if (lambda > 0.0) {
      for (i=0; i<n; i++) res[i] = 1.0/rgamma(lambda, 2.0/chi); 
    }
    else {
      for (i=0; i<n; i++) res[i] = rgamma(-lambda, 2.0/chi); 
    }    

  }

  else {
    double lambda_old = lambda;
    if (lambda < 0.) lambda = -lambda;
    alpha = sqrt(chi/psi);
    omega = sqrt(psi*chi);

    /* run generator */
    do {
      if (lambda > 2. || omega > 3.) {
        /* Ratio-of-uniforms with shift by 'mode', alternative implementation */
        _rgig_ROU_shift_alt(res, n, lambda, lambda_old, omega, alpha);
        break;
      }

      if (lambda >= 1.-2.25*omega*omega || omega > 0.2) {
        /* Ratio-of-uniforms without shift */
        _rgig_ROU_noshift(res, n, lambda, lambda_old, omega, alpha);
        break;
      }

      if (lambda >= 0. && omega > 0.) {
        /* New approach, constant hat in log-concave part. */
        _rgig_newapproach1(res, n, lambda, lambda_old, omega, alpha);
        break;
      }
      
      /* else */
      error("parameters must satisfy lambda>=0 and omega>0.");
      
    } while (0);
  }

  /* return result */
  UNPROTECT(1);
  return sexp_res;

} /* end of do_rgig() */

//[[Rcpp::export]]
double do_rgig1(double lambda, double chi, double psi)
{
  double omega, alpha;     /* parameters of standard distribution */
  SEXP sexp_res;           /* results */
  double *res;
  int i;

  double n = 1;

  /* check GIG parameters: */
  if ( !(R_FINITE(lambda) && R_FINITE(chi) && R_FINITE(psi)) ||
       (chi <  0. || psi < 0)      || 
       (chi == 0. && lambda <= 0.) ||
       (psi == 0. && lambda >= 0.) ) {
    error("invalid parameters for GIG distribution: lambda=%g, chi=%g, psi=%g",
    lambda, chi, psi);
  }

  /* allocate array for random sample */
  PROTECT(sexp_res = NEW_NUMERIC(n));
  res = REAL(sexp_res);

  if (chi < ZTOL) { 
    /* special cases which are basically Gamma and Inverse Gamma distribution */
    if (lambda > 0.0) {
      for (i=0; i<n; i++) res[i] = rgamma(lambda, 2.0/psi); 
    }
    else {
      for (i=0; i<n; i++) res[i] = 1.0/rgamma(-lambda, 2.0/psi); 
    }    
  }

  else if (psi < ZTOL) {
    /* special cases which are basically Gamma and Inverse Gamma distribution */
    if (lambda > 0.0) {
      for (i=0; i<n; i++) res[i] = 1.0/rgamma(lambda, 2.0/chi); 
    }
    else {
      for (i=0; i<n; i++) res[i] = rgamma(-lambda, 2.0/chi); 
    }    

  }

  else {
    double lambda_old = lambda;
    if (lambda < 0.) lambda = -lambda;
    alpha = sqrt(chi/psi);
    omega = sqrt(psi*chi);

    /* run generator */
    do {
      if (lambda > 2. || omega > 3.) {
        /* Ratio-of-uniforms with shift by 'mode', alternative implementation */
        _rgig_ROU_shift_alt(res, n, lambda, lambda_old, omega, alpha);
        break;
      }

      if (lambda >= 1.-2.25*omega*omega || omega > 0.2) {
        /* Ratio-of-uniforms without shift */
        _rgig_ROU_noshift(res, n, lambda, lambda_old, omega, alpha);
        break;
      }

      if (lambda >= 0. && omega > 0.) {
        /* New approach, constant hat in log-concave part. */
        _rgig_newapproach1(res, n, lambda, lambda_old, omega, alpha);
        break;
      }
      
      /* else */
      error("parameters must satisfy lambda>=0 and omega>0.");
      
    } while (0);
  }

  /* return result */
  UNPROTECT(1);
  return res[0];

} /* end of do_rgig() */

/*****************************************************************************/
/* Privat Functions                                                          */
/*****************************************************************************/

double _gig_mode(double lambda, double omega)
/*---------------------------------------------------------------------------*/
/* Compute mode of GIG distribution.                                         */
/*                                                                           */
/* Parameters:                                                               */
/*   lambda .. parameter for distribution                                    */
/*   omega ... parameter for distribution                                    */
/*                                                                           */
/* Return:                                                                   */
/*   mode                                                                    */
/*---------------------------------------------------------------------------*/
{
  if (lambda >= 1.)
    /* mode of fgig(x) */
    return (sqrt((lambda-1.)*(lambda-1.) + omega*omega)+(lambda-1.))/omega;
  else
    /* 0 <= lambda < 1: use mode of f(1/x) */
    return omega / (sqrt((1.-lambda)*(1.-lambda) + omega*omega)+(1.-lambda));
} /* end of _gig_mode() */

/*---------------------------------------------------------------------------*/

void _rgig_ROU_noshift (double *res, int n, double lambda, double lambda_old, double omega, double alpha)
/*---------------------------------------------------------------------------*/
/* Tpye 1:                                                                   */
/* Ratio-of-uniforms without shift.                                          */
/*   Dagpunar (1988), Sect.~4.6.2                                            */
/*   Lehner (1989)                                                           */
/*---------------------------------------------------------------------------*/
{
  double xm, nc;     /* location of mode; c=log(f(xm)) normalization constant */
  double ym, um;     /* location of maximum of x*sqrt(f(x)); umax of MBR */
  double s, t;       /* auxiliary variables */
  double U, V, X;    /* random variables */

  int i;             /* loop variable (number of generated random variables) */
  int count = 0;     /* counter for total number of iterations */

  /* -- Setup -------------------------------------------------------------- */

  /* shortcuts */
  t = 0.5 * (lambda-1.);
  s = 0.25 * omega;
  
  /* mode = location of maximum of sqrt(f(x)) */
  xm = _gig_mode(lambda, omega);

  /* normalization constant: c = log(sqrt(f(xm))) */
  nc = t*log(xm) - s*(xm + 1./xm);

  /* location of maximum of x*sqrt(f(x)):           */
  /* we need the positive root of                   */
  /*    omega/2*y^2 - (lambda+1)*y - omega/2 = 0    */
  ym = ((lambda+1.) + sqrt((lambda+1.)*(lambda+1.) + omega*omega))/omega;

  /* boundaries of minmal bounding rectangle:                   */
  /* we us the "normalized" density f(x) / f(xm). hence         */
  /* upper boundary: vmax = 1.                                  */
  /* left hand boundary: umin = 0.                              */
  /* right hand boundary: umax = ym * sqrt(f(ym)) / sqrt(f(xm)) */
  um = exp(0.5*(lambda+1.)*log(ym) - s*(ym + 1./ym) - nc);

  /* -- Generate sample ---------------------------------------------------- */

  for (i=0; i<n; i++) {
    do {
      ++count;
      U = um * unif_rand();        /* U(0,umax) */
      V = unif_rand();             /* U(0,vmax) */
      X = U/V;
    }                              /* Acceptance/Rejection */
    while (((log(V)) > (t*log(X) - s*(X + 1./X) - nc)));

    /* store random point */
    res[i] = (lambda_old < 0.) ? (alpha / X) : (alpha * X);
  }

  /* -- End ---------------------------------------------------------------- */

  return;
} /* end of _rgig_ROU_noshift() */


/*---------------------------------------------------------------------------*/

void _rgig_newapproach1 (double *res, int n, double lambda, double lambda_old, double omega, double alpha)
/*---------------------------------------------------------------------------*/
/* Type 4:                                                                   */
/* New approach, constant hat in log-concave part.                           */
/* Draw sample from GIG distribution.                                        */
/*                                                                           */
/* Case: 0 < lambda < 1, 0 < omega < 1                                       */
/*                                                                           */
/* Parameters:                                                               */
/*   n ....... sample size (positive integer)                                */
/*   lambda .. parameter for distribution                                    */
/*   omega ... parameter for distribution                                    */
/*                                                                           */
/* Return:                                                                   */
/*   random sample of size 'n'                                               */
/*---------------------------------------------------------------------------*/
{
  /* parameters for hat function */
  double A[3], Atot;  /* area below hat */
  double k0;          /* maximum of PDF */
  double k1, k2;      /* multiplicative constant */

  double xm;          /* location of mode */
  double x0;          /* splitting point T-concave / T-convex */
  double a;           /* auxiliary variable */

  double U, V, X;     /* random numbers */
  double hx;          /* hat at X */

  int i;              /* loop variable (number of generated random variables) */
  int count = 0;      /* counter for total number of iterations */

  /* -- Check arguments ---------------------------------------------------- */

  if (lambda >= 1. || omega >1.)
    error ("invalid parameters");

  /* -- Setup -------------------------------------------------------------- */

  /* mode = location of maximum of sqrt(f(x)) */
  xm = _gig_mode(lambda, omega);

  /* splitting point */
  x0 = omega/(1.-lambda);

  /* domain [0, x_0] */
  k0 = exp((lambda-1.)*log(xm) - 0.5*omega*(xm + 1./xm));     /* = f(xm) */
  A[0] = k0 * x0;

  /* domain [x_0, Infinity] */
  if (x0 >= 2./omega) {
    k1 = 0.;
    A[1] = 0.;
    k2 = pow(x0, lambda-1.);
    A[2] = k2 * 2. * exp(-omega*x0/2.)/omega;
  }
  
  else {
    /* domain [x_0, 2/omega] */
    k1 = exp(-omega);
    A[1] = (lambda == 0.) 
      ? k1 * log(2./(omega*omega))
      : k1 / lambda * ( pow(2./omega, lambda) - pow(x0, lambda) );

    /* domain [2/omega, Infinity] */
    k2 = pow(2/omega, lambda-1.);
    A[2] = k2 * 2 * exp(-1.)/omega;
  }

  /* total area */
  Atot = A[0] + A[1] + A[2];

  /* -- Generate sample ---------------------------------------------------- */

  for (i=0; i<n; i++) {
    do {
      ++count;

      /* get uniform random number */
      V = Atot * unif_rand();
      
      do {
	
	/* domain [0, x_0] */
	if (V <= A[0]) {
	  X = x0 * V / A[0];
	  hx = k0;
	  break;
	}
	
	/* domain [x_0, 2/omega] */
	V -= A[0];
	if (V <= A[1]) {
	  if (lambda == 0.) {
	    X = omega * exp(exp(omega)*V);
	    hx = k1 / X;
	  }
	  else {
	    X = pow(pow(x0, lambda) + (lambda / k1 * V), 1./lambda);
	    hx = k1 * pow(X, lambda-1.);
	  }
	  break;
	}
	
	/* domain [max(x0,2/omega), Infinity] */
	V -= A[1];
	a = (x0 > 2./omega) ? x0 : 2./omega;
	X = -2./omega * log(exp(-omega/2. * a) - omega/(2.*k2) * V);
	hx = k2 * exp(-omega/2. * X);
	break;
	
      } while(0);
      
      /* accept or reject */
      U = unif_rand() * hx;

      if (log(U) <= (lambda-1.) * log(X) - omega/2. * (X+1./X)) {
	/* store random point */
	res[i] = (lambda_old < 0.) ? (alpha / X) : (alpha * X);
	break;
      }
    } while(1);
    
  }

  /* -- End ---------------------------------------------------------------- */

  return;
} /* end of _rgig_newapproach1() */

/*---------------------------------------------------------------------------*/

void
_rgig_ROU_shift_alt (double *res, int n, double lambda, double lambda_old, double omega, double alpha)
/*---------------------------------------------------------------------------*/
/* Type 8:                                                                   */
/* Ratio-of-uniforms with shift by 'mode', alternative implementation.       */
/*   Dagpunar (1989)                                                         */
/*   Lehner (1989)                                                           */
/*---------------------------------------------------------------------------*/
{
  double xm, nc;     /* location of mode; c=log(f(xm)) normalization constant */
  double s, t;       /* auxiliary variables */
  double U, V, X;    /* random variables */

  int i;             /* loop variable (number of generated random variables) */
  int count = 0;     /* counter for total number of iterations */

  double a, b, c;    /* coefficent of cubic */
  double p, q;       /* coefficents of depressed cubic */
  double fi, fak;    /* auxiliary results for Cardano's rule */

  double y1, y2;     /* roots of (1/x)*sqrt(f((1/x)+m)) */

  double uplus, uminus;  /* maximum and minimum of x*sqrt(f(x+m)) */

  /* -- Setup -------------------------------------------------------------- */

  /* shortcuts */
  t = 0.5 * (lambda-1.);
  s = 0.25 * omega;

  /* mode = location of maximum of sqrt(f(x)) */
  xm = _gig_mode(lambda, omega);

  /* normalization constant: c = log(sqrt(f(xm))) */
  nc = t*log(xm) - s*(xm + 1./xm);

  /* location of minimum and maximum of (1/x)*sqrt(f(1/x+m)):  */

  /* compute coeffients of cubic equation y^3+a*y^2+b*y+c=0 */
  a = -(2.*(lambda+1.)/omega + xm);       /* < 0 */
  b = (2.*(lambda-1.)*xm/omega - 1.);
  c = xm;

  /* we need the roots in (0,xm) and (xm,inf) */

  /* substitute y=z-a/3 for depressed cubic equation z^3+p*z+q=0 */
  p = b - a*a/3.;
  q = (2.*a*a*a)/27. - (a*b)/3. + c;

  /* use Cardano's rule */
  fi = acos(-q/(2.*sqrt(-(p*p*p)/27.)));
  fak = 2.*sqrt(-p/3.);
  y1 = fak * cos(fi/3.) - a/3.;
  y2 = fak * cos(fi/3. + 4./3.*M_PI) - a/3.;

  /* boundaries of minmal bounding rectangle:                  */
  /* we us the "normalized" density f(x) / f(xm). hence        */
  /* upper boundary: vmax = 1.                                 */
  /* left hand boundary: uminus = (y2-xm) * sqrt(f(y2)) / sqrt(f(xm)) */
  /* right hand boundary: uplus = (y1-xm) * sqrt(f(y1)) / sqrt(f(xm)) */
  uplus  = (y1-xm) * exp(t*log(y1) - s*(y1 + 1./y1) - nc);
  uminus = (y2-xm) * exp(t*log(y2) - s*(y2 + 1./y2) - nc);

  /* -- Generate sample ---------------------------------------------------- */

  for (i=0; i<n; i++) {
    do {
      ++count;
      U = uminus + unif_rand() * (uplus - uminus);    /* U(u-,u+)  */
      V = unif_rand();                                /* U(0,vmax) */
      X = U/V + xm;
    }                                         /* Acceptance/Rejection */
    while ((X <= 0.) || ((log(V)) > (t*log(X) - s*(X + 1./X) - nc)));

    /* store random point */
    res[i] = (lambda_old < 0.) ? (alpha / X) : (alpha * X);
  }

  /* -- End ---------------------------------------------------------------- */

  return;
} /* end of _rgig_ROU_shift_alt() */

/*---------------------------------------------------------------------------*/

double
_unur_bessel_k_nuasympt (double x, double nu, int islog, int expon_scaled)
/*---------------------------------------------------------------------------*/
/* Asymptotic expansion of Bessel K_nu(x) function                           */
/* when BOTH  nu and x  are large.                                           */
/*                                                                           */
/* parameters:                                                               */
/*   x            ... argument for K_nu()                                    */
/*   nu           ... order or Bessel function                               */
/*   islog        ... return logarithm of result TRUE and result when FALSE  */
/*   expon_scaled ... return exp(-x)*K_nu(x) when TRUE and K_nu(x) when FALSE*/
/*                                                                           */
/*---------------------------------------------------------------------------*/
/*                                                                           */
/* references:                                                               */
/* ##  Abramowitz & Stegun , p.378, __ 9.7.8. __                             */
/*                                                                           */
/* ## K_nu(nu * z) ~ sqrt(pi/(2*nu)) * exp(-nu*eta)/(1+z^2)^(1/4)            */
/* ##                   * {1 - u_1(t)/nu + u_2(t)/nu^2 - ... }               */
/*                                                                           */
/* ## where   t = 1 / sqrt(1 + z^2),                                         */
/* ##       eta = sqrt(1 + z^2) + log(z / (1 + sqrt(1+z^2)))                 */
/* ##                                                                        */
/* ## and u_k(t)  from  p.366  __ 9.3.9 __                                   */
/*                                                                           */
/* ## u0(t) = 1                                                              */
/* ## u1(t) = (3*t - 5*t^3)/24                                               */
/* ## u2(t) = (81*t^2 - 462*t^4 + 385*t^6)/1152                              */
/* ## ...                                                                    */
/*                                                                           */
/* ## with recursion  9.3.10    for  k = 0, 1, .... :                        */
/* ##                                                                        */
/* ## u_{k+1}(t) = t^2/2 * (1 - t^2) * u'_k(t) +                             */
/* ##            1/8  \int_0^t (1 - 5*s^2)* u_k(s) ds                        */
/*---------------------------------------------------------------------------*/
/*                                                                           */
/* Original implementation in R code (R package "Bessel" v. 0.5-3) by        */
/*   Martin Maechler, Date: 23 Nov 2009, 13:39                               */
/*                                                                           */
/* Translated into C code by Kemal Dingic, Oct. 2011.                        */
/*                                                                           */
/* Modified by Josef Leydold on Tue Nov  1 13:22:09 CET 2011                 */
/*                                                                           */
/*---------------------------------------------------------------------------*/
{
#define M_LNPI     1.14472988584940017414342735135      /* ln(pi) */

  double z;                   /* rescaled argument for K_nu() */
  double sz, t, t2, eta;      /* auxiliary variables */
  double d, u1t,u2t,u3t,u4t;  /* (auxiliary) results for Debye polynomials */
  double res;                 /* value of log(K_nu(x)) [= result] */
  
  /* rescale: we comute K_nu(z * nu) */
  z = x / nu;

  /* auxiliary variables */
  sz = hypot(1,z);   /* = sqrt(1+z^2) */
  t = 1. / sz;
  t2 = t*t;

  eta = (expon_scaled) ? (1./(z + sz)) : sz;
  eta += log(z) - log1p(sz);                  /* = log(z/(1+sz)) */

  /* evaluate Debye polynomials u_j(t) */
  u1t = (t * (3. - 5.*t2))/24.;
  u2t = t2 * (81. + t2*(-462. + t2 * 385.))/1152.;
  u3t = t*t2 * (30375. + t2 * (-369603. + t2 * (765765. - t2 * 425425.)))/414720.;
  u4t = t2*t2 * (4465125. 
		 + t2 * (-94121676.
			 + t2 * (349922430. 
				 + t2 * (-446185740. 
					 + t2 * 185910725.)))) / 39813120.;
  d = (-u1t + (u2t + (-u3t + u4t/nu)/nu)/nu)/nu;

  /* log(K_nu(x)) */
  res = log(1.+d) - nu*eta - 0.5*(log(2.*nu*sz) - M_LNPI);

  return (islog ? res : exp(res));
} /* end of _unur_bessel_k_nuasympt() */

/*---------------------------------------------------------------------------*/


//lambda = p, psi=a, chi=b
// [[Rcpp::export]]
double gig_norm(double lambda, double chi, double psi) {
  double LOGNORMCONSTANT;
  /* compute normalization constant */
  if (psi == 0.) {
    /* case: Inverse Gamma distribution */
    LOGNORMCONSTANT = -lambda * log(0.5*chi) - lgammafn(-lambda);
  }
  else if (chi == 0.) {
    /* case: Gamma distribution */
    LOGNORMCONSTANT = lambda * log(0.5*psi) - lgammafn(lambda);
  }
  else {
    /* general case: */
    /*   pow(psi/chi, lambda/2.) / (2. * bessel_k(sqrt(psi*chi),lambda,1)); */
    double alambda = fabs(lambda);
    double beta = sqrt(psi*chi);
    LOGNORMCONSTANT = 0.5*lambda*log(psi/chi) - M_LN2;
    if (alambda < 50.) {
      /* threshold value 50 is selected by experiments */
      LOGNORMCONSTANT -= log(bessel_k(beta, alambda, 2)) - beta;
    }
    else {
      LOGNORMCONSTANT -= _unur_bessel_k_nuasympt(beta, alambda, TRUE, FALSE);
    }
  }
  return -LOGNORMCONSTANT;
}



