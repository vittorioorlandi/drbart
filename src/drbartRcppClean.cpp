#include <Rcpp.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <ctime>
#include <algorithm> 
#include <unistd.h>
#include <chrono>
#include <thread>
// #include <Rinterface.h>

#include "rng.h"
#include "tree.h"
#include "info.h"
#include "funs.h"
#include "bd.h"
#include "slice.h"


using namespace Rcpp;

// [[Rcpp::export]]
List drbartRcppClean(NumericVector y_, NumericVector x_, NumericVector xpred_, 
              List xinfo_list,
              int burn, int nd, int thin, int printevery, int m, 
              double alpha, double beta, double lambda, double nu, double kfac,
              IntegerVector trunc_below,
              CharacterVector treef_name_)
{
  
  std::string treef_name = as<std::string>(treef_name_); 
  std::ofstream treef(treef_name.c_str());
  
  RNGScope scope;  
  RNG gen; //this one random number generator is used in all draws
  
  /*****************************************************************************
   Read, format y
  *****************************************************************************/
  std::vector<double> y; //storage for y
  double miny = INFINITY, maxy = -INFINITY;
  sinfo allys;       //sufficient stats for all of y, use to initialize the bart trees.
  
  for (NumericVector::iterator it = y_.begin(); it != y_.end(); ++it) {
    y.push_back(*it);

    if (*it<miny) {
      miny=*it;
    }
    if (*it>maxy) {
      maxy=*it;
    }

    allys.sy += *it; // sum of y
    allys.sy2 += (*it) * (*it); // sum of y^2
  }
  size_t n = y.size();
  allys.n = n;
  
  double ybar = allys.sy/n; //sample mean
  double shat = sqrt((allys.sy2 - n * ybar * ybar) / (n - 1)); //sample standard deviation
  
  /*****************************************************************************
   Read, format X, Xpred
  *****************************************************************************/
  //read x   
  //the n*p numbers for x are stored as the p for first obs, then p for second, and so on.
  std::vector<double> x;
  for(NumericVector::iterator it=x_.begin(); it != x_.end(); ++it) {
    x.push_back(*it);
  }
  size_t p = x.size()/n;
  
  //x for predictions
  dinfo dip; //data information for prediction
  dip.n = 0;
  std::vector<double> xp;     //stored like x
  if (xpred_.size()) {
    for (NumericVector::iterator it = xpred_.begin(); it != xpred_.end(); ++it) {
       xp.push_back(*it);
    }
    size_t np = xp.size() / p;
    if (xp.size() != np * p) {
      Rcpp::stop("Wrong number of elements in prediction data set.");
    }
    if (np) {
      dip.n = np; 
      dip.p = p; 
      dip.x = &xp[0]; 
      dip.y = 0; //there are no y's!
    }
  }
  
  //x cutpoints
  xinfo xi;
  
  xi.resize(p);
  for (int i = 0; i < p; ++i) {
    NumericVector tmp = xinfo_list[i];
    std::vector<double> tmp2;
    for (size_t j = 0; j < tmp.size(); ++j) {
      tmp2.push_back(tmp[j]);
    }
    xi[i] = tmp2;
  }
  
  /*****************************************************************************
   Setup for MCMC
  *****************************************************************************/

  //trees
  std::vector<tree> t(m);
  for (size_t i = 0; i < m; i++) {
    t[i].setm(ybar / m); //if you sum the fit over the trees you get the fit.
  }
  
  //--------------------------------------------------
  //prior and mcmc
  pinfo pi;
  pi.pbd = 1.0; //prob of birth/death move
  pi.pb = .5; //prob of birth given  birth/death
  
  pi.alpha = alpha; //prior prob a bot node splits is alpha/(1+d)^beta, d is depth of node
  pi.beta = beta; //2 for bart means it is harder to build big trees.
  pi.tau = (maxy - miny) / (2 * kfac * sqrt((double)m)); //sigma_mu
  pi.sigma = shat;
    
  //--------------------------------------------------
  //dinfo
  double* allfit = new double[n]; //sum of fit of all trees
  for (size_t i = 0;i < n; i++) {
    allfit[i] = ybar;
  }
  double* r = new double[n]; //y-(allfit-ftemp) = y-allfit+ftemp
  double* ftemp = new double[n]; //fit of current tree
  dinfo di;
  di.n = n; di.p = p; di.x = &x[0]; di.y = r; //the y for each draw will be the residual 
  
  //--------------------------------------------------
  //storage for ouput
  
  //for sigma draw
  double rss, restemp;
  
  //save stuff to tree file
  treef << xi << endl; //cutpoints
  treef << m << endl;  //number of trees
  treef << p << endl;  //dimension of x's
  treef << nd << endl;

  NumericVector ssigma(nd);

  int niters = nd * thin + burn; 
  
  /*****************************************************************************
   MCMC
  *****************************************************************************/
  //begin dr bart
  
  tree::npv bnv;
  std::vector<tree::npv> bnvs;
  std::vector<std::vector<int> > leaf_counts(m);
  std::vector<double> lik(xi[0].size());
  std::vector<tree> using_u;
  std::vector<std::vector<double> > ucuts_post(nd);
  
  ld_bartU slice_density(0.0, 1.0);
  slice_density.xi = xi;
  slice_density.di = di;
  slice_density.i = 0;
  slice_density.using_u = using_u;
  //end dr bart
  
  for (size_t i = 0; i < niters; i++) {
    if (i % printevery == 0) {
      Rprintf("\r");
      Rprintf("Iteration %d / %d (%d%%)", i, niters, (int) 100 * i / niters);
      Rprintf("\r");
    }
    //draw trees
    for (size_t j = 0; j < m; j++) {
      fit(t[j], xi, di, ftemp);
      for (size_t k = 0; k < n; k++) {
        if (ftemp[k] != ftemp[k]) {
          Rcout << "tree " << j <<" obs "<< k<<" "<< endl;
          Rcout << t[j] << endl;
          stop("nan in ftemp");
        }
        allfit[k] = allfit[k] - ftemp[k];
        r[k] = y[k] - allfit[k];
      }

      bd(t[j], xi, di, pi, gen);
      drmu(t[j], xi, di, pi, gen);
      fit(t[j], xi, di, ftemp);

      for (size_t k = 0; k < n; k++) {
        allfit[k] += ftemp[k];
      }
    }
    //begin dr bart
    
    // impute censored values
    for (size_t k = 0; k < n; ++k) {
      if (trunc_below[k] > 0) {
        y[k] = rtnormlo(allfit[k], pi.sigma, y_[k]);// original y_ is obs value
      }
    }
    
    /*** sample u ***/
    using_u.clear();
    leaf_counts.clear();
    bnvs.clear();
    vector<std::map<tree::tree_cp,size_t> > bnmaps;
    std::set<size_t> ucuts; ucuts.insert(0); ucuts.insert(xi[0].size() - 1);
    
    //get trees splitting on u, the first variable
    for (int tt = 0; tt < m; ++tt) {
      if (t[tt].nuse(0)) {
        using_u.push_back(t[tt]);
      }
    }
    
    //update slice_density object
    slice_density.sigma = pi.sigma;
    slice_density.using_u = using_u;
    
    //get leaf counts for each tree splitting on u & also get partition of u
    for (size_t tt = 0; tt < using_u.size(); ++tt) {
      leaf_counts.push_back(counts(using_u[tt], xi, di, bnv)); //clears & populates bnv
      bnvs.push_back(bnv);
      using_u[tt].varsplits(ucuts, 0);
    }
    
    std::vector<size_t> ucutsv(ucuts.begin(), ucuts.end());
    std::vector<double> logpr(ucuts.size() - 1);

    //prebuild ix->bottom node maps for each tree splitting on u, big time saver.
    typedef tree::npv::size_type bvsz;
    for (size_t tt = 0; tt < using_u.size(); ++tt) {
      std::map<tree::tree_cp,size_t> bnmap;

      for (bvsz ii = 0; ii != bnvs[tt].size(); ii++) {
        bnmap[bnvs[tt][ii]] = ii; 
      }

      bnmaps.push_back(bnmap);
    }

    //loop over each observation
    std::vector<int> tmpcounts;
    std::vector<std::vector<int> > new_counts(using_u.size());
    size_t jj = 0; //<-- single latent variable for now.
    tree::tree_cp tbn;
    for (size_t k = 0; k < n; k++) {
      bool proceed = true;
//      
     int L, U;
     L = 0; U = xi[jj].size() - 1 ;
     
     //check that removing u won't result in bottom nodes 
     //todo: sample u uniformly from current partition? does that help?
     for (size_t tt = 0; tt < using_u.size(); ++tt) {
       tmpcounts = leaf_counts[tt];
       update_counts(k, tmpcounts, using_u[tt], xi, di, bnmaps[tt], -1, tbn); 
       new_counts[tt] = tmpcounts;
       if (*std::min_element(tmpcounts.begin(), tmpcounts.end()) < 5) {
         proceed = false;
         //Rcout << "skipping u because proceed is " << proceed << endl;
         break;
         
         //get the region where u has to live to avoid jumping bottom nodes
         // int L0, U0;
         // L0=0; U0 = xi[jj].size()-1;
         // tbn->rg(0, &L0, &U0);
         // L = std::max(L, L0);
         // U = std::min(U, U0);
       }
     }
      
      //resample u
      if (proceed) {
        int L,U;
        // L = 0; U = xi[jj].size() - 1;
        // //take out u, decrement counts - should use tmpcounts from above, right?? computing 2x?
        // for (size_t tt = 0; tt < using_u.size(); ++tt) {
        //   update_counts(k, leaf_counts[tt], using_u[tt], xi, di, bnmaps[tt], -1);
        //   if(*std::min_element(leaf_counts[tt].begin(), leaf_counts[tt].end()) < 5) {
        //     //proceed = false;
        //     //break;
            
        //     //get the region where u has to live to avoid jumping bottom nodes
        //     int L0, U0;
        //     L0 = 0; U0 = xi[jj].size() - 1;
        //     tbn->rg(0, &L0, &U0);
        //     L = std::max(L, L0);
        //     U = std::min(U, U0);  
        //   }
        // }

        leaf_counts = new_counts;
        
        double f = allfit[k] - fit_i(k, using_u, xi, di); //fit from trees that don't use u
        
        slice_density.i = k;
        slice_density.f = f;
        slice_density.yobs = y[k];
        double oldu = x[jj+k*p];
        double lo = xi[jj][L]; double hi = xi[jj][U];
        double newu = slice(oldu, &slice_density, 1.0, INFINITY, 0., 1.);
        x[jj + k * p] = newu;
        
        //update counts with new u
        for (size_t tt = 0; tt < using_u.size(); ++tt) {
          update_counts(k, leaf_counts[tt], using_u[tt], xi, di, bnmaps[tt], 1);
        }
        // add back the fit from trees splitting on u
        allfit[k] = f + fit_i(k, using_u, xi, di); //should save these in previous for loop?
      }
    }
    //end dr bart
    
    //draw sigma
    rss = 0.0;
    for (size_t k = 0; k < n; k++) {
      restemp = y[k] - allfit[k]; 
      rss += restemp * restemp;
    }

    pi.sigma = sqrt((nu * lambda + rss) / gen.chi_square(nu + n));

    if (i >= burn & i % thin == 0) {
      for (size_t uu = 0; uu < ucutsv.size(); ++uu) {
        ucuts_post[(i - burn) / thin].push_back(xi[jj][ucutsv[uu]]);
      }
      
      for (size_t j = 0; j < m; j++) {
        treef << t[j] << endl;
      }
      
      ssigma((i - burn) / thin) = pi.sigma;
    } 
  }

  t.clear();
  delete[] allfit;
  delete[] r;
  delete[] ftemp;
  
  treef.close();
  return(List::create(_["sigma"] = ssigma,
                      _["ucuts"] = ucuts_post));
}
