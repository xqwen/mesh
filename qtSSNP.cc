#include "qtSSNP.h"
#include <math.h>
#include <gsl/gsl_cdf.h>
#include <stdio.h>

void qtSSNP::load_data(int n_, double g_mean_, double gtg_, double y_mean_, double yty_, double gty_){
  
  n = n_;
  g_mean = g_mean_;
  gtg = gtg_;
  y_mean = y_mean_;
  yty = yty_;
  gty = gty_;


  rss0 = yty - n*y_mean*y_mean;
  vg = gtg - n*g_mean*g_mean;
  if(vg>1e-8){
    bhat = (gty-n*g_mean*y_mean)/vg;
    rss1 = rss0 - vg*bhat*bhat;
    sigmahat = sqrt(rss1/(n-2));
    sd_beta_orig = sd_beta = sigmahat/sqrt(vg);
    if(adjust)
      small_sample_correction();
    else
      chi2 = pow(bhat/sd_beta,2);


    ds2_ee = pow(sd_beta,2);
    ds2_es = 1/vg;

  }else{  
    vg = 0;
    rss1 = rss0;
    bhat = 0;
    sigmahat = sqrt(rss1/(n-2));
    ds2_ee = ds2_es = 0;
  }

}

void qtSSNP::load_data(int n_, double bhat_, double sd_beta_, double vx_){
  
  n = n_;
  g_mean = y_mean = 0;
  bhat = bhat_;
  sd_beta_orig=sd_beta = sd_beta_;
  vg = n*vx_;
  
  sigmahat = sqrt((n-1)*vx_*sd_beta*sd_beta);
  rss1 = (n-2)*sigmahat*sigmahat;
  yty = rss0 = rss1 + (n-1)*vx_*bhat*bhat;
  gtg = (n-1)*vx_;
  gty = gtg*bhat;
  
  
  
  if(adjust)
    small_sample_correction();
  else
    chi2 = pow(bhat/sd_beta,2);
  
  
  ds2_ee = pow(sd_beta,2);
  ds2_es = 1/vg;
    
}
  
  
void qtSSNP::small_sample_correction(){

  chi2 = pow(gsl_cdf_gaussian_Pinv(gsl_cdf_tdist_P(-fabs(bhat/sd_beta),n-2),1.0),2);
  if(chi2!=0){
    sd_beta =fabs(bhat)/sqrt(chi2);
    sigmahat = sd_beta*sqrt(vg);
  }
  return;
  
}


void qtSSNP::compute_rssb(double a2){
  // Bayesian RSS
  rssb = (vg*a2/(vg*a2+1))*rss1 + (1.0/(vg*a2+1))*rss0;
  eta = bhat*vg/(1+vg*a2);

}


double qtSSNP::compute_ABF(double a2, int option){

  if(vg == 0)
    return 0;
  double ds2 = ds2_es;

  if(option == 1){
    ds2 = ds2_ee;
  }

  double log10_ABF1 = .5*log10(ds2/(ds2+a2)) + (.5*chi2*a2/(ds2+a2))/log(10);
  return log10_ABF1;
  
}


