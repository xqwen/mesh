#include "qtMSNP_CEFEE.h"
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_cdf.h>


double cefee_target_f (double x, void *params){
  
  double *p = (double *)params;
  double val = 1;
  int size = int(p[0]);
  double k = p[1];
  double oa2 = p[2];
  //printf("init %f %f %f\n",double(size),k,oa2);
  for(int i=0;i<size;i++){
    double ds2 = pow(p[2*i+4],2.0);
    double beta = p[2*i+3];
    double t2 = pow(beta,2)/ds2;
    
    double dnum = ds2 + k*x*x;
    double fac = ds2/dnum;
    val *= sqrt(fac)*exp(-0.5*fac*t2-0.5*x*x/dnum+(beta/dnum)*x);
    //printf ("val = %e\n",val);
  }
  val *= exp(-.5*x*x/oa2);
  
  

  return val;
}
 

void qtMSNP_CEFEE::prepare_params(double k, double oa2){
  
  // find non-informative sets
  info_size = 0;
  
  for(int i=0;i<ssize;i++){
    if(sList[i].vg*indicator[i]!=0)
      info_size++;
  }
  
  
  if(info_size==0)
    return;



  param_list = new double[2*info_size+3];
  
  param_list[0] = info_size;
  param_list[1] = k;
  param_list[2] = oa2;
  
  double val = 0;
  int i=0;
  for(int j=0;j<ssize;j++){
    
    if(sList[j].vg*indicator[j] == 0)
      continue; 
    param_list[2*i+3] = sList[j].bhat;
    param_list[2*i+4] = sList[j].sd_beta;
    i++;
  }
  
}


double qtMSNP_CEFEE::compute_log10_BF(double k, double oa2){
  
  prepare_params(k,oa2);
  
  if(info_size==0)
    return 0.0;
  
  // numerical integration
  gsl_integration_workspace *w = gsl_integration_workspace_alloc (5000);
  double result, error;
  gsl_function F;
  F.function = &cefee_target_f;
  F.params = param_list;
    
  gsl_integration_qagi (&F, 0, 1e-7, 5000,w, &result, &error); 

  gsl_integration_workspace_free (w);
  delete[] param_list;



  double logBF = log(result)-.5*log(oa2)-.5*log(2*3.1415926536);
  
  for(int i=0;i<ssize;i++){
    if(sList[i].vg*indicator[i] == 0)
      continue;
    logBF += .5*sList[i].chi2;
  }
  
  return logBF/log(10.0);
  
  
}


double qtMSNP_CEFEE::compute_log10_ABF(double k, double oa2){
  return compute_log10_BF(k,oa2);
}



