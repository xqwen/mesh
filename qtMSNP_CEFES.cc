#include "qtMSNP_CEFES.h"
#include <math.h>
#include <gsl/gsl_integration.h>


double cefes_target_f (double x, void *params){
  
  double *p = (double *)params;
  double val = 1;
  int size = int(p[0]);
  double k = p[1];
  double oa2 = p[2];


  for(int i=0;i<size;i++){
    double ds2 = p[4*i+3];
    double t2 = p[4*i+4];
    double sn = p[4*i+5];
    double ss = p[4*i+6];
    

    double dnum = ds2 + k*x*x;
    double fac = ds2/dnum;
    val *= sqrt(fac)*exp(-0.5*fac*t2-0.5*x*x/dnum+(sn/dnum)*x);
    
    // correction term
    //double iss = 1.0/(ss-2);
    //double cor = (1 + iss/6.0 + .25*iss*pow( fac*t2-(sn/dnum)*x, 2)-iss*fac*t2+.75*iss*sn*x/dnum)/(1+iss/6.0);
    //val *= cor;
    
    //printf("correction = %f\n",cor);

    //printf ("val = %e\n",val);
  }
  val *= exp(-.5*x*x/oa2);
  
  return val;
}
 

void qtMSNP_CEFES::prepare_params(double k, double oa2){
  
  // find non-informative sets
  
  info_size = 0;
  for(int i=0;i<ssize;i++){
    if(sList[i].vg*indicator[i] != 0)
      info_size++;
  }
  
  
  if(info_size==0)
    return;

  param_list = new double[4*info_size+3];
  
  param_list[0] = info_size;
  param_list[1] = k;
  param_list[2] = oa2;
 
  int i=0;
  for(int j=0;j<ssize;j++){
    
    if(sList[j].vg*indicator[j] == 0)
      continue;
    
    param_list[4*i+3] = sList[j].ds2_es;
    param_list[4*i+4] = sList[j].chi2;
    param_list[4*i+5] = sList[j].bhat/sList[j].sigmahat;
    param_list[4*i+6] = sList[j].n;

    i++;
  }
  
}


double qtMSNP_CEFES::compute_log10_BF(double k, double oa2){
  
  
  // parameter setting
  prepare_params(k,oa2);
  //printf("info_size = %d\n",info_size);
  
  
  
  if(info_size==0)
    return 0.0;
  
  // numerical integration
  gsl_integration_workspace *w = gsl_integration_workspace_alloc (5000);
  double result, error;
  gsl_function F;
  F.function = &cefes_target_f;
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


double qtMSNP_CEFES::compute_log10_ABF(double k, double oa2){
  return compute_log10_BF(k,oa2);
}
