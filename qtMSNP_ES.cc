#include "qtMSNP_ES.h"
#include <math.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf.h>


double es_target_f (const gsl_vector *v, void *params){
  
  double *p = (double *)params;
  double val = 0;
  double evv = 0;
  
  int size = v->size;

  for(int i=0;i<size;i++){
    double tau = gsl_vector_get(v,i);
    val += (0.5*p[3*i]-1)*log(tau)- 0.5*tau*p[3*i+1];
    evv += sqrt(tau)*p[3*i+2];
  }
  
  val += 0.5*evv*evv/p[3*size];
  return -val;
}
 

void es_target_df (const gsl_vector *v, void *params, gsl_vector *df){
  
  double *p = (double *)params;
  vector<double> tempv1;
  vector<double> tempv2;
  double evv = 0;

  int size = v->size;

  for(int i=0;i<size;i++){
    double tau = gsl_vector_get(v,i);
    tempv1.push_back( (0.5*p[3*i]-1)/tau - .5*p[3*i+1]);
    tempv2.push_back( .5*p[3*i+2]/sqrt(tau) );
    evv += sqrt(tau)*p[3*i+2];
  }

  evv /= p[3*size];
  for(int i=0;i<size;i++)
    gsl_vector_set(df,i, -(tempv1[i]+evv*tempv2[i]) );
}
     

/* Compute both f and df together. */
void  es_target_fdf (const gsl_vector *x, void *params, double *f, gsl_vector *df){
  *f = es_target_f(x, params); 
  es_target_df(x, params, df);
} 




void qtMSNP_ES::prepare_params(double sa2, double oa2){
  
  // find non-informative sets
  
  info_size = 0;
  for(int i=0;i<ssize;i++){
    if(sList[i].vg*indicator[i] != 0)
      info_size++;
  }


  if(info_size==0)
    return;
  
  param_list = new double[3*info_size+1];
  double val = 0;
  int i=0;
  for(int j=0;j<ssize;j++){
    
    if(sList[j].vg*indicator[j] == 0)
      continue;
    
    sList[j].compute_rssb(sa2);

    param_list[3*i] = double(sList[j].n);
    param_list[3*i+1] = sList[j].rssb;
  
    double temp = sList[j].vg/(1.0+sa2*sList[j].vg);
    param_list[3*i+2] = sList[j].bhat*temp*sqrt(oa2);    
    
    val += temp;
    i++;
  }
  
  param_list[3*info_size] = oa2*val+1;  

}



double qtMSNP_ES::laplace_approx(double sa2, double oa2){

  size_t iter = 0;
  int status;
     
  const gsl_multimin_fdfminimizer_type *T;
  gsl_multimin_fdfminimizer *s;
     
  gsl_vector *x;
  gsl_multimin_function_fdf target_func;
  
  target_func.n  = info_size;
  target_func.f  = ::es_target_f;
  target_func.df = ::es_target_df;
  target_func.fdf = ::es_target_fdf;
  target_func.params = param_list;
     
  /* Starting point */
  x = gsl_vector_alloc(info_size);
  int j=0;
  for(int i=0;i<sList.size();i++){
    if(sList[i].vg*indicator[i] == 0)
      continue;
    double val = double(sList[i].n-2)/sList[i].rssb;
    gsl_vector_set (x, j++, val);
    
  }
  

  //T = gsl_multimin_fdfminimizer_conjugate_fr;
  //T = gsl_multimin_fdfminimizer_steepest_descent;
  T = gsl_multimin_fdfminimizer_vector_bfgs;
  s = gsl_multimin_fdfminimizer_alloc (T, info_size);
     
  gsl_multimin_fdfminimizer_set (s, &target_func, x, 0.01, 0.01);

  do{
    iter++;
    status = gsl_multimin_fdfminimizer_iterate (s);
    
    if (status)
      break;
  
    status = gsl_multimin_test_gradient (s->gradient, 1e-3);
    
    
    //printf ("%5d %7.5e  \n",iter,s->f);
    //for(int i=0;i<info_size;i++){
    //  double val = gsl_vector_get(s->x,i);
    //  printf(" %7.5e ",1.0/val);
    //}
    //printf("\n");
    

  }while (status == GSL_CONTINUE && iter < 5000);
  
  
  vector<double> tauv;
  for(int i=0;i<info_size;i++){
    double val = gsl_vector_get(s->x,i);  
    
    tauv.push_back(val);
  }
  
  // compute Hessian at the minimization point
  gsl_matrix *hessian = gsl_matrix_calloc(info_size,info_size);
  
  double evv = 0;
  for(int i=0;i<info_size;i++)
    evv += sqrt(tauv[i])*param_list[3*i+2];
  
  evv /= param_list[3*info_size];
  for(int i=0;i<info_size;i++){
    for(int j=i;j<info_size;j++){
      if(i==j){
	double val = -(0.5*param_list[3*i]-1)/(tauv[i]*tauv[i]) - 0.25*evv*param_list[3*i+2]*pow(tauv[i],-1.5)+0.25*param_list[3*i+2]*param_list[3*i+2]/(param_list[3*info_size]*tauv[i]);
	gsl_matrix_set(hessian,i,i,-val);
      }else{
	double val = 0.25*param_list[3*i+2]*param_list[3*j+2]/(param_list[3*info_size]*sqrt(tauv[i]*tauv[j]));
	gsl_matrix_set(hessian,i,j, -val);
	gsl_matrix_set(hessian,j,i, -val);
      }
    }
  }

    
  // LU decomp and det computation
  int ss;
  gsl_permutation * p = gsl_permutation_alloc (info_size);
     
  int flag = gsl_linalg_LU_decomp (hessian, p, &ss);
  double log_detVal = gsl_linalg_LU_lndet(hessian);
  
  double rst = 0.5*info_size*log(2*3.1415926536)-0.5*log_detVal - s->f;
  gsl_matrix_free(hessian);
  gsl_permutation_free (p);
  gsl_multimin_fdfminimizer_free (s);
  gsl_vector_free (x);
   
  return rst;
  
}
  


  

double qtMSNP_ES::compute_log10_BF(double sa2,double oa2){

  // parameter setting
  prepare_params(sa2,oa2);
  if(info_size==0)
    return 0.0;

  // laplace approximation of the integral
  double altVal = laplace_approx(sa2,oa2);
  
  // under the null
  double nulVal = 0;
  
  for(int i=0;i<ssize;i++){
    if(sList[i].vg*indicator[i]==0)
      continue;
    double z = 0.5*sList[i].n;
    nulVal +=  -z*log(sList[i].rss0/2.0);
    nulVal += 0.5*log(2*3.1415926535/z) + z*(log(z)-1);
    //nulVal += gsl_sf_lngamma(z);
    altVal += -.5*log(sa2*sList[i].vg+1); 
  }
  
  altVal += -0.5*log(param_list[3*info_size]);
  
  
  double log10_BF = (altVal - nulVal)/log(10); 
  delete[] param_list;

  return log10_BF;
  
}



// ========================================= ABF Computation ====================================== //


double qtMSNP_ES::compute_log10_ABF(double sa2,double oa2){

  double rst =0;
  double xi = 0;
  double Z = 0;
  for(int i=0;i<ssize;i++){
    if(indicator[i]*sList[i].vg==0)
      continue;
    rst += sList[i].compute_ABF(sa2,0);
    double vs = sList[i].vg;
    xi += vs/(1 + sa2*vs);
    Z  += (vs/(1 + sa2*vs))*sList[i].bhat/sList[i].sigmahat;
  }
  double Z2 = 0;
  if(xi>0)
    Z2 = Z*Z/xi;

  rst += -0.5*log10(1+xi*oa2) + 0.5*Z2*xi*oa2/(1+xi*oa2)/log(10) ;
  return rst;
}


