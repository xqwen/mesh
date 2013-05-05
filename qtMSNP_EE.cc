#include "qtMSNP_EE.h"
#include <math.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf.h>


double ee_target_f (const gsl_vector *v, void *params){
  
  double *p = (double *)params;
  int size = v->size;
  
  
  double sa2 = p[0];
  double oa2 = p[1];
  double eta2 = 0;
  
  for(int i=0;i<size;i++){
    double ds2= p[5*i+3];
    double t = fabs(gsl_vector_get(v,i));
    eta2 += t/(ds2+t*sa2);
  }
  
  eta2 = 1.0/eta2;
  
  double val = .5*log(eta2/(eta2+oa2));
  double val2 = 0;
  
  for(int i=0;i<size;i++){
    double t = fabs(gsl_vector_get(v,i));

    double n = p[5*i+2];
    double ds2 = p[5*i+3]; 
    double rss0 = p[5*i+4];
    double rss1 = p[5*i+5];
    

    double bhat = p[5*i+6];
    val += .5*log(ds2/(ds2 + t*sa2))+.5*(n-3)*log(t)-.5*t*(t*sa2*rss1+ds2*rss0)/(ds2+t*sa2); 
    
    val2 += t*bhat/(ds2+t*sa2);
    
  }

  val += 0.5*pow(val2,2)*eta2*oa2/(eta2 + oa2);
  
  return -val;
}
 

void ee_target_df (const gsl_vector *v, void *params, gsl_vector *df){
  
  double *p = (double *)params;
  
  int size = v->size;
  
  
  double sa2 = p[0];
  double oa2 = p[1];
  double eta2 = 0;
 
  double * DgDt = new double[size];
  double * Deta2Dt = new double[size];
  double * temp = new double[size];
  
  for(int i=0;i<size;i++){
    double t = gsl_vector_get(v,i);
    double ds2 = p[5*i+3];
    eta2 += t/(ds2+t*sa2);
    DgDt[i]=ds2*pow((ds2+t*sa2),-2); 
  }
  
  eta2 = 1.0/eta2;
  double val = 0;
  
  for(int i=0;i<size;i++){
    double t = gsl_vector_get(v,i);
    
    double n = p[5*i+2];
    double ds2 = p[5*i+3];
    double rss0 = p[5*i+4];
    double rss1 = p[5*i+5];
    double bhat = p[5*i+6];
    
    Deta2Dt[i] = -pow(eta2,2)*DgDt[i]; 
    temp[i] = .5*Deta2Dt[i]*(1/eta2 - 1/(eta2+oa2))-.5*sa2/(ds2+t*sa2)+.5*(n-3)/t -.5*(DgDt[i]*(t*sa2*rss1+ds2*rss0)+t*sa2*rss1/(ds2+t*sa2));
    val += t*bhat/(ds2+t*sa2);
  }
  
  for(int i=0;i<size;i++){
    double bhat = p[5*i+6];
    temp[i] += .5*oa2*( Deta2Dt[i]*pow(val,2)/(eta2+oa2) -pow(eta2+oa2,-2)*Deta2Dt[i]*eta2*pow(val,2) + 2*val*bhat*DgDt[i]*eta2/(eta2+oa2) );
    //temp[i] += .5*( pow(val,2)*(oa2*Deta2Dt[i]/(eta2+oa2) - pow(val,2)*eta2*oa2*Deta2Dt[i]/pow((eta2+oa2),2)) + 2*val*bhat*DgDt[i]*eta2*oa2/(eta2+oa2));   
    if(gsl_vector_get(v,i)<0)
      temp[i] *= -1;
    
    gsl_vector_set(df,i, -temp[i] );
  }
  
  delete[] temp;
  delete[] DgDt;
  delete[] Deta2Dt;
  
}
     

/* Compute both f and df together. */
void  ee_target_fdf (const gsl_vector *x, void *params, double *f, gsl_vector *df){
  *f = ee_target_f(x, params); 
  ee_target_df(x, params, df);
} 




void qtMSNP_EE::prepare_params(double sa2, double oa2){
  
  // find non-informative sets
 
  info_size = 0;
  for(int i=0;i<ssize;i++){
    if(sList[i].vg*indicator[i] != 0)
      info_size++;
  }
  
  
  if(info_size==0)
    return;

  param_list = new double[5*info_size+2];
  
  param_list[0] = sa2;
  param_list[1] = oa2;

  int i=0;
  for(int j=0;j<ssize;j++){
    
    if(sList[j].vg*indicator[j] == 0)
      continue;

    param_list[5*i+2] = sList[j].n;
    param_list[5*i+3] = 1.0/sList[j].vg;
    param_list[5*i+4] = sList[j].rss0;
    param_list[5*i+5] = sList[j].rss1;
    param_list[5*i+6] = sList[j].bhat;

    //printf("%f  %f  %f  %f  %f\n",param_list[5*i+2],param_list[5*i+3],param_list[5*i+4],param_list[5*i+5],param_list[5*i+6]);
    i++;
  }
  

}



double qtMSNP_EE::laplace_approx(double sa2, double oa2){

  size_t iter = 0;
  int status;
  

 
  const gsl_multimin_fdfminimizer_type *T;
  gsl_multimin_fdfminimizer *s;
     
  gsl_vector *x;
  gsl_multimin_function_fdf target_func;
  
  target_func.n  = info_size;
  target_func.f  = ::ee_target_f;
  target_func.df = ::ee_target_df;
  target_func.fdf = ::ee_target_fdf;
  target_func.params = param_list;

  
  /* Starting point */
  vector<double> tauv;
  x = gsl_vector_alloc(info_size);
  int j=0;
  for(int i=0;i<sList.size();i++){
    if(sList[i].vg*indicator[i] == 0)
      continue;
    double val = double(sList[i].n-2)/sList[i].rss1;
    gsl_vector_set (x, j++, val);  
    tauv.push_back(val);
  }

  double fval;
  
  //T = gsl_multimin_fdfminimizer_conjugate_fr;
  //T = gsl_multimin_fdfminimizer_steepest_descent;
  T = gsl_multimin_fdfminimizer_vector_bfgs2;
  s = gsl_multimin_fdfminimizer_alloc (T, info_size);
     
  gsl_multimin_fdfminimizer_set (s, &target_func, x, 0.01, 0.01);

  do{  
    iter++;  


    if(isnan(s->f)){
      break;
    }else{
      //book keeping
      fval=s->f;
      for(int i=0;i<info_size;i++)
    	tauv[i] = gsl_vector_get(s->x,i);
    }
    


    status = gsl_multimin_fdfminimizer_iterate (s);
    
    //printf("%d  %f\n",iter,s->f);

    if (status)
      break;
 
    status = gsl_multimin_test_gradient (s->gradient, 1e-3);

   
  }while (status == GSL_CONTINUE && iter < 5000);
  
  
  if(!isnan(s->f)){
    for(int i=0;i<info_size;i++){
      tauv[i] = gsl_vector_get(s->x,i);  
      //printf("max: %f  %f\n",val,sList[i].n/sList[i].rss1);
    }
    fval = s->f;
  }
  

  
  // compute Hessian at the minimization point
  
  gsl_matrix *hessian = gsl_matrix_calloc(info_size,info_size);
  
  int size = info_size;

  double * DgDt = new double[size];
  double * D2gDt2 = new double[size];
  
  double * Deta2Dt = new double[size];
  double ** D2eta2Dt2 = new double* [size];

  double * DthetaDt = new double[size];
  double * D2thetaDt2 = new double[size];
  
  double * DvDt = new double [size];
  double ** D2vDt2 = new double* [size];
  
  double *p = param_list;
  
  double eta2 = 0;
  double theta = 0;

  for(int i=0;i<size;i++){
    double t = tauv[i];
    double ds2 = p[5*i+3];
    double bhat = p[5*i+6];
    eta2 += t/(ds2+t*sa2);
    theta += bhat*t/(ds2+t*sa2);
    DgDt[i]=ds2*pow((ds2+t*sa2),-2);
    D2gDt2[i] = -2*ds2*sa2*pow((ds2+t*sa2),-3);
    DthetaDt[i] = bhat*DgDt[i];
    D2thetaDt2[i] = bhat*D2gDt2[i];
  }

  eta2 = 1.0/eta2;
  double v = eta2*oa2/(eta2+oa2);
  
  for(int i=0;i<size;i++){
    D2eta2Dt2[i] = new double[size];
    D2vDt2[i] = new double[size];
    Deta2Dt[i] = -pow(eta2,2)*DgDt[i];
    DvDt[i] = oa2*Deta2Dt[i]*(1.0/(eta2+oa2) - eta2*pow((eta2+oa2),-2));
  }
  
  // second derivatives
  for(int i=0;i<size;i++){
    for(int j=i;j<size;j++){
      if(i==j){
	D2eta2Dt2[i][i] = -2*eta2*Deta2Dt[i]*DgDt[i]-pow(eta2,2)*D2gDt2[i];
      }else{
	D2eta2Dt2[i][j] = D2eta2Dt2[j][i] = -2*eta2*Deta2Dt[j]*DgDt[i];
      }
    }
  }

  
  for(int i=0;i<size;i++){
    for(int j=i;j<size;j++){
      if(i==j){
        D2vDt2[i][i] = oa2*D2eta2Dt2[i][i]*(1.0/(eta2+oa2)-eta2*pow(eta2+oa2,-2))+2*oa2*pow(Deta2Dt[i],2)*(eta2*pow(eta2+oa2,-3)-pow(eta2+oa2,-2));
      }else{
	D2vDt2[i][j]= D2vDt2[j][i] = oa2*D2eta2Dt2[i][j]*(1.0/(eta2+oa2)-eta2*pow(eta2+oa2,-2))+2*oa2*Deta2Dt[i]*Deta2Dt[j]*(eta2*pow(eta2+oa2,-3)-pow(eta2+oa2,-2));
      }
    }
  }



  for(int i=0;i<size;i++){
    double t = tauv[i];
    double n = p[5*i+2];
    double ds2 = p[5*i+3];
    double rss0 = p[5*i+4];
    double rss1 = p[5*i+5];
    double bhat = p[5*i+6];
    for(int j=i;j<size;j++){
      double val = 0;
      if(i==j){
	val += 0.5*(D2eta2Dt2[i][i]*(1/eta2-1/(eta2+oa2)) + pow(Deta2Dt[i],2)*(pow(eta2+oa2,-2)-pow(eta2,-2)));
	val += .5*pow(sa2,2)*pow(ds2+t*sa2,-2);
	val += -.5*(n-3)*pow(t,-2);
	val += -.5*(2*DgDt[i]*sa2*rss1+(t*sa2*rss1+ds2*rss0)*D2gDt2[i]);
	val +=  .5*(pow(theta,2)*D2vDt2[i][i]+4*theta*DthetaDt[i]*DvDt[i]+2*v*pow(DthetaDt[i],2)+2*theta*v*D2thetaDt2[i]);
	gsl_matrix_set(hessian,i,i,-val);
      }else{
	val += 0.5*(D2eta2Dt2[i][j]*(1/eta2-1/(eta2+oa2)) + Deta2Dt[i]*Deta2Dt[j]*(pow(eta2+oa2,-2)-pow(eta2,-2)));
	val +=  .5*(pow(theta,2)*D2vDt2[i][j]+2*theta*(DthetaDt[i]*DvDt[j]+DthetaDt[j]*DvDt[i])+ 2*v*DthetaDt[i]*DthetaDt[j]);
	gsl_matrix_set(hessian,i,j,-val);
	gsl_matrix_set(hessian,j,i,-val);
      }
      
    }
  }
  
  
  delete[] DgDt;
  delete[] D2gDt2;
  delete[] Deta2Dt;
  delete[] DthetaDt;
  delete[] D2thetaDt2;

  delete[] DvDt;
  
  for(int i=0;i<size;i++){
    delete[] D2eta2Dt2[i];
    delete[] D2vDt2[i];
  }
  
  delete[] D2eta2Dt2;
  delete[] D2vDt2;


    
  // LU decomp and det computation
  int ss;
  gsl_permutation * pp = gsl_permutation_alloc (info_size);
     
  int flag = gsl_linalg_LU_decomp (hessian, pp, &ss);
  double log_detVal = gsl_linalg_LU_lndet(hessian);
  

  double rst = 0.5*info_size*log(2*3.1415926536)-0.5*log_detVal - fval;
  
  gsl_matrix_free(hessian);
  
  gsl_permutation_free (pp);
  gsl_multimin_fdfminimizer_free (s);
  gsl_vector_free (x);
  
  return rst;
}






double qtMSNP_EE::compute_log10_BF(double sa2,double oa2){

  // parameter setting
  prepare_params(sa2,oa2);
  if(info_size==0)
    return 0.0;

  // laplace approximation of the integral
  double altVal = laplace_approx(sa2,oa2);
  
  // under the null
  double nulVal = 0;
  
  for(int i=0;i<ssize;i++){
    if(sList[i].vg*indicator[i] == 0)
      continue;
    double z = 0.5*(sList[i].n-1);
    nulVal +=  -z*log(sList[i].rss0/2.0);
    nulVal += 0.5*log(2*3.1415926535/z) + z*(log(z)-1);
    //nulVal += gsl_sf_lngamma(z);
  }
  
  
  
  
  double log10_BF = (altVal - nulVal)/log(10); 
  delete[] param_list;

  return log10_BF;
  
}  




double qtMSNP_EE::compute_log10_ABF(double sa2,double oa2){

  double rst =0;
  double xi = 0;
  double Z = 0;
  for(int i=0;i<ssize;i++){
    if(indicator[i]*sList[i].vg==0)
      continue;
    rst += sList[i].compute_ABF(sa2,1);
    xi += 1.0/(sList[i].ds2_ee + sa2);
    Z  += sList[i].bhat/(sList[i].ds2_ee + sa2);;
  }
  double Z2 = 0;
  if(xi>0)
    Z2 = Z*Z/xi;

  rst += -0.5*log10(1+xi*oa2) + 0.5*Z2*xi*oa2/(1+xi*oa2)/log(10) ;
  return rst;
}
  


