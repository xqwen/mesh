#include "qtMSNP.h"
#include <math.h>

qtMSNP::qtMSNP(vector<qtSSNP> &list, string name){
  sList = list;
  ssize = sList.size();
  info_size = -1;
  snpName = name;

  for(int i=0;i<ssize;i++){
    indicator.push_back(1);
  }

}




double qtMSNP::log10_weighted_sum(double* vec, double *wts, int size){
  double max = vec[0];
  for(int i=0;i<size;i++){
    if(vec[i]>max)
      max = vec[i];
  }
  double sum = 0;
  for(int i=0;i<size;i++){
    sum += wts[i]*pow(10, (vec[i]-max));
  }

  return (max+log10(sum));
}


double qtMSNP::compute_log10_ABF(vector<double>& sa2_vec, vector<double>& oa2_vec){

  int size = int(sa2_vec.size());
  double *vec = new double[size];
  double *wts = new double[size];
  
  for(int i=0;i<size;i++){
    wts[i] = 1.0/size;
    vec[i] = compute_log10_ABF(sa2_vec[i],oa2_vec[i]);
  }
  
  double rst = log10_weighted_sum(vec,wts,size);
  delete[] vec;
  delete[] wts;
  
  return rst;
}

void qtMSNP::print_subgroup_info(FILE *fout){
  for(int i=0;i<sList.size();i++){
    fprintf(fout,"\t%7.3f\t%7.3f\t",sList[i].bhat, sList[i].sd_beta_orig);
  }
}



double qtMSNP::compute_log10_BF(vector<double>& sa2_vec, vector<double>& oa2_vec){

    int size = int(sa2_vec.size());
    double *vec = new double[size];
    double *wts = new double[size];

    for(int i=0;i<size;i++){
      wts[i] = 1.0/size;
      vec[i] = compute_log10_BF(sa2_vec[i],oa2_vec[i]);
    }

    double rst = log10_weighted_sum(vec,wts,size);
    delete[] vec;
    delete[] wts;
    return rst;
}

double qtMSNP::compute_log10_BF(double sa2, double oa2, int option){
 
 if(option==0)
    return compute_log10_BF(sa2,oa2);
  else
    return compute_log10_ABF(sa2,oa2);

}


double qtMSNP::compute_log10_BF(vector<double>& sa2_vec, vector<double> &oa2_vec, int option){
  
  if(option==0)
    return compute_log10_BF(sa2_vec,oa2_vec);
  else
    return compute_log10_ABF(sa2_vec,oa2_vec);
}



void qtMSNP::set_indicator(int mi){
  
  for(int i=0;i<ssize;i++){
    int r = mi%2;
    indicator[i] = r;
    mi = mi/2;
  }
}
 
  
