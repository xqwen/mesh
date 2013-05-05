#ifndef _QTMSNP_H_
#define _QTMSNP_H_

using namespace std;

#include "qtSSNP.h"
#include <vector>
#include <string>
#include <stdio.h>

class qtMSNP {

 protected:
  
  string  snpName;
  vector<qtSSNP> sList; 
  
  vector<int> indicator; // indicator for used subgroups in analysis
  
  int ssize;  // number of participants
  int info_size;  // number of informative participants: vg != 0   

  double *param_list;
  
  

  
 public:
  qtMSNP(vector<qtSSNP> &list, string name);

  virtual double compute_log10_BF(double sa2, double oa2) = 0;
  virtual double compute_log10_ABF(double sa2,double oa2) = 0;
 
  double compute_log10_BF(vector<double>& sa2_vec, vector<double> &oa2_vec); 
  double compute_log10_ABF(vector<double>& sa2_vec, vector<double> &oa2_vec); 

  double compute_log10_BF(double sa2, double oa2, int option);
  double compute_log10_BF(vector<double>& sa2_vec, vector<double> &oa2_vec, int option);
  void print_subgroup_info(FILE *fout);


  void set_indicator(vector<int> & indicator_){
    indicator = indicator_;
  }

  
  void set_indicator(int index);

 
  double log10_weighted_sum(double* vec, double *wts, int size);
 
  string get_name(){
    return snpName;
  }
  

  int get_ssize(){
    return ssize;
  }
};

#endif
