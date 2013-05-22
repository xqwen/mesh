using namespace std;

#include "qtMSNP.h"
#include "paramset.h"
#include <vector>
#include <string>
#include <stdio.h>

class controller {

 private:
  // storage
  vector<qtMSNP*> dataVec;
  
  // grid of omega^2 and phi^2
  vector<double> phi2_vec; 
  vector<double> omg2_vec; 
  
  // parameters need to be estimated
  int grid_size;  
  int msnp_option; // 1 -- ES; 2-- EE


  int abf_option;
  int output_subgrp;
  int use_config;
  int config_hm_output;
  int adjust_abf;
  int format;

  FILE *fout;
  
 private:
  void insert_msnp(vector<qtSSNP> & slit, string sid);

  void load_grid(char *gridfile);
  void load_data1(char *datafile);
  void load_data2(char *datafile);
  void compute_meta_BF();
  void compute_config_BF();

  
 public:
  
  void init(ParamSet &ps);
  void run();
  
 
};
