#ifndef __PARAMSET_H_
#define __PARAMSET_H_

using namespace std;

class ParamSet {

 public:

  // input files
  char gridfile[128];
  char datafile[128];

  // output file 
  char outfile[128]; // default (none specified) is stdout


  int format;
  // option
  // 1: (n-bhat-var-bhat-var-g format, default) 
  // 2: (n-ym-yy-gm-gg-yg format) 
  


  // model parameters 
  int model;
  // option
  // 1: ES (default)
  // 2: EE
  // 3: CEFES
  // 4: CEFEE

  int use_abf;
  // binary
  // 0: exact (default)
  

  
  int use_config;
  // binary
  // 0: default

  int adjust_abf;
  // binary: small sample correction for abf computing
  // 1: default
  

  // output options
  int output_subgroup;
  // NOT used with config
  // binary
  // 0: default
  
  int config_hm_output;
  // use with config
  // binary
  // 0: default


  int est_het;
  // estimate heterogeneity
  // binary
  // 0: default

  
  
  // constructor init everything
  ParamSet();
  void parse_params(char *paramfile);
  void print_settings();
  int check_params();

};

#endif

