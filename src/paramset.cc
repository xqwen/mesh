using namespace std;

#include "paramset.h"
#include <string>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


ParamSet::ParamSet(){
  
  // default value and initialization
  
  memset(outfile,0,128);
  memset(datafile,0,128);
  memset(gridfile,0,128);

  format = 1;
  
  model = 1;
  use_abf = 0;
  adjust_abf = 1;
  use_config = 0;

  
  output_subgroup = 0;
  config_hm_output = 0;
}


void ParamSet::parse_params(char *paramfile){
  
  ifstream pfs(paramfile);
  if(!pfs){
    fprintf(stderr,"Error: cannot open parameter file \"%s\".\n",paramfile);
    exit(1);
  }

  istringstream ins;
  string line;
  string param;
  string keyword;
  string value;

  while (getline(pfs,line) ){

    ins.clear();
    ins.str(line);

    ins>>ws;
    // empty line, skip                                                                                                                                                                                            
    if(ins.eof()){
      continue;
    }

    ins >> keyword;
    if(keyword!="#define")
      continue;

    ins >> param >> value;

    if(param == "OUTFILE"){
      if(strlen(outfile)==0){
        strcpy(outfile,value.c_str());
      }
      continue;
    }

    if(param == "DATAFILE"){
      if(strlen(datafile)==0){
        strcpy(datafile,value.c_str());
      }
      continue;
    }
    
    if(param == "GRIDFILE"){
      if(strlen(gridfile)==0){
        strcpy(gridfile,value.c_str());
      }
      continue;
    }

    if(param == "MODEL"){
      
      if(value == "ES" || value == "es" || atoi(value.c_str()) == 1)
	model = 1;
      else if(value == "EE" || value == "ee" || atoi(value.c_str()) == 2)
	model = 2;
      else if(value == "CEF_ES" || value == "cef_es" || value == "CEFES" || value == "cefes" || atoi(value.c_str()) == 3)
	model = 3;
      else if(value == "CEF_EE" || value == "cef_ee" || value == "CEFEE" || value == "cefee" || atoi(value.c_str()) == 4)
	model = 4;
      else
	fprintf(stderr, "Error: \"%s\" is not a valid option for MODEL\n",value.c_str());   
      continue;
    }

    
    if(param == "INPUT_FORMAT"){
      format = atoi(value.c_str());
      if(format != 1 && format != 2){
	fprintf(stderr, "Error: \"%s\" is not a valid option for INPUT_FORMAT\n",value.c_str());
      }
      continue;
    }



    if(param == "USE_ABF"){
      use_abf = atoi(value.c_str());
      continue;
    }

    if(param == "ADJUST_ABF"){
      adjust_abf = atoi(value.c_str());
    }



    if(param == "USE_CONFIG"){
      use_config = atoi(value.c_str()); 
      continue;
    }

    if(param == "OUTPUT_SUBGRP_DETAIL"){
      output_subgroup = atoi(value.c_str());
      continue;
    }

    if(param == "CONFIG_HM_OUTPUT"){
      config_hm_output = atoi(value.c_str());
      continue;
    }
    
    fprintf(stderr, "Warning: undefined parameter name \"%s\"\n", param.c_str());
    

  }

}



void ParamSet::print_settings(){
  fprintf(stderr,"\n");
  fprintf(stderr,"Input/Output Files\n");
  fprintf(stderr,"  Data File: \"%s\"\n",datafile);
  fprintf(stderr,"  Input Data File Format: ");
  if(format==1)
    fprintf(stderr,"N-bhat-se_bhat-var_g\n");
  if(format ==2)
    fprintf(stderr,"N-ym-yy-gm-gg-yg\n");
  if(format ==3)
    fprintf(stderr,"bhat-se_bhat\n");

  fprintf(stderr,"\n");
  fprintf(stderr,"  Grid File: \"%s\"\n",gridfile);
  fprintf(stderr,"  Output File: ");
  if(strlen(outfile)==0)
    fprintf(stderr," stdout (DEFAULT)");
  else
    fprintf(stderr," \"%s\"",outfile);
  fprintf(stderr,"\n\n");
  

  fprintf(stderr,"\n");
  
  fprintf(stderr,"Model Options\n");
  fprintf(stderr,"  Prior Choice: ");
  
  if(model==1)
    fprintf(stderr," ES\n");
  if(model==2)
    fprintf(stderr," EE\n");
  if(model==3)
    fprintf(stderr," CEF_ES\n");
  if(model==4)
    fprintf(stderr," CEF_EE\n");
  fprintf(stderr, "  Approximate Bayes Factor: %d\n",use_abf);
  if(use_abf){
    fprintf(stderr, "  Adjust ABF for Small Sample Size: %d\n", adjust_abf);
  }
  fprintf(stderr, "  Use Configurations: %d\n", use_config);
  fprintf(stderr,"\n");
  
  fprintf(stderr,"Output Options\n");
  fprintf(stderr,"  Print Subgroup Info: %d\n",output_subgroup);
  fprintf(stderr,"  Print Configuration Details: %d\n",config_hm_output);
  fprintf(stderr,"\n");


}


int ParamSet::check_params(){
  int status = 0;

  if(strlen(gridfile)==0){
    fprintf(stderr, "Error:  grid file is not specified.\n");
    status = -1;
  }else{
    if(!fopen(gridfile,"r")){
      fprintf(stderr, "Error: cannot open grid file \"%s\".\n", gridfile);
      status = -1;
    }
  }
  
  if(strlen(datafile)==0){
    fprintf(stderr, "Error: data file is not specified.\n");
    status = -1;
  }else{
    if(!fopen(datafile,"r")){
      fprintf(stderr, "Error: cannot open data file \"%s\".\n", datafile);
      status = -1;
    }
  }

  if(format == 3){
    if(model == 1 || model ==3){
      fprintf(stderr, "Error: Only EE and CEFEE models are applicable for minimum information format input.\n");
      status = -1;
    }
  }
  
  return status;


}
