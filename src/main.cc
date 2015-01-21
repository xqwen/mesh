#include "paramset.h"
#include "controller.h"
#include <fstream>
#include <sstream>
#include <string.h>
#include <stdlib.h>



void print_banner(){
}


void kill(){
  fprintf(stderr,"\nExiting due to error(s) listed above.\n\n");
  exit(1);
}



int main(int argc, char **argv){
  
  
  ParamSet ps;
  

  // parameters
  int format  = -1;
  int model = -1;
  int use_abf = -1; 
  int use_config = -1;
  int config_hm_output = -1;
  int output_subgroup = -1;
  int adjust_abf = -1;
  int verbose = 0;
  int est_het = 0;

  char gridfile[128];
  char datafile[128];
  char outfile[128];
  char paramfile[128];
  memset(gridfile,0,128);
  memset(datafile,0,128);
  memset(outfile,0,128);
  memset(paramfile,0,128);
  
  for(int i=1;i<argc;i++){
     
    if(strcmp(argv[i], "-g")==0 || strcmp(argv[i], "-grid")==0){
      strcpy(gridfile,argv[++i]);
      continue;
    } 
    
    if(strcmp(argv[i], "-d")==0 || strcmp(argv[i], "-data")==0){
      strcpy(datafile,argv[++i]);
      continue;
    }

    if(strcmp(argv[i], "-o")==0){
      strcpy(outfile,argv[++i]);
      continue;
    }

    
    if(strcmp(argv[i], "-format")==0){
      format = atoi(argv[++i]);
      continue;
    }

    
    if(strcmp(argv[i], "-es")==0){
      model = 1;
      continue;
    }

    
    
    if(strcmp(argv[i], "-ee")==0){
      model = 2;
      continue;
    }
    
    if(strcmp(argv[i], "-cefes")==0||strcmp(argv[i], "-cef_es")==0){
      model = 3;
      continue;
    }

    if(strcmp(argv[i], "-cefee")==0||strcmp(argv[i], "-cef_ee")==0){
      model = 4;
      continue;
    }

    
    if(strcmp(argv[i], "-abf")==0){
      use_abf = 1;
      continue;
    }


    if(strcmp(argv[i], "-est_het")==0){
      est_het = 1;
      continue;
    }


    if(strcmp(argv[i], "-use_config")==0){
      use_config = 1;
      continue;
    }

    if(strcmp(argv[i], "-print_config_grid")==0){
      output_subgroup = 1;
      continue;
    }

    
    if(strcmp(argv[i], "-prep_hm")==0){
      use_config = 1;
      config_hm_output=1;
      continue;
    }

    if(strcmp(argv[i], "-no_adjust")==0){
      adjust_abf = 0;
      continue;
    }

    if(strcmp(argv[i], "-min_info") ==0){
      format = 3;
      use_abf = 1;
      adjust_abf = 0;
      continue;
    }

   
    if(strcmp(argv[i], "-print_subgrp")==0){
      output_subgroup = 1;
      continue;
    }
    
    if(strcmp(argv[i], "-param") == 0){
      strcpy(paramfile,argv[++i]);
      continue;
    }

    if(strcmp(argv[i], "-v") == 0){
      verbose = 1;
      continue;
    }


    fprintf(stderr,"Error: unknown option \"%s\"\n",argv[i]);
    
    
  }

  
  if(strlen(paramfile)>0){
    ps.parse_params(paramfile);
  }

  // override param file settings if specified in command line
  
  if(format==3 && model ==-1){
    model =2;
  }
  
  if(format !=-1)
    ps.format = format;
  if(model != -1)
    ps.model = model;
  if(use_abf != -1)
    ps.use_abf = use_abf;
  if(use_config != -1)
    ps.use_config = use_config;
  if(config_hm_output!=-1)
    ps.config_hm_output = config_hm_output;
  if(output_subgroup != -1)
    ps.output_subgroup = output_subgroup;
  if(adjust_abf != -1)
    ps.adjust_abf = adjust_abf;
  
  ps.est_het = est_het;
  

  if(strlen(gridfile)!=0)
    strcpy(ps.gridfile,gridfile);
  if(strlen(datafile)!=0)
    strcpy(ps.datafile,datafile);
  if(strlen(outfile)!=0)
    strcpy(ps.outfile, outfile);
  
  if(ps.check_params() == -1)
    kill();

  if(verbose)
    ps.print_settings();


  
  controller con;
  con.init(ps);
  con.run();
}


