#include "classdef.h"
#include <fstream>
#include <sstream>



int main(int argc, char **argv){
  
  // creating the grid
  
  //olist.push_back(0.1);
  //phlist.push_back(0.05);

  char data_file[128];
  int csize = -1;
  int gsize = -1;
  int nthread = 1;
  memset(data_file,0,128);


  char init_file[128];
  int option = 0;

  char ci_file[128];
  memset(ci_file,0,128); 
  memset(data_file,0,128); 
  memset(init_file,0,128);
  
  
  double thresh = 0.05;
  

  for(int i=1;i<argc;i++){
    
    if(strcmp(argv[i], "-d")==0 || strcmp(argv[i], "-data")==0){
      strcpy(data_file,argv[++i]);
      continue;
    }
    
    if(strcmp(argv[i], "-s")==0){
      csize = atoi(argv[++i]);
      continue;	
    }
    
    if(strcmp(argv[i], "-g")==0 || strcmp(argv[i], "-grid")==0){
      gsize = atoi(argv[++i]);
      continue;
    }
    
    ////////////////////// optional ///////////////////////////

    if(strcmp(argv[i],"-i")==0 || strcmp(argv[i],"-init")==0){
      strcpy(init_file, argv[++i]);
      continue;
    }
    
   if(strcmp(argv[i], "-r")==0 || strcmp(argv[i],"-ran")==0){
     option = 1;
     continue;
   }
   
   if(strcmp(argv[i], "-t")==0 || strcmp(argv[i], "-thresh")==0){
     thresh = atof(argv[++i]);
     continue;
   }

   if(strcmp(argv[i], "-thread")==0){
     nthread = atof(argv[++i]);
     continue;
   }
   
   if(strcmp(argv[i], "-c")==0 || strcmp(argv[i],"-ci")==0){
     strcpy(ci_file, argv[++i]);
     continue;
   }
   

  }

  // checking mandatory arguments
  if(strlen(data_file)==0){
    fprintf(stderr,"Error: data file unspecified\n");
    exit(1);
  }

  if(csize==-1){
    fprintf(stderr,"Error: number of model types unspecified\n");
    exit(1);
  }
  
  if(gsize==-1){
    fprintf(stderr,"Error: number of model grids unspecified\n");
    exit(1);
  }

    
  // a global variable 
  eQTL_controller controller;
  controller.nthread = nthread;
  controller.output_option = 1;
  controller.load_data(data_file,csize,gsize);
  fprintf(stderr,"Finish loading ... \n");
  
  if(strlen(ci_file)>0){
    controller.init_params(ci_file);
    controller.estimate_profile_ci();
    return 0;
  }

  if(strlen(init_file)>0)
    controller.init_params(init_file);
  else
    controller.init_params(option);
  
  controller.run_EM(thresh);
  controller.compute_posterior();
  controller.print_result();
  controller.estimate_profile_ci();
}


