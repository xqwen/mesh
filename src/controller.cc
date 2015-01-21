#include "controller.h"
#include "qtMSNP_ES.h"
#include "qtMSNP_EE.h"
#include "qtMSNP_CEFES.h"
#include "qtMSNP_CEFEE.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


void controller::init(ParamSet &ps){

  est_het = ps.est_het;
  if(ps.est_het){
    load_het_grid(ps.gridfile);
  }else{
    load_grid(ps.gridfile);
  }
  
  msnp_option = ps.model;
  abf_option = ps.use_abf;
  adjust_abf = ps.adjust_abf;
  output_subgrp = ps.output_subgroup;
  use_config = ps.use_config;
  config_hm_output = ps.config_hm_output;
  format = ps.format;

  if(format == 3){
    abf_option = 1;
    adjust_abf = 0;// enforce
  }



  if(strlen(ps.outfile)!=0){
    fout = fopen(ps.outfile,"w");
    if(fout == 0){
      fprintf(stderr,"Error: cannot open FILE %s for writing\n",ps.outfile);
      fout = stdout;
    }
  }else
    fout = stdout;


  // important to load options first then parse data

  if(format == 1 || format == 3)
    load_data1(ps.datafile);
  if(format == 2)
    load_data2(ps.datafile);
  
  
    

}
  
void controller::run(){
  
  if(use_config){
    compute_config_BF();
  }else{
    if(est_het){
      compute_het_BF();
    }else{
      compute_meta_BF();
    }
  }
}

void controller::load_grid(char *gridfile){

  ifstream gfile(gridfile);
  string line;
  istringstream ins;
  while(getline(gfile,line)){
    
    ins.clear();
    ins.str(line);
    double phi;
    double omega;
    if(ins>>phi>>omega){
      double val = phi;
      if(msnp_option<=2)
	val = phi*phi;
      phi2_vec.push_back(val);
      omg2_vec.push_back(omega*omega);
    }
  }

  gfile.close();
  grid_size = omg2_vec.size();

}





void controller::load_het_grid(char *hgridfile){

  ifstream gfile(hgridfile);
  string line;
  istringstream ins;
  while(getline(gfile,line)){

    ins.clear();
    ins.str(line);
    string header;
    double val;
    ins>>header;
    if(header.compare("corr")==0){
      while(ins>>val){
	het_vec.push_back(val);
      }
    }

    if(header.compare("size")==0){
      while(ins>>val){
	size_vec.push_back(val);
      }
    }
    
  }

  gfile.close();
  
  std::sort(het_vec.begin(),het_vec.end());

  het_size = het_vec.size();
  for(int i=0;i<het_size;i++){
    vector<double> p2_vec;
    vector<double> o2_vec;
    double r = het_vec[i];
    for(int j=0;j<size_vec.size();j++){
      double s2 = size_vec[j]*size_vec[j];
      double omg2 = r*s2;
      double phi2 = s2 - omg2;
      p2_vec.push_back(phi2);
      o2_vec.push_back(omg2);
    }
    het_phi2_vec.push_back(p2_vec);
    het_omg2_vec.push_back(o2_vec);
  }


}





void controller::load_data1(char *datafile){
  
  
  // reading data file
  vector<qtSSNP> slist;
  double size,bhat,sd_bhat,varx;
  
  string sid;
  string stype;
  
  string curr_sid;  

  ifstream infile(datafile);
  istringstream ins;
  istringstream ins_counter;
  string line;
  int line_count=0;

  while(getline(infile,line)){
    
    
    line_count++;
    
    ins_counter.clear();
    ins_counter.str(line);
    int word_count = 0;
    string word;

    while(ins_counter >> word) { word_count++;}


    if(format==1 && word_count!= 6){
      fprintf(stderr,"Error: unexpected input data format in file %s at line %d (6 entries expected, %d entries observed)\n",datafile, line_count, word_count);
      exit(1);
    }
    

    if(format==3 && word_count!= 4){
      fprintf(stderr,"Error: unexpected input data format in file %s at line %d (4 entries expected, %d entries observed)\n",datafile, line_count, word_count);
      exit(1);
    }




    
    ins.clear();
    ins.str(line);
    if(format == 1)
      ins>>sid>>stype>>size>>bhat>>sd_bhat>>varx;
    
    if(format == 3){
      ins>>sid>>stype>>bhat>>sd_bhat;
      size = 100; // irrelevant set to arbitary value
      varx = 1;   // irrelevant
    }
    
    qtSSNP ssnp(adjust_abf);
    ssnp.load_data(int(size),bhat,sd_bhat,varx);
    

    if(strcmp(sid.c_str(), curr_sid.c_str())!=0){
      
      if(slist.size()!=0){
	insert_msnp(slist,curr_sid);
	slist.clear();
      }
      curr_sid = sid;
    }
    
    slist.push_back(ssnp);    

  }
  
  if(slist.size()!=0){
    insert_msnp(slist,curr_sid);
  }


}









void controller::load_data2(char *datafile){
  
  
  // reading data file
  vector<qtSSNP> slist;
  double size,ym,yty,gm,gtg,gty;
  
  string sid;
  string stype;
  
  string curr_sid; 
  

  ifstream infile(datafile);
  istringstream ins;
  string line;

  int line_count=0;
  istringstream ins_counter;
    
  while(getline(infile,line)){


    line_count++;
    
    ins_counter.clear();
    ins_counter.str(line);
    int word_count = 0;
    string word;
    while(ins_counter >> word) { word_count++;}
    if(word_count!= 8){
      fprintf(stderr,"Error: unexpected input data format in file %s at line %d (8 entries expected, %d entries observed)\n",datafile, line_count, word_count);
      exit(1);
    }
    
    
    ins.clear();
    ins.str(line);
    ins>>sid>>stype>>size>>ym>>yty>>gm>>gtg>>gty;
    
        
    qtSSNP ssnp(adjust_abf);
    ssnp.load_data(int(size),gm,gtg,ym,yty,gty);
    

    if(strcmp(sid.c_str(), curr_sid.c_str())!=0){
      if(slist.size()!=0){
	insert_msnp(slist,curr_sid);
	slist.clear();
      }
      curr_sid = sid;
    }
    
    slist.push_back(ssnp);    

  }
  
  if(slist.size()!=0){
    insert_msnp(slist,sid);
  }
}




void controller::insert_msnp(vector<qtSSNP>& slist, string sid){
  
    if(msnp_option == 1){
      qtMSNP_ES *qm = new qtMSNP_ES(slist,sid);
      dataVec.push_back(qm);
    }

    if(msnp_option == 2){
      qtMSNP_EE *qm = new qtMSNP_EE(slist,sid);
      dataVec.push_back(qm);
    }
    
    if(msnp_option == 3){
      qtMSNP_CEFES *qm = new qtMSNP_CEFES(slist,sid);
      dataVec.push_back(qm);
    }

    if(msnp_option == 4){
      qtMSNP_CEFEE *qm = new qtMSNP_CEFEE(slist,sid);
      dataVec.push_back(qm);
    }



}



void controller::compute_config_BF(){
  if(config_hm_output){
    for(int i=0;i<dataVec.size();i++){    
      int size = (1<<dataVec[i]->get_ssize())-1;
      for(int mi = 1; mi<=size; mi++){
	fprintf(fout,"%s\t%d",(dataVec[i]->get_name()).c_str(),mi);
	dataVec[i]->set_indicator(mi);
	for(int k=0;k<phi2_vec.size();k++){
	  fprintf(fout,"\t%8.3f", dataVec[i]->compute_log10_BF(phi2_vec[k], omg2_vec[k], abf_option));
	}
	fprintf(fout,"\n");
      }
    }
  }else{
    for(int i=0;i<dataVec.size();i++){    
      int size = (1<<dataVec[i]->get_ssize())-1;
      fprintf(fout,"%s",(dataVec[i]->get_name()).c_str());
      for(int mi = 1; mi<=size; mi++){
	dataVec[i]->set_indicator(mi);
	
	fprintf(fout,"\t%8.3f", dataVec[i]->compute_log10_BF(phi2_vec, omg2_vec, abf_option));
	
      }
      fprintf(fout,"\n");
    }
  }
    
}
  



void controller::compute_meta_BF(){
   
  for(int i=0;i<dataVec.size();i++){
    double rst = 0;

    rst = dataVec[i]->compute_log10_BF(phi2_vec,omg2_vec,abf_option);
    
    fprintf(fout,"%20s\t%9.4f",(dataVec[i]->get_name()).c_str(), rst);
    
    if(output_subgrp){
      fprintf(fout,"\t\t");
      dataVec[i]->print_subgroup_info(fout);
    }
    
    fprintf(fout,"\n");
    
  }

}




void controller::compute_het_BF(){
 
  for(int i=0;i<dataVec.size();i++){
    double rst = 0;
    vector<double> hrst_vec;
    for(int j=0;j<het_size;j++){      
      double hrst = 0;
      hrst = dataVec[i]->compute_log10_BF(het_phi2_vec[j],het_omg2_vec[j],abf_option);
      hrst_vec.push_back(hrst);
    }
    fprintf(fout,"%20s\t%9.4f\t%9.4f\t%7.3f\n",(dataVec[i]->get_name()).c_str(), hrst_vec[0],hrst_vec[het_size-1], compute_het_weight(hrst_vec));
   
  }
  

}




double controller::compute_het_weight(vector<double> &vec){
  
  double max = -99999;
  for(int i=0;i<vec.size();i++){
    if(vec[i]>=max){
      max = vec[i];
    }
  }

  vector<double> vec2;
  double sum = 0;
  for(int i=0;i<vec.size();i++){
    double val = pow(10,vec[i]-max);
    vec2.push_back(val);
    sum += val;
  }

  double rst = 0;
  for(int i=0;i<het_vec.size();i++){
    rst += het_vec[i]*vec2[i]/sum;
  }
  /*
  for(int i=0;i<vec.size();i++){
    printf(" %4.1f ",vec[i]);
  }
  */
  return rst;
  



}  
  
