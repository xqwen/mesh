using namespace std;

#include "qtMSNP.h"

class qtMSNP_EE : public qtMSNP {

 public:
   
 qtMSNP_EE(vector<qtSSNP> &list, string name) : qtMSNP(list, name){};
  
  double compute_log10_BF(double sa2, double oa2);
  double compute_log10_ABF(double sa2,double oa2);

 private:
 
  void prepare_params(double sa2, double oa2);
  double laplace_approx(double sa2, double oa2);
  

};
  
