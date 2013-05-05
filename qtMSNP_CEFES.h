using namespace std;

#include "qtMSNP.h"

class qtMSNP_CEFES : public qtMSNP {

 public:
   
 qtMSNP_CEFES(vector<qtSSNP> &list, string name) : qtMSNP(list, name){};
  
  double compute_log10_BF(double k, double oa2);
  double compute_log10_ABF(double k,double oa2);

 private:
 
  void prepare_params(double k, double oa2);
  double laplace_approx(double k, double oa2);
  

};
  
