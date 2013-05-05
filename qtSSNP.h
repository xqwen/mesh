using namespace std;



class qtSSNP {

 public:

  int typed;

  int adjust;


  // sample size                                                                    
  int n;

  // genotype info                                                                  
  double gtg; // G'G                                                                
  double g_mean; // n \bar G^2                                                      

  // phenotype info                                                                 
  double yty;  // Y'Y                                                               
  double y_mean;  // n \bar Y^2                                                     

  // correlation                                                                    
  double gty; // G'Y                                                                


  // useful summary info for meta-analysis                                          
  double vg;
  double rss0; // total variance                                                    
  double rss1; // rss by ols                                                        
  double chi2;

  double bhat;

  double sigmahat;
  double sd_beta;
  double sd_beta_orig;

  double ds2_es;
  double ds2_ee;

  // used by meta analysis                                                          

  double rssb; // rss bayesian                                                      
  double eta;

  double log10_BF;  // single study Bayes Factor                                    
  double log10_ABF;

  qtSSNP(){
    adjust = 1;
  }

  qtSSNP(int flag){
    adjust = flag;
  }


  void compute_rssb(double a2);

  void load_data(int n_, double g_mean_, double gtg_, double y_mean_, double yty_, double gty_);
  void load_data(int n_, double bhat_, double sd_beta_, double vx_);
  void small_sample_correction();

  void set_adjust_flag(int flag){
    adjust = flag;
  }
  double compute_ABF(double a2, int option=0);

};
