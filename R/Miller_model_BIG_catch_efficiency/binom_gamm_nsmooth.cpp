#include <TMB.hpp>
#include <iostream>

#define see(object) std::cout << #object ":\n" << object << "\n";

template<class Type>
Type objective_function<Type>::operator() ()
{

  //init_int n_sta
  //DATA_IVECTOR(n_per_sta); //n_sta
  DATA_IVECTOR(id); //nobs
  DATA_VECTOR(n); //nobs 
  DATA_VECTOR(y); //nobs
  DATA_VECTOR(offst); //nobs
  DATA_MATRIX(X); //nobs x (n_beta_fixed + n_station_smooth_re) covariates for "population level" model.
  DATA_IVECTOR(beta_ns_ind); //which of the full beta vector are fixed effects.
  DATA_IVECTOR(fp_mean_smooth_re); //(nsmu) first column of X for each smoothers random effects
  DATA_IVECTOR(lp_mean_smooth_re); //(nsmu) last column of X for each smoothers random effects
  DATA_MATRIX(Z); //(nobs x (n_beta_re + n_station_smooth_re) covariates for any random effects other than smoothers
  DATA_IVECTOR(beta_re_ns_ind); //which columns of the full beta_re vector are not smoother random effects.
  DATA_IVECTOR(fp_station_smooth_re); //(nsmu_re) first column of Z for each smoothers random effects
  DATA_IVECTOR(lp_station_smooth_re); //(nsmu_re) last column of Z for each smoothers random effects
  //DATA_MATRIX(X_lambda_station);
  DATA_VECTOR(Dplus_diag); // diagonal of 1/variance components for random effects for all smoothers
  DATA_VECTOR(Dplus_diag_station); // diagonal of 1/variance components for station-specfic random effects for all smoothers
  DATA_IMATRIX(id_Dplus_diag_station) // n_station x ncol(station_smooth_re) indicate which component of Dplus_diag_station to apply to each station_smooth_re
  //DATA_IMATRIX(id_lambda_station) // n_station x ncol(station_smooth_re) indicate which lambda_par_station to apply to each station_smooth_re
  DATA_INTEGER(use_beta_re); 
  DATA_INTEGER(use_beta_re_cor); //have beta_re correlated
  DATA_INTEGER(use_mean_smooth_re);
  DATA_INTEGER(use_station_smooth_re);
  DATA_INTEGER(do_pred); //predict and report linear predictions of new data provided in Xpred and Zpred
  DATA_MATRIX(Xpred);
  DATA_MATRIX(Zpred);
 
  PARAMETER_VECTOR(beta); // n_beta_fixed //fixed effects
  PARAMETER_VECTOR(mean_smooth_re); //n_smooth // population-level average random effects portion of smoothers
  PARAMETER_MATRIX(station_smooth_re); //n_sta x n_station_smooth_re //random effects part of station-specific smooths
  PARAMETER_VECTOR(lambda_par); //log(mean) smoothing parameters across stations
  PARAMETER_MATRIX(lambda_par_station);  //nstation x n  random smoothers for each station
  PARAMETER_MATRIX(beta_re); //n_sta x n_beta_re //random effects for fixed effects part of station-specific smooth
  PARAMETER_MATRIX(beta_re_var_pars); //nstation x ncol(beta_re), //variance of random effects for fixed effects part of station-specific smooth use chol decomp to ensure pos-def
  PARAMETER_MATRIX(beta_re_cor_pars); //nstation x ncol(beta_re)*(ncol(beta_re)-1)/2, //correlation of random effects for fixed effects part of station-specific smooth use chol decomp to ensure pos-def
  using namespace density;
  Type nll = 0.0;
  int n_stations = lambda_par_station.rows();
  //vector<Type> log_lambda_station = X_lambda_station * lambda_par_station;
  int N = n.size();
  vector<Type> mu(N);
  
  vector<Type> beta_use(X.cols());
  beta_use.fill(0);
  matrix<Type> beta_re_use(beta_re.rows(), Z.cols());
  beta_re_use.fill(0);
  int nsmu = fp_mean_smooth_re.size();
  int nsmu_re = fp_station_smooth_re.size();
  REPORT(nsmu);
  REPORT(nsmu_re);
  
  matrix<Type> nll_station_smooth_re(station_smooth_re.rows(),station_smooth_re.cols());
  nll_station_smooth_re.fill(0);
  if(use_station_smooth_re == 1)
  {
    for(int i = 0; i < n_stations; i++) 
    {
      //distribution of random effects portion of smoother given lambda at the station
      //vector<Type> sd_station_smooth_re_i = exp(Type(-0.5) * (log(Dplus_diag) + lambda_par_stations)); 
      //for(int j = 0; j < station_smooth_re.cols(); j++) 
      int pp = 0;
      for(int s = 0; s < nsmu_re; s++)  for(int k = fp_station_smooth_re(s); k <= lp_station_smooth_re(s); k++) //
      {
        Type sd = exp(Type(-0.5) * (log(Dplus_diag_station(id_Dplus_diag_station(i,pp))) + lambda_par_station(i,s)));
        nll_station_smooth_re(i,pp) -= dnorm(station_smooth_re(i,pp), Type(0.0), sd,1);
        SIMULATE station_smooth_re(i,pp) = rnorm(Type(0.0), sd);
        beta_re_use(i,k) = station_smooth_re(i,pp);
        pp++;
      }
    }
    nll += nll_station_smooth_re.sum();
    REPORT(nll_station_smooth_re);
  }
  //using namespace density;
  
	//distribution of population-level random effects for smoother portion
  vector<Type> nll_mean_smooth(mean_smooth_re.size());
  nll_mean_smooth.fill(0);
  for(int i = 0; i < beta_ns_ind.size(); i++) beta_use(beta_ns_ind(i)) = beta(i);
	if(use_mean_smooth_re == 1) 
  {
    int pp = 0;
	  if(nsmu > 0) for(int s = 0; s < nsmu; s++) for(int i = fp_mean_smooth_re(s); i <= lp_mean_smooth_re(s); i++)//for(int j = 0; j < mean_smooth_re.size(); j++) 
	  {
	    nll_mean_smooth(pp) -= dnorm(mean_smooth_re(pp), Type(0.0), exp(-0.5*(log(Dplus_diag(pp)) + lambda_par(s))), 1);
	    SIMULATE mean_smooth_re(pp) = rnorm(Type(0.0), exp(-0.5*(log(Dplus_diag(pp)) + lambda_par(s))));
      beta_use(i) = mean_smooth_re(pp);
      pp++;
	  }
    nll += nll_mean_smooth.sum();
    REPORT(nll_mean_smooth);
	}
	//distribution of random effects on fixed effects portion of mean (intercept and some coefficients of smoother)
	//series of conditionals.
  vector<Type> nll_beta_re(beta_re.rows());
  nll_beta_re.fill(0);
	if(use_beta_re == 1) 
  {
    if(use_beta_re_cor == 0) for(int i = 0; i <beta_re.rows(); i++) for(int j = 0; j < beta_re.cols(); j++) 
    {
      nll_beta_re(i) -= dnorm(beta_re(i,j),Type(0), exp(beta_re_var_pars(i,j)),1);
      SIMULATE beta_re(i,j) = rnorm(Type(0), exp(beta_re_var_pars(i,j)));
    }
    else
    {
      for(int i = 0; i <beta_re.rows(); i++)
      {
        vector<Type> beta_re_i_cor_par = beta_re_cor_pars.row(i);
        UNSTRUCTURED_CORR_t<Type> nll_beta_i(beta_re_i_cor_par); //diag = 1
        vector<Type> beta_re_i_log_sds = beta_re_var_pars.row(i);
        vector<Type> beta_re_i = beta_re.row(i);
        nll_beta_re(i) += VECSCALE(nll_beta_i,exp(beta_re_i_log_sds))(beta_re_i); //diag = sig^2
        SIMULATE
        {
          beta_re.row(i) = nll_beta_i.simulate() * exp(beta_re_i_log_sds);
        }
        matrix<Type> beta_re_cor = nll_beta_i.cov();
        REPORT(beta_re_cor);
      }
    }
    nll += nll_beta_re.sum();
    SIMULATE REPORT(beta_re);
    REPORT(nll_beta_re);
    for(int i = 0; i < beta_re_ns_ind.size(); i++) beta_re_use.col(beta_re_ns_ind(i)) = beta_re.col(i);
  }
   
  
  vector<Type> eta = offst + X * beta_use;
  vector<Type> nll_obs(N);
  nll_obs.fill(0);
  for(int j = 0; j < N; j++) 
  {
		vector<Type> beta_re_use_i = beta_re_use.row(id(j));
    vector<Type> Zj = Z.row(j);
    eta(j) +=  (Zj * beta_re_use_i).sum();
    mu(j) = Type(1.0)/(Type(1.0) + exp(-eta(j)));
    nll_obs(j) -= dbinom(y(j),n(j),mu(j),1);
    SIMULATE y(j) = rbinom(n(j), mu(j));
  }
  nll += nll_obs.sum();
  REPORT(nll_obs);
  SIMULATE 
  {
    REPORT(n);
    if(use_station_smooth_re == 1) REPORT(station_smooth_re);
    if(use_beta_re == 1) REPORT(beta_re);
    if(use_mean_smooth_re == 1) REPORT(mean_smooth_re);
  }
	REPORT(eta);
	REPORT(mu);
	if(do_pred)
	{
	  int n_pred = Xpred.rows();// = n_per_sta.sum();
    vector<Type> mean_pred_eta = Xpred * beta_use;  
	  matrix<Type> station_pred_eta(beta_re.rows(),n_pred);
	  for(int i = 0; i < beta_re.rows(); i++) 
    {
      vector<Type> beta_re_use_pred = beta_re_use.row(i);
      station_pred_eta.row(i) = mean_pred_eta;
      for(int j = 0; j < n_pred; j++)
      {
        vector<Type> Zj = Zpred.row(j);
		    station_pred_eta(i,j) += (Zj * beta_re_use_pred).sum();
      }
    }
    REPORT(station_pred_eta);
    REPORT(mean_pred_eta);
    REPORT(beta_re_use);
    REPORT(beta_use);
    ADREPORT(mean_pred_eta);
  }
	return(nll);
}

