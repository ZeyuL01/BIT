//#include "META.h"
#include "RcppArmadillo.h"
//#include "pgdraw.h"
#include "polyagamma_wrapper.h"
#include "PolyaGamma.h"
#include <Rcpp.h>
#include <mvnorm.h>
#include <truncnorm.h>
#include <RcppArmadilloExtensions/sample.h>
#include "PolyaGammaApproxSP.h"
#include <map>

using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppDist)]]

//Claim used matrix
arma::vec xct(int N);
arma::vec pct(int N);
arma::vec nct(int N);

//labels have to be imported
arma::vec tf_labels(int N);

//functions to assign M/Mc/Ji/label indicators
List auxiliary_list_generator(arma::vec tf_labels);
arma::mat label_mat(arma::vec tf_labels);


//functions to draw parameters
double draw_tau0S(int I,double mu0,arma::vec theta_i,const double a0,const double b0);
double draw_mu0(int I,double tau0S,arma::vec theta_i,double upperli,double lowerli);
double draw_tau1S(arma::vec Mc_indicator,arma::vec theta_ij,arma::vec theta_i,const double a1=0.01,const double b1=0.01);
arma::vec draw_sigmaS(arma::vec M_indicator,arma::vec theta_ij,arma::vec theta_i,const double a2,const double b2);
arma::vec draw_theta_ij(arma::vec xct,arma::vec nct,arma::rowvec M_indicator,arma::rowvec Mc_indicator,
                        arma::vec sigmaS,double tau1S,arma::vec lambda_ij,arma::vec theta_i);
arma::vec draw_theta_i(arma::vec theta_ij,arma::rowvec M_indicator,arma::rowvec Mc_indicator,
                       double tau1S,arma::vec sigmaS,double mu0,double tau0S);
arma::vec draw_lambda_ij(arma::vec nct, arma::vec theta_ij);

//initialize parameters
double tau0S_0(arma::vec xct,arma::vec nct);
double mu0_0(arma::vec xct,arma::vec nct);
double tau1S_0(arma::vec xct,arma::vec nct);
arma::vec theta_ij_0(arma::vec xct,arma::vec nct);
arma::vec theta_i_0(arma::vec xct,arma::vec nct);
arma::vec lambda_ij_0(arma::vec xct,arma::vec nct);
arma::vec sigmaS_0(arma::vec xct,arma::vec nct);


//main sampling function
// [[Rcpp::export]]
List Main_Sampling(int M,arma::vec xct,arma::vec nct,arma::vec tf_labels){
    int i;
    int N = xct.n_rows;
    arma::vec kappa_ij = xct - nct / 2;

    //generate auxiliary lists for the ease of computation.
    List auxiliary_lists = auxiliary_list_generator(tf_labels);
    arma::rowvec M_indicator = auxiliary_lists["M_indicator"];
    arma::rowvec Mc_indicator = auxiliary_lists["Mc_indicator"];
    arma::rowvec Ji_counts = auxiliary_lists["Ji_counts"];
    arma::mat label_mat = auxiliary_lists["Label_mat"];

    //initialize parameters



    return();
};


//functions to update parameters
double draw_tau0S(int I,double mu0,arma::vec theta_i,const double a0=0.01,const double b0=0.01){
  double rate = a0 + I / 2;
  double scale = b0 + sum(pow((theta_i - mu0),2)) / 2;

  double tau0S_new = R::rgamma(rate,scale);

  return(1/tau0S_new);
};


double draw_mu0(int I,double tau0S,arma::vec theta_i,double upperli,double lowerli){
  double mean = sum(theta_i) / I;
  double var = tau0S / I;

  double mu0_new = r_truncnorm(mean,var,lowerli,upperli);

  return(mu0_new);
};

double draw_tau1S(arma::rowvec Mc_indicator,arma::vec theta_ij,arma::vec theta_i,const double a1=0.01,const double b1=0.01){
  double rate = a1 + sum(Mc_indicator) / 2;
  double scale = b1 + sum(Mc_indicator%pow(theta_ij-theta_i,2)) / 2;

  double tau1S_new = R::rgamma(rate,scale);

  return(1/tau1S_new);
};

arma::rowvec draw_sigmaS(arma::rowvec Ji,arma::mat label_mat,arma::vec theta_ij,
                      arma::vec theta_i,const double a2,const double b2){
  int i;

  arma::rowvec sigmaS_new(label_mat.n_cols,fill::zeros);
  arma::rowvec part = pow(theta_ij-theta_i,2);

  for(i=0;i<sigmaS_new.n_cols;i++){
    double a2_new = a2 + Ji(i) / 2;
    double b2_new = b2 + sum(label_mat.row(i) % part) / 2;
    sigmaS_new(i) = 1/(R::rgamma(a2_new,b2_new));
  }

  return(sigmaS_new);
};

arma::rowvec draw_theta_ij(arma::vec xct,arma::vec nct,arma::vec kappa_ij,arma::rowvec M_indicator,arma::rowvec Mc_indicator,
                        arma::vec sigmaS,arma::vec tf_labels,double tau1S,arma::vec lambda_ij,arma::vec theta_i){
  int i;
  double V1 = 0;
  double m1 = 0;
  double V2 = 0;
  double m2 = 0;

  double sigmaS_used = 0;
  arma::rowvec theta_ij_new(xct.n_rows);
  arma::rowvec theta_i1_new(xct.n_rows);

  for(i=0;i<theta_ij_new.n_cols;i++){
    sigmaS_used = sigmaS(tf_labels(i)-1);
    V1 = 1 / (lambda_ij(i) + 1 / sigmaS_used);
    m1 = V1 * (kappa_ij(i) + theta_i(i)/sigmaS_used);
    theta_ij_new(i) = R::rnorm(m1,V1);

    V2 = 1 / (lambda_ij(i) + 1 / tau1S);
    m2 = V2 * (kappa_ij(i) + theta_i(i)/tau1S);
    theta_i1_new(i) = R::rnorm(m1,V1);
  }

  return((M_indicator % theta_ij_new)+(Mc_indicator % theta_i1_new));
}

arma::vec draw_theta_i(arma::vec theta_ij,arma::rowvec M_indicator,arma::rowvec Mc_indicator,
                       double tau1S,arma::vec sigmaS,double mu0,double tau0S,arma::mat label_mat,arma::vec Ji_counts){
  int i;
  arma::vec theta_i_new(theta_ij.n_cols,fill::zeros);

  for(i=0;i<label_mat.n_rows;i++){
    double mu_star;
    double tau_star;

    if(Ji_counts(i)==1){
      arma::mat mu_star_mat = (mu0 * tau1S + label_mat.row(i) * theta_ij * tau0S) / (tau0S + tau1S);
      double mu_star = mu_star_mat(0,0);
      double tau_star = 1 / (1 / tau1S + 1 / tau0S);
    }else{
      arma::mat mu_star_mat = (mu0 * sigmaS(i) + label_mat.row(i) * theta_ij * tau0S) / (sigmaS(i) + Ji_counts(i) * tau0S);
      double mu_star = mu_star_mat(0,0);
      double tau_star = 1 / (1 / tau0S + Ji_counts(i) / sigmaS(i));
    }
    double theta_i_rd = R::rnorm(mu_star,tau_star);
    theta_i_new = theta_i_new + label_mat.row(i) * theta_i_rd;
  }
  return(theta_i_new);
};

arma::vec draw_lambda_ij(arma::vec nct, arma::vec theta_ij){
  arma::vec lambda_ij_new(theta_ij.n_cols);
  int i;

  for(i=0;i<lambda_ij_new.n_rows;i++){
    double theta_ij_val = theta_ij(i);
    double nct_val = nct(i);
    lambda_ij_new(i) = rpg_hybrid(nct_val,theta_ij_val);
  }

  return(lambda_ij_new);
};

//function to indicate whether the TF has repeated ChIP-seq datasets.
//also counts the occurrence time of ChIP-seq datasets of same TF.
//generate a label mat for the ease of computation later.

List auxiliary_list_generator(arma::vec tf_labels){
  arma::rowvec M(tf_labels.n_rows,fill::zeros);
  arma::rowvec Mc(tf_labels.n_rows,fill::zeros);

  int i;
  std::map<int,int> label_counts;
  std::map<int,int>::iterator iter;

  for(i=0;i<tf_labels.n_rows;i++){
    ++label_counts[tf_labels[i]];
  }

  for(i=0;i<tf_labels.n_rows;i++){
    if(label_counts[tf_labels(i)]>1){
      M(i) = 1;
    }else{
      Mc(i) = 1;
    }
  }

  arma::rowvec Ji(label_counts.size(),fill::zeros);
  for(iter=label_counts.begin();iter!=label_counts.end();iter++){
    Ji((iter -> first) - 1) = iter -> second;
  }

  arma::mat label_m = label_mat(tf_labels);

  return(List::create(Rcpp::Named("M_indicator") = M,
                      Rcpp::Named("Mc_indicator") = Mc,
                      Rcpp::Named("Ji_counts") = Ji,
                      Rcpp::Named("Label_mat") = label_m));
};

//function transform labels to a matrix of label indicator, for the ease of computation.
//with index of row corresponds to TF.
//and column index corresponds to chip-seq dataset.
arma::mat label_mat(arma::vec tf_labels){
  int i;
  std::map<int,int> label_counts;
  std::map<int,int>::iterator iter;

  for(i=0;i<tf_labels.n_rows;i++){
    ++label_counts[tf_labels[i]];
  }

  arma::mat label_mat(label_counts.size(),tf_labels.n_rows,fill::zeros);

  for(i=0;i<tf_labels.n_rows;i++){
    label_mat(tf_labels(i)-1,i)=1;
  }

  return(label_mat);
}

//functions to initialize parameters.
double tau0S_0(arma::vec xct,arma::vec nct){
  return(arma::var(log((xct+0.5)/(nct-xct+0.5))));
};

double mu0_0(arma::vec xct,arma::vec nct){
  return(mean(log((xct+0.5)/(nct-xct+0.5))));
};

double tau1S_0(arma::vec xct, arma::vec nct, arma::vec Mc_indicator){
  int i;
  arma::uvec theta_i1_indices = arma::find(Mc_indicator);
  arma::vec theta_i1_0(sum(Mc_indicator));
  double x;
  double n;

  for(i=0;i<theta_i1_0.n_rows;i++){
    x = xct(theta_i1_indices(i));
    n = nct(theta_i1_indices(i));
    theta_i1_0(i) = log((x+0.5)/(n-x+0.5));
  }

  return(arma::var(theta_i1_0));
};

arma::vec sigmaS_0(arma::vec xct,arma::vec nct, arma::mat label_mat){
  int i;
  double x;
  double n;
  arma::vec sigmaS_0(label_mat.n_cols,fill::zeros);

  for(i=0;i<sigmaS_0.n_rows;i++){
    int j;
    arma::uvec indices = find(label_mat.row(i));
    arma::vec temp_vals(indices.n_rows);

    for(j=0;j<indices.n_rows;j++){
      x = xct(indices(j));
      n = nct(indices(j));
      temp_vals(j) = log((x+0.5)/(n-x+0.5));
    }

    sigmaS_0(i) = arma::var(temp_vals);
  }

  return(sigmaS_0);
};

arma::vec theta_ij_0(arma::vec xct,arma::vec nct){
  return(log((xct+0.5)/(nct-xct+0.5)));
};

arma::vec theta_i_0(arma::vec theta_ij_0,arma::mat label_mat){
  int i;
  arma::vec theta_i_0(theta_ij_0.n_rows,fill::zeros);

  for(i=0;i<label_mat.n_cols;i++){
    arma::mat theta_i_val = (theta_ij_0 * label_mat.row(i)) / sum(label_mat.row(i));
    theta_i_0 = theta_i_0 + trans(label_mat.row(i) * theta_i_val(0,0));
  }
  return(theta_i_0);
};

arma::vec lambda_ij_0(arma::vec xct,arma::vec nct);

//For test only, delete later.
// [[Rcpp::export]]



