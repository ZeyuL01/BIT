// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppDist)]]
// [[Rcpp::depends(RcppProgress)]]

#include <RcppArmadillo.h>
#include "RcppArmadillo.h"
#include "polyagamma_wrapper.h"
#include <Rcpp.h>
#include <mvnorm.h>
#include <RcppArmadilloExtensions/sample.h>
#include "PolyaGammaApproxSP.h"
#include <progress.hpp>
#include <progress_bar.hpp>


using namespace Rcpp;
using namespace arma;

//Claim used matrix
arma::vec xct(int N);
arma::vec pct(int N);
arma::vec nct(int N);

//functions to assign M/Mc/Ji/label indicators
List auxiliary_list_generator(arma::vec tf_labels);
arma::mat label_mat(arma::vec tf_labels);
// After compile, this function will be immediately called using
// the below snippet and results will be sent to the R console.

//functions to draw parameters
double draw_tau0S(int I, double mu0, arma::vec theta_i, arma::rowvec unique_theta_i, const double a0 = 0.01, const double b0 = 0.01);
double draw_mu0(int I, double tau0S, arma::vec theta_i, arma::rowvec unique_theta_i);
double draw_tau1S(arma::rowvec Mc_indicator,arma::vec theta_ij,arma::vec theta_i,const double a1=0.01,const double b1=0.01);
arma::vec draw_sigmaS(double tau1S, arma::rowvec Ji,arma::mat label_mat,arma::vec theta_ij,
                      arma::vec theta_i,const double a2 = 0.01,const double b2 = 0.01);
arma::vec draw_theta_ij(arma::vec xct, arma::vec nct, arma::vec kappa_ij, arma::vec sigmaS, arma::vec tf_labels,
                        double tau1S, arma::vec lambda_ij, arma::vec theta_i, arma::rowvec Ji_counts);
arma::vec draw_theta_i(arma::vec theta_ij, double tau1S, arma::vec sigmaS, double mu0, double tau0S,
                       arma::mat label_mat,arma::rowvec Ji_counts);
arma::vec draw_lambda_ij(arma::vec nct, arma::vec theta_ij);

//initialize parameters
double tau0S_0(arma::vec theta_i, arma::rowvec unique_theta_i);
double mu0_0(arma::vec xct,arma::vec nct);
double tau1S_0(arma::vec xct, arma::vec nct, arma::rowvec Mc_indicator);
arma::vec theta_ij_0(arma::vec xct,arma::vec nct);
arma::vec theta_i_0(arma::vec theta_ij_0,arma::mat label_mat);
arma::vec lambda_ij_0(arma::vec theta_ij_0,arma::vec nct);
arma::vec sigmaS_0(arma::vec xct,arma::vec nct, arma::mat label_mat);

// [[Rcpp::export]]
List Main_Sampling(int N,arma::vec xct,arma::vec nct,arma::vec tf_labels, bool display_progress=true){
  int i;
  int I;
  arma::vec kappa_ij = xct - nct / 2;
  Progress p(N, display_progress);

  //generate auxiliary lists for the ease of computation.
  List auxiliary_lists = auxiliary_list_generator(tf_labels);
  arma::rowvec M_indicator = auxiliary_lists["M_indicator"];
  arma::rowvec Mc_indicator = auxiliary_lists["Mc_indicator"];
  arma::rowvec Ji_counts = auxiliary_lists["Ji_counts"];
  arma::mat label_mat = auxiliary_lists["Label_mat"];
  arma::rowvec unique_theta_i = auxiliary_lists["unique_theta_i"];

  //initialize parameters
  I = label_mat.n_rows;

  arma::vec tau1S_vec(N, fill::zeros);
  tau1S_vec(0) = tau1S_0(xct, nct, Mc_indicator);

  arma::mat theta_ij_mat(xct.n_rows,N,fill::zeros);
  theta_ij_mat.col(0) = theta_ij_0(xct, nct);

  arma::mat theta_i_mat(xct.n_rows,N,fill::zeros);
  theta_i_mat.col(0) = theta_i_0(theta_ij_mat.col(0), label_mat);

  arma::vec mu0_vec(N, fill::zeros);
  mu0_vec(0) = mu0_0(xct, nct);

  arma::vec tau0S_vec(N, fill::zeros);
  tau0S_vec(0) = tau0S_0(theta_i_mat.col(0), unique_theta_i);

  arma::mat sigmaS_mat(label_mat.n_rows,N,fill::zeros);
  sigmaS_mat.col(0) = sigmaS_0(xct, nct, label_mat);

  arma::mat lambda_ij_mat(xct.n_rows,N,fill::zeros);
  lambda_ij_mat.col(0) = lambda_ij_0(theta_ij_mat.col(0), nct);

  for(i=0;i<N-1;i++){
    p.increment();

    mu0_vec(i+1) = draw_mu0(I, tau0S_vec(i), theta_i_mat.col(i), unique_theta_i);

    tau0S_vec(i+1) = draw_tau0S(I, mu0_vec(i+1), theta_i_mat.col(i), unique_theta_i);

    tau1S_vec(i+1) = draw_tau1S(Mc_indicator, theta_ij_mat.col(i), theta_i_mat.col(i));

    sigmaS_mat.col(i+1) = draw_sigmaS(tau1S_vec(i+1), Ji_counts, label_mat, theta_ij_mat.col(i+1), theta_i_mat.col(i+1));

    theta_ij_mat.col(i+1) = draw_theta_ij(xct, nct, kappa_ij, sigmaS_mat.col(i), tf_labels,
                     tau1S_vec(i+1), lambda_ij_mat.col(i), theta_i_mat.col(i),Ji_counts);

    theta_i_mat.col(i+1) = draw_theta_i(theta_ij_mat.col(i+1), tau1S_vec(i+1), sigmaS_mat.col(i),
                    mu0_vec(i+1), tau0S_vec(i+1), label_mat, Ji_counts);

    lambda_ij_mat.col(i+1) = draw_lambda_ij(nct, theta_ij_mat.col(i+1));

  }

  return(List::create(Rcpp::Named("mu0") = mu0_vec,
                      Rcpp::Named("tau0S") = tau0S_vec,
                      Rcpp::Named("tau1S") = tau1S_vec,
                      Rcpp::Named("theta_ij") = theta_ij_mat,
                      Rcpp::Named("theta_i") = theta_i_mat,
                      Rcpp::Named("sigmaS") = sigmaS_mat,
                      Rcpp::Named("lambda_ij") = lambda_ij_mat,
                      Rcpp::Named("label_mat") = label_mat));
};

//functions to update parameters
double draw_tau0S(int I, double mu0, arma::vec theta_i, arma::rowvec unique_theta_i, const double a0, const double b0){
  double shape = a0 + I / 2;
  double scale = b0 + sum(unique_theta_i * pow((theta_i - mu0),2)) / 2;

  double tau0S_new = R::rgamma(shape,1 / scale);
  return(1 / tau0S_new);
};


double draw_mu0(int I, double tau0S, arma::vec theta_i, arma::rowvec unique_theta_i){

  double mean_new = sum(unique_theta_i * theta_i) / I;
  double sd_new = sqrt(tau0S / I);
  double mu0_new = R::rnorm(mean_new,sd_new);

  return(mu0_new);
};

double draw_tau1S(arma::rowvec Mc_indicator,arma::vec theta_ij,arma::vec theta_i, const double a1,const double b1){
  double shape = a1 + sum(Mc_indicator) / 2;
  double scale = b1 + sum(Mc_indicator * pow(theta_ij-theta_i,2)) / 2;

  double tau1S_new = R::rgamma(shape,1 / scale);

  return(1/tau1S_new);
};

arma::vec draw_sigmaS(double tau1S,arma::rowvec Ji, arma::mat label_mat, arma::vec theta_ij,
                      arma::vec theta_i, const double a2, const double b2){

  double a2_new = 0;
  double b2_new = 0;
  int i;

  arma::vec sigmaS_new(label_mat.n_rows,fill::zeros);
  arma::vec part = pow(theta_ij-theta_i,2);

  for(i=0;i<sigmaS_new.n_rows;i++){
    if(Ji(i)==1){
      sigmaS_new(i) = tau1S;
      continue;
    }
    a2_new = a2 + Ji(i) / 2;
    b2_new = b2 + sum(label_mat.row(i) * part) / 2;
    sigmaS_new(i) = 1/(R::rgamma(a2_new,1/b2_new));
  }

  return(sigmaS_new);
};

arma::vec draw_theta_ij(arma::vec xct, arma::vec nct, arma::vec kappa_ij, arma::vec sigmaS, arma::vec tf_labels,
                        double tau1S, arma::vec lambda_ij, arma::vec theta_i, arma::rowvec Ji_counts){
  int i;
  double V1 = 0;
  double m1 = 0;
  arma::vec theta_ij_new(xct.n_rows);

  for(i=0;i<theta_ij_new.n_rows;i++){
    V1 = 1 / (lambda_ij(i) + 1 / sigmaS(tf_labels(i)-1));
    m1 = V1 * (kappa_ij(i) + theta_i(i) / sigmaS(tf_labels(i)-1));
    theta_ij_new(i) = R::rnorm(m1,sqrt(V1));
  }

  return(theta_ij_new);
}

arma::vec draw_theta_i(arma::vec theta_ij, double tau1S, arma::vec sigmaS, double mu0, double tau0S,
                       arma::mat label_mat, arma::rowvec Ji_counts){
  int i;
  double mu_star;
  double tau_star;
  int Ji;
  arma::rowvec theta_i_new(theta_ij.n_rows,fill::zeros);

  for(i=0;i<label_mat.n_rows;i++){
    Ji = Ji_counts(i);
    mu_star = sum((mu0 * sigmaS(i) + sum(label_mat.row(i) * theta_ij * tau0S)) / (sigmaS(i) + Ji * tau0S));
    tau_star = 1 / (1 / tau0S + Ji / sigmaS(i));
    double theta_i_rd = R::rnorm(mu_star,sqrt(tau_star));
    theta_i_new = theta_i_new + label_mat.row(i) * theta_i_rd;
  }
  return(trans(theta_i_new));
};

arma::vec draw_lambda_ij(arma::vec nct, arma::vec theta_ij){
  arma::vec lambda_ij_new(theta_ij.n_rows);
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
  arma::rowvec unique_theta_i(tf_labels.n_rows,fill::zeros);

  int i;
  std::map<int,int> label_counts;
  std::map<int,int>::iterator iter;

  for(i=0;i<tf_labels.n_rows;i++){
    ++label_counts[tf_labels[i]];
    if(label_counts[tf_labels[i]]==1){
      unique_theta_i(i) = 1;
    }
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
                      Rcpp::Named("Label_mat") = label_m,
                      Rcpp::Named("unique_theta_i") = unique_theta_i));
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
double tau0S_0(arma::vec theta_i, arma::rowvec unique_theta_i){
  arma::uvec theta_i_indices = arma::find(unique_theta_i);
  return(arma::var(theta_i.elem(theta_i_indices)));
};

double mu0_0(arma::vec xct,arma::vec nct){
  return(mean(log((xct+0.1)/(nct-xct+0.1))));
};

double tau1S_0(arma::vec xct, arma::vec nct, arma::rowvec Mc_indicator){
  int i;
  arma::uvec theta_i1_indices = arma::find(Mc_indicator);
  arma::vec theta_i1_0(sum(Mc_indicator));
  double x;
  double n;

  for(i=0;i<theta_i1_0.n_rows;i++){
    x = xct(theta_i1_indices(i));
    n = nct(theta_i1_indices(i));
    theta_i1_0(i) = log((x+0.1)/(n-x+0.1));
  }

  return(arma::var(theta_i1_0));
};

arma::vec sigmaS_0(arma::vec xct,arma::vec nct, arma::mat label_mat){
  int i;
  int j;
  double x;
  double n;
  arma::vec sigmaS_0(label_mat.n_rows,fill::zeros);

  for(i=0;i<sigmaS_0.n_rows;i++){
    arma::uvec indices = find(label_mat.row(i));
    arma::vec temp_vals(indices.n_rows);

    for(j=0;j<indices.n_rows;j++){
      x = xct(indices(j));
      n = nct(indices(j));
      temp_vals(j) = log((x+0.1)/(n-x+0.1));
    }

    sigmaS_0(i) = arma::var(temp_vals);
    if(sigmaS_0(i)==0){
      sigmaS_0(i)=0.01;
    }
  }

  return(sigmaS_0);
};

arma::vec theta_ij_0(arma::vec xct,arma::vec nct){
  return(log((xct+0.1)/(nct-xct+0.1)));
};

arma::vec theta_i_0(arma::vec theta_ij_0,arma::mat label_mat){
  int i;
  arma::vec theta_i_0(theta_ij_0.n_rows,fill::zeros);

  for(i=0;i<label_mat.n_rows;i++){
    arma::mat theta_i_val = (label_mat.row(i) * theta_ij_0 ) / sum(label_mat.row(i));
    theta_i_0 = theta_i_0 + trans(label_mat.row(i) * theta_i_val(0,0));
  }
  return(theta_i_0);
};

arma::vec lambda_ij_0(arma::vec theta_ij_0,arma::vec nct){
  int i;
  arma::vec lambda_ij_0(theta_ij_0.n_rows,fill::zeros);

  for(i=0;i<lambda_ij_0.n_rows;i++){
    lambda_ij_0(i) = rpg_hybrid(nct(i),theta_ij_0(i));
  }

  return(lambda_ij_0);
};



