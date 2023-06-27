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

using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppDist)]]

//Generate data
arma::mat nct(int k);
arma::mat pct(double mu,int k,double theta,double tauS,double w);
arma::mat xct(arma::mat nct, arma::mat pct);

//functions to update parameters
List draw_phi(double mu_0,double sigmaS,double theta_0,double tauS,double w,arma::mat phi_iter,arma::mat nct,arma::mat xct);
double draw_sigma2(double mu_0,arma::mat phi_ct,double w,double a,double b);
double draw_tau2(double theta_0,arma::mat phi_ct,double a,double b);
double draw_theta_0(double sigmaS,double mu_0,double tauS,arma::mat phi_ct,double w,double lowerli,double upperli);
double draw_mu_0(double sigmaS,double theta_0,double tauS,arma::mat phi_ct,double w,double lowerli,double upperli);
double draw_w(double mu_0,double sigmaS,double theta_0,double tauS,arma::mat phi_ct);
double draw_w_continous(double mu_0,double sigmaS,arma::mat phi_ct);
List bnv(double mu_0,double sigmaS,double theta_0,double tauS,double w);

//initial parameters
rowvec SAA(arma::mat,arma::mat);
double SA(arma::mat,arma::mat);
double tauS0(arma::mat,arma::mat);
double sigmaS0(arma::mat,arma::mat);
double mu0(arma::mat,arma::mat);
rowvec mu00(arma::mat,arma::mat);

// [[Rcpp::export]]
List main_draw(int M, arma::mat xct, arma::mat nct){
  //double init_mu
  int i;
  int k = xct.n_cols;
  //set the starting points for all parameters.
  vec mu_0(M, fill::zeros);
  mu_0(0) = mu0(xct,nct);
  //
  vec sigmaS(M, fill::zeros);
  sigmaS(0) = sigmaS0(xct,nct);
  //
  vec theta_0(M, fill::zeros);
  theta_0(0) = SA(xct,nct);
  //theta_0(0) = init_mu;
  //
  vec tauS(M, fill::zeros);
  tauS(0) = tauS0(xct,nct);
  //
  vec w(M, fill::zeros);
  w(0) = 0;

  cube phi_test(2,k,M,fill::randu);

  List ss = bnv(mu_0(0),sigmaS(0),theta_0(0),tauS(0),w(0));
  arma::mat phi_iterst = rmvnorm(k,ss["mu_phi"],ss["Sigma_phi"]);
  phi_test.slice(0) = phi_iterst.t();

  //arma::mat phi_check(2,M, fill::zeros);
  //arma::mat lambda_check(2,M, fill::zeros);

  double Lu = min(mu00(xct,nct))-5;
  double Uu = max(mu00(xct,nct))+5;
  double Ltheta = min(SAA(xct,nct))-5;
  double Utheta = max(SAA(xct,nct))+5;

  for (i=0; i<M-1; i++){
    List phi_update = draw_phi(mu_0(i),sigmaS(i),
                               theta_0(i),tauS(i),
                               w(i),phi_test.slice(i),
                               nct,xct);
    arma::mat lambda = phi_update["lambda_ct"];
    arma::mat phi_new =  phi_update["phi_ct"];
    phi_test.slice(i+1) = phi_new;

    //phi_check(0,i+1) = lambda(0,1);
    //phi_check(1,i+1) = lambda(1,1);
    //lambda_check(0,i+1) = phi_new(0,1);
    //lambda_check(1,i+1) = phi_new(1,1);

    double mu0_update = draw_mu_0(sigmaS(i),theta_0(i),tauS(i),phi_new,w(i),Lu,Uu);
    mu_0(i+1) = mu0_update;

    double theta0_update = draw_theta_0(sigmaS(i),mu0_update,tauS(i),phi_new,w(i),Ltheta,Utheta);
    theta_0(i+1) = theta0_update;

    double sigmaS_update = draw_sigma2(mu0_update,phi_new,w(i),0.01,0.01);
    sigmaS(i+1) = sigmaS_update;

    double tauS_update = draw_tau2(theta0_update,phi_new,0.01,0.01);
    tauS(i+1) = tauS_update;

    double w_update = draw_w(mu0_update,sigmaS_update,theta0_update,tauS_update,phi_new);
    //double w_update = draw_w_continous(mu0_update,sigmaS_update,phi_new);
    w(i+1) = w_update;
  }


  return List::create(Rcpp::Named("theta0") = theta_0,
                      Rcpp::Named("mu0") = mu_0,
                      Rcpp::Named("sigmaS") = sigmaS,
                      Rcpp::Named("tauS") = tauS,
                      Rcpp::Named("w")=w,
                      Rcpp::Named("phi")=phi_test);

}



List draw_phi(double mu_0,double sigmaS,double theta_0,double tauS,double w,arma::mat phi_iter,arma::mat nct,arma::mat xct){
  double var_phi_1 = sigmaS+w*w*tauS;
  double var_phi_2 = sigmaS+(1-w)*(1-w)*tauS;
  double cov_phi_1_2 = sigmaS-w*(1-w)*tauS;
  double mu_phi_1 = mu_0-w*theta_0;
  double mu_phi_2 = mu_0+(1-w)*theta_0;

  arma::mat Sigma_phi = {{var_phi_1,cov_phi_1_2},{cov_phi_1_2,var_phi_2}};
  arma::mat mu_phi = {mu_phi_1,mu_phi_2};

  int k = nct.n_cols;
  int i;
  //generate two matrices
  arma::mat phi_ct(2,k);
  arma::mat lambda_ct(2,k);
  for(i=0; i<k; i++){

    //arma::vec phi_1_vec = {phi_iter(0,i)};
    //arma::vec phi_2_vec = {phi_iter(1,i)};
    double phi_1_vec = phi_iter(0,i);
    double phi_2_vec = phi_iter(1,i);
    // PolyaGamma obj(1000);

    //arma::vec lambda_1 = rcpp_pgdraw(nct(0,i),phi_1_vec);
    //arma::vec lambda_2 = rcpp_pgdraw(nct(1,i),phi_2_vec);
    // double lambda_1 = obj.draw_sum_of_gammas(nct(0,i),phi_1_vec);
    // double lambda_2 = obj.draw_sum_of_gammas(nct(1,i),phi_2_vec);
    double lambda_1 = rpg_hybrid(nct(0,i),phi_1_vec);
    double lambda_2 = rpg_hybrid(nct(1,i),phi_2_vec);

    arma::mat Lambda = {{lambda_1,0},{0,lambda_2}};
    arma::mat V_lambda_i = inv(Lambda+inv(Sigma_phi));

    arma::mat kappa_i  = {xct(0,i)-nct(0,i)*0.5,xct(1,i)-nct(1,i)*0.5};
    arma::mat m_i = V_lambda_i * (kappa_i.t() + (inv(Sigma_phi) * mu_phi.t()));
    arma::mat phi_update = rmvnorm(1,m_i,V_lambda_i);
    phi_ct.col(i) = phi_update.t();
    arma::mat lambda_temp ={lambda_1,lambda_2};
    lambda_ct.col(i) = lambda_temp.t();

  }
  return List::create(Rcpp::Named("phi_ct") = phi_ct,
                      Rcpp::Named("lambda_ct") = lambda_ct);

}
double draw_sigma2(double mu_0,arma::mat phi_ct,double w,double a=0.01,double b=0.01){
  rowvec phi_ct_1_mean = phi_ct.row(0);
  rowvec phi_ct_2_mean = phi_ct.row(1);
  int k = phi_ct.n_cols;

  double rate = b+sum(pow(w*(phi_ct_2_mean-phi_ct_1_mean)+phi_ct_1_mean-mu_0,2))*0.5;
  double scale = 1/rate;

  double sigma_temp = R::rgamma(a+k*0.5,scale);
  return(1/sigma_temp);
}
double draw_tau2(double theta_0,arma::mat phi_ct,double a=0.01,double b=0.01){
  rowvec phi_ct_1_mean = phi_ct.row(0);
  rowvec phi_ct_2_mean = phi_ct.row(1);
  int k = phi_ct.n_cols;

  double rate = b+sum(pow(phi_ct_1_mean-phi_ct_2_mean+theta_0,2))*0.5;
  double scale = 1/rate;

  double tau_temp =R::rgamma(a+k*0.5,scale);
  return(1/tau_temp);
}
double draw_theta_0(double sigmaS,double mu_0,double tauS,arma::mat phi_ct,double w,double lowerli,double upperli){

  double var_phi_1 = sigmaS+w*w*tauS;
  double var_phi_2 = sigmaS+(1-w)*(1-w)*tauS;
  double cov_phi_1_2 = sigmaS-w*(1-w)*tauS;

  int k = phi_ct.n_cols;

  rowvec phi_ct_1_mean = phi_ct.row(0);
  rowvec phi_ct_2_mean = phi_ct.row(1);

  double mean_theta_0 = (var_phi_1*(1-w)*(mean(phi_ct_2_mean)-mu_0)+
                         var_phi_2*w*(mu_0-mean(phi_ct_1_mean))+
                         cov_phi_1_2*(w*(mean(phi_ct_2_mean)-mu_0)+
                         (1-w)*(-mean(phi_ct_1_mean)+mu_0)))/sigmaS;
  double var_theta0 = tauS/k;
  double theta0_update = r_truncnorm(mean_theta_0,sqrt(var_theta0),lowerli,upperli);
  return(theta0_update);

}
double draw_mu_0(double sigmaS,double theta_0,double tauS,arma::mat phi_ct,double w,double lowerli,double upperli){
  double var_phi_1 = sigmaS+w*w*tauS;
  double var_phi_2 = sigmaS+(1-w)*(1-w)*tauS;
  double cov_phi_1_2 = sigmaS-w*(1-w)*tauS;

  int k = phi_ct.n_cols;

  rowvec phi_ct_1_mean = phi_ct.row(0);
  rowvec phi_ct_2_mean = phi_ct.row(1);

  double mean_mu_0 = (var_phi_2*(mean(phi_ct_1_mean)+w*theta_0)+
                      var_phi_1*(mean(phi_ct_2_mean)-(1-w)*theta_0)-
                      cov_phi_1_2*(mean(phi_ct_1_mean)+
                      mean(phi_ct_2_mean )+w*theta_0-(1-w)*theta_0))/(var_phi_2+var_phi_1-2*cov_phi_1_2);

  double var_mu0 = sigmaS/(k*(2*w*(1-w)+1));
  double mu0_update = r_truncnorm(mean_mu_0,sqrt(var_mu0),lowerli,upperli);

  return(mu0_update);

}
double draw_w_continous(double mu_0,double sigmaS,arma::mat phi_ct){
  rowvec phi_ct_1_mean = phi_ct.row(0);
  rowvec phi_ct_2_mean = phi_ct.row(1);

  double mean_w = sum((phi_ct_2_mean-phi_ct_1_mean)%(mu_0-phi_ct_1_mean))/sum(pow(phi_ct_2_mean-phi_ct_1_mean,2));
  double var_w =  sigmaS/sum(pow(phi_ct_2_mean-phi_ct_1_mean,2));
  double w_update = r_truncnorm(mean_w,sqrt(var_w),0,1);
  return w_update;
}
double draw_w(double mu_0,double sigmaS,double theta_0,double tauS,arma::mat phi_ct){
  int k = phi_ct.n_cols;

  List w0 = bnv(mu_0,sigmaS,theta_0,tauS,0);
  List wh = bnv(mu_0,sigmaS,theta_0,tauS,0.5);
  List w1 = bnv(mu_0,sigmaS,theta_0,tauS,1);

  vec w0v = dmvnorm(phi_ct.t(),w0["mu_phi"],w0["Sigma_phi"]);
  vec whv = dmvnorm(phi_ct.t(),wh["mu_phi"],wh["Sigma_phi"]);
  vec w1v = dmvnorm(phi_ct.t(),w1["mu_phi"],w1["Sigma_phi"]);

  vec w0v_temp = cumprod(500*w0v);
  vec whv_temp =cumprod(500*whv);
  vec w1v_temp = cumprod(500*w1v);

  double sum_weight = w0v_temp(k-1)+whv_temp(k-1)+w1v_temp(k-1);
  vec x = {0,0.5,1};
  vec prob0 = {w0v_temp(k-1)/sum_weight,
               whv_temp(k-1)/sum_weight,w1v_temp(k-1)/sum_weight};
  vec w_new = Rcpp::RcppArmadillo::sample(x,1,false,prob0);

  double w_update = w_new(0);
  return w_update;
}
List bnv(double mu_0,double sigmaS,double theta_0,double tauS,double w){
  double var_phi_1 = sigmaS+w*w*tauS;
  double var_phi_2 = sigmaS+(1-w)*(1-w)*tauS;
  double cov_phi_1_2 = sigmaS-w*(1-w)*tauS;
  double mu_phi_1 = mu_0-w*theta_0;
  double mu_phi_2 = mu_0+(1-w)*theta_0;

  arma::mat Sigma_phi ={{var_phi_1,cov_phi_1_2},{cov_phi_1_2,var_phi_2}};
  vec mu_phi = {mu_phi_1,mu_phi_2 };

  return List::create(Rcpp::Named("Sigma_phi") = Sigma_phi,
                      Rcpp::Named("mu_phi") = mu_phi);
}
arma::mat xct(arma::mat nct, arma::mat pct){
  int k,i;
  double xc,xt;
  k = nct.n_cols;
  arma::mat mat0(2,k);
  for(i=0; i<k; i++){
    xc = R::rbinom(nct(0,i),pct(0,i));
    xt = R::rbinom(nct(1,i),pct(1,i));
    mat0(0,i) = xc;
    mat0(1,i) = xt;
  }
  return mat0;
}
arma::mat pct(double mu,int k,double theta,double tauS,double w){
  double epsi1,epsi2,pc,pt;
  int i;
  arma::mat mat0(2,k);
  for(i=0; i<k; i++){
    epsi1 = R::rnorm(0,sqrt(0.5));
    epsi2 = R::rnorm(0,sqrt(tauS));
    pc = exp(mu+epsi1-w*(theta+epsi2))/(1+exp(mu+epsi1-w*(theta+epsi2)));
    pt = exp(mu+epsi1+(1-w)*(theta+epsi2))/(1+exp(mu+epsi1+(1-w)*(theta+epsi2)));
    mat0(0,i) = pc;
    mat0(1,i) = pt;
  }
  return mat0;
}
arma::mat nct(int k){
  double x,y;
  int i;
  //generate a matrix
  arma::mat mat0(2, k, fill::randu);
  for (i=0; i<k; i++){
    x = R::runif(50,1000);
    y = R::runif(50,1000);
    mat0(0,i) = round(x);
    mat0(1,i) = round(y);

  }
  return mat0;
}

// starting points for parameters.
rowvec SAA(arma::mat xct, arma::mat nct){
  rowvec xc =xct.row(0);
  rowvec nc = nct.row(0);
  rowvec xt =xct.row(1);
  rowvec nt = nct.row(1);

  return(log((xt+0.5)/(nt-xt+0.5))-log((xc+0.5)/(nc-xc+0.5)));
}


double SA(arma::mat xct, arma::mat nct){
  rowvec xc =xct.row(0);
  rowvec nc = nct.row(0);
  rowvec xt =xct.row(1);
  rowvec nt = nct.row(1);
  return (mean(log((xt+0.5)/(nt-xt+0.5))-log((xc+0.5)/(nc-xc+0.5))));
}


double tauS0(arma::mat xct, arma::mat nct){
  rowvec xc =xct.row(0);
  rowvec nc = nct.row(0);
  rowvec xt =xct.row(1);
  rowvec nt = nct.row(1);
  return(arma::var(log((xt+0.5)/(nt-xt+0.5))-log((xc+0.5)/(nc-xc+0.5))));
}

double sigmaS0(arma::mat xct, arma::mat nct){
  rowvec xc =xct.row(0);
  rowvec nc = nct.row(0);
  return(arma::stddev((log((xc+0.5)/(nc-xc+0.5)))));
}


rowvec mu00(arma::mat xct, arma::mat nct){
  rowvec xc =xct.row(0);
  rowvec nc = nct.row(0);
  return(log((xc+0.5)/(nc-xc+0.5)));
}

double mu0(arma::mat xct, arma::mat nct)
{
  rowvec xc =xct.row(0);
  rowvec nc = nct.row(0);
  return(mean(log((xc+0.5)/(nc-xc+0.5))));
}



