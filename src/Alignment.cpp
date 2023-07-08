//For quick alignment between input peaks and reference datasets.

#include <iostream>
#include <cmath>
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
List Alignment(NumericVector input_vec, NumericVector ref_vec){
  int i = 0;
  int j = 0;
  int m = 0;

  int Xi = 0;
  int Ni = 0;

  while(i<input_vec.length() & j<ref_vec.length()){
    if(m==0){
      if(ref_vec[j]==input_vec[i]){
        ++Xi;
      }else if(ref_vec[j]>input_vec[i]){
        m = 1;
        ++i;
        continue;
      }
      ++j;
    }

    if(m==1){
      if(input_vec[i]==ref_vec[j]){
        ++Xi;
      }else if(input_vec[i]>ref_vec[j]){
        m = 0;
        ++j;
        continue;
      }
      ++i;
    }
  }

  Ni = input_vec.length() + ref_vec.length() - Xi;

  return(List::create(Rcpp::Named("Xi") = Xi,
                      Rcpp::Named("Ni") = Ni));
}
