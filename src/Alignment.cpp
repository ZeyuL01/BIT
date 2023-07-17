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

  int Xi_GOOD = 0;
  int Ni_TOTAL = 0;
  int Xi_BAD;

  while(i<input_vec.length() & j<ref_vec.length()){
    if(m==0){
      if(ref_vec[j]==input_vec[i]){
        ++Xi_GOOD;
      }else if(ref_vec[j] > input_vec[i]){
        m = 1;
        ++i;
        continue;
      }
      ++j;
    }

    if(m==1){
      if(input_vec[i]==ref_vec[j]){
        ++Xi_GOOD;
      }else if(input_vec[i] > ref_vec[j]){
        m = 0;
        ++j;
        continue;
      }
      ++i;
    }
  }

  Ni_TOTAL = input_vec.length() + ref_vec.length() - Xi_GOOD;
  Xi_BAD = input_vec.length() - Xi_GOOD;

  return(List::create(Rcpp::Named("Xi_GOOD") = Xi_GOOD,
                      Rcpp::Named("Xi_BAD") = Xi_BAD,
                      Rcpp::Named("Ni_TOTAL") = Ni_TOTAL));
}
