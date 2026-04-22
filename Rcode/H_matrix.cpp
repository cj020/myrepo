#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix propagator_1st_Rcpp(NumericMatrix neighbor,
                                  NumericVector diffu_rate,
                                  NumericVector growth_rate,
                                  NumericVector growth_rate2,
                                  double dx, double dy, double dt) {
  int n = neighbor.nrow();
  NumericMatrix H1(n, n);
  
  for(int i = 0; i < n; i++){
    if(sum(neighbor(i,_)>0)==4){
      H1(i,i) = 1-2*diffu_rate(i)*dt/pow(dx, 2.0)-2*diffu_rate(i)*dt/pow(dy, 2.0)+
        growth_rate(i)*dt;
      H1(i,neighbor(i,0)-1) = diffu_rate(i)*dt/pow(dx, 2.0);
      H1(i,neighbor(i,1)-1) = diffu_rate(i)*dt/pow(dx, 2.0);
      H1(i,neighbor(i,2)-1) = diffu_rate(i)*dt/pow(dy, 2.0);
      H1(i,neighbor(i,3)-1) = diffu_rate(i)*dt/pow(dy, 2.0);
    }
  }
  
  return H1;
}

// [[Rcpp::export]]
NumericMatrix propagator_2nd_Rcpp(NumericMatrix neighbor,
                                  NumericVector diffu_rate,
                                  NumericVector growth_rate,
                                  NumericVector growth_rate2,
                                  double dx, double dy, double dt){
  int n = neighbor.nrow();
  NumericMatrix H2(n, n);
  
  for(int i = 0; i < n; i++){
    if(sum(neighbor(i,_)>0)==4){
      H2(i,i) = -growth_rate2(i)*dt;
    }
  }
  
  return H2;
}

// [[Rcpp::export]]
NumericVector propagator_2nd_vec_Rcpp(NumericMatrix neighbor,
                                      NumericVector diffu_rate,
                                      NumericVector growth_rate,
                                      NumericVector growth_rate2,
                                      double dx, double dy, double dt){
  int n = neighbor.nrow();
  NumericVector H2_vec(n);
  
  for(int i = 0; i < n; i++){
    if(sum(neighbor(i,_)>0)==4){
      H2_vec(i) = -growth_rate2(i)*dt;
    }
  }
  
  return H2_vec;
}