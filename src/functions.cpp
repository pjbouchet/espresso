#include <algorithm>
#include <Rcpp.h>
#include <iostream>
#include <vector>
#include <math.h>       /* isinf, sqrt */
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector dt_norm(NumericVector x,
                      NumericVector location,
                      NumericVector scale,
                      NumericVector L,
                      NumericVector U,
                      bool do_log) {

  int n = std::max({x.size(), location.size(), scale.size(), L.size(), U.size()});
  NumericVector X = rep_len(x, n);
  NumericVector LOC = rep_len(location, n);
  NumericVector SCALE = rep_len(scale, n);
  NumericVector LOW = rep_len(L, n);
  NumericVector UP = rep_len(U, n);

  NumericVector out(n);
 
  for (int i=0; i<n; i++) {
    if (X[i] < LOW[i]) out[i] = -100000;
    if (X[i] > UP[i]) out[i] = -100000;
    if (X[i] > LOW[i] & X[i] < UP[i]) out[i] = R::dnorm(X[i], LOC[i], SCALE[i], TRUE) -
      (log(R::pnorm(UP[i], LOC[i], SCALE[i], TRUE, FALSE) - R::pnorm(LOW[i], LOC[i], SCALE[i], TRUE, FALSE)));
    if(isinf(out[i])) out[i] = -100000;
  }
  if(!do_log){out = exp(out);}
  return out;
}


// [[Rcpp::export]]
NumericVector rt_norm(int n,
                      NumericVector location,
                      NumericVector scale,
                      NumericVector L,
                      NumericVector U) {
  
  int N = std::max({location.size(), scale.size(), L.size(), U.size()});
  if (n > N) N = n;
  
  NumericVector LOC = rep_len(location, N);
  NumericVector SCALE = rep_len(scale, N);
  NumericVector LOW = rep_len(L, N);
  NumericVector UP = rep_len(U, N);
  
  NumericVector out(N);
  
  for (int i=0; i<N; i++) {
    auto tot = R::pnorm(LOW[i], LOC[i], SCALE[i], TRUE, FALSE) + R::runif(0,1) * (R::pnorm(UP[i], LOC[i], SCALE[i], TRUE, FALSE) - R::pnorm(LOW[i], LOC[i], SCALE[i], TRUE, FALSE));
    out[i] = LOC[i] + SCALE[i] * R::qnorm(tot, 0, 1, TRUE, FALSE);
  }
  
  return out;
}

// [[Rcpp::export]]
NumericVector d_Binom(NumericVector x,
                      int size,
                      NumericVector prob,
                      bool do_log) {
  
  int n = std::max({x.size(), prob.size()});
  NumericVector X = rep_len(x, n);
  NumericVector PROB = rep_len(prob, n);
  
  NumericVector out(n);
  
  for (int i=0; i<n; i++) {
    out[i] = R::dbinom(X[i], size, PROB[i], do_log);
    if(isinf(out[i])) out[i] = -100000;
  }
  return out;
}

// Function unique_cpp("unique");
// 
// // [[Rcpp::export]]
// List propose_alpha(List rjobj, double sd = 1) {
//   
//   // Extract elements (Reminder: indices start at 0)
//   Rcpp::List MLIST = rjobj["mlist"];
//   String CURMOD = rjobj["current.model"];
//   IntegerVector ITER = rjobj["iter"];
//   NumericMatrix ALPHA = rjobj["alpha"];
//   
//   IntegerVector mod_index = MLIST[CURMOD];
//   int i = ITER("alpha") - 1; 
// 
//   // Unique values of alpha
//   NumericVector m = unique_cpp(ALPHA(i, _));
//   
//   NumericVector out(m.size());
// 
//   // Propose values
//   for (int j=0; j<m.size(); j++) {
//   out[j] = R::rnorm(m[j], sd);
//   }
// 
//   return List::create(Named("alpha") = out[mod_index - 1]);
// }