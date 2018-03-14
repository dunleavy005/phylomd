#ifndef PHYLOMD_SUBSTITUTION_MODELS_H_
#define PHYLOMD_SUBSTITUTION_MODELS_H_

#include <RcppArmadillo.h>


//
// This header file declares the CTMC substitution model functions.
//


Rcpp::List JC69(double mu, bool scale);


Rcpp::List K80(double alpha, double beta, bool scale);


Rcpp::List F81(double mu, const arma::vec& pi, bool scale);


Rcpp::List HKY85(double alpha, double beta, const arma::vec& pi, bool scale);


Rcpp::List GTR(double rAC, double rAG, double rAT, double rCG, double rCT,
               double rGT, const arma::vec& pi, bool scale);


#endif  // PHYLOMD_SUBSTITUTION_MODELS_H_