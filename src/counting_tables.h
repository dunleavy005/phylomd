#ifndef PHYLOMD_COUNTING_TABLES_H_
#define PHYLOMD_COUNTING_TABLES_H_

#include <RcppArmadillo.h>


//
// This header file declares the counting table functions.
//


arma::mat choose_table(int n);


arma::vec factorial_table(int n);


arma::mat stirling_num_table(int n);


#endif  // PHYLOMD_COUNTING_TABLES_H_