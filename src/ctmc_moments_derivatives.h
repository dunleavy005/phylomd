#ifndef PHYLOMD_CTMC_MOMENTS_DERIVATIVES_H_
#define PHYLOMD_CTMC_MOMENTS_DERIVATIVES_H_

#include <string>

#include <RcppArmadillo.h>


//
// This header file declares the CTMC moment/derivative functions.
//


enum class Mode { NSUBS_MOMENTS, REWARD_MOMENTS, Q_DERIVATIVES, T_DERIVATIVES };


arma::cube ctmc_moments_Q_derivatives_aux(double t, const arma::mat& Q,
                                          const arma::mat& B, int max_order,
                                          Mode mode);


arma::cube ctmc_t_derivatives_aux(double t, const arma::mat& Q, int max_order);


arma::cube ctmc_moments_derivatives(double t, const arma::mat& Q,
                                    const arma::mat& B, int max_order,
                                    Mode mode);


arma::cube ctmc_nsubs_moments(double t, const Rcpp::List& subst_mod,
                              const arma::mat& L, int max_order);


arma::cube ctmc_reward_moments(double t, const Rcpp::List& subst_mod,
                               const arma::vec& w, int max_order);


arma::cube ctmc_Q_derivatives(double t, const Rcpp::List& subst_mod,
                              const std::string& param_name, int max_order);


arma::cube ctmc_t_derivatives(double t, const Rcpp::List& subst_mod,
                              int max_order);


#endif  // PHYLOMD_CTMC_MOMENTS_DERIVATIVES_H_