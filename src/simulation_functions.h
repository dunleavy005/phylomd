#ifndef PHYLOMD_SIMULATION_FUNCTIONS_H_
#define PHYLOMD_SIMULATION_FUNCTIONS_H_

#include <string>
#include <utility>

#include <RcppArmadillo.h>
#include "phylomd_types.h"


//
// This header file declares the simulation functions.
//


Rcpp::DataFrame ctmc_sim_aux(double t, const Vector<std::string>& states,
                             const Vector<int>& state_inds, const arma::mat& Q,
                             int init_state_ind);


std::pair<Vector<std::string>, Vector<int>> asr_sim_aux(
    const arma::imat& edge, int num_edges, int num_term_nodes,
    int root_node_ind, int num_int_nodes, const arma::vec& edge_lengths,
    const Vector<std::string>& states, const Vector<int>& state_inds,
    const arma::mat& Q, const arma::vec& pi, const arma::ivec& tip_data);


Rcpp::DataFrame ctmc_sim(double t, const Rcpp::List& subst_mod,
                         const std::string& init_state);


Vector<std::string> asr_sim(const Rcpp::List& tree, const Rcpp::List& subst_mod,
                            const Vector<std::string>& tip_states);


Vector<Rcpp::DataFrame> smap_sim(const Rcpp::List& tree,
                                 const Rcpp::List& subst_mod,
                                 const Vector<std::string>& tip_states);


Vector<std::string> tips_sim(const Rcpp::List& tree,
                             const Rcpp::List& subst_mod);


#endif  // PHYLOMD_SIMULATION_FUNCTIONS_H_