#ifndef PHYLOMD_PHYLO_MOMENTS_DERIVATIVES_H_
#define PHYLOMD_PHYLO_MOMENTS_DERIVATIVES_H_

#include <string>

#include <RcppArmadillo.h>
#include "EdgeSet.h"
#include "ListID.h"
#include "MomentDerivativeID.h"
#include "ctmc_moments_derivatives.h"
#include "phylomd_types.h"


//
// This header file declares the phylogenetic moment/derivative functions.
//


void find_connected_moment_derivative_ids_aux(
    const FlatListID& flat_list_id, const Vector<EdgeSet>& esets,
    FlatListID::const_iterator curr_it, const Map<int, int>& eset_labels_orders,
    Vector<MomentDerivativeID>& md_ids);


Vector<MomentDerivativeID> find_connected_moment_derivative_ids(
    const FlatListID& flat_list_id, const Vector<EdgeSet>& esets);


int get_connected_counting_coef(const FlatMomentDerivativeID& flat_md_id,
                                FlatListID& flat_list_id);


Map<std::string, double> phylo_moments_derivatives(
    arma::imat& edge, const Vector<std::string>& tip_labels, int num_int_nodes,
    const arma::vec& edge_lengths, const arma::mat& Q, const arma::mat& B,
    const arma::vec& pi, VectorVector<int>& esets_inp, int max_order, Mode mode,
    const arma::ivec& tip_data);


double phylo_likelihood(const Rcpp::List& tree, const Rcpp::List& subst_mod,
                        const Vector<std::string>& tip_states);


Map<std::string, double> phylo_nsubs_moments(
    const Rcpp::List& tree, const Rcpp::List& subst_mod, const arma::mat& L,
    VectorVector<int> edge_sets, int max_order,
    const Vector<std::string>& tip_states);


Map<std::string, double> phylo_reward_moments(
    const Rcpp::List& tree, const Rcpp::List& subst_mod, const arma::vec& w,
    VectorVector<int> edge_sets, int max_order,
    const Vector<std::string>& tip_states);


Map<std::string, double> phylo_Q_derivatives(
    const Rcpp::List& tree, const Rcpp::List& subst_mod,
    const std::string& param_name, int max_order,
    const Vector<std::string>& tip_states);


Map<std::string, double> phylo_t_derivatives(
    const Rcpp::List& tree, const Rcpp::List& subst_mod, int max_order,
    const Vector<std::string>& tip_states);


#endif  // PHYLOMD_PHYLO_MOMENTS_DERIVATIVES_H_