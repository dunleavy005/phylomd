// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "phylomd_types.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// ctmc_tpm
arma::mat ctmc_tpm(double t, const Rcpp::List& subst_mod);
RcppExport SEXP _phylomd_ctmc_tpm(SEXP tSEXP, SEXP subst_modSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type t(tSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type subst_mod(subst_modSEXP);
    rcpp_result_gen = Rcpp::wrap(ctmc_tpm(t, subst_mod));
    return rcpp_result_gen;
END_RCPP
}
// ctmc_nsubs_moments
arma::cube ctmc_nsubs_moments(double t, const Rcpp::List& subst_mod, const arma::mat& L, int max_order);
RcppExport SEXP _phylomd_ctmc_nsubs_moments(SEXP tSEXP, SEXP subst_modSEXP, SEXP LSEXP, SEXP max_orderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type t(tSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type subst_mod(subst_modSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type L(LSEXP);
    Rcpp::traits::input_parameter< int >::type max_order(max_orderSEXP);
    rcpp_result_gen = Rcpp::wrap(ctmc_nsubs_moments(t, subst_mod, L, max_order));
    return rcpp_result_gen;
END_RCPP
}
// ctmc_reward_moments
arma::cube ctmc_reward_moments(double t, const Rcpp::List& subst_mod, const arma::vec& w, int max_order);
RcppExport SEXP _phylomd_ctmc_reward_moments(SEXP tSEXP, SEXP subst_modSEXP, SEXP wSEXP, SEXP max_orderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type t(tSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type subst_mod(subst_modSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< int >::type max_order(max_orderSEXP);
    rcpp_result_gen = Rcpp::wrap(ctmc_reward_moments(t, subst_mod, w, max_order));
    return rcpp_result_gen;
END_RCPP
}
// ctmc_Q_derivatives
arma::cube ctmc_Q_derivatives(double t, const Rcpp::List& subst_mod, const std::string& param_name, int max_order);
RcppExport SEXP _phylomd_ctmc_Q_derivatives(SEXP tSEXP, SEXP subst_modSEXP, SEXP param_nameSEXP, SEXP max_orderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type t(tSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type subst_mod(subst_modSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type param_name(param_nameSEXP);
    Rcpp::traits::input_parameter< int >::type max_order(max_orderSEXP);
    rcpp_result_gen = Rcpp::wrap(ctmc_Q_derivatives(t, subst_mod, param_name, max_order));
    return rcpp_result_gen;
END_RCPP
}
// ctmc_t_derivatives
arma::cube ctmc_t_derivatives(double t, const Rcpp::List& subst_mod, int max_order);
RcppExport SEXP _phylomd_ctmc_t_derivatives(SEXP tSEXP, SEXP subst_modSEXP, SEXP max_orderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type t(tSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type subst_mod(subst_modSEXP);
    Rcpp::traits::input_parameter< int >::type max_order(max_orderSEXP);
    rcpp_result_gen = Rcpp::wrap(ctmc_t_derivatives(t, subst_mod, max_order));
    return rcpp_result_gen;
END_RCPP
}
// phylo_likelihood
double phylo_likelihood(const Rcpp::List& tree, const Rcpp::List& subst_mod, const std::vector<std::string>& tip_states);
RcppExport SEXP _phylomd_phylo_likelihood(SEXP treeSEXP, SEXP subst_modSEXP, SEXP tip_statesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type tree(treeSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type subst_mod(subst_modSEXP);
    Rcpp::traits::input_parameter< const std::vector<std::string>& >::type tip_states(tip_statesSEXP);
    rcpp_result_gen = Rcpp::wrap(phylo_likelihood(tree, subst_mod, tip_states));
    return rcpp_result_gen;
END_RCPP
}
// phylo_nsubs_moments
Map<std::string, double> phylo_nsubs_moments(const Rcpp::List& tree, const Rcpp::List& subst_mod, const arma::mat& L, VectorVector<int> edge_sets, int max_order, const std::vector<std::string>& tip_states);
RcppExport SEXP _phylomd_phylo_nsubs_moments(SEXP treeSEXP, SEXP subst_modSEXP, SEXP LSEXP, SEXP edge_setsSEXP, SEXP max_orderSEXP, SEXP tip_statesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type tree(treeSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type subst_mod(subst_modSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type L(LSEXP);
    Rcpp::traits::input_parameter< VectorVector<int> >::type edge_sets(edge_setsSEXP);
    Rcpp::traits::input_parameter< int >::type max_order(max_orderSEXP);
    Rcpp::traits::input_parameter< const std::vector<std::string>& >::type tip_states(tip_statesSEXP);
    rcpp_result_gen = Rcpp::wrap(phylo_nsubs_moments(tree, subst_mod, L, edge_sets, max_order, tip_states));
    return rcpp_result_gen;
END_RCPP
}
// phylo_reward_moments
Map<std::string, double> phylo_reward_moments(const Rcpp::List& tree, const Rcpp::List& subst_mod, const arma::vec& w, VectorVector<int> edge_sets, int max_order, const std::vector<std::string>& tip_states);
RcppExport SEXP _phylomd_phylo_reward_moments(SEXP treeSEXP, SEXP subst_modSEXP, SEXP wSEXP, SEXP edge_setsSEXP, SEXP max_orderSEXP, SEXP tip_statesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type tree(treeSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type subst_mod(subst_modSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< VectorVector<int> >::type edge_sets(edge_setsSEXP);
    Rcpp::traits::input_parameter< int >::type max_order(max_orderSEXP);
    Rcpp::traits::input_parameter< const std::vector<std::string>& >::type tip_states(tip_statesSEXP);
    rcpp_result_gen = Rcpp::wrap(phylo_reward_moments(tree, subst_mod, w, edge_sets, max_order, tip_states));
    return rcpp_result_gen;
END_RCPP
}
// phylo_Q_derivatives
Map<std::string, double> phylo_Q_derivatives(const Rcpp::List& tree, const Rcpp::List& subst_mod, const std::string& param_name, int max_order, const std::vector<std::string>& tip_states);
RcppExport SEXP _phylomd_phylo_Q_derivatives(SEXP treeSEXP, SEXP subst_modSEXP, SEXP param_nameSEXP, SEXP max_orderSEXP, SEXP tip_statesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type tree(treeSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type subst_mod(subst_modSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type param_name(param_nameSEXP);
    Rcpp::traits::input_parameter< int >::type max_order(max_orderSEXP);
    Rcpp::traits::input_parameter< const std::vector<std::string>& >::type tip_states(tip_statesSEXP);
    rcpp_result_gen = Rcpp::wrap(phylo_Q_derivatives(tree, subst_mod, param_name, max_order, tip_states));
    return rcpp_result_gen;
END_RCPP
}
// phylo_t_derivatives
Map<std::string, double> phylo_t_derivatives(const Rcpp::List& tree, const Rcpp::List& subst_mod, int max_order, const std::vector<std::string>& tip_states);
RcppExport SEXP _phylomd_phylo_t_derivatives(SEXP treeSEXP, SEXP subst_modSEXP, SEXP max_orderSEXP, SEXP tip_statesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type tree(treeSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type subst_mod(subst_modSEXP);
    Rcpp::traits::input_parameter< int >::type max_order(max_orderSEXP);
    Rcpp::traits::input_parameter< const std::vector<std::string>& >::type tip_states(tip_statesSEXP);
    rcpp_result_gen = Rcpp::wrap(phylo_t_derivatives(tree, subst_mod, max_order, tip_states));
    return rcpp_result_gen;
END_RCPP
}
// ctmc_sim
Rcpp::DataFrame ctmc_sim(double t, const Rcpp::List& subst_mod, const std::string& init_state);
RcppExport SEXP _phylomd_ctmc_sim(SEXP tSEXP, SEXP subst_modSEXP, SEXP init_stateSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type t(tSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type subst_mod(subst_modSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type init_state(init_stateSEXP);
    rcpp_result_gen = Rcpp::wrap(ctmc_sim(t, subst_mod, init_state));
    return rcpp_result_gen;
END_RCPP
}
// asr_sim
std::vector<std::string> asr_sim(const Rcpp::List& tree, const Rcpp::List& subst_mod, const std::vector<std::string>& tip_states);
RcppExport SEXP _phylomd_asr_sim(SEXP treeSEXP, SEXP subst_modSEXP, SEXP tip_statesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type tree(treeSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type subst_mod(subst_modSEXP);
    Rcpp::traits::input_parameter< const std::vector<std::string>& >::type tip_states(tip_statesSEXP);
    rcpp_result_gen = Rcpp::wrap(asr_sim(tree, subst_mod, tip_states));
    return rcpp_result_gen;
END_RCPP
}
// smap_sim
std::vector<Rcpp::DataFrame> smap_sim(const Rcpp::List& tree, const Rcpp::List& subst_mod, const std::vector<std::string>& tip_states);
RcppExport SEXP _phylomd_smap_sim(SEXP treeSEXP, SEXP subst_modSEXP, SEXP tip_statesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type tree(treeSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type subst_mod(subst_modSEXP);
    Rcpp::traits::input_parameter< const std::vector<std::string>& >::type tip_states(tip_statesSEXP);
    rcpp_result_gen = Rcpp::wrap(smap_sim(tree, subst_mod, tip_states));
    return rcpp_result_gen;
END_RCPP
}
// tips_sim
std::vector<std::string> tips_sim(const Rcpp::List& tree, const Rcpp::List& subst_mod);
RcppExport SEXP _phylomd_tips_sim(SEXP treeSEXP, SEXP subst_modSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type tree(treeSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type subst_mod(subst_modSEXP);
    rcpp_result_gen = Rcpp::wrap(tips_sim(tree, subst_mod));
    return rcpp_result_gen;
END_RCPP
}
// JC69
Rcpp::List JC69(double mu, bool scale);
RcppExport SEXP _phylomd_JC69(SEXP muSEXP, SEXP scaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< bool >::type scale(scaleSEXP);
    rcpp_result_gen = Rcpp::wrap(JC69(mu, scale));
    return rcpp_result_gen;
END_RCPP
}
// K80
Rcpp::List K80(double alpha, double beta, bool scale);
RcppExport SEXP _phylomd_K80(SEXP alphaSEXP, SEXP betaSEXP, SEXP scaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< bool >::type scale(scaleSEXP);
    rcpp_result_gen = Rcpp::wrap(K80(alpha, beta, scale));
    return rcpp_result_gen;
END_RCPP
}
// F81
Rcpp::List F81(double mu, const arma::vec& pi, bool scale);
RcppExport SEXP _phylomd_F81(SEXP muSEXP, SEXP piSEXP, SEXP scaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type pi(piSEXP);
    Rcpp::traits::input_parameter< bool >::type scale(scaleSEXP);
    rcpp_result_gen = Rcpp::wrap(F81(mu, pi, scale));
    return rcpp_result_gen;
END_RCPP
}
// HKY85
Rcpp::List HKY85(double alpha, double beta, const arma::vec& pi, bool scale);
RcppExport SEXP _phylomd_HKY85(SEXP alphaSEXP, SEXP betaSEXP, SEXP piSEXP, SEXP scaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type pi(piSEXP);
    Rcpp::traits::input_parameter< bool >::type scale(scaleSEXP);
    rcpp_result_gen = Rcpp::wrap(HKY85(alpha, beta, pi, scale));
    return rcpp_result_gen;
END_RCPP
}
// GTR
Rcpp::List GTR(double rAC, double rAG, double rAT, double rCG, double rCT, double rGT, const arma::vec& pi, bool scale);
RcppExport SEXP _phylomd_GTR(SEXP rACSEXP, SEXP rAGSEXP, SEXP rATSEXP, SEXP rCGSEXP, SEXP rCTSEXP, SEXP rGTSEXP, SEXP piSEXP, SEXP scaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type rAC(rACSEXP);
    Rcpp::traits::input_parameter< double >::type rAG(rAGSEXP);
    Rcpp::traits::input_parameter< double >::type rAT(rATSEXP);
    Rcpp::traits::input_parameter< double >::type rCG(rCGSEXP);
    Rcpp::traits::input_parameter< double >::type rCT(rCTSEXP);
    Rcpp::traits::input_parameter< double >::type rGT(rGTSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type pi(piSEXP);
    Rcpp::traits::input_parameter< bool >::type scale(scaleSEXP);
    rcpp_result_gen = Rcpp::wrap(GTR(rAC, rAG, rAT, rCG, rCT, rGT, pi, scale));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_phylomd_ctmc_tpm", (DL_FUNC) &_phylomd_ctmc_tpm, 2},
    {"_phylomd_ctmc_nsubs_moments", (DL_FUNC) &_phylomd_ctmc_nsubs_moments, 4},
    {"_phylomd_ctmc_reward_moments", (DL_FUNC) &_phylomd_ctmc_reward_moments, 4},
    {"_phylomd_ctmc_Q_derivatives", (DL_FUNC) &_phylomd_ctmc_Q_derivatives, 4},
    {"_phylomd_ctmc_t_derivatives", (DL_FUNC) &_phylomd_ctmc_t_derivatives, 3},
    {"_phylomd_phylo_likelihood", (DL_FUNC) &_phylomd_phylo_likelihood, 3},
    {"_phylomd_phylo_nsubs_moments", (DL_FUNC) &_phylomd_phylo_nsubs_moments, 6},
    {"_phylomd_phylo_reward_moments", (DL_FUNC) &_phylomd_phylo_reward_moments, 6},
    {"_phylomd_phylo_Q_derivatives", (DL_FUNC) &_phylomd_phylo_Q_derivatives, 5},
    {"_phylomd_phylo_t_derivatives", (DL_FUNC) &_phylomd_phylo_t_derivatives, 4},
    {"_phylomd_ctmc_sim", (DL_FUNC) &_phylomd_ctmc_sim, 3},
    {"_phylomd_asr_sim", (DL_FUNC) &_phylomd_asr_sim, 3},
    {"_phylomd_smap_sim", (DL_FUNC) &_phylomd_smap_sim, 3},
    {"_phylomd_tips_sim", (DL_FUNC) &_phylomd_tips_sim, 2},
    {"_phylomd_JC69", (DL_FUNC) &_phylomd_JC69, 2},
    {"_phylomd_K80", (DL_FUNC) &_phylomd_K80, 3},
    {"_phylomd_F81", (DL_FUNC) &_phylomd_F81, 3},
    {"_phylomd_HKY85", (DL_FUNC) &_phylomd_HKY85, 4},
    {"_phylomd_GTR", (DL_FUNC) &_phylomd_GTR, 8},
    {NULL, NULL, 0}
};

RcppExport void R_init_phylomd(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
