#include "simulation_functions.h"

#include <algorithm>
#include <cstddef>
#include <tuple>

#include <RcppArmadilloExtensions/sample.h>
#include "EdgeList.h"
#include "ListID.h"
#include "NodeList.h"
#include "PartitionSet.h"
#include "counting_tables.h"
#include "ctmc_moments_derivatives.h"
#include "postorder_traverse.h"


//
// This source file defines the simulation functions.
//


Rcpp::DataFrame ctmc_sim_aux(double t, const Vector<std::string>& states,
                             const arma::mat& Q, const Vector<int>& state_inds,
                             int init_state_ind) {
  // Initialize the CTMC simulation process.
  Vector<double> sim_times = {0.0};
  Vector<std::string> sim_states = {states[init_state_ind]};
  int curr_state_ind = init_state_ind;

  // Run the CTMC simulation process.
  while (sim_times.back() <= t) {
    // Sample the exponentially distributed waiting time.
    double wait_time = Rcpp::rexp(1, -Q(curr_state_ind, curr_state_ind))[0];

    // Sample the CTMC state.
    arma::vec state_probs = Q.row(curr_state_ind).t();
    state_probs(curr_state_ind) = 0;
    int next_state_ind =
        Rcpp::RcppArmadillo::sample(state_inds, 1, true, state_probs)[0];

    // Cache the current sample.
    sim_times.push_back(sim_times.back() + wait_time);
    sim_states.push_back(states[next_state_ind]);
    curr_state_ind = next_state_ind;
  }

  // Discard the final sample.
  sim_times.back() = t;
  sim_states.back() = sim_states.end()[-2];

  return Rcpp::DataFrame::create(Rcpp::Named("time") = sim_times,
                                 Rcpp::Named("state") = sim_states,
                                 Rcpp::Named("stringsAsFactors") = false);
}


std::pair<Vector<std::string>, Vector<int>> asr_sim_aux(
    arma::imat& edge, const Vector<std::string>& tip_labels, int num_int_nodes,
    const arma::vec& edge_lengths, const Vector<std::string>& states,
    const arma::mat& Q, const arma::vec& pi, const Vector<int>& state_inds,
    const arma::ivec& tip_data) {
  // Modify the edge matrix and create some useful variables.
  edge -= 1;
  int num_edges = edge.n_rows;
  int num_term_nodes = tip_labels.size();
  int root_node_ind = num_term_nodes;

  // Initialize the storage of the CTMC transition probability matrix at each
  // edge.
  Vector<arma::cube> edges_tpm(num_edges);

  // Create the "empty" (i.e. partial likelihood) node list and edge list.
  arma::mat choose = choose_table(0);
  Vector<ListID> list_ids = {
      ListID(Vector<ListIDElement>(), Map<Ref<const PartitionSet>, int>())};
  Vector<NodeList> nlists = {NodeList(list_ids[0], list_ids, choose, Q.n_rows,
                                      num_int_nodes + num_term_nodes)};
  Vector<EdgeList> elists = {
      EdgeList(list_ids[0], list_ids, choose, Q.n_rows, num_edges)};

  // Perform the postorder recursion.
  postorder_traverse(edge, edge_lengths, Q, arma::mat(), 0, Mode::T_DERIVATIVES,
                     tip_data, Map<int, const PartitionSet&>(), num_edges,
                     num_term_nodes, root_node_ind, edges_tpm, nlists, elists);

  // Initialize the storage of the ancestral state samples.
  Vector<std::string> asr_states(num_int_nodes + num_term_nodes);
  Vector<int> asr_state_inds(num_int_nodes + num_term_nodes);

  // Sample the root node ancestral state.
  arma::vec state_probs = nlists[0].elems().col(root_node_ind) % pi;
  asr_state_inds[root_node_ind] =
      Rcpp::RcppArmadillo::sample(state_inds, 1, true, state_probs)[0];
  asr_states[root_node_ind] = states[asr_state_inds[root_node_ind]];

  // Sample the remaining ancestral states.
  for (int edge_ind = 0; edge_ind < num_edges; ++edge_ind) {
    int parent_node_ind = edge(edge_ind, 0);
    int child_node_ind = edge(edge_ind, 1);

    // Is the current child node an internal node?
    if (child_node_ind >= root_node_ind) {
      // If so, sample the ancestral state at the current child node.
      arma::vec state_probs =
          nlists[0].elems().col(child_node_ind) %
          edges_tpm[edge_ind].slice(0).row(asr_state_inds[parent_node_ind]).t();
      asr_state_inds[child_node_ind] =
          Rcpp::RcppArmadillo::sample(state_inds, 1, true, state_probs)[0];
    } else {
      // Otherwise, fix the "ancestral state" to be the observed tip state.
      asr_state_inds[child_node_ind] = tip_data(child_node_ind);
    }

    asr_states[child_node_ind] = (asr_state_inds[child_node_ind] != -1)
                                     ? states[asr_state_inds[child_node_ind]]
                                     : "N";
  }

  return {asr_states, asr_state_inds};
}


//' CTMC path simulation
//' 
//' Simulates a CTMC sample path.
//' 
//' @param t A nonnegative numeric scalar representing the CTMC time interval 
//'   length.
//' @param subst_mod An S3 object of class \code{"substitution.model"}.
//' @param init_state A string that specifies the initial state of the CTMC.
//'   
//' @return A data frame with two columns.  The \code{"time"} column holds the 
//'   times at which each of the simulated CTMC states is entered and the 
//'   \code{"state"} column stores the corresponding CTMC states.
//'   
//' @references Nielsen R (2002) \dQuote{Mapping Mutations on Phylogenies}, 
//'   \emph{Systematic Biology}, 51(5):729-739.
//'   
//' @seealso \code{\link{asr.sim}}, \code{\link{smap.sim}}, 
//'   \code{\link{tips.sim}}
//'   
//' @export
// [[Rcpp::export(name = "ctmc.sim")]]
Rcpp::DataFrame ctmc_sim(double t, const Rcpp::List& subst_mod,
                         const std::string& init_state) {
  if (t < 0.0) Rcpp::stop("'t' cannot be less than 0.");
  if (!subst_mod.inherits("substitution.model"))
    Rcpp::stop("'subst_mod' must be an object of class 'substitution.model'.");

  const Vector<std::string>& states = subst_mod["states"];
  const arma::mat& Q = subst_mod["Q"];
  Vector<int> state_inds;
  state_inds.reserve(states.size());
  for (std::size_t state_ind = 0; state_ind < states.size(); ++state_ind) {
    state_inds.push_back(state_ind);
  }

  int init_state_ind;
  auto find_it = std::find(states.begin(), states.end(), init_state);
  if (find_it != states.end()) {
    init_state_ind = find_it - states.begin();
  } else {
    Rcpp::stop("'init_state' must be a valid state.");
  }

  return ctmc_sim_aux(t, states, Q, state_inds, init_state_ind);
}


//' Phylogenetic ancestral state reconstruction
//' 
//' Samples ancestral states at the internal nodes of a phylogeny conditional on
//' observed tip states.
//' 
//' More information on S3 objects of class \code{"phylo"} is found at 
//' \url{http://ape-package.ird.fr/misc/FormatTreeR_24Oct2012.pdf}.
//' 
//' If an observed tip state is invalid, then it is treated as an ambiguous 
//' character state.
//' 
//' @param tree An S3 object of class \code{"phylo"}.
//' @param subst_mod An S3 object of class \code{"substitution.model"}.
//' @param tip_states A character vector of observed tip states.
//'   
//' @return A character vector that holds the sampled ancestral states and 
//'   observed tip states.
//'   
//' @references Nielsen R (2002) \dQuote{Mapping Mutations on Phylogenies}, 
//'   \emph{Systematic Biology}, 51(5):729-739.
//'   
//' @seealso \code{\link{smap.sim}}, \code{\link{tips.sim}}, 
//'   \code{\link{ctmc.sim}}
//'   
//' @export
// [[Rcpp::export(name = "asr.sim")]]
std::vector<std::string> asr_sim(const Rcpp::List& tree,
                                 const Rcpp::List& subst_mod,
                                 const std::vector<std::string>& tip_states) {
  if (!tree.inherits("phylo"))
    Rcpp::stop("'tree' must be an object of class 'phylo'.");
  if (!tree.hasAttribute("order") ||
      Rcpp::as<std::string>(tree.attr("order")) != "cladewise")
    Rcpp::stop("The edge matrix must be in 'cladewise' order.");
  if (!tree.containsElementNamed("edge.length"))
    Rcpp::stop("'tree' must contain a vector of edge lengths.");
  if (!subst_mod.inherits("substitution.model"))
    Rcpp::stop("'subst_mod' must be an object of class 'substitution.model'.");

  arma::imat edge = tree["edge"];
  const Vector<std::string>& tip_labels = tree["tip.label"];
  int num_int_nodes = tree["Nnode"];
  const arma::vec& edge_lengths = tree["edge.length"];
  const Vector<std::string>& states = subst_mod["states"];
  const arma::mat& Q = subst_mod["Q"];
  const arma::vec& pi = subst_mod["pi"];
  Vector<int> state_inds;
  state_inds.reserve(states.size());
  for (std::size_t state_ind = 0; state_ind < states.size(); ++state_ind) {
    state_inds.push_back(state_ind);
  }

  if (arma::find(edge.col(0) == tip_labels.size() + 1).eval().n_elem > 2)
    Rcpp::stop("'tree' must be a rooted tree.");
  if (tip_states.size() != tip_labels.size())
    Rcpp::stop("'tip_states' must be compatible with 'tree'.");

  arma::ivec tip_data(tip_states.size(), arma::fill::none);
  for (std::size_t i = 0; i < tip_states.size(); ++i) {
    auto find_it = std::find(states.begin(), states.end(), tip_states[i]);
    tip_data(i) = (find_it != states.end()) ? find_it - states.begin() : -1;
  }

  return asr_sim_aux(edge, tip_labels, num_int_nodes, edge_lengths, states, Q,
                     pi, state_inds, tip_data)
      .first;
}


//' Phylogenetic stochastic mapping simulation
//' 
//' Samples a stochastic mapping on a phylogeny conditional on observed tip 
//' states.
//' 
//' More information on S3 objects of class \code{"phylo"} is found at 
//' \url{http://ape-package.ird.fr/misc/FormatTreeR_24Oct2012.pdf}.
//' 
//' If an observed tip state is invalid, then it is treated as an ambiguous 
//' character state.
//' 
//' @param tree An S3 object of class \code{"phylo"}.
//' @param subst_mod An S3 object of class \code{"substitution.model"}.
//' @param tip_states A character vector of observed tip states.
//'   
//' @return A list that stores the simulated CTMC trajectory at each edge.
//'   
//' @references Nielsen R (2002) \dQuote{Mapping Mutations on Phylogenies}, 
//'   \emph{Systematic Biology}, 51(5):729-739.
//'   
//' @seealso \code{\link{asr.sim}}, \code{\link{tips.sim}}, 
//'   \code{\link{ctmc.sim}}
//'   
//' @export
// [[Rcpp::export(name = "smap.sim")]]
std::vector<Rcpp::DataFrame> smap_sim(
    const Rcpp::List& tree, const Rcpp::List& subst_mod,
    const std::vector<std::string>& tip_states) {
  if (!tree.inherits("phylo"))
    Rcpp::stop("'tree' must be an object of class 'phylo'.");
  if (!tree.hasAttribute("order") ||
      Rcpp::as<std::string>(tree.attr("order")) != "cladewise")
    Rcpp::stop("The edge matrix must be in 'cladewise' order.");
  if (!tree.containsElementNamed("edge.length"))
    Rcpp::stop("'tree' must contain a vector of edge lengths.");
  if (!subst_mod.inherits("substitution.model"))
    Rcpp::stop("'subst_mod' must be an object of class 'substitution.model'.");

  arma::imat edge = tree["edge"];
  const Vector<std::string>& tip_labels = tree["tip.label"];
  int num_int_nodes = tree["Nnode"];
  const arma::vec& edge_lengths = tree["edge.length"];
  int num_edges = edge.n_rows;
  const Vector<std::string>& states = subst_mod["states"];
  const arma::mat& Q = subst_mod["Q"];
  const arma::vec& pi = subst_mod["pi"];
  Vector<int> state_inds;
  state_inds.reserve(states.size());
  for (std::size_t state_ind = 0; state_ind < states.size(); ++state_ind) {
    state_inds.push_back(state_ind);
  }

  if (arma::find(edge.col(0) == tip_labels.size() + 1).eval().n_elem > 2)
    Rcpp::stop("'tree' must be a rooted tree.");
  if (tip_states.size() != tip_labels.size())
    Rcpp::stop("'tip_states' must be compatible with 'tree'.");

  arma::ivec tip_data(tip_states.size(), arma::fill::none);
  for (std::size_t i = 0; i < tip_states.size(); ++i) {
    auto find_it = std::find(states.begin(), states.end(), tip_states[i]);
    tip_data(i) = (find_it != states.end()) ? find_it - states.begin() : -1;
  }

  // Sample the ancestral states.
  Vector<std::string> asr_states;
  Vector<int> asr_state_inds;
  std::tie(asr_states, asr_state_inds) =
      asr_sim_aux(edge, tip_labels, num_int_nodes, edge_lengths, states, Q, pi,
                  state_inds, tip_data);

  // Initialize the storage of the stochastic mapping sample.
  Vector<Rcpp::DataFrame> smap(num_edges);

  // Run the CTMC simulation process at each edge.
  for (int edge_ind = 0; edge_ind < num_edges; ++edge_ind) {
    int parent_node_ind = edge(edge_ind, 0);
    int child_node_ind = edge(edge_ind, 1);

    do {
      smap[edge_ind] =
          ctmc_sim_aux(edge_lengths(edge_ind), states, Q, state_inds,
                       asr_state_inds[parent_node_ind]);
    } while (Rcpp::as<Vector<std::string>>(smap[edge_ind]["state"]).back() !=
                 asr_states[child_node_ind] &&
             asr_states[child_node_ind] != "N");
  }

  return smap;
}


//' Phylogenetic tip state simulation
//' 
//' Simulates the tip states of a phylogeny.
//' 
//' More information on S3 objects of class \code{"phylo"} is found at 
//' \url{http://ape-package.ird.fr/misc/FormatTreeR_24Oct2012.pdf}.
//' 
//' @param tree An S3 object of class \code{"phylo"}.
//' @param subst_mod An S3 object of class \code{"substitution.model"}.
//'   
//' @return A character vector of simulated tip states.
//'   
//' @seealso \code{\link{asr.sim}}, \code{\link{smap.sim}}, 
//'   \code{\link{ctmc.sim}}
//'   
//' @export
// [[Rcpp::export(name = "tips.sim")]]
std::vector<std::string> tips_sim(const Rcpp::List& tree,
                                  const Rcpp::List& subst_mod) {
  if (!tree.inherits("phylo"))
    Rcpp::stop("'tree' must be an object of class 'phylo'.");
  if (!tree.hasAttribute("order") ||
      Rcpp::as<std::string>(tree.attr("order")) != "cladewise")
    Rcpp::stop("The edge matrix must be in 'cladewise' order.");
  if (!tree.containsElementNamed("edge.length"))
    Rcpp::stop("'tree' must contain a vector of edge lengths.");
  if (!subst_mod.inherits("substitution.model"))
    Rcpp::stop("'subst_mod' must be an object of class 'substitution.model'.");

  arma::imat edge = tree["edge"];
  const Vector<std::string>& tip_labels = tree["tip.label"];
  int num_int_nodes = tree["Nnode"];
  const arma::vec& edge_lengths = tree["edge.length"];
  edge -= 1;
  int num_edges = edge.n_rows;
  int num_term_nodes = tip_labels.size();
  int root_node_ind = num_term_nodes;
  const Vector<std::string>& states = subst_mod["states"];
  const arma::mat& Q = subst_mod["Q"];
  arma::vec pi = subst_mod["pi"];
  Vector<int> state_inds;
  state_inds.reserve(states.size());
  for (std::size_t state_ind = 0; state_ind < states.size(); ++state_ind) {
    state_inds.push_back(state_ind);
  }

  if (arma::find(edge.col(0) == tip_labels.size() + 1).eval().n_elem > 2)
    Rcpp::stop("'tree' must be a rooted tree.");

  // Initialize the storage of the observed tip states.
  Vector<std::string> tip_states(num_term_nodes);
  Vector<int> node_state_inds(num_int_nodes + num_term_nodes);

  // Sample the root node state.
  node_state_inds[root_node_ind] =
      Rcpp::RcppArmadillo::sample(state_inds, 1, true, pi)[0];

  // Sample the remaining node states.
  for (int edge_ind = 0; edge_ind < num_edges; ++edge_ind) {
    int parent_node_ind = edge(edge_ind, 0);
    int child_node_ind = edge(edge_ind, 1);

    arma::cube edge_tpm = ctmc_moments_derivatives(
        edge_lengths(edge_ind), Q, arma::mat(), 0, Mode::T_DERIVATIVES);
    arma::vec state_probs =
        edge_tpm.slice(0).row(node_state_inds[parent_node_ind]).t();
    node_state_inds[child_node_ind] =
        Rcpp::RcppArmadillo::sample(state_inds, 1, true, state_probs)[0];

    // Is the current child node a terminal node?
    if (child_node_ind < root_node_ind) {
      // If so, cache the observed tip state.
      tip_states[child_node_ind] = states[node_state_inds[child_node_ind]];
    }
  }

  return tip_states;
}