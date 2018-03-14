#include "phylo_moments_derivatives.h"

#include <algorithm>
#include <cstddef>
#include <tuple>
#include <utility>

#include "counting_tables.h"


//
// This source file defines the phylogenetic moment/derivative functions.
//


void update_node_lists(int child_node_ind, const arma::uvec& child_edge_inds,
                       const Vector<EdgeList>& elists,
                       Vector<NodeList>& nlists) {
  for (std::size_t nlist_ind = 0; nlist_ind < nlists.size(); ++nlist_ind) {
    for (const auto& ri_tup : nlists[nlist_ind].recursion_info()) {
      // Extract the recursion information 3-tuple elements.
      int left_elist_ind = std::get<0>(ri_tup);
      int right_elist_ind = std::get<1>(ri_tup);
      int choose_coef = std::get<2>(ri_tup);

      // Update the current node list.
      nlists[nlist_ind].elems().col(child_node_ind) +=
          elists[left_elist_ind].elems().col(child_edge_inds(0)) %
          elists[right_elist_ind].elems().col(child_edge_inds(1)) * choose_coef;
    }
  }
}


void update_edge_lists(
    int edge_ind, int child_node_ind, const arma::cube& edge_mds,
    Map<int, const PartitionSet&>::const_iterator edge_pset_it,
    Map<int, const PartitionSet&>::const_iterator edge_psets_end_it,
    const Vector<NodeList>& nlists, Vector<EdgeList>& elists) {
  for (std::size_t elist_ind = 0; elist_ind < elists.size(); ++elist_ind) {
    // Is the current edge in any of the partition sets?
    if (edge_pset_it != edge_psets_end_it) {
      // If so, is this partition set in the current edge list ID?
      auto pset_inds_it =
          elists[elist_ind].recursion_info_inds().find(edge_pset_it->second);
      if (pset_inds_it != elists[elist_ind].recursion_info_inds().end()) {
        // If so, we must account for the possible assignments of ID elements to
        // the current edge.
        int begin_ind, end_ind;
        std::tie(begin_ind, end_ind) = pset_inds_it->second;

        for (int ri_ind = begin_ind; ri_ind < end_ind; ++ri_ind) {
          // Extract the recursion information 3-tuple elements.
          const auto& ri_tup = elists[elist_ind].recursion_info()[ri_ind];
          const ListIDElement& curr_id_elem = std::get<0>(ri_tup);
          int nlist_ind = std::get<1>(ri_tup);
          int choose_coef = std::get<2>(ri_tup);

          // Update the current edge list.
          elists[elist_ind].elems().col(edge_ind) +=
              edge_mds.slice(curr_id_elem.order()) *
              nlists[nlist_ind].elems().col(child_node_ind) * choose_coef;
        }
      }
    }

    // Update the current edge list.
    // (Note: this update accounts for the case where no ID element is assigned
    // to the current edge.)
    elists[elist_ind].elems().col(edge_ind) +=
        edge_mds.slice(0) * nlists[elist_ind].elems().col(child_node_ind);
  }
}


void find_connected_moment_derivative_ids_aux(
    const FlatListID& flat_list_id, const Vector<EdgeSet>& esets,
    FlatListID::const_iterator curr_it, const Map<int, int>& eset_labels_orders,
    Vector<MomentDerivativeID>& md_ids) {
  // If we have traversed over all the flattened list ID elements, then cache
  // the current moment/derivative ID and exit the function.
  if (curr_it == flat_list_id.end()) {
    MomentDerivativeID curr_md_id;
    curr_md_id.reserve(eset_labels_orders.size());

    for (const auto& eset_label_order : eset_labels_orders) {
      curr_md_id.emplace_back(esets[eset_label_order.first],
                              eset_label_order.second);
    }

    md_ids.push_back(std::move(curr_md_id));
    return;
  }

  // Otherwise, loop over the possible edge set labels in the current flattened
  // list ID element.
  for (auto eset_label : curr_it->set().label()) {
    // Construct the next (`eset_label`, `order`) map by inserting the current
    // (`eset_label`, 1) pair into the current map.
    Map<int, int> next_eset_labels_orders = eset_labels_orders;
    auto insert_results = next_eset_labels_orders.emplace(eset_label, 1);
    if (!insert_results.second) insert_results.first->second += 1;

    // Recurse over the next flattened list ID element.
    find_connected_moment_derivative_ids_aux(flat_list_id, esets, curr_it + 1,
                                             next_eset_labels_orders, md_ids);
  }
}


Vector<MomentDerivativeID> find_connected_moment_derivative_ids(
    const FlatListID& flat_list_id, const Vector<EdgeSet>& esets) {
  Vector<MomentDerivativeID> md_ids;
  find_connected_moment_derivative_ids_aux(flat_list_id, esets,
                                           flat_list_id.begin(), {}, md_ids);

  // Keep only the unique moment/derivative IDs.
  std::sort(md_ids.begin(), md_ids.end());
  auto end_it = std::unique(md_ids.begin(), md_ids.end());
  md_ids.resize(end_it - md_ids.begin());

  return md_ids;
}


int get_connected_counting_coef(const FlatMomentDerivativeID& flat_md_id,
                                FlatListID& flat_list_id) {
  int counting_coef = 0;

  // Loop through the permutations of the flattened list ID and count those that
  // are compatible with the flattened moment/derivative ID.
  do {
    bool is_compatible = std::equal(
        flat_md_id.begin(), flat_md_id.end(), flat_list_id.begin(),
        [](const FlatMomentDerivativeIDElement& flat_md_id_elem,
           const FlatListIDElement& flat_list_id_elem) {
          return std::binary_search(flat_list_id_elem.set().label().begin(),
                                    flat_list_id_elem.set().label().end(),
                                    flat_md_id_elem.set().label());
        });
    if (is_compatible) counting_coef += 1;
  } while (std::next_permutation(flat_list_id.begin(), flat_list_id.end()));

  return counting_coef;
}


Map<std::string, double> phylo_moments_derivatives(
    arma::imat& edge, const Vector<std::string>& tip_labels, int num_int_nodes,
    const arma::vec& edge_lengths, const arma::mat& Q, const arma::mat& B,
    const arma::vec& pi, VectorVector<int>& esets_inp, int max_order, Mode mode,
    const arma::ivec& tip_data) {
  // Modify the edge matrix and create some useful variables.
  edge -= 1;
  int num_edges = edge.n_rows;
  int num_term_nodes = tip_labels.size();
  int root_node_ind = num_term_nodes;

  // Initialize the counting table.
  arma::mat choose = choose_table(max_order);

  // Construct the edge sets and partition sets.
  Vector<EdgeSet> esets = create_edge_sets(esets_inp);
  Vector<PartitionSet> psets = partition_edge_sets(esets);

  // Form the (edge index, partition set) inverse map.
  Map<int, const PartitionSet&> edge_psets = create_edge_pset_map(psets);

  // Generate the moment/derivative IDs and list IDs.
  Vector<MomentDerivativeID> md_ids =
      get_moment_derivative_ids(esets, max_order);
  Vector<ListID> list_ids = get_list_ids(psets, max_order);

  // Flatten the moment/derivative IDs and list IDs.
  // Create the node lists and edge lists.
  Vector<FlatMomentDerivativeID> flat_md_ids(md_ids.size());
  Vector<FlatListID> flat_list_ids(list_ids.size());

  Vector<NodeList> nlists;
  nlists.reserve(list_ids.size());
  Vector<EdgeList> elists;
  elists.reserve(list_ids.size());

  for (std::size_t i = 0; i < md_ids.size(); ++i) {
    for (const auto& id_elem : md_ids[i]) {
      flat_md_ids[i].insert(flat_md_ids[i].end(), id_elem.order(),
                            id_elem.set());
    }
  }

  for (std::size_t i = 0; i < list_ids.size(); ++i) {
    for (const auto& id_elem : list_ids[i]) {
      flat_list_ids[i].insert(flat_list_ids[i].end(), id_elem.order(),
                              id_elem.set());
    }

    nlists.emplace_back(list_ids[i], list_ids, choose, Q.n_rows,
                        num_int_nodes + num_term_nodes);
    elists.emplace_back(list_ids[i], list_ids, choose, Q.n_rows, num_edges);
  }

  // Perform the postorder recursion.
  // (Note: we iterate up the edge matrix because, by default, it is arranged in
  // a preorder format.)
  for (int edge_ind = num_edges - 1; edge_ind >= 0; --edge_ind) {
    int child_node_ind = edge(edge_ind, 1);

    ////
    //// First, we consider the node lists at the current child node.
    ////

    // Is the current child node an internal node?
    if (child_node_ind >= root_node_ind) {
      // If so, update each node list using the appropriate edge lists at the
      // two branches below the current child node.
      arma::uvec child_edge_inds = arma::find(edge.col(0) == child_node_ind);
      update_node_lists(child_node_ind, child_edge_inds, elists, nlists);
    } else {
      // Otherwise, initialize the "empty" (i.e. partial likelihood) node list
      // at the given terminal node.
      // (Note: the "empty" node list is always the first element in the node
      // list vector.)

      // Do we have observed data at the terminal nodes?
      if ((int)tip_data.n_elem == num_term_nodes) {
        nlists[0].elems()(tip_data(child_node_ind), child_node_ind) = 1;
      } else {
        nlists[0].elems().col(child_node_ind).fill(1);
      }
    }

    ////
    //// Second, we consider the edge lists at the current edge.
    ////

    // Is the current edge in any of the partition sets?
    auto edge_pset_it = edge_psets.find(edge_ind);

    // If so, compute the necessary CTMC moments/derivatives at the current
    // edge.
    // Otherwise, compute only the zeroth CTMC moment/derivative (i.e. the
    // transition probability matrix) at the current edge.
    int edge_max_order = (edge_pset_it != edge_psets.end()) ? max_order : 0;
    arma::cube edge_mds = ctmc_moments_derivatives(edge_lengths(edge_ind), Q, B,
                                                   edge_max_order, mode);

    // Update the edge lists at the current edge.
    update_edge_lists(edge_ind, child_node_ind, edge_mds, edge_pset_it,
                      edge_psets.end(), nlists, elists);
  }

  // Update the node lists at the root node.
  arma::uvec root_edge_inds = arma::find(edge.col(0) == root_node_ind);
  update_node_lists(root_node_ind, root_edge_inds, elists, nlists);

  // Calculate the moments/derivatives.
  arma::vec phylo_mds(md_ids.size(), arma::fill::zeros);
  for (std::size_t i = 0; i < nlists.size(); ++i) {
    // Compute the contribution from the current node list.
    double md_part = arma::dot(pi, nlists[i].elems().col(root_node_ind));

    // Identify the moment/derivative IDs connected to the current node list ID.
    Vector<MomentDerivativeID> conn_md_ids =
        find_connected_moment_derivative_ids(flat_list_ids[i], esets);

    // Add the current node list contribution to the connected
    // moments/derivatives.
    for (const auto& conn_md_id : conn_md_ids) {
      int conn_md_id_ind = find_moment_derivative_id_index(md_ids, conn_md_id);

      // Determine the connected counting coefficient.
      int counting_coef = get_connected_counting_coef(
          flat_md_ids[conn_md_id_ind], flat_list_ids[i]);

      // Update the connected moment/derivative.
      phylo_mds(conn_md_id_ind) += md_part * counting_coef;
    }
  }

  // Prepare the output map.
  Map<std::string, double> phylo_mds_map;
  for (std::size_t i = 0; i < md_ids.size(); ++i) {
    // Create the moment/derivative ID label.
    std::string md_id_label = create_moment_derivative_id_label(md_ids[i]);

    // Cache the moment/derivative.
    phylo_mds_map.emplace(std::move(md_id_label), phylo_mds(i));
  }

  return phylo_mds_map;
}


//' Phylogenetic stochastic mapping moments of labeled substitution counts
//' 
//' Computes the phylogenetic stochastic mapping moments of labeled substitution
//' counts for all orders less than or equal to the user-specified maximum order
//' of interest.
//' 
//' More information on S3 objects of class \code{"phylo"} is found at 
//' \url{http://ape-package.ird.fr/misc/FormatTreeR_24Oct2012.pdf}.
//' 
//' The zeroth-order stochastic mapping moment is defined to be the phylogenetic
//' likelihood.
//' 
//' If the number of observed tip states (\code{length(tip_states)}) is not 
//' equal to the number of tips on the phylogeny 
//' (\code{length(tree$tip.label)}), then the stochastic mapping moments are 
//' computed using ambiguous character states.
//' 
//' @param tree An S3 object of class \code{"phylo"}.
//' @param subst_mod An S3 object of class \code{"substitution.model"}.
//' @param L An indicator matrix that labels the substitutions of interest.
//' @param edge_sets A list of integer vectors that defines the edge sets over 
//'   which the moments are calculated.
//' @param max_order An integer specifying the maximum moment order of interest.
//' @param tip_states A character vector of observed tip states.
//'   
//' @return A named numeric vector that holds the computed stochastic mapping 
//'   moments.  Each element name is formed by concatenating strings of the form
//'   \code{"(E):O"} with \code{"--"} separators, where \code{E} denotes an edge
//'   set index, \code{O} represents the exponent order associated with the 
//'   additive mapping summary defined over the \code{E}th edge set, and 
//'   \code{"--"} symbolizes the product operator between two mapping summaries.
//'   However, the empty string (\code{""}) is the element name of the 
//'   zeroth-order stochastic mapping moment.
//'   
//' @references Minin VN and Suchard MA (2008) \dQuote{Fast, Accurate and 
//'   Simulation-Free Stochastic Mapping}, \emph{Philosophical Transactions of 
//'   the Royal Society B: Biological Sciences}, 363(1512):3985-3995.
//'   
//'   Dhar A and Minin VN (2017) \dQuote{Calculating Higher-Order Moments of 
//'   Phylogenetic Stochastic Mapping Summaries in Linear Time}, \emph{Journal 
//'   of Computational Biology}, 24(5):377-399.
//'   
//'   OUR PAPER!!
//'   
//' @seealso \code{\link{phylo.reward.moments}}, 
//'   \code{\link{ctmc.nsubs.moments}}
//'   
//' @export
// [[Rcpp::export(name = "phylo.nsubs.moments")]]
Map<std::string, double> phylo_nsubs_moments(
    const Rcpp::List& tree, const Rcpp::List& subst_mod, const arma::mat& L,
    VectorVector<int> edge_sets, int max_order,
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
  if (arma::any(arma::abs(arma::vectorise(L) - 0) >= arma::datum::eps &&
                arma::abs(arma::vectorise(L) - 1) >= arma::datum::eps))
    Rcpp::stop("'L' must be an indicator matrix.");
  if (max_order < 0) Rcpp::stop("'max_order' cannot be less than 0.");

  arma::imat edge = tree["edge"];
  const Vector<std::string>& tip_labels = tree["tip.label"];
  int num_int_nodes = tree["Nnode"];
  const arma::vec& edge_lengths = tree["edge.length"];
  const Vector<std::string>& states = subst_mod["states"];
  const arma::mat& Q = subst_mod["Q"];
  const arma::vec& pi = subst_mod["pi"];

  if (arma::find(edge.col(0) == tip_labels.size() + 1).eval().n_elem > 2)
    Rcpp::stop("'tree' must be a rooted tree.");
  if (arma::size(Q) != arma::size(L))
    Rcpp::stop("The rate matrix and 'L' must have the same dimensions.");
  for (const auto& edge_set : edge_sets) {
    for (auto edge_label : edge_set) {
      if (edge_label < 1 || edge_label > (int)edge.n_rows)
        Rcpp::stop("'edge_sets' must contain valid edges.");
    }
  }

  arma::ivec tip_data(tip_states.size(), arma::fill::none);
  for (std::size_t i = 0; i < tip_states.size(); ++i) {
    auto find_it = std::find(states.begin(), states.end(), tip_states[i]);
    if (find_it != states.end()) {
      tip_data(i) = find_it - states.begin();
    } else {
      Rcpp::stop("'tip_states' must contain valid states.");
    }
  }

  // Divide the computed moments by the likelihood.
  Map<std::string, double> moments = phylo_moments_derivatives(
      edge, tip_labels, num_int_nodes, edge_lengths, Q, Q % L, pi, edge_sets,
      max_order, Mode::NSUBS_MOMENTS, tip_data);
  double likelihood = moments.at("");

  for (auto it = moments.begin(); it != moments.end(); ++it) {
    if (it->first != "") it->second /= likelihood;
  }

  return moments;
}


//' Phylogenetic stochastic mapping moments of total rewards
//' 
//' Computes the phylogenetic stochastic mapping moments of total rewards for 
//' all orders less than or equal to the user-specified maximum order of 
//' interest.
//' 
//' More information on S3 objects of class \code{"phylo"} is found at 
//' \url{http://ape-package.ird.fr/misc/FormatTreeR_24Oct2012.pdf}.
//' 
//' The zeroth-order stochastic mapping moment is defined to be the phylogenetic
//' likelihood.
//' 
//' If the number of observed tip states (\code{length(tip_states)}) is not 
//' equal to the number of tips on the phylogeny 
//' (\code{length(tree$tip.label)}), then the stochastic mapping moments are 
//' computed using ambiguous character states.
//' 
//' @param tree An S3 object of class \code{"phylo"}.
//' @param subst_mod An S3 object of class \code{"substitution.model"}.
//' @param w A numeric vector that defines the reward per unit time (reward 
//'   rate) for each CTMC state.
//' @param edge_sets A list of integer vectors that defines the edge sets over 
//'   which the moments are calculated.
//' @param max_order An integer specifying the maximum moment order of interest.
//' @param tip_states A character vector of observed tip states.
//'   
//' @return A named numeric vector that holds the computed stochastic mapping 
//'   moments.  Each element name is formed by concatenating strings of the form
//'   \code{"(E):O"} with \code{"--"} separators, where \code{E} denotes an edge
//'   set index, \code{O} represents the exponent order associated with the 
//'   additive mapping summary defined over the \code{E}th edge set, and 
//'   \code{"--"} symbolizes the product operator between two mapping summaries.
//'   However, the empty string (\code{""}) is the element name of the 
//'   zeroth-order stochastic mapping moment.
//'   
//' @references Minin VN and Suchard MA (2008) \dQuote{Fast, Accurate and 
//'   Simulation-Free Stochastic Mapping}, \emph{Philosophical Transactions of 
//'   the Royal Society B: Biological Sciences}, 363(1512):3985-3995.
//'   
//'   Dhar A and Minin VN (2017) \dQuote{Calculating Higher-Order Moments of 
//'   Phylogenetic Stochastic Mapping Summaries in Linear Time}, \emph{Journal 
//'   of Computational Biology}, 24(5):377-399.
//'   
//'   OUR PAPER!!
//'   
//' @seealso \code{\link{phylo.nsubs.moments}}, 
//'   \code{\link{ctmc.reward.moments}}
//'   
//' @export
// [[Rcpp::export(name = "phylo.reward.moments")]]
Map<std::string, double> phylo_reward_moments(
    const Rcpp::List& tree, const Rcpp::List& subst_mod, const arma::vec& w,
    VectorVector<int> edge_sets, int max_order,
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
  if (max_order < 0) Rcpp::stop("'max_order' cannot be less than 0.");

  arma::imat edge = tree["edge"];
  const Vector<std::string>& tip_labels = tree["tip.label"];
  int num_int_nodes = tree["Nnode"];
  const arma::vec& edge_lengths = tree["edge.length"];
  const Vector<std::string>& states = subst_mod["states"];
  const arma::mat& Q = subst_mod["Q"];
  const arma::vec& pi = subst_mod["pi"];

  if (arma::find(edge.col(0) == tip_labels.size() + 1).eval().n_elem > 2)
    Rcpp::stop("'tree' must be a rooted tree.");
  if (Q.n_rows != w.n_elem)
    Rcpp::stop("The rate matrix and 'w' must have compatible dimensions.");
  for (const auto& edge_set : edge_sets) {
    for (auto edge_label : edge_set) {
      if (edge_label < 1 || edge_label > (int)edge.n_rows)
        Rcpp::stop("'edge_sets' must contain valid edges.");
    }
  }

  arma::ivec tip_data(tip_states.size(), arma::fill::none);
  for (std::size_t i = 0; i < tip_states.size(); ++i) {
    auto find_it = std::find(states.begin(), states.end(), tip_states[i]);
    if (find_it != states.end()) {
      tip_data(i) = find_it - states.begin();
    } else {
      Rcpp::stop("'tip_states' must contain valid states.");
    }
  }

  // Divide the computed moments by the likelihood.
  Map<std::string, double> moments = phylo_moments_derivatives(
      edge, tip_labels, num_int_nodes, edge_lengths, Q, arma::diagmat(w), pi,
      edge_sets, max_order, Mode::REWARD_MOMENTS, tip_data);
  double likelihood = moments.at("");

  for (auto it = moments.begin(); it != moments.end(); ++it) {
    if (it->first != "") it->second /= likelihood;
  }

  return moments;
}


//' Phylogenetic likelihood derivatives with respect to rate matrix parameters
//' 
//' Computes the phylogenetic likelihood derivatives with respect to a given 
//' rate matrix parameter for all orders less than or equal to the 
//' user-specified maximum order of interest.
//' 
//' More information on S3 objects of class \code{"phylo"} is found at 
//' \url{http://ape-package.ird.fr/misc/FormatTreeR_24Oct2012.pdf}.
//' 
//' The zeroth-order likelihood derivative is defined to be the phylogenetic 
//' likelihood.
//' 
//' @param tree An S3 object of class \code{"phylo"}.
//' @param subst_mod An S3 object of class \code{"substitution.model"}.
//' @param param_name The rate matrix parameter name of interest.
//' @param max_order An integer specifying the maximum derivative order of 
//'   interest.
//' @param tip_states A character vector of observed tip states.
//'   
//' @return A named numeric vector that holds the computed likelihood 
//'   derivatives.  Each element name is a string of the form \code{"(1):O"}, 
//'   where \code{O} represents the derivative order of the corresponding 
//'   element.  However, the empty string (\code{""}) is the element name of the
//'   zeroth-order likelihood derivative.
//'   
//' @references Kenney T and Gu H (2012) \dQuote{Hessian Calculation for 
//'   Phylogenetic Likelihood based on the Pruning Algorithm and its 
//'   Applications}, \emph{Statistical Applications in Genetics and Molecular 
//'   Biology}, 11(4).
//'   
//'   OUR PAPER!!
//'   
//' @seealso \code{\link{phylo.t.derivatives}}, \code{\link{ctmc.Q.derivatives}}
//'   
//' @export
// [[Rcpp::export(name = "phylo.Q.derivatives")]]
Map<std::string, double> phylo_Q_derivatives(
    const Rcpp::List& tree, const Rcpp::List& subst_mod,
    const std::string& param_name, int max_order,
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
  if (max_order < 0) Rcpp::stop("'max_order' cannot be less than 0.");

  arma::imat edge = tree["edge"];
  const Vector<std::string>& tip_labels = tree["tip.label"];
  int num_int_nodes = tree["Nnode"];
  const arma::vec& edge_lengths = tree["edge.length"];
  const Vector<std::string>& states = subst_mod["states"];
  const arma::mat& Q = subst_mod["Q"];
  const arma::vec& pi = subst_mod["pi"];
  std::string d_param_name = "d_" + param_name;

  if (arma::find(edge.col(0) == tip_labels.size() + 1).eval().n_elem > 2)
    Rcpp::stop("'tree' must be a rooted tree.");
  if (!subst_mod.containsElementNamed(d_param_name.c_str()))
    Rcpp::stop("'param_name' is not a valid 'subst_mod' parameter name.");
  if (tip_states.size() != tip_labels.size())
    Rcpp::stop("'tip_states' must be compatible with 'tree'.");

  const arma::mat& dQ = subst_mod[d_param_name];
  VectorVector<int> edge_sets(1);
  edge_sets[0].reserve(edge.n_rows);
  for (std::size_t edge_label = 1; edge_label <= edge.n_rows; ++edge_label) {
    edge_sets[0].push_back(edge_label);
  }

  arma::ivec tip_data(tip_states.size(), arma::fill::none);
  for (std::size_t i = 0; i < tip_states.size(); ++i) {
    auto find_it = std::find(states.begin(), states.end(), tip_states[i]);
    if (find_it != states.end()) {
      tip_data(i) = find_it - states.begin();
    } else {
      Rcpp::stop("'tip_states' must contain valid states.");
    }
  }

  return phylo_moments_derivatives(edge, tip_labels, num_int_nodes,
                                   edge_lengths, Q, dQ, pi, edge_sets,
                                   max_order, Mode::Q_DERIVATIVES, tip_data);
}


//' Phylogenetic likelihood derivatives with respect to branch lengths
//' 
//' Computes the phylogenetic likelihood derivatives with respect to the branch 
//' lengths for all orders less than or equal to the user-specified maximum 
//' order of interest.
//' 
//' More information on S3 objects of class \code{"phylo"} is found at 
//' \url{http://ape-package.ird.fr/misc/FormatTreeR_24Oct2012.pdf}.
//' 
//' The zeroth-order likelihood derivative is defined to be the phylogenetic 
//' likelihood.
//' 
//' @param tree An S3 object of class \code{"phylo"}.
//' @param subst_mod An S3 object of class \code{"substitution.model"}.
//' @param max_order An integer specifying the maximum derivative order of 
//'   interest.
//' @param tip_states A character vector of observed tip states.
//'   
//' @return A named numeric vector that holds the computed likelihood 
//'   derivatives.  Each element name is formed by concatenating strings of the 
//'   form \code{"(E):O"} with \code{"--"} separators, where \code{E} denotes an
//'   edge index, \code{O} represents the derivative order associated with the 
//'   edge index \code{E}, and \code{"--"} symbolizes the partial 
//'   differentiation operator with respect to two branch lengths.  However, the
//'   empty string (\code{""}) is the element name of the zeroth-order 
//'   likelihood derivative.
//'   
//' @references Kenney T and Gu H (2012) \dQuote{Hessian Calculation for 
//'   Phylogenetic Likelihood based on the Pruning Algorithm and its 
//'   Applications}, \emph{Statistical Applications in Genetics and Molecular 
//'   Biology}, 11(4).
//'   
//'   OUR PAPER!!
//'   
//' @seealso \code{\link{phylo.Q.derivatives}}, \code{\link{ctmc.t.derivatives}}
//'   
//' @export
// [[Rcpp::export(name = "phylo.t.derivatives")]]
Map<std::string, double> phylo_t_derivatives(
    const Rcpp::List& tree, const Rcpp::List& subst_mod, int max_order,
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
  if (max_order < 0) Rcpp::stop("'max_order' cannot be less than 0.");

  arma::imat edge = tree["edge"];
  const Vector<std::string>& tip_labels = tree["tip.label"];
  int num_int_nodes = tree["Nnode"];
  const arma::vec& edge_lengths = tree["edge.length"];
  const Vector<std::string>& states = subst_mod["states"];
  const arma::mat& Q = subst_mod["Q"];
  const arma::vec& pi = subst_mod["pi"];

  if (arma::find(edge.col(0) == tip_labels.size() + 1).eval().n_elem > 2)
    Rcpp::stop("'tree' must be a rooted tree.");
  if (tip_states.size() != tip_labels.size())
    Rcpp::stop("'tip_states' must be compatible with 'tree'.");

  VectorVector<int> edge_sets(edge.n_rows);
  for (std::size_t edge_label = 1; edge_label <= edge.n_rows; ++edge_label) {
    edge_sets[edge_label - 1].reserve(1);
    edge_sets[edge_label - 1].push_back(edge_label);
  }

  arma::ivec tip_data(tip_states.size(), arma::fill::none);
  for (std::size_t i = 0; i < tip_states.size(); ++i) {
    auto find_it = std::find(states.begin(), states.end(), tip_states[i]);
    if (find_it != states.end()) {
      tip_data(i) = find_it - states.begin();
    } else {
      Rcpp::stop("'tip_states' must contain valid states.");
    }
  }

  return phylo_moments_derivatives(edge, tip_labels, num_int_nodes,
                                   edge_lengths, Q, arma::mat(), pi, edge_sets,
                                   max_order, Mode::T_DERIVATIVES, tip_data);
}