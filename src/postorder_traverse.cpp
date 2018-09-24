#include "postorder_traverse.h"

#include <cstddef>
#include <tuple>
#include <utility>

#include "ListID.h"


//
// This source file defines the post-order tree traversal functions.
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
    int edge_ind, int child_node_ind,
    Map<int, const PartitionSet&>::const_iterator edge_pset_it,
    Map<int, const PartitionSet&>::const_iterator edge_psets_end_it,
    const arma::cube& edge_mds, const Vector<NodeList>& nlists,
    Vector<EdgeList>& elists) {
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


void postorder_traverse(const arma::imat& edge, int num_edges,
                        int num_term_nodes, int root_node_ind,
                        const arma::vec& edge_lengths, const arma::mat& Q,
                        const arma::mat& B,
                        const Map<int, const PartitionSet&>& edge_psets,
                        int max_order, Mode mode, const arma::ivec& tip_data,
                        Vector<arma::cube>& edges_mds, Vector<NodeList>& nlists,
                        Vector<EdgeList>& elists) {
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

      // Do we have observed data at this terminal node?
      if (tip_data(child_node_ind) != -1) {
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
    edges_mds[edge_ind] = ctmc_moments_derivatives(edge_lengths(edge_ind), Q, B,
                                                   edge_max_order, mode);

    // Update the edge lists at the current edge.
    update_edge_lists(edge_ind, child_node_ind, edge_pset_it, edge_psets.end(),
                      edges_mds[edge_ind], nlists, elists);
  }

  // Update the node lists at the root node.
  arma::uvec root_edge_inds = arma::find(edge.col(0) == root_node_ind);
  update_node_lists(root_node_ind, root_edge_inds, elists, nlists);
}