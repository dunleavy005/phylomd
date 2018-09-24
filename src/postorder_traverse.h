#ifndef PHYLOMD_POSTORDER_TRAVERSE_H_
#define PHYLOMD_POSTORDER_TRAVERSE_H_

#include <RcppArmadillo.h>
#include "EdgeList.h"
#include "NodeList.h"
#include "PartitionSet.h"
#include "ctmc_moments_derivatives.h"
#include "phylomd_types.h"


//
// This header file declares the post-order tree traversal functions.
//


void update_node_lists(int child_node_ind, const arma::uvec& child_edge_inds,
                       const Vector<EdgeList>& elists,
                       Vector<NodeList>& nlists);


void update_edge_lists(
    int edge_ind, int child_node_ind,
    Map<int, const PartitionSet&>::const_iterator edge_pset_it,
    Map<int, const PartitionSet&>::const_iterator edge_psets_end_it,
    const arma::cube& edge_mds, const Vector<NodeList>& nlists,
    Vector<EdgeList>& elists);


void postorder_traverse(const arma::imat& edge, int num_edges,
                        int num_term_nodes, int root_node_ind,
                        const arma::vec& edge_lengths, const arma::mat& Q,
                        const arma::mat& B,
                        const Map<int, const PartitionSet&>& edge_psets,
                        int max_order, Mode mode, const arma::ivec& tip_data,
                        Vector<arma::cube>& edges_mds, Vector<NodeList>& nlists,
                        Vector<EdgeList>& elists);


#endif  // PHYLOMD_POSTORDER_TRAVERSE_H_