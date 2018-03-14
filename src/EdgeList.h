#ifndef PHYLOMD_EDGELIST_H_
#define PHYLOMD_EDGELIST_H_

#include <tuple>
#include <utility>

#include <RcppArmadillo.h>
#include "ListID.h"
#include "PartitionSet.h"
#include "phylomd_types.h"


//
// This header file defines the EdgeList class and declares the EdgeList-related
// member functions.
//


class EdgeList {
 private:
  const ListID& id_;
  // This is a vector of 3-tuples that summarizes the recurrence relation
  // between this edge list and the node lists.
  //
  // For any given edge, the first tuple element represents a possible
  // assignment of an ID element, the second element contains the corresponding
  // node list index associated with the child node, and the third element
  // stores the related counting coefficient.
  Vector<std::tuple<const ListIDElement&, int, int>> recursion_info_;
  Map<Ref<const PartitionSet>, std::pair<int, int>> recursion_info_inds_;
  arma::mat elems_;

  void init_recursion_info_aux(
      const Vector<ListID>& ids, const arma::mat& choose,
      Vector<std::tuple<const ListIDElement&, int, int>>& recursion_info) const;
  Vector<std::tuple<const ListIDElement&, int, int>> init_recursion_info(
      const Vector<ListID>& ids, const arma::mat& choose) const;
  Map<Ref<const PartitionSet>, std::pair<int, int>> init_recursion_info_inds()
      const;

 public:
  EdgeList(const ListID& id, const Vector<ListID>& ids, const arma::mat& choose,
           int nrow, int ncol)
      : id_(id),
        recursion_info_(init_recursion_info(ids, choose)),
        recursion_info_inds_(init_recursion_info_inds()),
        elems_(arma::mat(nrow, ncol, arma::fill::zeros)) {}

  arma::mat& elems() { return elems_; }
  const arma::mat& elems() const { return elems_; }

  const ListID& id() const { return id_; }
  const Vector<std::tuple<const ListIDElement&, int, int>>& recursion_info()
      const {
    return recursion_info_;
  }
  const Map<Ref<const PartitionSet>, std::pair<int, int>>& recursion_info_inds()
      const {
    return recursion_info_inds_;
  }
};


#endif  // PHYLOMD_EDGELIST_H_