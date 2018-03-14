#ifndef PHYLOMD_NODELIST_H_
#define PHYLOMD_NODELIST_H_

#include <tuple>

#include <RcppArmadillo.h>
#include "ListID.h"
#include "phylomd_types.h"


//
// This header file defines the NodeList class and declares the NodeList-related
// member functions.
//


class NodeList {
 private:
  const ListID& id_;
  // This is a vector of 3-tuples that summarizes the recurrence relation
  // between this node list and the edge lists.
  //
  // For any given node, the first tuple element contains the edge list index
  // associated with the left child edge, the second element represents the edge
  // list index associated with the right child edge, and the third element
  // stores the related counting coefficient.
  Vector<std::tuple<int, int, int>> recursion_info_;
  arma::mat elems_;

  void init_recursion_info_aux(
      const Vector<ListID>& ids, const arma::mat& choose,
      ListID::const_iterator curr_it,
      const Vector<ListIDElement>& left_elist_id_elems,
      Vector<std::tuple<int, int, int>>& recursion_info) const;
  Vector<std::tuple<int, int, int>> init_recursion_info(
      const Vector<ListID>& ids, const arma::mat& choose) const;

 public:
  NodeList(const ListID& id, const Vector<ListID>& ids, const arma::mat& choose,
           int nrow, int ncol)
      : id_(id),
        recursion_info_(init_recursion_info(ids, choose)),
        elems_(arma::mat(nrow, ncol, arma::fill::zeros)) {}

  arma::mat& elems() { return elems_; }
  const arma::mat& elems() const { return elems_; }

  const ListID& id() const { return id_; }
  const Vector<std::tuple<int, int, int>>& recursion_info() const {
    return recursion_info_;
  }
};


#endif  // PHYLOMD_NODELIST_H_